import os
import subprocess
from pathlib import Path
from typing import Union, Optional, Literal
import shutil

from pdbfixer import PDBFixer
from openmm.app import PDBFile

import uuid

from .utils.viewer import view_molecule


# Constants
MGLTOOLS_PATH = (Path(__file__).parent / "bin" / "mgltools").resolve()
MGL_PYTHON_EXE = (MGLTOOLS_PATH / "python.exe").resolve()
PREPARE_RECEPTOR_SCRIPT = (
    MGLTOOLS_PATH / "Lib" / "site-packages" /
    "AutoDockTools" / "Utilities24" / "prepare_receptor4.py"
).resolve()

OBABEL_EXE = (Path(__file__).parent / "bin" /
              "obabel" / "obabel.exe").resolve()
TEMP_DIR = (Path(__file__).parent / "temp").resolve()
TEMP_DIR.mkdir(parents=True, exist_ok=True)


class Protein:
    """Handles protein preparation for docking using PDBFixer, Open Babel, and AutoDockTools.

    This class automates the preprocessing of protein structures by:
    - Converting various input formats to PDB using Open Babel
    - Optionally fixing missing residues and atoms using PDBFixer
    - Optionally removing water molecules, adding hydrogens and charges
    - Converting the fixed PDB file to PDBQT format using AutoDockTools
    """

    SUPPORTED_INPUTS = {".pdb", ".mol2", ".sdf",
                        ".pdbqt", ".cif", ".ent", ".xyz"}

    def __init__(self, file_path: Union[str, Path]):
        """Initialize a Protein object with a given file path.

        Args:
            file_path (str | Path): Path to the protein input file.

        Raises:
            FileNotFoundError: If the provided file does not exist.
            ValueError: If the file format is not supported.
        """
        self.file_path = Path(file_path).resolve()
        self.pdb_path: Optional[Path] = None
        self.pdbqt_path: Optional[Path] = None

        if not self.file_path.is_file():
            raise FileNotFoundError(
                f"❌ Protein file not found: {self.file_path}")

        self.ext = self.file_path.suffix.lower()
        if self.ext not in self.SUPPORTED_INPUTS:
            raise ValueError(
                f"❌ Unsupported file format '{self.ext}'. Supported formats: {self.SUPPORTED_INPUTS}")

    def prepare(
        self,
        fix_pdb: bool = True,
        remove_heterogens: bool = True,
        add_hydrogens: bool = True,
        remove_water: bool = True,
        add_charges: bool = True,
        preserve_charge_types: Optional[list[str]] = None,

    ) -> None:
        """
        Handles protein preparation for docking using PDBFixer, Open Babel, and AutoDockTools (ADT).

        Args:
            fix_pdb (bool, optional): Fix missing residues, atoms, replace non-standard residues using PDBFixer. Defaults to True.
            remove_heterogens (bool, optional): Remove ligands/heterogens during fixing. Defaults to True.
            add_hydrogens (bool, optional): Add hydrogens in PDBQT preparation. Defaults to True.
            remove_water (bool, optional): Remove water molecules in PDBQT preparation. Defaults to True.
            add_charges (bool, optional): Assign Gasteiger charges during PDBQT generation. If False, all input charges are preserved. Default is True.
            preserve_charge_types (list[str], optional): Atom types (e.g., ["Zn", "Fe"])
                whose charges should be preserved. Ignored if `add_charges=False`. Defaults to None.

        Raises:
            RuntimeError: If Open Babel or AutoDockTools commands fail.
        """
        # Create a unique temp directory per object
        self.temp_dir = TEMP_DIR / "Proteins" / f"{self.file_path.stem}_{uuid.uuid4().hex[:8]}"
        self.temp_dir.mkdir(parents=True, exist_ok=True)

        # Convert to .pdb if needed
        if self.ext != ".pdb":
            self.pdb_path = self.temp_dir / f"{self.file_path.stem}.pdb"
            cmd = [str(OBABEL_EXE), str(self.file_path),
                   "-O", str(self.pdb_path), "--gen3d"]
            subprocess.run(cmd, check=True)
        else:
            self.pdb_path = self.file_path

        # Fix structure using PDBFixer
        fixer = PDBFixer(filename=str(self.pdb_path))
        if fix_pdb:
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
        if remove_heterogens:
            fixer.removeHeterogens(keepWater=True)


        # Save fixed PDB
        fixed_pdb_path = self.temp_dir / f"{self.file_path.stem}_fixed.pdb"
        with open(fixed_pdb_path, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        self.pdb_path = fixed_pdb_path

        # Convert to PDBQT using AutoDockTools
        output_pdbqt = self.temp_dir / f"{self.file_path.stem}.pdbqt"
        U_flag = "nphs_lps_waters" if remove_water else "nphs_lps"
        cmd = [
            str(MGL_PYTHON_EXE),
            str(PREPARE_RECEPTOR_SCRIPT),
            "-r", str(self.pdb_path),
            "-o", str(output_pdbqt),
            "-U", U_flag
        ]
        if add_hydrogens:
            cmd += ["-A", "hydrogens"]

        # Control charges
        if not add_charges:
            cmd += ["-C"]  # disable Gasteiger charges
        elif preserve_charge_types:
            for atom in preserve_charge_types:
                cmd += ["-p", atom]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"❌ Error preparing PDBQT:\n{result.stderr}")

        self.pdbqt_path = output_pdbqt

    



    def view_molecule(self):
        """Visualize the ligand structure in a Jupyter notebook.

        Uses nglview to render either the prepared or original file.

        Returns:
            object: An nglview.NGLWidget object for rendering.

        Raises:
            FileNotFoundError: If neither prepared nor input file exists.
        """
        path = Path(self.pdbqt_path if self.pdbqt_path else self.file_path).resolve()
        return view_molecule(file_path=path)





    def save_pdbqt(self, save_to: Union[str, Path] = ".") -> Path:
        """Save the prepared PDBQT file to a given location.

        Args:
            save_to (str | Path, optional): Destination path.
                - If a directory: file will be saved with original name.
                - If a file path: saved with the given name.

        Returns:
            Path: Final saved file path.

        Raises:
            RuntimeError: If `prepare()` has not been run or PDBQT file is missing.
        """
        if self.pdbqt_path is None or not self.pdbqt_path.exists():
            raise RuntimeError("❌ Protein not prepared. Run prepare() first.")

        save_to = Path(save_to).expanduser().resolve()

        # treat as file only if it has a suffix (e.g., .pdbqt)
        if not save_to.suffix:
            save_to = save_to / self.pdbqt_path.name

        save_to.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(self.pdbqt_path, save_to)
        return save_to

