import os
import subprocess
from pathlib import Path
from typing import Union, Optional, Literal
import shutil

from pdbfixer import PDBFixer
from openmm.app import PDBFile

from .utils.viewer import view_molecule


# Constants - Paths to bundled executables
MGLTOOLS_PATH = (Path(__file__).parent / "bin" / "mgltools").resolve()
MGL_PYTHON_EXE = (MGLTOOLS_PATH / "python.exe").resolve()
PREPARE_RECEPTOR_SCRIPT = (MGLTOOLS_PATH / "Lib" / "site-packages" /"AutoDockTools" / "Utilities24" / "prepare_receptor4.py").resolve()

OBABEL_EXE = (Path(__file__).parent / "bin" /"obabel" / "obabel.exe").resolve()



class Protein:
    """Protein structure preparation for molecular docking.

        This module provides automated protein preparation using a combination of
        PDBFixer, Open Babel, and AutoDockTools (MGLTools). It handles format conversion,
        structure fixing, and PDBQT generation required for molecular docking.

        The preparation workflow:
            1. Converting various input formats to PDB (using Open Babel)
            2. Optionally fixing missing residues/atoms and nonstandard residues (using PDBFixer)
            3. Optionally removing heterogens(using PDBFixer) and removing water molecules (using AutoDockTools).
            4. Optionally adding hydrogens and gasteiger charges (using AutoDockTools).
            5. Converting the PDB file to PDBQT format (using AutoDockTools).

        Note:
            Intermediate files (converted PDB, fixed PDB) are saved in an
            ``intermediate_proteins/`` subfolder within the output directory.
    """
    SUPPORTED_INPUTS = {".pdb", ".mol2", ".sdf",
                        ".pdbqt", ".cif", ".ent", ".xyz"}

    def __init__(
        self,
        input: Union[str, Path],
        fix_pdb: bool = True,
        remove_heterogens: bool = True,
        remove_water: bool = True,
        add_hydrogens: bool = True,
        add_charges: bool = True,
        preserve_charge_types: Optional[list[str]] = None,
    ):

        """   
        Initialize a Protein object with a given file path and preparation parameters.

        Args:
            input (str | Path): Path to the protein input file.
            fix_pdb (bool, optional): Fix missing residues/atoms. Defaults to True.
            remove_heterogens (bool, optional): Remove ligands/heterogens. Defaults to True.
            remove_water (bool, optional): Remove water molecules. Defaults to True.
            add_hydrogens (bool, optional): Add hydrogens. Defaults to True.
            add_charges (bool, optional): Assign Gasteiger charges. Defaults to True.
            preserve_charge_types (list[str], optional): Atom types (e.g.,["Zn", "Fe"]) whose charges are preserved; 
                others get Gasteiger charges; ignored if add_charges=False. Defaults to None.

        Raises:
            FileNotFoundError: If the provided file does not exist.
            ValueError: If the file format is not supported.

        """
        self.file_path = Path(input).resolve()
        self.pdb_path: Optional[Path] = None
        self.fixed_pdb_path: Optional[Path] = None
        self.pdbqt_path: Optional[Path] = None

        if not self.file_path.is_file():
            raise FileNotFoundError(
                f"❌ Protein file not found: {self.file_path}")

        self.ext = self.file_path.suffix.lower()
        if self.ext not in self.SUPPORTED_INPUTS:
            raise ValueError(
                f"❌ Unsupported file format '{self.ext}'. Supported formats: {self.SUPPORTED_INPUTS}")

        # Store preparation parameters
        self.fix_pdb = fix_pdb
        self.remove_heterogens = remove_heterogens
        self.add_hydrogens = add_hydrogens
        self.remove_water = remove_water
        self.add_charges = add_charges
        self.preserve_charge_types = preserve_charge_types

    def prepare(self, save_to: Union[str, Path] = ".") -> Path:
        """
        Handles protein preparation for docking using PDBFixer, Open Babel, and AutoDockTools (ADT).
        Saves the prepared PDBQT file to the specified location.

        Args:
            save_to (str | Path, optional): Destination path for the PDBQT file.
                - If a directory: file will be saved with original name.
                - If a file path: saved with the given name.
                Defaults to current directory.

        Returns:
            Path: Path to the saved PDBQT file.

        Raises:
            RuntimeError: If Open Babel or AutoDockTools commands fail.
        """
        # Resolve output directory early so intermediate files go there
        save_to = Path(save_to).expanduser().resolve()
        if not save_to.suffix:
            save_to = save_to / f"{self.file_path.stem}.pdbqt"
        save_to.parent.mkdir(parents=True, exist_ok=True)
        output_dir = save_to.parent

        # Intermediate files go in a subfolder of the output directory
        intermediates_dir = output_dir / "intermediate_proteins"
        intermediates_dir.mkdir(parents=True, exist_ok=True)

        # Convert to .pdb if needed
        if self.ext != ".pdb":
            self.pdb_path = intermediates_dir / f"{self.file_path.stem}.pdb"
            cmd = [str(OBABEL_EXE), str(self.file_path),"-O", str(self.pdb_path), "--gen3d"]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"❌ Open Babel conversion failed:\n{result.stderr}")

        else:
            self.pdb_path = self.file_path

        # Fix structure using PDBFixer
        fixer = PDBFixer(filename=str(self.pdb_path))
        if self.fix_pdb:
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
        if self.remove_heterogens:
            fixer.removeHeterogens(keepWater=True)


        # Save fixed PDB
        fixed_pdb_path = intermediates_dir / f"{self.file_path.stem}_fixed.pdb"
        with open(fixed_pdb_path, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        self.fixed_pdb_path = fixed_pdb_path



        # Convert to PDBQT using AutoDockTools

        U_flag = "nphs_lps_waters" if self.remove_water else "nphs_lps"
        cmd = [
            str(MGL_PYTHON_EXE),
            str(PREPARE_RECEPTOR_SCRIPT),
            "-r", str(self.fixed_pdb_path),
            "-o", str(save_to),
            "-U", U_flag
        ]
        if self.add_hydrogens:
            cmd += ["-A", "hydrogens"]

        # Control charges
        if not self.add_charges:
            cmd += ["-C"]  # disable Gasteiger charges
        elif self.preserve_charge_types:
            for atom in self.preserve_charge_types:
                cmd += ["-p", atom]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"❌ Error preparing PDBQT:\n{result.stderr}")

        self.pdbqt_path = save_to

        print(f"✅ Protein prepared successfully: {self.pdbqt_path}")
        return self.pdbqt_path

    def view_molecule(self):
        """Visualize the protein structure in a Jupyter notebook.

        Uses NGLView to render either the prepared PDBQT file or the original
        input file in an interactive 3D viewer.

        Returns:
            nglview.NGLWidget: Interactive 3D molecular viewer widget.

        Raises:
            FileNotFoundError: If neither prepared nor input file exists.
            ImportError: If nglview is not installed.

        Note:
            This method requires a Jupyter Notebook/Lab environment and the
            nglview package.

        """
        path = Path(self.pdbqt_path if self.pdbqt_path else self.file_path).resolve()
        return view_molecule(file_path=path)


