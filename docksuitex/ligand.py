"""Ligand structure preparation for molecular docking.

This module provides automated ligand preparation using Open Babel for format
conversion and energy minimization, combined with AutoDockTools (MGLTools) for
PDBQT generation required for molecular docking.

The preparation workflow:
    1. Format conversion to MOL2 using Open Babel
    2. 3D coordinate generation (if needed)
    3. Optional energy minimization with forcefields (MMFF94, UFF, GAFF)
    4. Water molecule removal
    5. Hydrogen addition and charge assignment
    6. PDBQT format conversion using AutoDockTools

Example:
    Basic ligand preparation::

        from docksuitex import Ligand

        # Prepare ligand from SDF file
        ligand = Ligand(
            input="ligand.sdf",
            minimize="mmff94",
            remove_water=True
        )
        
        # Generate PDBQT file
        pdbqt_path = ligand.prepare(save_to="prepared_ligand.pdbqt")
        
        # Visualize in Jupyter
        ligand.view_molecule()

    Preparation from SMILES string::

        ligand = Ligand(
            input="aspirin.smi",
            minimize="uff",
            add_hydrogens=True
        )
        pdbqt_path = ligand.prepare()

"""

import shutil
import subprocess
from pathlib import Path
from typing import Optional, Union, Literal
import uuid

from docksuitex.utils.viewer import view_molecule

# Configuration - Paths to bundled executables
MGLTOOLS_PATH = (Path(__file__).parent / "bin" / "mgltools").resolve()
MGL_PYTHON_EXE = (MGLTOOLS_PATH / "python.exe").resolve()
PREPARE_LIGAND_SCRIPT = (MGLTOOLS_PATH / "Lib" / "site-packages" /
                         "AutoDockTools" / "Utilities24" / "prepare_ligand4.py").resolve()
OBABEL_EXE = (Path(__file__).parent / "bin" /
              "obabel" / "obabel.exe").resolve()

TEMP_DIR = (Path.cwd() / "temp").resolve()
TEMP_DIR.mkdir(exist_ok=True)



class Ligand:
    """Automated ligand structure preparation for molecular docking.

    This class provides a complete workflow for preparing ligand structures
    for docking simulations. It integrates Open Babel for format conversion
    and energy minimization, and AutoDockTools for PDBQT generation.

    The preparation process handles:
        - Format conversion from various structure formats to MOL2
        - 3D coordinate generation for 2D structures or SMILES
        - Energy minimization with multiple forcefields
        - Water molecule removal
        - Hydrogen addition (polar or all)
        - Gasteiger charge assignment
        - PDBQT format conversion for docking


    Example:
        Standard ligand preparation::

            from docksuitex import Ligand

            ligand = Ligand("compound.sdf", minimize="mmff94")
            pdbqt_file = ligand.prepare(save_to="ligand.pdbqt")

        Preparation from SMILES::

            ligand = Ligand("molecule.smi", minimize="uff")
            pdbqt_file = ligand.prepare()

    Note:
        Temporary files are created in a `temp/Ligands/` directory and are
        not automatically cleaned up. The final PDBQT file is saved to the
        specified location.
    """

    SUPPORTED_INPUTS = {"mol2", "sdf", "pdb", "mol", "smi"}
    SUPPORTED_FORCEFIELDS = {"mmff94", "mmff94s", "uff", "gaff"}

    def __init__(
        self,
        input: Union[str, Path],
        minimize: Optional[str] = None,
        remove_water: bool = True,
        add_hydrogens: bool = True,
        add_charges: bool = True,
        preserve_charge_types: Optional[list[str]] = None,
    ):
        """Initialize a Ligand object with a given input file and preparation parameters.

        Args:
            input (str | Path): Path to the ligand input file.
            minimize (str, optional): Forcefield for energy minimization ("mmff94", "mmff94s", "uff", "gaff").
            remove_water (bool, optional): Remove water molecules. Defaults to True.
            add_hydrogens (bool, optional): Add polar hydrogens. Defaults to True.
            add_charges (bool, optional): Assign Gasteiger charges. Defaults to True.
            preserve_charge_types (list[str], optional): Atom types to preserve charges for. Defaults to None.

        Raises:
            FileNotFoundError: If the input file does not exist.
            ValueError: If the file extension is unsupported.
        """
        self.file_path = Path(input).resolve()
        self.mol2_path: Optional[Path] = None
        self.pdbqt_path: Optional[Path] = None

        if not self.file_path.is_file():
            raise FileNotFoundError(
                f"❌ Ligand file not found: {self.file_path}")

        ext = self.file_path.suffix.lower().lstrip(".")
        if ext not in self.SUPPORTED_INPUTS:
            raise ValueError(
                f"❌ Unsupported file format '.{ext}'. Supported formats: {self.SUPPORTED_INPUTS}")
        self.input_format = ext

        # Store preparation parameters
        self.minimize = minimize
        self.remove_water = remove_water
        self.add_hydrogens = add_hydrogens
        self.add_charges = add_charges
        self.preserve_charge_types = preserve_charge_types

    def prepare(self, save_to: Union[str, Path] = ".") -> Path:
        """
        Prepare the ligand by converting to MOL2, optionally minimizing energy, 
        and generating a final PDBQT file using AutoDockTools (from MGLTools).
        Saves the prepared PDBQT file to the specified location.

        Args:
            save_to (str | Path, optional): Destination file or directory.
                - If directory: saves with the original filename.
                - If file path: saves with the given name.
                Defaults to current directory.

        Returns:
            Path: Path to the saved PDBQT file.

        Raises:
            ValueError: If an unsupported forcefield or input format is provided.
            RuntimeError: If AutoDockTools fails to generate the PDBQT file.
        """
        # Create a unique temp directory per object
        self.temp_dir = TEMP_DIR / "Ligands" / f"{self.file_path.stem}_{uuid.uuid4().hex[:8]}"
        self.temp_dir.mkdir(parents=True, exist_ok=True)

        # === Step 1: Convert + Gen3D + Minimize to MOL2 ===
        self.mol2_path = self.temp_dir / f"{self.file_path.stem}.mol2"
        cmd = [
            str(OBABEL_EXE), "-i", self.input_format, str(self.file_path),
            "-o", "mol2", "-O", str(self.mol2_path),
            "--gen3d"
        ]

        # Universal water removal: works for PDB (HOH) + all other formats ([#8H2])
        if self.remove_water:
            cmd += ["--delete", "HOH", "--delete", "[#8H2]"]

        if self.minimize:
            forcefield = self.minimize.lower()
            if forcefield not in self.SUPPORTED_FORCEFIELDS:
                raise ValueError(
                    f"❌ Unsupported forcefield '{forcefield}'. Supported: {self.SUPPORTED_FORCEFIELDS}")
            cmd += ["--minimize", "--ff", forcefield]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"❌ OpenBabel failed:\n{result.stderr}")

        # === Step 2: MGLTools to PDBQT ===

        save_to = Path(save_to).expanduser().resolve()
        
        # treat as file only if it has a suffix (e.g., .pdbqt)
        if not save_to.suffix:
            save_to = save_to / f"{self.file_path.stem}.pdbqt"
        
        save_to.parent.mkdir(parents=True, exist_ok=True)

        mgl_cmd = [
            str(MGL_PYTHON_EXE), str(PREPARE_LIGAND_SCRIPT),
            "-l", self.mol2_path, "-o", save_to,
            "-U", "nphs_lps"
        ]
        # ADT prepare_ligand4.py doesn't have -U waters flag, remove water is handled by obabel

        if self.add_hydrogens:
            mgl_cmd += ["-A", "hydrogens"]
        else:
            mgl_cmd += ["-A", "None"]

        # Charge options
        if not self.add_charges:
            mgl_cmd += ["-C"]  # preserve all charges
        elif self.preserve_charge_types:
            for atom_type in self.preserve_charge_types:
                mgl_cmd += ["-p", atom_type]

        result = subprocess.run(
            mgl_cmd,
            text=True,
            capture_output=True,
            cwd=self.temp_dir
        )

        if result.returncode != 0:
            raise RuntimeError(f"❌ MGLTools ligand preparation failed:\n{result.stderr}")

        self.pdbqt_path = save_to

        # # Save logic (merged from save_pdbqt)
        # save_to = Path(save_to).expanduser().resolve()

        # # treat as file only if it has a suffix (e.g., .pdbqt)
        # if not save_to.suffix:
        #     save_to = save_to / self.pdbqt_path.name

        # save_to.parent.mkdir(parents=True, exist_ok=True)
        # shutil.copy2(self.pdbqt_path, save_to)
        print(f"✅ Ligand prepared successfully: {self.pdbqt_path}")
        return self.pdbqt_path

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
