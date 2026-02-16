import shutil
import subprocess
from pathlib import Path
from typing import Optional, Union, Literal
import os

from docksuitex.utils.viewer import view_molecule

MGLTOOLS_PATH = (Path(__file__).parent / "bin" / "mgltools").resolve()
MGL_PYTHON_EXE = (MGLTOOLS_PATH / "python.exe").resolve()
PREPARE_LIGAND_SCRIPT = (MGLTOOLS_PATH / "Lib" / "site-packages" / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py").resolve()
OBABEL_EXE = (Path(__file__).parent / "bin" / "obabel" / "obabel.exe").resolve()



class Ligand:
    """Ligand structure preparation for molecular docking.

    This module provides automated ligand preparation using Open Babel for format
    conversion and energy minimization, combined with AutoDockTools (MGLTools) for
    PDBQT generation required for molecular docking.

        The preparation workflow:
            1. Converting various input formats to MOL2, 3D coordinate generation if needed (using Open Babel)
            2. Optional energy minimization with forcefields: MMFF94, MMFF94S UFF, GAFF (using Open Babel)
            3. Optional water molecule removal (using Open Babel)
            4. Optional hydrogen addition and gasteiger charge assignment (using AutoDockTools)
            5. Converting the MOL2 file to PDBQT format (using AutoDockTools)

        Supported Input Formats:
            DockSuiteX supports a variety of ligand input formats including Tripos MOL2 (.mol2), 
            Structure Data File (.sdf), PDB (.pdb), MDL Molfile (.mol), and SMILES (.smi).

    Note:
        Intermediate files (e.g., MOL2) are saved in an
        ``intermediate_ligands/`` subfolder within the output directory.
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
            preserve_charge_types (list[str], optional): Atom types (e.g.,["Zn", "Fe"]) whose charges are preserved; 
                others get Gasteiger charges; ignored if add_charges=False. Defaults to None.

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
        Handles ligand preparation for docking using Open Babel and AutoDockTools (ADT). 
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
        # Resolve output directory early so intermediate files go there
        save_to = Path(save_to).expanduser().resolve()
        if not save_to.suffix:
            save_to = save_to / f"{self.file_path.stem}.pdbqt"
        save_to.parent.mkdir(parents=True, exist_ok=True)
        output_dir = save_to.parent

        # Intermediate files go in a subfolder of the output directory
        intermediates_dir = output_dir / "intermediate_ligands"
        intermediates_dir.mkdir(parents=True, exist_ok=True)

        # === Step 1: Convert + Gen3D + remove water + Minimize to MOL2 ===
        self.mol2_path = intermediates_dir / f"{self.file_path.stem}.mol2"
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



        # Fix for missing OpenBabel data in bundled app: use MGLTools' OB data
        env = os.environ.copy()
        mgl_ob_data = MGLTOOLS_PATH / "OpenBabel-2.3.2" / "data"
        if mgl_ob_data.exists():
            env["BABEL_DATADIR"] = str(mgl_ob_data)

        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if result.returncode != 0:
            raise RuntimeError(f"❌ OpenBabel failed:\n{result.stderr}")

        # === Step 2: MGLTools to PDBQT ===

        mgl_cmd = [
            str(MGL_PYTHON_EXE), str(PREPARE_LIGAND_SCRIPT),
            "-l", str(self.mol2_path), "-o", str(save_to),
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
            cwd=str(self.mol2_path.parent),
        )

        if result.returncode != 0:
            raise RuntimeError(f"❌ MGLTools ligand preparation failed:\n{result.stderr}")

        self.pdbqt_path = save_to

        print(f"✅ Ligand prepared successfully: {self.pdbqt_path}")
        return self.pdbqt_path

    def view_molecule(self):
        """Visualize the ligand structure in a Jupyter notebook.

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
