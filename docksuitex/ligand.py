import shutil
import subprocess
from pathlib import Path
from typing import Optional, Union, Literal
import uuid

from docksuitex.utils.viewer import view_molecule

# === CONFIGURATION ===
MGLTOOLS_PATH = (Path(__file__).parent / "bin" / "mgltools").resolve()
MGL_PYTHON_EXE = (MGLTOOLS_PATH / "python.exe").resolve()
PREPARE_LIGAND_SCRIPT = (MGLTOOLS_PATH / "Lib" / "site-packages" /
                         "AutoDockTools" / "Utilities24" / "prepare_ligand4.py").resolve()
OBABEL_EXE = (Path(__file__).parent / "bin" /
              "obabel" / "obabel.exe").resolve()

TEMP_DIR = (Path(__file__).parent / "temp").resolve()
TEMP_DIR.mkdir(exist_ok=True)


class Ligand:
    """Ligand preparation pipeline using Open Babel and MGLTools.

    This class automates preprocessing of ligands by:
    - Converting input formats to MOL2
    - Optionally minimizing energy with forcefields
    - Converting MOL2 to PDBQT using AutoDockTools
    """

    SUPPORTED_INPUTS = {"mol2", "sdf", "pdb", "mol", "smi"}
    SUPPORTED_FORCEFIELDS = {"mmff94", "mmff94s", "uff", "gaff"}

    def __init__(self, file_path: Union[str, Path]):
        """Initialize a Ligand object with a given input file.

        Args:
            file_path (str | Path): Path to the ligand input file.

        Raises:
            FileNotFoundError: If the input file does not exist.
            ValueError: If the file extension is unsupported.
        """
        self.file_path = Path(file_path).resolve()
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

    def prepare(
        self,
        minimize: Optional[str] = None,
        remove_water: bool = True,
        add_hydrogens: bool = True,
        add_charges: bool = True,
        preserve_charge_types: Optional[list[str]] = None,
    ) -> None:
        """
        Prepare the ligand by converting to MOL2, optionally minimizing energy, 
        and generating a final PDBQT file using AutoDockTools (from MGLTools).

        Args:
            minimize (str, optional): Forcefield to use for energy minimization 
                ("mmff94", "mmff94s", "uff", or "gaff"). If None, no minimization 
                is performed. Default is None.
            remove_water (bool, optional): If True, remove water molecules during Open Babel preprocessing. Default is True.
            add_hydrogens (bool, optional): Add polar hydrogens during PDBQT preparation. Default is True.
            add_charges (bool, optional): Assign Gasteiger charges during PDBQT preparation. If False, all input charges are preserved. Default is True.
            preserve_charge_types (list[str], optional): Atom types (e.g., ["Zn", "Fe"]) whose
                charges should be preserved. Ignored if `add_charges=False`. Defaults to None.

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
        if remove_water:
            cmd += ["--delete", "HOH", "--delete", "[#8H2]"]

        if minimize:
            forcefield = minimize.lower()
            if forcefield not in self.SUPPORTED_FORCEFIELDS:
                raise ValueError(
                    f"❌ Unsupported forcefield '{forcefield}'. Supported: {self.SUPPORTED_FORCEFIELDS}")
            cmd += ["--minimize", "--ff", forcefield]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"❌ OpenBabel failed:\n{result.stderr}")

        # === Step 2: MGLTools to PDBQT ===
        pdbqt_filename = f"{self.mol2_path.stem}.pdbqt"
        mgl_cmd = [
            str(MGL_PYTHON_EXE), str(PREPARE_LIGAND_SCRIPT),
            "-l", self.mol2_path.name, "-o", pdbqt_filename,
            "-U", "nphs_lps"
        ]
        # ADT prepare_ligand4.py doesn't have -U waters flag, remove water is handled by obabel

        if add_hydrogens:
            mgl_cmd += ["-A", "hydrogens"]
        else:
            mgl_cmd += ["-A", "None"]

        # Charge options
        if not add_charges:
            mgl_cmd += ["-C"]  # preserve all charges
        elif preserve_charge_types:
            for atom_type in preserve_charge_types:
                mgl_cmd += ["-p", atom_type]

        result = subprocess.run(
            mgl_cmd,
            text=True,
            capture_output=True,
            cwd=self.temp_dir
        )

        if result.returncode != 0:
            raise RuntimeError(f"❌ MGLTools ligand preparation failed:\n{result.stderr}")

        self.pdbqt_path = self.temp_dir / pdbqt_filename



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
        """Save the prepared PDBQT file to the specified location.

        Args:
            save_to (str | Path, optional): Destination file or directory.
                - If directory: saves with the original filename.
                - If file path: saves with the given name.

        Returns:
            Path: Path to the saved PDBQT file.

        Raises:
            RuntimeError: If `prepare()` has not been run or the PDBQT file is missing.
        """
        if self.pdbqt_path is None or not self.pdbqt_path.exists():
            raise RuntimeError("❌ Ligand not prepared. Run prepare() first.")

        save_to = Path(save_to).expanduser().resolve()

        # treat as file only if it has a suffix (e.g., .pdbqt)
        if not save_to.suffix:
            save_to = save_to / self.pdbqt_path.name

        save_to.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(self.pdbqt_path, save_to)
        return save_to
