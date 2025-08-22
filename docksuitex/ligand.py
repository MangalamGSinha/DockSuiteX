import shutil
import subprocess
from pathlib import Path
from typing import Optional, Union
import uuid

# === CONFIGURATION ===
MGLTOOLS_PATH = (Path(__file__).parent / "bin" / "mgltools").resolve()
MGL_PYTHON_EXE = (MGLTOOLS_PATH / "python.exe").resolve()
PREPARE_LIGAND_SCRIPT = (MGLTOOLS_PATH / "Lib" / "site-packages" /
                         "AutoDockTools" / "Utilities24" / "prepare_ligand4.py").resolve()
OBABEL_EXE = (Path(__file__).parent / "bin" /
              "OpenBabel-3.1.1" / "obabel.exe").resolve()

TEMP_DIR = (Path(__file__).parent / "temp").resolve()
TEMP_DIR.mkdir(exist_ok=True)


class Ligand:
    """
    Ligand preparation pipeline using Open Babel and MGLTools.

    This class supports automatic conversion, 3D generation, optional energy minimization, 
    and conversion to the PDBQT format required by AutoDock Vina.

    Supported input formats: mol2, sdf, pdb, mol, smi  
    Supported forcefields: mmff94, mmff94s, uff, gaff
    """

    SUPPORTED_INPUTS = {"mol2", "sdf", "pdb", "mol", "smi"}
    SUPPORTED_FORCEFIELDS = {"mmff94", "mmff94s", "uff", "gaff"}

    def __init__(self, file_path: Union[str, Path]):
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
            remove_water (bool): If True, remove water molecules during Open Babel 
                preprocessing (HOH residues and [#8H2] pattern). Default is True.
            add_hydrogens (bool): Whether to explicitly add hydrogens in ADT. 
                Default is True.
            add_charges (bool): Whether to assign Gasteiger charges in ADT. 
                If False, all input charges are preserved. Default is True.
            preserve_charge_types (list[str], optional): Atom types (e.g., ["Zn", "Fe"]) 
                whose charges should be preserved. Other atoms get Gasteiger charges. 
                Default is None.

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

    def save_pdbqt(self, save_path: Union[str, Path] = ".") -> Path:
        """
        Save the prepared PDBQT file to the specified location.

        Args:
            save_path (str | Path): Destination file or directory where the PDBQT file should be saved. 
                If a directory is given, the original filename is preserved.

        Raises:
            RuntimeError: If prepare() has not been called or the PDBQT file does not exist.
        """
        if self.pdbqt_path is None or not self.pdbqt_path.exists():
            raise RuntimeError("❌ Ligand not prepared. Run prepare() first.")
        
        save_path = Path(save_path).resolve()
        if save_path.is_dir():
            save_path = save_path / self.pdbqt_path.name

        save_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(self.pdbqt_path, save_path)
        return save_path
