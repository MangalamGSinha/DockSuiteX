"""Ligand structure preparation module for DockSuiteX."""

import subprocess
import shutil
from pathlib import Path
from typing import Optional, Union
import os

from docksuitex.utils.viewer import view_molecule
from .platform_config import (
    MGLTOOLS_PATH, MGL_PYTHON_EXE, PREPARE_LIGAND_SCRIPT,
    get_mgltools_env,
)

# Locate obabel from PATH (installed via: pip install openbabel-wheel)
OBABEL_EXE = shutil.which("obabel")



class Ligand:
    """Ligand structure preparation for molecular docking.

    This module provides automated ligand preparation using Open Babel for
    optional energy minimization and water removal, combined with
    AutoDockTools (MGLTools) for PDBQT generation.

    The preparation workflow:
        1. Optional energy minimization with forcefields: MMFF94, MMFF94S, UFF, GAFF (using Open Babel).
        2. Optional water molecule removal (using Open Babel).
        3. Optional pH-aware hydrogen addition (using Open Babel) and Gasteiger charge assignment (using AutoDockTools).
        4. Removing non-polar hydrogens and converting to PDBQT format (using AutoDockTools).

    Supported Input Formats:
        SDF, MOL2 (Single Molecule only).

    Note:
        If you have files in other formats (PDB, SMILES, etc.) or
        multi-molecule libraries (multi-SDF), please use the
        **Format Converter** utility to prepare a single-molecule
        SDF or MOL2 file before using this module.

        Intermediate processed files (prepared MOL2) are
        saved in an ``intermediate_ligands/`` subfolder within the
        output directory.
    """

    SUPPORTED_INPUTS = {"sdf", "mol2"}
    SUPPORTED_FORCEFIELDS = {"mmff94", "mmff94s", "uff", "gaff"}

    def __init__(
        self,
        input: Union[str, Path],
        minimize: Optional[str] = None,
        remove_water: bool = True,
        add_hydrogens: bool = True,
        ph: float = 7.4,
        add_charges: bool = True,
        preserve_charge_types: Optional[list[str]] = None,
    ):
        """Initialize a Ligand object with a given input file and preparation parameters.

        Args:
            input (str | Path): Path to the single-molecule ligand file (SDF or MOL2).
            minimize (str, optional): Forcefield for energy minimization ("mmff94", "mmff94s", "uff", "gaff").
            remove_water (bool, optional): Remove water molecules. Defaults to True.
            add_hydrogens (bool, optional): Add polar hydrogens via Open Babel at specified pH. Defaults to True.
            ph (float, optional): Custom pH for Open Babel hydrogen addition. Defaults to 7.4.
            add_charges (bool, optional): Assign Gasteiger charges. Defaults to True.
            preserve_charge_types (list[str], optional): Atom types (e.g.,["Zn", "Fe"]) whose charges are preserved; 
                others get Gasteiger charges; ignored if add_charges=False. Defaults to None.

        Raises:
            FileNotFoundError: If the input file does not exist.
            ValueError: If the file extension is not a supported format.
        """
        self.file_path = Path(input).resolve()
        self.prepared_mol2_path: Optional[Path] = None
        self.pdbqt_path: Optional[Path] = None

        if not self.file_path.is_file():
            raise FileNotFoundError(
                f"❌ Ligand file not found: {self.file_path}")

        ext = self.file_path.suffix.lower().lstrip(".")
        if ext not in self.SUPPORTED_INPUTS:
            raise ValueError(
                f"❌ Unsupported file format '.{ext}'. Ligand preparation only accepts .sdf and .mol2 files. "
                f"Please use the Format Converter utility to convert other formats.")
        self.input_format = ext

        # Single molecule check
        with open(self.file_path, "r", errors="ignore") as f:
            content = f.read()
            if self.input_format == "sdf":
                mol_count = content.count("$$$$")
                if mol_count > 1:
                    raise ValueError(
                        f"❌ Ligand file contains {mol_count} molecules. Multi-molecule SDF libraries are not supported. "
                        f"Please use the Format Converter utility to split the file before preparation.")
            elif self.input_format == "mol2":
                mol_count = content.count("@<TRIPOS>MOLECULE")
                if mol_count > 1:
                    raise ValueError(
                        f"❌ Ligand file contains {mol_count} molecules. Multi-molecule MOL2 files are not supported. "
                        f"Please use the Format Converter utility to split the file before preparation.")

        # Store preparation parameters
        self.minimize = minimize
        self.remove_water = remove_water
        self.add_hydrogens = add_hydrogens
        self.ph = ph
        self.add_charges = add_charges
        self.preserve_charge_types = preserve_charge_types

    def _is_input_3d(self) -> bool:
        """Check if the input molecule has 3D coordinates.

        Optimized check:
        1. SMILES/CDX are always 2D.
        2. SDF/MOL are parsed to check Z-coordinates.
        3. PDB/MOL2/PDBQT/XYZ are assumed 3D.
        """
        if self.input_format in ["smi", "smiles", "cdx"]:
            return False

        if self.input_format in ["sdf", "mol"]:
            try:
                # Simple parser for SDF/MOL atom block Z-coordinates
                with open(self.file_path, "r", encoding="utf-8", errors="ignore") as f:
                    # Move to counts line (line 4)
                    for _ in range(3):
                        f.readline()
                    counts_line = f.readline()
                    if not counts_line:
                        return False
                    
                    try:
                        num_atoms = int(counts_line[:3].strip())
                    except (ValueError, IndexError):
                        return False
                    
                    # Check first atoms in the block
                    for _ in range(num_atoms):
                        line = f.readline()
                        if not line or len(line) < 30:
                            break
                        try:
                            z_val = float(line[20:30].strip())
                            if abs(z_val) > 1e-5:
                                return True
                        except (ValueError, IndexError):
                            continue
                return False
            except Exception:
                return False

        # Structural formats are assumed 3D
        return True

    def prepare(self, save_to: Union[str, Path] = ".") -> Path:
        """Handle ligand preparation for docking using Open Babel and AutoDockTools (ADT).

        Saves the prepared PDBQT file to the specified location.

        Args:
            save_to (str | Path, optional): Destination file or directory.
                - If directory: saves with the original filename.
                - If file path: saves with the given name.
                Defaults to current directory.

        Returns:
            Path: Path to the saved PDBQT file.

        Raises:
            ValueError: If an unsupported forcefield is provided.
            RuntimeError: If Open Babel or AutoDockTools fails.
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

        self.prepared_mol2_path = intermediates_dir / f"{self.file_path.stem}_prepared.mol2"

        # Convert input → MOL2 via Open Babel (with optional water removal & minimization)
        # MGLTools' prepare_ligand4.py requires MOL2 (or PDB) format as input.
        # We always run Open Babel to ensure the `_prepared.mol2` file is created.
        cmd = [
            OBABEL_EXE,
            "-i", self.input_format, str(self.file_path),
            "-o", "mol2", "-O", str(self.prepared_mol2_path),
        ]

        if self.add_hydrogens:
            cmd += ["-p", str(self.ph)]

        if self.remove_water:
            cmd += ["--delete", "HOH", "--delete", "[#8H2]"]

        # Conditional 3D generation
        if not self._is_input_3d():
            cmd.append("--gen3d")

        if self.minimize:
            forcefield = self.minimize.lower()
            if forcefield not in self.SUPPORTED_FORCEFIELDS:
                raise ValueError(
                    f"❌ Unsupported forcefield '{forcefield}'. Supported: {self.SUPPORTED_FORCEFIELDS}")
                
            cmd += ["--minimize", "--ff", forcefield, "--steps", "2500"]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0 or not self.prepared_mol2_path.exists():
            raise RuntimeError(
                f"❌ OpenBabel ligand processing failed:\n"
                f"--- STDERR ---\n{result.stderr}\n"
                f"--- STDOUT ---\n{result.stdout}"
            )

        # Validate that Open Babel actually wrote a non-empty molecule file
        if not self.prepared_mol2_path.exists() or self.prepared_mol2_path.stat().st_size == 0:
            raise RuntimeError(
                "❌ OpenBabel produced an empty MOL2 file. "
                "The input molecule may be invalid or incompatible with "
                f"the '{self.minimize}' forcefield."
            )

        # === MGLTools: Convert MOL2 to PDBQT ===
        mgl_cmd = [
            str(MGL_PYTHON_EXE), str(PREPARE_LIGAND_SCRIPT),
            "-l", str(self.prepared_mol2_path), "-o", str(save_to),
            "-U", "nphs_lps"
        ]

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
            cwd=str(self.prepared_mol2_path.parent),
            env=get_mgltools_env(),
        )

        # MGLTools may exit 0 even on failure; check stderr and output file
        if result.returncode != 0 or not save_to.exists() or save_to.stat().st_size == 0 or "Traceback" in result.stderr:
            raise RuntimeError(
                f"❌ MGLTools ligand preparation failed:\n"
                f"--- STDERR ---\n{result.stderr}\n"
                f"--- STDOUT ---\n{result.stdout}"
            )

        self.pdbqt_path = save_to

        print(f"✅ Ligand prepared successfully: {self.pdbqt_path}")
        return self.pdbqt_path

    def view_molecule(self) -> "nv.NGLWidget":
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
