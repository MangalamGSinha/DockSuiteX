"""Protein structure preparation module for DockSuiteX."""

import subprocess
import shutil
import os
from pathlib import Path
from typing import Union, Optional

from pdbfixer import PDBFixer
from openmm.app import PDBFile

from .utils.viewer import view_molecule
from .platform_config import (
    MGLTOOLS_PATH, MGL_PYTHON_EXE, PREPARE_RECEPTOR_SCRIPT,
    get_mgltools_env,
)

# Locate obabel from PATH (installed via: pip install openbabel-wheel)
OBABEL_EXE = shutil.which("obabel")


class Protein:
    """Protein structure preparation for molecular docking.

    This module provides automated protein preparation using PDBFixer for
    fixing structures and AutoDockTools (MGLTools) for PDBQT generation.

    The preparation workflow:
        1. Optionally fixing missing residues/atoms and nonstandard residues (using PDBFixer).
        2. Optionally removing heterogens (using PDBFixer) and removing water molecules (using AutoDockTools).
        3. Optionally adding pH-aware hydrogens (using PDBFixer) and Gasteiger charges (using AutoDockTools).
        4. Removing non-polar hydrogens and converting to PDBQT format (using AutoDockTools).

    Supported Input Formats:
        PDB (Single Model only).

    Note:
        If you have files in other formats (CIF, MOL2, etc.) or multi-model
        files (NMR ensembles), use the **Format Converter** utility
        to prepare a single-model PDB file before using this module.

        Intermediate files (fixed PDB) are saved in an
        ``intermediate_proteins/`` subfolder within the output directory.
    """

    SUPPORTED_INPUTS = {".pdb"}

    def __init__(
        self,
        input: Union[str, Path],
        fix_pdb: bool = True,
        remove_heterogens: bool = True,
        remove_water: bool = True,
        add_hydrogens: bool = True,
        ph: float = 7.4,
        add_charges: bool = True,
        preserve_charge_types: Optional[list[str]] = None,
    ):
        """Initialize a Protein object with a given file path and preparation parameters.

        Args:
            input (str | Path): Path to the single-model PDB file.
            fix_pdb (bool, optional): Fix missing residues/atoms. Defaults to True.
            remove_heterogens (bool, optional): Remove ligands/heterogens. Defaults to True.
            remove_water (bool, optional): Remove water molecules. Defaults to True.
            add_hydrogens (bool, optional): Add hydrogens using PDBFixer at the specified pH. Defaults to True.
            ph (float, optional): Custom pH for adding missing hydrogens. Defaults to 7.4.
            add_charges (bool, optional): Assign Gasteiger charges. Defaults to True.
            preserve_charge_types (list[str], optional): Atom types (e.g.,["Zn", "Fe"]) whose charges are preserved; 
                others get Gasteiger charges; ignored if add_charges=False. Defaults to None.

        Raises:
            FileNotFoundError: If the provided file does not exist.
            ValueError: If the file extension is not a supported format.
        """
        self.file_path = Path(input).resolve()
        self.pdb_path: Optional[Path] = None
        self.prepared_pdb_path: Optional[Path] = None
        self.pdbqt_path: Optional[Path] = None

        if not self.file_path.is_file():
            raise FileNotFoundError(
                f"❌ Protein file not found: {self.file_path}")

        self.ext = self.file_path.suffix.lower()
        if self.ext not in self.SUPPORTED_INPUTS:
            raise ValueError(
                f"❌ Unsupported file format '{self.ext}'. Protein preparation only accepts .pdb files. "
                f"Please use the Format Converter utility to convert other formats.")

        # Single model check
        with open(self.file_path, "r", errors="ignore") as f:
            model_count = sum(1 for line in f if line.startswith("MODEL"))
            if model_count > 1:
                raise ValueError(
                    f"❌ Protein file contains {model_count} models. Multi-model files are not supported. "
                    f"Please use the Format Converter utility to split the file before preparation.")

        # Store preparation parameters
        self.fix_pdb = fix_pdb
        self.remove_heterogens = remove_heterogens
        self.add_hydrogens = add_hydrogens
        self.ph = ph
        self.remove_water = remove_water
        self.add_charges = add_charges
        self.preserve_charge_types = preserve_charge_types

    def prepare(self, save_to: Union[str, Path] = ".") -> Path:
        """Handle protein preparation for docking using PDBFixer, Open Babel, and AutoDockTools (ADT).

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

        # Convert to PDB via Open Babel if input is not already PDB
        if self.ext != ".pdb":
            input_format = self.ext.lstrip(".")
            converted_pdb_path = intermediates_dir / f"{self.file_path.stem}.pdb"

            cmd = [
                OBABEL_EXE,
                "-i", input_format, str(self.file_path),
                "-o", "pdb", "-O", str(converted_pdb_path),
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0 or not converted_pdb_path.exists():
                raise RuntimeError(
                    f"❌ OpenBabel conversion to PDB failed:\n"
                    f"--- STDERR ---\n{result.stderr}\n"
                    f"--- STDOUT ---\n{result.stdout}"
                )

            self.pdb_path = converted_pdb_path
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
        if self.add_hydrogens:
            fixer.addMissingHydrogens(self.ph)

        # Save prepared PDB
        self.prepared_pdb_path = intermediates_dir / f"{self.file_path.stem}_prepared.pdb"
        with open(self.prepared_pdb_path, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)


        # Convert to PDBQT using AutoDockTools

        U_flag = "nphs_lps_waters" if self.remove_water else "nphs_lps"
        cmd = [
            str(MGL_PYTHON_EXE),
            str(PREPARE_RECEPTOR_SCRIPT),
            "-r", str(self.prepared_pdb_path),
            "-o", str(save_to),
            "-U", U_flag
        ]

        # Control charges
        if not self.add_charges:
            cmd += ["-C"]  # disable Gasteiger charges
        elif self.preserve_charge_types:
            for atom in self.preserve_charge_types:
                cmd += ["-p", atom]

        result = subprocess.run(cmd, capture_output=True, text=True, env=get_mgltools_env())
        if result.returncode != 0 or not save_to.exists():
            raise RuntimeError(
                f"❌ Error preparing PDBQT:\n"
                f"--- STDERR ---\n{result.stderr}\n"
                f"--- STDOUT ---\n{result.stdout}"
            )

        self.pdbqt_path = save_to

        print(f"✅ Protein prepared successfully: {self.pdbqt_path}")
        return self.pdbqt_path

    def view_molecule(self) -> "nv.NGLWidget":
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


