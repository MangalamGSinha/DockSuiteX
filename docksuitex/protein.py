"""Protein structure preparation for molecular docking.

This module provides automated protein preparation using a combination of
PDBFixer, Open Babel, and AutoDockTools (MGLTools). It handles format conversion,
structure fixing, and PDBQT generation required for molecular docking.

The preparation workflow:
    1. Format conversion (if needed) using Open Babel
    2. Structure fixing with PDBFixer (missing residues/atoms, nonstandard residues)
    3. Heteroatom and water removal (optional)
    4. Hydrogen addition and charge assignment
    5. PDBQT format conversion using AutoDockTools

Example:
    Basic protein preparation::

        from docksuitex import Protein

        # Prepare protein from PDB file
        protein = Protein(
            input="protein.pdb",
            fix_pdb=True,
            remove_water=True,
            add_hydrogens=True
        )
        
        # Generate PDBQT file
        pdbqt_path = protein.prepare(save_to="prepared_protein.pdbqt")
        
        # Visualize in Jupyter
        protein.view_molecule()

    Advanced preparation with charge preservation::

        protein = Protein(
            input="protein_with_metal.pdb",
            fix_pdb=True,
            preserve_charge_types=["Zn", "Fe", "Mg"]
        )
        pdbqt_path = protein.prepare()

Attributes:
    MGLTOOLS_PATH (Path): Path to bundled MGLTools installation.
    MGL_PYTHON_EXE (Path): Path to MGLTools Python interpreter.
    PREPARE_RECEPTOR_SCRIPT (Path): Path to prepare_receptor4.py script.
    OBABEL_EXE (Path): Path to Open Babel executable.
    TEMP_DIR (Path): Temporary directory for intermediate files.
"""

import os
import subprocess
from pathlib import Path
from typing import Union, Optional, Literal
import shutil

from pdbfixer import PDBFixer
from openmm.app import PDBFile

import uuid

from .utils.viewer import view_molecule


# Constants - Paths to bundled executables
MGLTOOLS_PATH = (Path(__file__).parent / "bin" / "mgltools").resolve()
MGL_PYTHON_EXE = (MGLTOOLS_PATH / "python.exe").resolve()
PREPARE_RECEPTOR_SCRIPT = (
    MGLTOOLS_PATH / "Lib" / "site-packages" /
    "AutoDockTools" / "Utilities24" / "prepare_receptor4.py"
).resolve()

OBABEL_EXE = (Path(__file__).parent / "bin" /
              "obabel" / "obabel.exe").resolve()
TEMP_DIR = (Path.cwd() / "temp").resolve()
TEMP_DIR.mkdir(parents=True, exist_ok=True)



class Protein:
    """Automated protein structure preparation for molecular docking.

    This class provides a complete workflow for preparing protein structures
    for docking simulations. It integrates PDBFixer for structure repair,
    Open Babel for format conversion, and AutoDockTools for PDBQT generation.

    The preparation process handles:
        - Format conversion from various structure formats to PDB
        - Missing residue and atom reconstruction
        - Nonstandard residue replacement
        - Heteroatom and water molecule removal
        - Hydrogen addition (polar or all)
        - Gasteiger charge assignment
        - PDBQT format conversion for docking

    Example:
        Standard protein preparation::

            from docksuitex import Protein

            protein = Protein("1abc.pdb", fix_pdb=True)
            pdbqt_file = protein.prepare(save_to="receptor.pdbqt")

        Preparation with metal ion charge preservation::

            protein = Protein(
                "metalloprotein.pdb",
                preserve_charge_types=["Zn", "Fe", "Ca"]
            )
            pdbqt_file = protein.prepare()

    Note:
        Temporary files are created in a `temp/Proteins/` directory and are
        not automatically cleaned up. The final PDBQT file is saved to the
        specified location.
    """

    SUPPORTED_INPUTS = {".pdb", ".mol2", ".sdf",
                        ".pdbqt", ".cif", ".ent", ".xyz"}

    def __init__(
        self,
        input: Union[str, Path],
        fix_pdb: bool = True,
        remove_heterogens: bool = True,
        add_hydrogens: bool = True,
        remove_water: bool = True,
        add_charges: bool = True,
        preserve_charge_types: Optional[list[str]] = None,
    ):
        """Initialize a Protein object with a given file path and preparation parameters.

        Args:
            input (str | Path): Path to the protein input file.
            fix_pdb (bool, optional): Fix missing residues/atoms using PDBFixer. Defaults to True.
            remove_heterogens (bool, optional): Remove ligands/heterogens. Defaults to True.
            add_hydrogens (bool, optional): Add hydrogens. Defaults to True.
            remove_water (bool, optional): Remove water molecules. Defaults to True.
            add_charges (bool, optional): Assign Gasteiger charges. Defaults to True.
            preserve_charge_types (list[str], optional): Atom types to preserve charges for. Defaults to None.

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
        # Create a unique temp directory per object
        self.temp_dir = TEMP_DIR / "Proteins" / f"{self.file_path.stem}_{uuid.uuid4().hex[:8]}"
        self.temp_dir.mkdir(parents=True, exist_ok=True)

        # Convert to .pdb if needed
        if self.ext != ".pdb":
            self.pdb_path = self.temp_dir / f"{self.file_path.stem}.pdb"
            cmd = [str(OBABEL_EXE), str(self.file_path),"-O", str(self.pdb_path), "--gen3d"]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                raise PreparationError(f"Open Babel conversion failed:\n{result.stderr}")

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
        fixed_pdb_path = self.temp_dir / f"{self.file_path.stem}_fixed.pdb"
        with open(fixed_pdb_path, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        self.fixed_pdb_path = fixed_pdb_path



        # Convert to PDBQT using AutoDockTools
        save_to = Path(save_to).expanduser().resolve()
        
        # treat as file only if it has a suffix (e.g., .pdbqt)
        if not save_to.suffix:
            save_to = save_to / f"{self.file_path.stem}.pdbqt"
        
        save_to.parent.mkdir(parents=True, exist_ok=True)

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

        # # Save logic (merged from save_pdbqt)
        # save_to = Path(save_to).expanduser().resolve()
        
        # # treat as file only if it has a suffix (e.g., .pdbqt)
        # if not save_to.suffix:
        #     save_to = save_to / self.pdbqt_path.name

        # save_to.parent.mkdir(parents=True, exist_ok=True)
        # shutil.copy2(self.pdbqt_path, save_to)
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

        Example:
            ::

                protein = Protein("protein.pdb")
                protein.prepare()
                protein.view_molecule()  # Opens 3D viewer
        """
        path = Path(self.pdbqt_path if self.pdbqt_path else self.file_path).resolve()
        return view_molecule(file_path=path)


