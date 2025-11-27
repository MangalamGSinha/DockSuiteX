"""Utility functions module initialization.

This module provides utility functions for DockSuiteX including molecular visualization,
structure fetching from online databases, output file parsing, and cleanup operations.

The utilities support the main docking workflow by providing tools for:
    - 3D molecular visualization in Jupyter notebooks (NGLView)
    - Fetching protein structures from RCSB PDB
    - Fetching ligand structures from PubChem
    - Parsing docking results into CSV format
    - Managing temporary files and binaries

Example:
    Using utility functions::

        from docksuitex.utils import (
            fetch_pdb,
            fetch_sdf,
            view_molecule,
            parse_vina_log_to_csv
        )

        # Fetch structures
        protein_path = fetch_pdb("1UBQ", save_to="structures")
        ligand_path = fetch_sdf("2244", save_to="structures")  # Aspirin

        # Visualize in Jupyter
        view_molecule(protein_path)

        # Parse docking results
        parse_vina_log_to_csv(
            results_dir="vina_results",
            output_csv="summary.csv"
        )

Modules:
    viewer: 3D molecular visualization with NGLView
    cleaner: Temporary file and binary cleanup utilities
    fetcher: Structure fetching from online databases (PDB, PubChem)
    parser: Docking output file parsing to CSV format
"""

from .viewer import view_molecule, view_results
from .fetcher import fetch_pdb, fetch_sdf
from .parser import parse_vina_log_to_csv, parse_ad4_dlg_to_csv
from .cleaner import delete_binaries

__all__ = [
    "view_molecule",
    "view_results",
    "fetch_pdb",
    "fetch_sdf",
    "parse_vina_log_to_csv",
    "parse_ad4_dlg_to_csv",
    "delete_binaries"
]
