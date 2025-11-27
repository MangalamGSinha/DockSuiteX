"""Batch docking module initialization.

This module provides batch processing capabilities for high-throughput molecular
docking workflows. It includes classes for parallel preparation of proteins and
ligands, batch pocket finding, and batch docking with AutoDock Vina and AutoDock4.

The batch processing classes use Python's ProcessPoolExecutor for parallel execution,
significantly speeding up workflows when processing multiple structures.

Example:
    Batch protein preparation and docking workflow::

        from docksuitex.batch_docking import (
            BatchProtein,
            BatchLigand,
            BatchPocketFinder,
            BatchVinaDocking
        )

        # Prepare multiple proteins
        batch_proteins = BatchProtein(
            inputs="proteins_folder",
            fix_pdb=True
        )
        batch_proteins.prepare_all(save_to="prepared_proteins")

        # Prepare multiple ligands
        batch_ligands = BatchLigand(
            inputs="ligands_folder",
            minimize="mmff94"
        )
        batch_ligands.prepare_all(save_to="prepared_ligands")

        # Find pockets for all proteins
        batch_pockets = BatchPocketFinder(inputs="prepared_proteins")
        pockets_dict = batch_pockets.run_all(save_to="pocket_results")

        # Run batch docking
        batch_docking = BatchVinaDocking(
            receptors_with_centers=pockets_dict,
            ligands="prepared_ligands"
        )
        results = batch_docking.run_all(save_to="docking_results")

Modules:
    batch_autodock4: Batch AutoDock4 docking
    batch_vina: Batch AutoDock Vina docking
    batch_protein: Batch protein preparation
    batch_ligand: Batch ligand preparation
    batch_pocket_finder: Batch pocket prediction
"""

from .batch_autodock4 import BatchAD4Docking
from .batch_vina import BatchVinaDocking
from .batch_protein import BatchProtein
from .batch_ligand import BatchLigand
from .batch_pocket_finder import BatchPocketFinder

__all__ = [
    "BatchAD4Docking",
    "BatchVinaDocking",
    "BatchProtein",
    "BatchLigand",
    "BatchPocketFinder" 
]
