"""Batch ligand preparation with parallel processing.

This module provides high-throughput ligand preparation capabilities using
parallel execution. It automates the preparation of multiple ligand structures
simultaneously for docking workflows, including energy minimization.

Example:
    Batch ligand preparation::

        from docksuitex.batch_docking import BatchLigand

        # Prepare all ligands with energy minimization
        batch = BatchLigand(
            inputs="ligands_folder",
            minimize="mmff94",
            remove_water=True
        )

        # Process all ligands in parallel
        results = batch.prepare_all(save_to="prepared_ligands")

        # Check results
        successful = [r for r in results if r["status"] == "success"]
        print(f"Prepared {len(successful)} ligands")
"""

import os
from pathlib import Path
from typing import List, Union, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback

from ..ligand import Ligand



class BatchLigand:
    """
    Handles batch ligand preparation in parallel.

    Example:
        batch = BatchLigand(
            inputs="ligands_folder",
            minimize="mmff94",
            remove_water=True,
            add_hydrogens=True,
            add_charges=True,
        )
        batch.prepare_all(save_to="prepared_ligands")
    """

    def __init__(
        self,
        inputs: Union[List[Union[str, Path]], str, Path],
        minimize: Optional[str] = None,
        remove_water: bool = True,
        add_hydrogens: bool = True,
        add_charges: bool = True,
        preserve_charge_types: Optional[list[str]] = None,
    ):
        """
        Initialize the batch processor and preparation parameters.

        Args:
            inputs (list[str | Path] | str | Path): 
                List of ligand files or directory containing ligand files.
            minimize, remove_water, add_hydrogens, add_charges: 
                Preparation parameters (same as Ligand.prepare()).
            preserve_charge_types (list[str], optional): Atom types to preserve charges for.
        """
        # Handle input paths
        if isinstance(inputs, (str, Path)):
            path = Path(inputs).resolve()
            if path.is_dir():
                self.files = [
                    f.resolve() for f in path.glob("*")
                    if f.suffix.lower().lstrip(".") in Ligand.SUPPORTED_INPUTS
                ]
            elif path.is_file():
                if path.suffix.lower().lstrip(".") in Ligand.SUPPORTED_INPUTS:
                    self.files = [path]
                else:
                    raise ValueError(f"❌ Invalid file type: {path.suffix}. Supported: {Ligand.SUPPORTED_INPUTS}")
            else:
                 raise ValueError(f"❌ Input path does not exist: {inputs}")
        elif isinstance(inputs, list):
            self.files = [Path(f).resolve() for f in inputs]
        else:
            raise ValueError("❌ Invalid input. Provide a list of files, a directory path, or a single file path.")

        if not self.files:
            raise ValueError("❌ No valid ligand files found.")

        # Store preparation parameters
        self.minimize = minimize
        self.remove_water = remove_water
        self.add_hydrogens = add_hydrogens
        self.add_charges = add_charges
        self.preserve_charge_types = preserve_charge_types

        self.results: List[Dict[str, Union[str, Path, bool]]] = []

    @staticmethod
    def _process_one(
        file_path: Union[str, Path],
        minimize: Optional[str],
        remove_water: bool,
        add_hydrogens: bool,
        add_charges: bool,
        preserve_charge_types: Optional[list[str]],
        save_to: Path,
    ) -> Dict[str, Union[str, Path, bool]]:
        """
        Internal worker to process one ligand file.
        Returns a result dictionary.
        """
        try:
            ligand = Ligand(
                input=file_path,
                minimize=minimize,
                remove_water=remove_water,
                add_hydrogens=add_hydrogens,
                add_charges=add_charges,
                preserve_charge_types=preserve_charge_types,
            )
            
            # Prepare and save in one step
            save_path = ligand.prepare(save_to=save_to)

            return {
                "file": str(file_path),
                "status": "success",
                "pdbqt_path": str(save_path),
            }
        except Exception as e:
            return {
                "file": str(file_path),
                "status": "error",
                "error": str(e),
                "traceback": traceback.format_exc(),
            }

    def prepare_all(
        self,
        save_to: Union[str, Path],
        cpu: int = os.cpu_count() or 1,
    ) -> List[Dict[str, Union[str, Path, bool]]]:
        """
        Prepare all ligands in batch and save PDBQT files to the specified folder.

        Args:
            save_to (str | Path): Directory to save all prepared PDBQT files.
            cpu (int, optional): Total number of CPU cores to use.
                Defaults to all available cores. Each worker uses 1 CPU.

        Returns:
            list[dict]: Result dictionary for each ligand file.
        """
        save_to = Path(save_to).resolve()
        save_to.mkdir(parents=True, exist_ok=True)

        # Simple strategy: Divide total CPUs among workers
        n_files = len(self.files)
        max_workers = min(cpu, n_files)

        print(f"Starting ligand preparation for {n_files} files...")
        print(f"Using {max_workers} parallel workers")
        print(f"Total CPU allocation: {max_workers}")
        print(f"Output directory: {save_to}")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    self._process_one,
                    file,
                    self.minimize,
                    self.remove_water,
                    self.add_hydrogens,
                    self.add_charges,
                    self.preserve_charge_types,
                    save_to,
                ): file for file in self.files
            }

            for future in as_completed(futures):
                result = future.result()
                self.results.append(result)
                if result["status"] == "success":
                    print(f"✅ {Path(result['file']).name} → saved to {Path(result['pdbqt_path']).name}")
                else:
                    print(f"❌ {Path(result['file']).name} failed. Error: {result.get('error', 'Unknown')}")

        print("Batch processing completed!")
        return self.results
