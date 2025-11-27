"""Batch protein preparation with parallel processing.

This module provides high-throughput protein preparation capabilities using
parallel execution. It automates the preparation of multiple protein structures
simultaneously for docking workflows.

Example:
    Batch protein preparation::

        from docksuitex.batch_docking import BatchProtein

        # Prepare all proteins in a folder
        batch = BatchProtein(
            inputs="proteins_folder",
            fix_pdb=True,
            remove_water=True
        )

        # Process all proteins in parallel
        results = batch.prepare_all(save_to="prepared_proteins")

        # Check results
        for result in results:
            if result["status"] == "success":
                print(f"Prepared: {result['pdbqt_path']}")
"""

import os
from pathlib import Path
from typing import List, Union, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback

from ..protein import Protein



class BatchProtein:

    def __init__(
        self,
        inputs: Union[List[Union[str, Path]], str, Path],
        fix_pdb: bool = True,
        remove_heterogens: bool = True,
        add_hydrogens: bool = True,
        remove_water: bool = True,
        add_charges: bool = True,
        preserve_charge_types: Optional[list[str]] = None,
    ):
        """
        Initialize the batch processor and preparation parameters.

        Args:
            inputs (list[str | Path] | str | Path): 
                List of protein files or directory containing protein files.
            fix_pdb, remove_heterogens, add_hydrogens, remove_water, add_charges: 
                Preparation parameters (same as Protein.prepare()).
            preserve_charge_types (list[str], optional): Atom types to preserve charges for.
        """
        # Handle input paths
        if isinstance(inputs, (str, Path)):
            path = Path(inputs).resolve()
            if path.is_dir():
                self.files = [
                    f.resolve() for f in path.glob("*")
                    if f.suffix.lower() in Protein.SUPPORTED_INPUTS
                ]
            elif path.is_file():
                if path.suffix.lower() in Protein.SUPPORTED_INPUTS:
                    self.files = [path]
                else:
                    raise ValueError(f"❌ Invalid file type: {path.suffix}. Supported: {Protein.SUPPORTED_INPUTS}")
            else:
                 raise ValueError(f"❌ Input path does not exist: {inputs}")
        elif isinstance(inputs, list):
            self.files = [Path(f).resolve() for f in inputs]
        else:
            raise ValueError("❌ Invalid input. Provide a list of files, a directory path, or a single file path.")

        if not self.files:
            raise ValueError("❌ No valid protein files found.")

        # Store preparation parameters
        self.fix_pdb = fix_pdb
        self.remove_heterogens = remove_heterogens
        self.add_hydrogens = add_hydrogens
        self.remove_water = remove_water
        self.add_charges = add_charges
        self.preserve_charge_types = preserve_charge_types

        self.results: List[Dict[str, Union[str, Path, bool]]] = []

    @staticmethod
    def _process_one(
        file_path: Union[str, Path],
        fix_pdb: bool,
        remove_heterogens: bool,
        add_hydrogens: bool,
        remove_water: bool,
        add_charges: bool,
        preserve_charge_types: Optional[list[str]],
        save_to: Path,
    ) -> Dict[str, Union[str, Path, bool]]:
        """
        Internal worker to process one protein file.
        Returns a result dictionary.
        """
        try:
            protein = Protein(
                input=file_path,
                fix_pdb=fix_pdb,
                remove_heterogens=remove_heterogens,
                add_hydrogens=add_hydrogens,
                remove_water=remove_water,
                add_charges=add_charges,
                preserve_charge_types=preserve_charge_types,
            )
            
            # Prepare and save in one step
            save_path = protein.prepare(save_to=save_to)

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
        Prepare all proteins in batch and save PDBQT files to the specified folder.

        Args:
            save_to (str | Path): Directory to save all prepared PDBQT files.
            cpu (int, optional): Total number of CPU cores to use.
                Defaults to all available cores. Each worker uses 1 CPU.

        Returns:
            list[dict]: Result dictionary for each protein file.
        """
        save_to = Path(save_to).resolve()
        save_to.mkdir(parents=True, exist_ok=True)

        # Simple strategy: Divide total CPUs among workers
        n_files = len(self.files)
        max_workers = min(cpu, n_files)

        print(f"Starting protein preparation for {n_files} files...")
        print(f"Using {max_workers} parallel workers")
        print(f"Total CPU allocation: {max_workers}")
        print(f"Output directory: {save_to}")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    self._process_one,
                    file,
                    self.fix_pdb,
                    self.remove_heterogens,
                    self.add_hydrogens,
                    self.remove_water,
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
