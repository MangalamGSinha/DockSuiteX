"""Batch pocket finding with parallel processing.

This module provides high-throughput binding pocket prediction using P2Rank
with parallel execution. It automates pocket finding for multiple protein
structures simultaneously.

Example:
    Batch pocket finding::

        from docksuitex.batch_docking import BatchPocketFinder

        # Find pockets for all proteins
        batch = BatchPocketFinder(
            inputs="prepared_proteins",
            max_workers=4
        )

        # Run pocket finding in parallel
        pockets_dict = batch.run_all(save_to="pocket_results")

        # pockets_dict maps protein paths to pocket centers
        for protein_path, centers in pockets_dict.items():
            print(f"{protein_path}: {len(centers)} pockets found")
"""

import os
from pathlib import Path
from typing import List, Union, Dict, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback

from ..pocket_finder import PocketFinder
from ..protein import Protein



class BatchPocketFinder:

    def __init__(
        self,
        inputs: Union[List[Union[str, Path]], str, Path],
    ):
        """
        Initialize the batch pocket finder.

        Args:
            inputs (list[str | Path] | str | Path): 
                List of protein files (pdb/pdbqt) or directory containing them.
        """
        # Handle input paths
        # Handle input paths
        if isinstance(inputs, (str, Path)):
            path = Path(inputs).resolve()
            if path.is_dir():
                # We accept .pdb and .pdbqt for pocket finding
                self.files = [
                    f.resolve() for f in path.glob("*")
                    if f.suffix.lower() in [".pdb", ".pdbqt"]
                ]
            elif path.is_file():
                if path.suffix.lower() in [".pdb", ".pdbqt"]:
                    self.files = [path]
                else:
                    raise ValueError(f"❌ Invalid file type: {path.suffix}. Supported: .pdb, .pdbqt")
            else:
                 raise ValueError(f"❌ Input path does not exist: {inputs}")
        elif isinstance(inputs, list):
            self.files = [Path(f).resolve() for f in inputs]
        else:
            raise ValueError("❌ Invalid input. Provide a list of files, a directory path, or a single file path.")

        if not self.files:
            raise ValueError("❌ No valid protein files found.")
        
        # Sort files to ensure deterministic order (though we return a dict now)
        self.files.sort()

        self.results: Dict[str, List[Tuple[float, float, float]]] = {}

    @staticmethod
    def _process_one(
        file_path: Union[str, Path],
        save_to: Path,
        cpu_per_worker: int,
    ) -> Dict[str, Union[str, List, bool]]:
        """
        Internal worker to find pockets for one protein file.
        
        Args:
            file_path (Union[str, Path]): Path to protein file.
            save_to (Path): Directory to save results.
            cpu_per_worker (int): Number of CPU cores allocated to this worker.
            
        Returns:
            Dict: Result dictionary with status and pockets information.
        """
        try:
            # If it's a pdbqt, we can pass it directly. If it's pdb, PocketFinder handles it.
            
            # Save results to a specific folder for this protein
            protein_name = Path(file_path).stem
            protein_output_dir = save_to / f"{protein_name}_pockets"
            
            # Create PocketFinder with allocated CPU cores
            finder = PocketFinder(file_path, _cpu=cpu_per_worker)
            # Run and save directly to output_dir
            pockets = finder.run(save_to=protein_output_dir)

            return {
                "file": str(file_path),
                "status": "success",
                "pockets": pockets,
                "output_dir": str(protein_output_dir)
            }
        except Exception as e:
            return {
                "file": str(file_path),
                "status": "error",
                "error": str(e),
                "traceback": traceback.format_exc(),
            }

    def run_all(
        self,
        save_to: Union[str, Path],
        cpu: int = os.cpu_count() or 1,
    ) -> Dict[str, List[Tuple[float, float, float]]]:
        """
        Run pocket finding for all proteins in batch.

        Args:
            save_to (str | Path): Directory to save all results.
            cpu (int, optional): Total number of CPU cores to use.
                Defaults to all available cores. CPUs are divided among workers,
                with each worker receiving multiple CPUs for P2Rank.

        Returns:
            Dict[str, List[Tuple[float, float, float]]]: 
                Dictionary mapping absolute protein file path (str) to a list of center coordinates.
                Example: {"/path/to/prot1.pdbqt": [(1.0, 2.0, 3.0), ...]}
        """
        save_to = Path(save_to).resolve()
        save_to.mkdir(parents=True, exist_ok=True)

        # Simple strategy: Divide total CPUs among workers
        # Each worker gets at least 1 CPU
        n_files = len(self.files)
        
        # Calculate number of workers and CPUs per worker
        # Start with number of files, but limit by available CPUs
        max_workers = min(cpu, n_files)  # Can't have more workers than CPUs or files
        cpu_per_worker = max(1, cpu // max_workers)  # At least 1 CPU per worker
        
        print(f"Starting pocket finding for {n_files} proteins...")
        print(f"Using {max_workers} parallel workers, {cpu_per_worker} CPUs per worker")
        print(f"Total CPU allocation: {max_workers * cpu_per_worker}")
        print(f"Output directory: {save_to}")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    self._process_one,
                    file,
                    save_to,
                    cpu_per_worker,
                ): file for file in self.files
            }

            for future in as_completed(futures):
                result = future.result()
                file_path = result["file"]
                
                if result["status"] == "success":
                    pockets = result.get("pockets", [])
                    # Extract centers from pockets list (assuming pocket is a dict with "center" key)
                    # Based on usage in notebook: centers = [item["center"] for item in pockets]
                    centers = [p["center"] for p in pockets]
                    self.results[file_path] = centers
                    
                    count = len(centers)
                    print(f"✅ {Path(file_path).name} → found {count} pockets, saved to {Path(result['output_dir']).name}")
                else:
                    print(f"❌ {Path(file_path).name} failed. Error: {result.get('error', 'Unknown')}")
                    # We don't add failed files to results, or we could add empty list
                    # self.results[file_path] = [] 

        print("Batch processing completed!")
        return self.results
