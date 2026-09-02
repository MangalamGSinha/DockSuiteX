"""Batch protein preparation module for DockSuiteX."""

import os
from pathlib import Path
from typing import List, Union, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback

from ..protein import Protein



class BatchProtein:
    """Handles batch protein preparation with parallel processing."""

    def __init__(
        self,
        inputs: Union[List[Union[str, Path]], str, Path],
        fix_pdb: bool = True,
        remove_heterogens: bool = True,
        remove_water: bool = True,
        add_hydrogens: bool = True,
        ph: float = 7.4,
        add_charges: bool = True,
        preserve_charge_types: Optional[list[str]] = None,
    ):
        """Initialize the batch processor and preparation parameters.

        Args:
            inputs (Union[List[Union[str, Path]], str, Path]): List or directory of single-model PDB files.
            fix_pdb (bool, optional): Fix missing residues/atoms. Defaults to True.
            remove_heterogens (bool, optional): Remove ligands/heterogens. Defaults to True.
            remove_water (bool, optional): Remove water molecules. Defaults to True.
            add_hydrogens (bool, optional): Add hydrogens using PDBFixer. Defaults to True.
            ph (float, optional): Custom pH for adding missing hydrogens via PDBFixer. Defaults to 7.4.
            add_charges (bool, optional): Assign Gasteiger charges. Defaults to True.
            preserve_charge_types (Optional[list[str]], optional): Atom types to preserve charges for. Defaults to None.
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
        self.ph = ph
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
        ph: float,
        remove_water: bool,
        add_charges: bool,
        preserve_charge_types: Optional[list[str]],
        save_to: Path,
    ) -> Dict[str, Union[str, Path, bool]]:
        """Process one protein file.

        This is an internal worker used for parallel processing.
        """
        try:
            protein = Protein(
                input=file_path,
                fix_pdb=fix_pdb,
                remove_heterogens=remove_heterogens,
                add_hydrogens=add_hydrogens,
                ph=ph,
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
        cpu: int = (os.cpu_count() or 2) - 1,
    ) -> List[Dict[str, Union[str, Path, bool]]]:
        """Prepare all proteins in batch and save PDBQT files to the specified folder.

        Args:
            save_to (str | Path): Directory to save all prepared PDBQT files.
            cpu (int, optional): Total number of CPU cores to use.
                Defaults to ``os.cpu_count() - 1``. Each worker uses 1 CPU.

        Returns:
            list[dict]: Result dictionary for each protein file.
        """
        save_to = Path(save_to).resolve()
        save_to.mkdir(parents=True, exist_ok=True)

        # Simple strategy: Divide total CPUs among workers
        n_files = len(self.files)
        max_workers = min(cpu, n_files)

        print(f"Starting protein preparation for {n_files} files...")
        print(f"Using {max_workers} parallel workers, 1 CPUs per worker")
        print(f"Output directory: {save_to}")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    self._process_one,
                    file,
                    self.fix_pdb,
                    self.remove_heterogens,
                    self.add_hydrogens,
                    self.ph,
                    self.remove_water,
                    self.add_charges,
                    self.preserve_charge_types,
                    save_to,
                ): file for file in self.files
            }

            for future in as_completed(futures):
                result = future.result()
                self.results.append(result)

        print("✅ Batch processing completed!")
        return self.results
