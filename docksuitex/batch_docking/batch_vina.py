"""Batch AutoDock Vina docking module for DockSuiteX."""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Union, Sequence, Dict, List
import os
from ..vina import VinaDocking



class BatchVinaDocking:
    """Batch docking manager for AutoDock Vina.

    Runs AutoDock Vina docking for multiple ligands against multiple proteins,
    each with its own set of binding pocket centers, in parallel using a process pool.
    """

    def __init__(
        self,
        receptors_with_pockets: Dict[Union[str, Path], List[Dict]],
        ligands: Union[Sequence[Union[str, Path]], str, Path],
        exhaustiveness: int = 8,
        num_modes: int = 9,
        verbosity: int = 1,
        seed: int | None = 42,
    ):
        """Initialize a batch Vina docking job.

        Args:
            receptors_with_pockets (Dict[str | Path, List[Dict]]):
                Dictionary mapping receptor PDBQT files to their list of binding
                pockets (each pocket MUST be a dict with 'center' and 'grid_size').
            ligands (Sequence[str | Path]): List of ligand PDBQT files.
            exhaustiveness (int, optional):
                Sampling exhaustiveness. Higher values increase accuracy but
                also computation time. Defaults to 8.
            num_modes (int, optional):
                Maximum number of binding modes. Defaults to 9.
            verbosity (int, optional):
                Verbosity level (0 = quiet, 1 = normal, 2 = verbose).
                Defaults to 1.
            seed (int, optional):
                Random seed. If None, Vina selects automatically. Defaults to 42.
        """
        self.receptors = receptors_with_pockets
        
        # Handle ligands input
        if isinstance(ligands, (str, Path)):
            path = Path(ligands).resolve()
            if path.is_dir():
                self.ligands = list(path.glob("*.pdbqt"))
            elif path.is_file():
                if path.suffix.lower() == ".pdbqt":
                    self.ligands = [path]
                else:
                    raise ValueError(f"❌ Invalid ligand file type: {path.suffix}. Must be .pdbqt")
            else:
                 raise ValueError(f"❌ Ligand input path does not exist: {ligands}")
        elif isinstance(ligands, list):
            self.ligands = [Path(l).expanduser().resolve() for l in ligands]
        else:
            raise ValueError("❌ Invalid ligands input. Provide a list of files, a directory path, or a single file path.")

        if not self.ligands:
            raise ValueError("❌ No valid ligand PDBQT files found.")

        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.seed = seed
        self.verbosity = verbosity
        self.results: dict[tuple[str, str, tuple[float, float, float]], Union[Path, str]] = {}
        self.output_dir: Path | None = None

    def _dock_one(
        self,
        save_to: Union[str, Path],
        receptor: Path,
        ligand: Path,
        center: tuple[float, float, float],
        grid_size: tuple[float, float, float],
        vina_cpu: int,
    ) -> tuple[str, str, tuple[float, float, float], Path]:
        """Dock a single ligand at a single pocket center of a receptor.

        Args:
            save_to (str | Path): Directory to save results.
            receptor (Path): Receptor PDBQT file.
            ligand (Path): Ligand PDBQT file.
            center (tuple[float, float, float]): Grid center coordinates.
            vina_cpu (int): Number of CPUs assigned to this docking job.

        Returns:
            tuple: (receptor_name, ligand_name, center, result_path)
        """
        vina = VinaDocking(
            receptor=receptor,
            ligand=ligand,
            grid_center=center,
            grid_size=grid_size,
            exhaustiveness=self.exhaustiveness,
            num_modes=self.num_modes,
            verbosity=self.verbosity,
            seed=self.seed,
            _cpu=vina_cpu,  # CPU cores per worker
        )

        center_str = "_".join(f"{c:.2f}" for c in center)
        vina.run(save_to=Path(save_to) / f"{receptor.stem}_{ligand.stem}_center_{center_str}")
        
        # vina.run() now returns a DataFrame, but we need the output path here
        result_path = vina.output_dir

        return receptor.name, ligand.name, center, result_path

    def run_all(
        self,
        cpu: int = (os.cpu_count() or 2) - 1,
        save_to: Union[str, Path] = "./batch_vina_results",
    ) -> dict[tuple[str, str, tuple[float, float, float]], Union[Path, str]]:
        """Run AutoDock Vina docking for all ligands × all centers × all receptors in parallel.

        Args:
            cpu (int, optional): Total number of CPU cores to use.
                Defaults to ``os.cpu_count() - 1``. CPUs are divided among workers,
                with each worker receiving multiple CPUs for Vina.
            save_to (str | Path, optional): Directory where docking
                results will be stored. Defaults to "./batch_vina_results".

        Returns:
            dict[tuple[str, str, tuple[float, float, float]], Path | str]:
                Mapping from (receptor_name, ligand_name, center) to:
                - Path: Path to the docking result file, if successful.
                - str: Error message if the docking failed.
        """
        save_to = Path(save_to).expanduser().resolve()
        save_to.mkdir(parents=True, exist_ok=True)
        self.output_dir = save_to

        # Simple strategy: Divide total CPUs among workers
        # Each worker gets at least 1 CPU
        total_tasks = sum(len(self.ligands) * len(pockets) for pockets in self.receptors.values())
        
        # Calculate number of workers and CPUs per worker
        max_workers = min(cpu, total_tasks)  # Can't have more workers than CPUs or tasks
        vina_cpu = max(1, cpu // max_workers)  # At least 1 CPU per worker

        print(f"Starting AutoDock Vina docking for {total_tasks} tasks...")
        print(f"Using {max_workers} parallel workers, {vina_cpu} CPUs per worker")
        print(f"Output directory: {save_to}")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for receptor_path, pockets in self.receptors.items():
                receptor = Path(receptor_path).expanduser().resolve()
                for p in pockets:
                    center = p["center"]
                    grid_size = p["grid_size"]
                    for lig in self.ligands:
                        future = executor.submit(
                            self._dock_one, save_to, receptor, lig, center, grid_size, vina_cpu
                        )
                        futures[future] = (receptor, lig, center)

            for future in as_completed(futures):
                receptor, lig, center = futures[future]
                try:
                    rec_name, lig_name, ctr, path = future.result()
                    self.results[(rec_name, lig_name, ctr)] = path
                except Exception as e:
                    self.results[(receptor.name, lig.name, center)] = f"❌ Failed: {e}"

        print("✅ Batch processing completed!\n")
        return self.results

    def parse_results(self, save_to: Union[str, Path, None] = None) -> "pd.DataFrame":
        """Parse all Vina docking log files into a single CSV summary.

        Args:
            save_to (str | Path, optional): Path to save the CSV file.
                Defaults to `<output_dir>/vina_summary.csv`.

        Returns:
            pd.DataFrame: DataFrame containing parsed docking results.

        Raises:
            FileNotFoundError: If run_all() has not been called yet.
        """
        if not self.output_dir:
            raise FileNotFoundError("No output directory found. Call run_all() first.")

        from ..utils.parser import parse_vina_log

        if save_to is None:
            save_to = self.output_dir / "vina_summary.csv"

        log_files = [
            path / "log.txt"
            for path in self.results.values()
            if isinstance(path, Path)
        ]

        return parse_vina_log(log_file=log_files, save_to=save_to)

    def interaction_profile(
        self,
        cpu: int = (os.cpu_count() or 2) - 1,
        save_to: Union[str, Path, None] = None,
    ) -> "pd.DataFrame":
        """Compute ProLIF interaction fingerprints for all successful docking results.

        Iterates over each successful docking run and executes
        :class:`~docksuitex.interaction_profiler.InteractionProfiler`.
        The protein PDB is automatically resolved from each receptor PDBQT
        using the ``intermediate_proteins/<stem>_prepared.pdb`` convention.

        Must be called after :meth:`run_all`.

        Args:
            cpu (int, optional): Number of CPUs for ProLIF
                per docking result. Defaults to ``os.cpu_count() - 1``.
            save_to (str | Path, optional): Base directory where ProLIF
                results will be saved. Each result gets a ``prolif_results``
                subfolder inside its docking output directory. If provided,
                all results are saved under this directory instead.
                Defaults to None (saves inside each docking output directory).

        Returns:
            pd.DataFrame: Combined DataFrame of ProLIF interaction fingerprints
                from all successful docking results. Returns an empty DataFrame
                if no results could be profiled.

        Raises:
            FileNotFoundError: If run_all() has not been called yet.
        """
        if not self.output_dir:
            raise FileNotFoundError(
                "No output directory found. Call run_all() first.")

        from .batch_interaction_profiler import batch_interaction_profile

        return batch_interaction_profile(
            results=self.results,
            receptors=self.receptors,
            output_dir=self.output_dir,
            cpu=cpu,
            save_to=save_to,
        )
