"""Batch AutoDock4 docking module for DockSuiteX."""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Union, Sequence, Dict, List
import os
from ..autodock4 import AD4Docking



class BatchAD4Docking:
    """Batch docking manager for AutoDock4.

    Runs AutoDock4 docking for multiple ligands against multiple proteins,
    each with its own set of binding pocket centers, in parallel using a process pool.
    """

    def __init__(
        self,
        receptors_with_pockets: Dict[Union[str, Path], List[Dict]],
        ligands: Union[Sequence[Union[str, Path]], str, Path],
        spacing: float = 0.375,
        dielectric: float = -0.1465,
        smooth: float = 0.5,
        ga_pop_size: int = 150,
        ga_num_evals: int = 2_500_000,
        ga_num_generations: int = 27_000,
        ga_elitism: int = 1,
        ga_mutation_rate: float = 0.02,
        ga_crossover_rate: float = 0.8,
        ga_run: int = 10,
        rmstol: float = 2.0,
        seed: tuple[Union[int, str], Union[int, str]] = (27, 6)
    ):
        """Initialize a batch docking job.

        Args:
            receptors_with_pockets (Dict[str | Path, List[Dict]]):
                Dictionary mapping receptor PDBQT files to their list of binding
                pockets (each pocket MUST be a dict with 'center' and 'grid_size').
            ligands (Sequence[str | Path]): List of ligand PDBQT files.
            spacing (float, optional): Grid spacing in Å. Defaults to 0.375.
            dielectric (float, optional): Dielectric constant. Defaults to -0.1465.
            smooth (float, optional): Smoothing factor for potential maps. Defaults to 0.5.
            ga_pop_size (int, optional): Genetic algorithm population size. Defaults to 150.
            ga_num_evals (int, optional): Max energy evaluations. Defaults to 2,500,000.
            ga_num_generations (int, optional): Max generations. Defaults to 27,000.
            ga_elitism (int, optional): Elite individuals preserved. Defaults to 1.
            ga_mutation_rate (float, optional): GA mutation rate. Defaults to 0.02.
            ga_crossover_rate (float, optional): GA crossover rate. Defaults to 0.8.
            ga_run (int, optional): Independent GA runs. Defaults to 10.
            rmstol (float, optional): RMSD tolerance for clustering. Defaults to 2.0.
            seed (tuple[int | str, int | str], optional): Random seed for docking. 
                Each element can be an integer or the keywords "pid" or "time". Defaults to (27, 6).

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
        
        self.spacing = spacing
        self.dielectric = dielectric
        self.smooth = smooth
        self.ga_pop_size = ga_pop_size
        self.ga_num_evals = ga_num_evals
        self.ga_num_generations = ga_num_generations
        self.ga_elitism = ga_elitism
        self.ga_mutation_rate = ga_mutation_rate
        self.ga_crossover_rate = ga_crossover_rate
        self.ga_run = ga_run
        self.rmstol = rmstol
        self.seed = seed
        self.results: dict[tuple[str, str, tuple[float, float, float]], Union[Path, str]] = {}
        self.output_dir: Path | None = None

    def _dock_one(
        self,
        save_to: Union[str, Path],
        receptor: Path,
        ligand: Path,
        center: tuple[float, float, float],
        grid_size: tuple[float, float, float],
    ) -> tuple[str, str, tuple[float, float, float], Path]:
        """Dock a single ligand at a single pocket center of a receptor.

        Args:
            save_to (str | Path): Directory to save results.
            receptor (Path): Receptor PDBQT file.
            ligand (Path): Ligand PDBQT file.
            center (tuple[float, float, float]): Grid center coordinates.

        Returns:
            tuple: (receptor_name, ligand_name, center, result_path)
        """
        ad4 = AD4Docking(
            receptor=receptor,
            ligand=ligand,
            grid_center=center,
            grid_size=grid_size,
            spacing=self.spacing,
            dielectric=self.dielectric,
            smooth=self.smooth,
            ga_pop_size=self.ga_pop_size,
            ga_num_evals=self.ga_num_evals,
            ga_num_generations=self.ga_num_generations,
            ga_elitism=self.ga_elitism,
            ga_mutation_rate=self.ga_mutation_rate,
            ga_crossover_rate=self.ga_crossover_rate,
            ga_run=self.ga_run,
            rmstol=self.rmstol,
            seed=self.seed
        )
        
        center_str = "_".join(f"{c:.2f}" for c in center)
        ad4.run(save_to=Path(save_to) / f"{receptor.stem}_{ligand.stem}_center_{center_str}", verbose=False)
        
        # ad4.run() now returns a DataFrame, so we use ad4.output_dir for the path
        result_path = ad4.output_dir

        return receptor.name, ligand.name, center, result_path

    def run_all(
        self,
        cpu: int = (os.cpu_count() or 2) - 1,
        save_to: Union[str, Path] = "./batch_ad4_results",
    ) -> dict[tuple[str, str, tuple[float, float, float]], Union[Path, str]]:
        """Run AutoDock4 docking for all ligands × all centers × all receptors in parallel.

        Args:
            cpu (int, optional): Total number of CPU cores to use.
                Defaults to ``os.cpu_count() - 1``. Each worker uses 1 CPU
                (AutoDock4 is single-threaded).
            save_to (str | Path, optional): Directory where results are stored. 
                Defaults to "./batch_ad4_results".

        Returns:
            dict[tuple[str, str, tuple[float, float, float]], Path | str]:  
                Mapping from (receptor_name, ligand_name, center) to:
                - Path: Path to the docking result directory, if successful.
                - str: Error message if the docking failed.
        """
        save_to = Path(save_to).expanduser().resolve()
        save_to.mkdir(parents=True, exist_ok=True)
        self.output_dir = save_to

        total_tasks = sum(len(self.ligands) * len(pockets) for pockets in self.receptors.values())
        max_workers = min(cpu, total_tasks)

        print(f"Starting AutoDock4 docking for {total_tasks} tasks...")
        print(f"Using {max_workers} parallel workers, 1 CPUs per worker")
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
                            self._dock_one, save_to, receptor, lig, center, grid_size
                        )
                        futures[future] = (receptor, lig, center)

            for future in as_completed(futures):
                receptor, lig, center = futures[future]
                try:
                    rec_name, lig_name, ctr, path = future.result()
                    self.results[(rec_name, lig_name, ctr)] = path
                except Exception as e:
                    self.results[(receptor.name, lig.name, center)] = f"❌ Failed: {e}"

        print("✅ Batch processing completed!")
        return self.results

    def parse_results(self, save_to: Union[str, Path, None] = None) -> "pd.DataFrame":
        """Parse all AutoDock4 DLG result files into a single CSV summary.

        Args:
            save_to (str | Path, optional): Path to save the CSV file.
                Defaults to `<output_dir>/ad4_summary.csv`.

        Returns:
            pd.DataFrame: DataFrame containing parsed docking results.

        Raises:
            FileNotFoundError: If run_all() has not been called yet.
        """
        if not self.output_dir:
            raise FileNotFoundError("No output directory found. Call run_all() first.")

        from ..utils.parser import parse_ad4_dlg

        if save_to is None:
            save_to = self.output_dir / "ad4_summary.csv"

        dlg_files = [
            path / "results.dlg"
            for path in self.results.values()
            if isinstance(path, Path)
        ]

        return parse_ad4_dlg(dlg_file=dlg_files, save_to=save_to)

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
