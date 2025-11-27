"""Batch AutoDock4 docking with parallel processing.

This module provides high-throughput docking capabilities using AutoDock4
with parallel execution. It handles multiple ligands, multiple receptors, and
multiple binding pockets simultaneously using process pools.

Example:
    Batch AutoDock4 docking::

        from docksuitex.batch_docking import BatchAD4Docking

        # Define receptors with their pocket centers
        receptors_with_centers = {
            "receptor1.pdbqt": [(10.0, 15.0, 20.0)],
            "receptor2.pdbqt": [(8.0, 12.0, 16.0)]
        }

        # Initialize batch docking
        batch = BatchAD4Docking(
            receptors_with_centers=receptors_with_centers,
            ligands=["lig1.pdbqt", "lig2.pdbqt"],
            ga_run=50
        )

        # Run all docking jobs in parallel
        results = batch.run_all(cpu=16, save_to="ad4_batch_results")
"""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Union, Sequence, Dict
import os
from ..autodock4 import AD4Docking



class BatchAD4Docking:
    """Batch docking manager for AutoDock4.

    Runs AutoDock4 docking for multiple ligands against multiple proteins,
    each with its own set of binding pocket centers, in parallel using a process pool.
    """

    def __init__(
        self,
        receptors_with_centers: Dict[Union[str, Path], Sequence[tuple[float, float, float]]],
        ligands: Union[Sequence[Union[str, Path]], str, Path],
        grid_size: tuple[int, int, int] = (60, 60, 60),
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
        seed: tuple[Union[int, str], Union[int, str]] = ("pid", "time")
    ):
        """Initialize a batch docking job.

        Args:
            receptors_with_centers (Dict[str | Path, Sequence[tuple[float, float, float]]]): 
                Dictionary mapping receptor PDBQT files to their list of binding pocket centers.
            ligands (Sequence[str | Path]): List of ligand PDBQT files.
            grid_size (tuple[int, int, int], optional): Number of grid points per axis. Defaults to (60, 60, 60).
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
                Each element can be an integer or the keywords "pid" or "time".
        """
        self.receptors = receptors_with_centers
        
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
        
        self.grid_size = grid_size
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

    def _dock_one(
        self,
        save_to: Union[str, Path],
        receptor: Path,
        ligand: Path,
        center: tuple[float, float, float],
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
            grid_size=self.grid_size,
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
        result_path = ad4.run(save_to=Path(save_to) / f"{receptor.stem}_{ligand.stem}_center_{center_str}")

        return receptor.name, ligand.name, center, result_path

    def run_all(
        self,
        cpu: int = os.cpu_count() or 1,
        save_to: Union[str, Path] = "./batch_ad4_results",
    ) -> dict[tuple[str, str, tuple[float, float, float]], Union[Path, str]]:
        """Run AutoDock4 docking for all ligands × all centers × all receptors in parallel.

        Args:
            cpu (int, optional): Total number of CPU cores to use.
                Defaults to all available cores. Each worker uses 1 CPU
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

        total_tasks = sum(len(self.ligands) * len(centers) for centers in self.receptors.values())
        max_workers = min(cpu, total_tasks)

        print(f"Starting AutoDock4 docking for {total_tasks} tasks...")
        print(f"Using {max_workers} parallel workers")
        print(f"Total CPU allocation: {max_workers}")
        print(f"Output directory: {save_to}")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for receptor_path, centers in self.receptors.items():
                receptor = Path(receptor_path).expanduser().resolve()
                for center in centers:
                    for lig in self.ligands:
                        future = executor.submit(self._dock_one, save_to, receptor, lig, center)
                        futures[future] = (receptor, lig, center)

            for future in as_completed(futures):
                receptor, lig, center = futures[future]
                try:
                    rec_name, lig_name, ctr, path = future.result()
                    self.results[(rec_name, lig_name, ctr)] = path
                    center_str = "_".join(f"{c:.2f}" for c in ctr)
                    print(f"✅ {rec_name} + {lig_name} @ center {center_str} → saved to {Path(path).name}")
                except Exception as e:
                    self.results[(receptor.name, lig.name, center)] = f"❌ Failed: {e}"
                    center_str = "_".join(f"{c:.2f}" for c in center)
                    print(f"❌ {receptor.name} + {lig.name} @ center {center_str} failed. Error: {e}")

        print("Batch processing completed!")
        return self.results

