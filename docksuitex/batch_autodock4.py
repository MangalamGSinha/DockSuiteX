from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Union, Sequence
import os
from .autodock4 import AD4Docking


class BatchAD4Docking:
    """Batch docking manager for AutoDock4.

    Runs AutoDock4 docking for multiple ligands and multiple binding pocket
    centers in parallel, using a process pool.
    """

    def __init__(
        self,
        receptor: Union[str, Path],
        ligand_list: Sequence[Union[str, Path]],
        center_list: Sequence[tuple[float, float, float]],
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
    ):
        """Initialize a batch docking job.

        Args:
            receptor (str | Path): Path to receptor PDBQT file.
            ligand_list (Sequence[str | Path]): List of ligand PDBQT files.
            center_list (Sequence[tuple[float, float, float]]): 
                List of docking box centers.
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
        """
        self.receptor = Path(receptor).expanduser().resolve()
        self.ligand_list = [Path(l).expanduser().resolve() for l in ligand_list]
        self.center_list = center_list
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

    def _dock_one(
        self,
        save_to: Union[str, Path],
        ligand: Path,
        center: tuple[float, float, float],
    ) -> tuple[str, tuple[float, float, float], Path]:
        """Dock a single ligand at a single pocket center.

        Args:
            save_to (str | Path): Directory to save results.
            ligand (Path): Ligand PDBQT file.
            center (tuple[float, float, float]): Grid center coordinates.

        Returns:
            tuple: (ligand_name, center, result_path)
        """
        ad4 = AD4Docking(
            receptor=self.receptor,
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
        )
        ad4.run()
        
        center_str = "_".join(f"{c:.2f}" for c in center)
        result_path = ad4.save_results(save_to=Path(save_to) / f"{self.receptor.stem}_{ligand.stem}_center_{center_str}")

        return ligand.name, center, result_path

    def run_all(
        self,
        cpu: int = os.cpu_count() or 1,
        save_to: Union[str, Path] = "./batch_ad4_results",
    ) -> dict[tuple[str, tuple[float, float, float]], Union[Path, str]]:
        """Run AutoDock4 docking for all ligands × all centers in parallel.

        Args:
            cpu (int, optional): Number of CPU cores to use. Defaults to all available cores.
            save_to (str | Path, optional): Directory where results are stored. 
                Defaults to "./batch_ad4_results".

        Returns:
            dict[tuple[str, tuple[float, float, float]], Path | str]:  
                Mapping from (ligand_name, center) →  
                - Path to result file if successful,  
                - Error message if docking failed.
        """
        save_to = Path(save_to).expanduser().resolve()
        save_to.mkdir(parents=True, exist_ok=True)

        max_workers = min(cpu, len(self.ligand_list) * len(self.center_list))
        results = {}

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(self._dock_one, save_to, lig, center): (lig, center)
                for lig in self.ligand_list
                for center in self.center_list
            }

            for future in as_completed(futures):
                lig, center = futures[future]
                try:
                    lig_name, ctr, path = future.result()
                    results[(lig_name, ctr)] = path
                except Exception as e:
                    results[(lig.name, center)] = f"❌ Failed: {e}"

        return results
