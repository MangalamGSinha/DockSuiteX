from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Union, Sequence
import os
from .vina import VinaDocking


class BatchVinaDocking:
    """Batch docking manager for AutoDock Vina.

    Runs AutoDock Vina docking for multiple ligands and multiple binding pocket
    centers in parallel using a process pool.
    """

    def __init__(
        self,
        receptor: Union[str, Path],
        ligand_list: Sequence[Union[str, Path]],
        center_list: Sequence[tuple[float, float, float]],
        grid_size: tuple[int, int, int] = (20, 20, 20),
        exhaustiveness: int = 8,
        num_modes: int = 9,
        verbosity: int = 1,
        seed: int | None = None,
    ):
        """Initialize a batch Vina docking job.

        Args:
            receptor (str | Path): Path to receptor PDBQT file.
            ligand_list (Sequence[str | Path]): List of ligand PDBQT files.
            center_list (Sequence[tuple[float, float, float]]): 
                List of docking box centers.
            grid_size (tuple[int, int, int], optional):
                Dimensions of the search box in Å. Defaults to (20, 20, 20).
            exhaustiveness (int, optional):
                Sampling exhaustiveness. Higher values increase accuracy but
                also computation time. Defaults to 8.
            num_modes (int, optional):
                Maximum number of binding modes. Defaults to 9.
            verbosity (int, optional):
                Verbosity level (0 = quiet, 1 = normal, 2 = verbose).
                Defaults to 1.
            seed (int, optional):
                Random seed. If None, Vina selects automatically.
        """
        self.receptor = Path(receptor).expanduser().resolve()
        self.ligand_list = [Path(l).expanduser().resolve() for l in ligand_list]
        self.center_list = center_list
        self.grid_size = grid_size
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.seed = seed
        self.verbosity = verbosity

    def _dock_one(
        self,
        save_to: Union[str, Path],
        ligand: Path,
        center: tuple[float, float, float],
        vina_cpu: int,
    ) -> tuple[str, tuple[float, float, float], Path]:
        """Dock a single ligand at a single pocket center.

        Args:
            save_to (str | Path): Directory to save results.
            ligand (Path): Ligand PDBQT file.
            center (tuple[float, float, float]): Grid center coordinates.
            vina_cpu (int): Number of CPUs assigned to this docking job.

        Returns:
            tuple: (ligand_name, center, result_path)
        """
        vina = VinaDocking(
            receptor=self.receptor,
            ligand=ligand,
            grid_center=center,
            grid_size=self.grid_size,
            exhaustiveness=self.exhaustiveness,
            num_modes=self.num_modes,
            cpu=vina_cpu,
            verbosity=self.verbosity,
            seed=self.seed,
        )
        vina.run()

        center_str = "_".join(f"{c:.2f}" for c in center)
        result_path = vina.save_results(save_to=Path(save_to) / f"{self.receptor.stem}_{ligand.stem}_center_{center_str}")

        return ligand.name, center, result_path

    def run_all(
        self,
        cpu: int = os.cpu_count() or 1,
        save_to: Union[str, Path] = "./batch_vina_results",
    ) -> dict[tuple[str, tuple[float, float, float]], Union[Path, str]]:
        """Run AutoDock Vina docking for all ligands × all centers in parallel.

        Args:
            cpu (int, optional): Number of CPU cores to use.
                Defaults to all available cores.
            save_to (str | Path, optional): Directory where docking
                results will be stored. Defaults to "./batch_vina_results".

        Returns:
            dict[tuple[str, tuple[float, float, float]], Path | str]:
                Mapping from (ligand_name, center) to:
                - Path: Path to the docking result file, if successful.
                - str: Error message if the docking failed.
        """
        save_to = Path(save_to).expanduser().resolve()
        save_to.mkdir(parents=True, exist_ok=True)

        total_tasks = len(self.ligand_list) * len(self.center_list)
        n_jobs = min(cpu, total_tasks)
        vina_cpu = max(1, cpu // n_jobs)

        results: dict[tuple[str, tuple[float, float, float]], Union[Path, str]] = {}

        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = {
                executor.submit(self._dock_one, save_to, lig, center, vina_cpu): (lig, center)
                for lig in self.ligand_list
                for center in self.center_list
            }
            for future in as_completed(futures):
                lig, center = futures[future]
                try:
                    lig_name, ctr, path = future.result()
                    results[(lig_name, ctr)] = path
                except Exception as e:
                    results[(lig.name, center)] = f"Failed: {e}"

        return results
