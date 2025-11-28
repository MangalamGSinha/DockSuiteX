from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Union, Sequence, Dict
import os
from ..vina import VinaDocking



class BatchVinaDocking:
    """Batch docking manager for AutoDock Vina.

    Runs AutoDock Vina docking for multiple ligands against multiple proteins,
    each with its own set of binding pocket centers, in parallel using a process pool.
    """

    def __init__(
        self,
        receptors_with_centers: Dict[Union[str, Path], Sequence[tuple[float, float, float]]],
        ligands: Union[Sequence[Union[str, Path]], str, Path],
        grid_size: tuple[int, int, int] = (20, 20, 20),
        exhaustiveness: int = 8,
        num_modes: int = 9,
        verbosity: int = 1,
        seed: int | None = None,
    ):
        """Initialize a batch Vina docking job.

        Args:
            receptors_with_centers (Dict[str | Path, Sequence[tuple[float, float, float]]]): 
                Dictionary mapping receptor PDBQT files to their list of binding pocket centers.
            ligands (Sequence[str | Path]): List of ligand PDBQT files.
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
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.seed = seed
        self.verbosity = verbosity
        self.results: dict[tuple[str, str, tuple[float, float, float]], Union[Path, str]] = {}

    def _dock_one(
        self,
        save_to: Union[str, Path],
        receptor: Path,
        ligand: Path,
        center: tuple[float, float, float],
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
            grid_size=self.grid_size,
            exhaustiveness=self.exhaustiveness,
            num_modes=self.num_modes,
            verbosity=self.verbosity,
            seed=self.seed,
            _cpu=vina_cpu,  # Internal parameter for CPU control
        )

        center_str = "_".join(f"{c:.2f}" for c in center)
        result_path = vina.run(save_to=Path(save_to) / f"{receptor.stem}_{ligand.stem}_center_{center_str}")


        return receptor.name, ligand.name, center, result_path

    def run_all(
        self,
        cpu: int = os.cpu_count() or 1,
        save_to: Union[str, Path] = "./batch_vina_results",
    ) -> dict[tuple[str, str, tuple[float, float, float]], Union[Path, str]]:
        """Run AutoDock Vina docking for all ligands × all centers × all receptors in parallel.

        Args:
            cpu (int, optional): Total number of CPU cores to use.
                Defaults to all available cores. CPUs are divided among workers,
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

        # Simple strategy: Divide total CPUs among workers
        # Each worker gets at least 1 CPU
        total_tasks = sum(len(self.ligands) * len(centers) for centers in self.receptors.values())
        
        # Calculate number of workers and CPUs per worker
        max_workers = min(cpu, total_tasks)  # Can't have more workers than CPUs or tasks
        vina_cpu = max(1, cpu // max_workers)  # At least 1 CPU per worker

        print(f"Starting AutoDock Vina docking for {total_tasks} tasks...")
        print(f"Using {max_workers} parallel workers, {vina_cpu} CPUs per worker")
        print(f"Total CPU allocation: {max_workers * vina_cpu}")
        print(f"Output directory: {save_to}")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for receptor_path, centers in self.receptors.items():
                receptor = Path(receptor_path).expanduser().resolve()
                for center in centers:
                    for lig in self.ligands:
                        future = executor.submit(self._dock_one, save_to, receptor, lig, center, vina_cpu)
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
