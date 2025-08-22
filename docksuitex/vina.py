import subprocess
from pathlib import Path
from typing import Optional, Union
import shutil
import os
import uuid

# Path to AutoDock Vina executable
VINA_PATH = (Path(__file__).parent / "bin" / "Vina" / "vina.exe").resolve()
TEMP_DIR = (Path(__file__).parent / "temp").resolve()
TEMP_DIR.mkdir(parents=True, exist_ok=True)


class VinaDocking:
    """
    A wrapper class for running AutoDock Vina to perform molecular docking and extract results.
    """

    def __init__(
        self,
        receptor: Union[str, Path],
        ligand: Union[str, Path],
        grid_center: tuple[float, float, float],
        grid_size: tuple[int, int, int] = (20, 20, 20),
        exhaustiveness: int = 8,
        num_modes: int = 9,
        cpu: int = os.cpu_count() or 1,
        verbosity: int = 1,
        seed: Optional[int] = None,
    ):
        """
        Initialize a Vina docking run.

        Parameters
        ----------
        receptor : str or Path
            Path to the receptor PDBQT file.
        ligand : str or Path
            Path to the ligand PDBQT file.
        grid_center : tuple of float
            Coordinates (x, y, z) of the grid box grid_center in Å.
        grid_size : tuple of int, default=(20, 20, 20)
            grid_size of the search box in Å along each dimension (x, y, z).
        exhaustiveness : int, default=8
            How exhaustively the conformational space is sampled. Higher values
            increase accuracy but also computation time.
        num_modes : int, default=9
            Maximum number of binding modes to output.
        cpu : int, default=os.cpu_count()
            Number of CPU cores to use for docking.
        verbosity : int, default=1
            Verbosity level of Vina output (0 = quiet, 1 = normal, 2 = verbose).
        seed : int, optional
            Random seed for reproducibility. If None, Vina chooses automatically.
        """

        self.receptor = Path(receptor).resolve()
        self.ligand = Path(ligand).resolve()
        self._validate_inputs(grid_center, grid_size)

        if not self.receptor.is_file():
            raise FileNotFoundError(
                f"❌ Receptor file not found: {self.receptor}")
        if not self.ligand.is_file():
            raise FileNotFoundError(f"❌ Ligand file not found: {self.ligand}")
        
        if self.receptor.suffix.lower() != ".pdbqt":
            raise ValueError("⚠️ Receptor must be a .pdbqt file.")
        if self.ligand.suffix.lower() != ".pdbqt":
            raise ValueError("⚠️ Ligand must be a .pdbqt file.")


        self.grid_center = grid_center
        self.grid_size = grid_size
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.cpu = cpu
        self.seed = seed
        self.verbosity = verbosity


        # Temp directories
        # Use receptor/ligand names + timestamp for uniqueness
        self.temp_dir = TEMP_DIR / "vina_results" / f"{self.receptor.stem}_{self.ligand.stem}_docked_ad4_{uuid.uuid4().hex[:8]}"
        self.temp_dir.mkdir(parents=True, exist_ok=True)

        shutil.copy2(self.receptor, self.temp_dir / self.receptor.name)
        shutil.copy2(self.ligand, self.temp_dir / self.ligand.name)

        self.receptor = self.temp_dir / self.receptor.name
        self.ligand = self.temp_dir / self.ligand.name

        # Output files
        self.output_pdbqt = self.temp_dir / f"output.pdbqt"
        self.output_log = self.temp_dir / f"log.txt"

        self._vina_output: Optional[str] = None

    def _validate_inputs(self, grid_center, grid_size):
        if not (isinstance(grid_center, tuple) and len(grid_center) == 3):
            raise ValueError("⚠️ 'grid_center' must be a 3-tuple of floats.")
        if not (isinstance(grid_size, tuple) and len(grid_size) == 3):
            raise ValueError("⚠️ 'grid_size' must be a 3-tuple of floats.")
        if any(not isinstance(v, (float, int)) for v in grid_center + grid_size):
            raise TypeError(
                "⚠️ Grid grid_center and grid_size values must be float or int.")

    def run(self):
        """
        Executes AutoDock Vina with the specified parameters.

        Raises:
            RuntimeError: If Vina execution fails or produces no output.
        """
        cmd = [
            str(VINA_PATH),
            "--receptor", str(self.receptor),
            "--ligand", str(self.ligand),
            "--center_x", str(self.grid_center[0]),
            "--center_y", str(self.grid_center[1]),
            "--center_z", str(self.grid_center[2]),
            "--size_x", str(self.grid_size[0]),
            "--size_y", str(self.grid_size[1]),
            "--size_z", str(self.grid_size[2]),
            "--out", str(self.output_pdbqt),
            "--exhaustiveness", str(self.exhaustiveness),
            "--num_modes", str(self.num_modes),
            "--cpu", str(self.cpu),
            "--verbosity", str(self.verbosity),
        ]

        if self.seed is not None:
            cmd += ["--seed", str(self.seed)]


        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"❌ Error running AutoDock Vina:\n{result.stderr}")

        self._vina_output = result.stdout

        if self._vina_output:
            with open(self.output_log, "w") as log_file:
                log_file.write(self._vina_output)

    def save_results(self, save_dir: Union[str, Path] = Path("./vina_results")) -> Path:
        """
        Copies the docking result files (.pdbqt, .log) to the specified directory.

        Args:
            save_dir (Union[str, Path], optional): Destination directory to save the output files.
                Defaults to './vina_results'.

        Returns:
            Path: The full resolved path of the target output folder.

        Raises:
            RuntimeError: If any result file is missing (docking not run or failed).
        """
        # Check if results exist before proceeding
        if not self.output_pdbqt.exists() or not self.output_log.exists():
            raise RuntimeError(
                "❌ Docking results are missing. Please run docking before saving results."
            )
        save_dir = Path(save_dir).resolve()

        shutil.copytree(self.temp_dir, save_dir, dirs_exist_ok=True)

        return save_dir
