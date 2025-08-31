import subprocess
from pathlib import Path
from typing import Optional, Union
import shutil
import os
import uuid
from .protein import Protein
from .ligand import Ligand
from .utils.viewer import view_results

# Path to AutoDock Vina executable
VINA_PATH = (Path(__file__).parent / "bin" / "vina" / "vina.exe").resolve()
TEMP_DIR = (Path(__file__).parent / "temp").resolve()
TEMP_DIR.mkdir(parents=True, exist_ok=True)


class VinaDocking:
    """A wrapper for running AutoDock Vina to perform molecular docking.

    This class automates running Vina, and saving/visualizing results.
    """

    def __init__(
        self,
        receptor: Union[str, Path, "Protein"],
        ligand: Union[str, Path, "Ligand"],
        grid_center: tuple[float, float, float],
        grid_size: tuple[int, int, int] = (20, 20, 20),
        exhaustiveness: int = 8,
        num_modes: int = 9,
        cpu: int = os.cpu_count() or 1,
        verbosity: int = 1,
        seed: Optional[int] = None,
    ):
        """Initialize a Vina docking job.

        Args:
            receptor (Union[str, Path, Protein]): Path to receptor `.pdbqt` file or a `Protein` object.
            ligand (Union[str, Path, Ligand]): Path to ligand `.pdbqt` file or a `Ligand` object.
            grid_center (Tuple[float, float, float]): Grid box center coordinates (x, y, z) in Å.
            grid_size (Tuple[int, int, int], optional): Grid box size along each axis in Å. Defaults to (20, 20, 20).
            exhaustiveness (int, optional): Search exhaustiveness. Higher = slower but more accurate. Defaults to 8.
            num_modes (int, optional): Maximum number of binding modes to output. Defaults to 9.
            cpu (int, optional): Number of CPU cores to use. Defaults to all available cores.
            verbosity (int, optional): Verbosity level (0=quiet, 1=normal, 2=verbose). Defaults to 1.
            seed (Optional[int], optional): Random seed for reproducibility. Defaults to None.

        Raises:
            FileNotFoundError: If receptor or ligand file is missing.
            ValueError: If input files are not `.pdbqt` or grid parameters are invalid.
            RuntimeError: If Protein or Ligand objects are not prepared.
        """


        # normalize receptor
        if isinstance(receptor, Protein):
            if receptor.pdbqt_path is None or not Path(receptor.pdbqt_path).exists():
                raise RuntimeError("❌ Protein not prepared. Run Protein.prepare() first.")
            self.receptor = Path(receptor.pdbqt_path)
        else:
            self.receptor = Path(receptor).resolve()

        # normalize ligand
        if isinstance(ligand, Ligand):
            if ligand.pdbqt_path is None or not Path(ligand.pdbqt_path).exists():
                raise RuntimeError("❌ Ligand not prepared. Run Ligand.prepare() first.")
            self.ligand = Path(ligand.pdbqt_path)
        else:
            self.ligand = Path(ligand).resolve()


        if not (isinstance(grid_center, tuple) and len(grid_center) == 3):
            raise ValueError("⚠️ 'grid_center' must be a 3-tuple of floats.")
        if not (isinstance(grid_size, tuple) and len(grid_size) == 3):
            raise ValueError("⚠️ 'grid_size' must be a 3-tuple of floats.")
        if any(not isinstance(v, (float, int)) for v in grid_center + grid_size):
            raise TypeError(
                "⚠️ Grid grid_center and grid_size values must be float or int.")


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
        self.temp_dir = TEMP_DIR / "vina_results" / f"{self.receptor.stem}_{self.ligand.stem}_docked_vina_{uuid.uuid4().hex[:8]}"
        self.temp_dir.mkdir(parents=True, exist_ok=True)

        shutil.copy2(self.receptor, self.temp_dir / self.receptor.name)
        shutil.copy2(self.ligand, self.temp_dir / self.ligand.name)

        self.receptor = self.temp_dir / self.receptor.name
        self.ligand = self.temp_dir / self.ligand.name

        # Output files
        self.output_pdbqt = self.temp_dir / f"output.pdbqt"
        self.output_log = self.temp_dir / f"log.txt"

        self._vina_output: Optional[str] = None



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


    def view_results(self):
        """Visualize docking results using NGLView.

        Returns:
            nglview.NGLWidget: Interactive 3D viewer for receptor-ligand docking.
        """
        view_results(protein_file=self.receptor, ligand_file=self.output_pdbqt)



    def save_results(self, save_to: Union[str, Path] = Path("./vina_results")) -> Path:
        """Save docking results to a specified directory.

        Args:
            save_to (Union[str, Path], optional): Destination folder. Defaults to './vina_results'.

        Returns:
            Path: Absolute path to the saved results folder.

        Raises:
            RuntimeError: If result files are missing (e.g., `run()` not called or failed).
        """
        # Check if results exist before proceeding
        if not self.output_pdbqt.exists() or not self.output_log.exists():
            raise RuntimeError(
                "❌ Docking results are missing. Please run docking before saving results."
            )
        save_to = Path(save_to).resolve()

        shutil.copytree(self.temp_dir, save_to, dirs_exist_ok=True)

        return save_to