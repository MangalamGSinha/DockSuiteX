import subprocess
from pathlib import Path
from typing import Optional, Union
import shutil
import os
from .utils.viewer import view_results

# Path to AutoDock Vina executable bundled with DockSuiteX
VINA_PATH = (Path(__file__).parent / "bin" / "vina" / "vina.exe").resolve()


class VinaDocking:
    """Python wrapper for AutoDock Vina molecular docking.

    This class provides a high-level interface to AutoDock Vina, handling
    parameter validation, command execution, result management, and visualization.
    It automates the docking workflow from input validation to result generation.

    The docking workflow:
        1. Directory setup: Creates output directory and copies input files.
        2. Docking (Vina):
            - Constructs command line arguments including receptor, ligand, and grid parameters.
            - Runs AutoDock Vina to perform docking.
        3. Result Processing:
            - Captures and saves the execution log.
            - Validates output files (pdbqt poses and log).

    Note:
        Grid center coordinates should be determined from binding pocket
        prediction (e.g., using PocketFinder) or known binding site information.
    """

    def __init__(
        self,
        receptor: Union[str, Path],
        ligand: Union[str, Path],
        grid_center: tuple[float, float, float],
        grid_size: tuple[int, int, int] = (20, 20, 20),
        exhaustiveness: int = 8,
        num_modes: int = 9,
        verbosity: int = 1,
        seed: Optional[int] = None,

        _cpu: int = os.cpu_count() or 1,
    ):
        """Initialize a Vina docking run.

        Args:
            receptor (Union[str, Path): Path to receptor PDBQT file.
                Must be a prepared protein structure in PDBQT format.
            ligand (Union[str, Path): Path to ligand PDBQT file.
                Must be a prepared ligand structure in PDBQT format.
            grid_center (tuple[float, float, float]): Grid box center coordinates
                (x, y, z) in Angstroms. Defines the center of the search space.
            grid_size (tuple[int, int, int], optional): Grid box dimensions
                (x, y, z) in Angstroms. Defines the size of the search space.
                Defaults to (20, 20, 20).
            exhaustiveness (int, optional): Search exhaustiveness parameter.
                Higher values increase accuracy but also computation time.
                Typical range: 1-32. Defaults to 8.
            num_modes (int, optional): Maximum number of binding modes to generate.
                Vina will output up to this many poses ranked by predicted affinity.
                Defaults to 9.
            verbosity (int, optional): Output verbosity level.
                0 = quiet, 1 = normal, 2 = verbose. Defaults to 1.
            seed (Optional[int], optional): Random seed for reproducibility.
                If None, Vina uses a random seed. Defaults to None.

        Raises:
            FileNotFoundError: If receptor or ligand file does not exist.
            ValueError: If input files are not PDBQT format or grid parameters
                are invalid (wrong tuple size or non-numeric values).
            TypeError: If grid_center or grid_size contain non-numeric values.
        """

        self.receptor = Path(receptor).resolve()
        self.ligand = Path(ligand).resolve()

        if not self.receptor.is_file():
            raise FileNotFoundError(
                f"❌ Receptor file not found: {self.receptor}")
        if not self.ligand.is_file():
            raise FileNotFoundError(f"❌ Ligand file not found: {self.ligand}")
        
        if self.receptor.suffix.lower() != ".pdbqt":
            raise ValueError("⚠️ Receptor must be a .pdbqt file.")
        if self.ligand.suffix.lower() != ".pdbqt":
            raise ValueError("⚠️ Ligand must be a .pdbqt file.")

        if not (isinstance(grid_center, tuple) and len(grid_center) == 3):
            raise ValueError("⚠️ 'grid_center' must be a 3-tuple of floats.")
        if not (isinstance(grid_size, tuple) and len(grid_size) == 3):
            raise ValueError("⚠️ 'grid_size' must be a 3-tuple of floats.")
        if any(not isinstance(v, (float, int)) for v in grid_center + grid_size):
            raise TypeError(
                "⚠️ Grid grid_center and grid_size values must be float or int.")

        self.grid_center = grid_center
        self.grid_size = grid_size
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.cpu = _cpu
        self.seed = seed
        self.verbosity = verbosity


    def run(self, save_to: Union[str, Path] = None) -> Path:
        """Execute AutoDock Vina docking simulation.

        Runs the Vina docking calculation with the configured parameters,
        saves results to the specified directory, and generates output files
        including docked poses and a log file with binding energies.

        Args:
            save_to (Union[str, Path], optional): Directory path where docking
                results will be saved. If None, creates a directory named
                "vina_docked_{receptor}_{ligand}_center_{x}_{y}_{z}" in the current directory.
                Defaults to None.

        Returns:
            Path: Absolute path to the output directory containing:
                - Receptor and ligand PDBQT files (copies of inputs)
                - output.pdbqt: Docked ligand poses ranked by affinity
                - log.txt: Vina output log with binding energies and RMSD values

        Raises:
            RuntimeError: If Vina execution fails (non-zero return code) or
                if expected output files (output.pdbqt, log.txt) are not created.

        """

        if save_to is None:
            center_str = "_".join(f"{c:.2f}" for c in self.grid_center)
            save_to = f"vina_docked_{self.receptor.stem}_{self.ligand.stem}_center_{center_str}"
        self.output_dir = Path(save_to).resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)

        shutil.copy2(self.receptor, self.output_dir / self.receptor.name)
        shutil.copy2(self.ligand, self.output_dir / self.ligand.name)

        self.receptor = self.output_dir / self.receptor.name
        self.ligand = self.output_dir / self.ligand.name

        # Output files
        self.output_pdbqt = self.output_dir / f"output.pdbqt"
        self.output_log = self.output_dir / f"log.txt"

        self._vina_output: Optional[str] = None

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

        # Check if results exist before proceeding
        if not self.output_pdbqt.exists() or not self.output_log.exists():
            raise RuntimeError(
                "❌ Docking results are missing. Check log.txt for details."
            )

        print(f"✅ Vina docking completed. Results saved to: {self.output_dir}")
        return self.output_dir


    def view_results(self):
        """Visualize docking results using NGLView in Jupyter Notebook.

        Creates an interactive 3D visualization of the receptor-ligand complex
        with controls to browse through different docking poses. This method
        requires a Jupyter Notebook environment and the nglview package.

        The visualization includes:
            - Receptor structure displayed as cartoon representation
            - Ligand poses displayed as ball-and-stick models
            - Interactive controls to step through poses
            - Play/pause animation of poses
            - Speed control slider

        Returns:
            None: Displays the interactive widget directly in the notebook.

        Raises:
            AttributeError: If run() has not been called yet (output files not set).
            FileNotFoundError: If output files have been deleted or moved.
            ImportError: If nglview is not installed.

        Note:
            This method must be called after run() to ensure output files exist.
            It only works in Jupyter Notebook/Lab environments.

        """
        view_results(protein_file=self.receptor, ligand_file=self.output_pdbqt)