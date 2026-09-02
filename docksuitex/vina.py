"""AutoDock Vina molecular docking module for DockSuiteX."""

import subprocess
from pathlib import Path
from typing import Optional, Union
import shutil
import os
import time
from .utils.viewer import view_docked_poses
from .platform_config import VINA_PATH


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
        prediction (e.g., using GridCalculator) or known binding site information.
    """

    def __init__(
        self,
        receptor: Union[str, Path],
        ligand: Union[str, Path],
        grid_center: Union[tuple[float, float, float], str] = "blind_docking",
        grid_size: Union[tuple[int, int, int], str] = "blind_docking",
        exhaustiveness: int = 8,
        num_modes: int = 9,
        verbosity: int = 1,
        seed: Optional[int] = 42,

        _cpu: int = (os.cpu_count() or 2) - 1,
    ):
        """Initialize a Vina docking run.

        Args:
            receptor (Union[str, Path): Path to receptor PDBQT file.
                Must be a prepared protein structure in PDBQT format.
            ligand (Union[str, Path): Path to ligand PDBQT file.
                Must be a prepared ligand structure in PDBQT format.
            grid_center (Union[tuple[float, float, float], str]): Grid box center
                coordinates (x, y, z) in Angstroms, or ``"blind_docking"`` to
                auto-compute from receptor heavy atoms. Defaults to ``"blind_docking"``.
            grid_size (Union[tuple[int, int, int], str], optional): Grid box
                dimensions (x, y, z) in Angstroms, or ``"blind_docking"`` to
                auto-compute from receptor heavy atoms.
                Defaults to ``"blind_docking"``.
            exhaustiveness (int, optional): Search exhaustiveness parameter.
                Higher values increase accuracy but also computation time.
                Typical range: 1-32. Defaults to 8.
            num_modes (int, optional): Maximum number of binding modes to generate.
                Vina will output up to this many poses ranked by predicted affinity.
                Defaults to 9.
            verbosity (int, optional): Output verbosity level.
                0 = quiet, 1 = normal, 2 = verbose. Defaults to 1.
            seed (Optional[int], optional): Random seed for reproducibility.
                If None, Vina uses a random seed. Defaults to 42.
            _cpu (int, optional): Number of CPU cores for Vina.
                Defaults to ``os.cpu_count() - 1``.

        Raises:
            FileNotFoundError: If receptor or ligand file does not exist.
            ValueError: If input files are not PDBQT format or grid parameters
                are invalid (wrong tuple size or non-numeric values).
            TypeError: If grid_center or grid_size contain non-numeric values.
        """
        self.receptor = Path(receptor).resolve()
        self._original_receptor = self.receptor
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

        # --- Blind docking: auto-compute grid from receptor heavy atoms ---
        _center_blind = isinstance(grid_center, str) and grid_center.lower() == "blind_docking"
        _size_blind   = isinstance(grid_size, str)   and grid_size.lower()   == "blind_docking"

        if _center_blind != _size_blind:
            raise ValueError(
                "⚠️ 'grid_center' and 'grid_size' must both be 'blind_docking' "
                "or both be numeric tuples.")

        if _center_blind:  # both are "blind_docking"
            from .grid_calculator import compute_blind_box
            grid_center, grid_size = compute_blind_box(self.receptor)
            print(
                f"Blind docking mode — auto-computed grid box:\n"
                f"   Center : {grid_center}\n"
                f"   Size   : {tuple(round(s, 2) for s in grid_size)}\n"
            )
        else:
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


    @staticmethod
    def _safe_copy(src: Path, dst: Path) -> Path:
        """Copy a file, handling Windows file-locking errors gracefully.

        Skips the copy if src and dst resolve to the same path. Falls back
        to shutil.copy (content-only) if shutil.copy2 raises PermissionError,
        and retries once after a brief delay if the fallback also fails.
        """
        src, dst = Path(src).resolve(), Path(dst).resolve()
        if src == dst:
            return dst
        try:
            shutil.copy2(src, dst)
        except PermissionError:
            try:
                shutil.copy(src, dst)
            except PermissionError:
                time.sleep(0.5)
                shutil.copy(src, dst)
        return dst

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

        self.receptor = self._safe_copy(self.receptor, self.output_dir / self.receptor.name)
        self.ligand = self._safe_copy(self.ligand, self.output_dir / self.ligand.name)

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
            "--energy_range", "1000",
            "--cpu", str(self.cpu),
            "--verbosity", str(self.verbosity),
        ]

        if self.seed is not None:
            cmd += ["--seed", str(self.seed)]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"❌ Error running AutoDock Vina: {result.stderr}")

        self._vina_output = result.stdout

        if self._vina_output:
            with open(self.output_log, "w") as log_file:
                log_file.write(self._vina_output)

        # Check if results exist before proceeding
        if not self.output_pdbqt.exists() or not self.output_log.exists():
            raise RuntimeError(
                "❌ Docking results are missing. Check log.txt for details."
            )

        # Remove blank lines from output.pdbqt
        lines = self.output_pdbqt.read_text().splitlines()
        cleaned = [line for line in lines if line.strip()]
        self.output_pdbqt.write_text("\n".join(cleaned) + "\n")

        print(f"✅ Vina docking completed. Results saved to: {self.output_dir}\n")
        return self.output_dir

    def parse_results(self, save_to: Union[str, Path, None] = None) -> "pd.DataFrame":
        """Parse the Vina log file and return a DataFrame of docking results.

        Delegates to :func:`docksuitex.utils.parser.parse_vina_log`
        to extract receptor/ligand names, grid parameters, and binding mode
        details (affinity, RMSD) from the log file generated by run().

        Args:
            save_to (str | Path, optional): Path to save the CSV file.
                Defaults to `<output_dir>/vina_summary.csv`.

        Returns:
            pd.DataFrame: DataFrame with columns including
                Receptor, Ligand, Mode, Affinity (kcal/mol), RMSD LB, RMSD UB, etc.

        Raises:
            FileNotFoundError: If run() has not been called yet or the
                output directory does not exist.
        """
        if not hasattr(self, 'output_dir') or not self.output_dir.exists():
            raise FileNotFoundError("No output directory found. Call run() first.")

        from .utils.parser import parse_vina_log

        if save_to is None:
            save_to = self.output_dir / "vina_summary.csv"

        return parse_vina_log(log_file=self.output_log, save_to=save_to)


    def view_grid_box(self) -> "nv.NGLWidget":
        """Visualize the docking grid box on the receptor using NGLView.

        Draws a red wireframe box representing the search space, overlaid
        on the receptor structure, with a sphere marking the grid center.
        Can be called before or after :meth:`run`.

        Returns:
            nv.NGLWidget: Interactive 3D widget (Jupyter Notebook only).
        """
        from .utils.viewer import view_grid_box
        return view_grid_box(
            protein_file=self.receptor,
            grid_center=self.grid_center,
            grid_size=self.grid_size,
        )

    def view_docked_poses(self) -> None:
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
        view_docked_poses(protein_file=self.receptor, ligand_file=self.output_pdbqt)

    def interaction_profile(
        self,
        cpu: int = (os.cpu_count() or 2) - 1,
        save_to: Union[str, Path, None] = None,
    ) -> "pd.DataFrame":
        """Compute protein–ligand interaction fingerprints with ProLIF.

        Convenience wrapper around
        :class:`~docksuitex.interaction_profiler.InteractionProfiler`.
        Must be called after :meth:`run`.

        The protein PDB is automatically resolved from the receptor PDBQT
        path using the ``intermediate_proteins/<stem>_prepared.pdb`` convention
        established by :meth:`Protein.prepare`.

        Args:
            cpu (int, optional): Number of CPUs for ProLIF.
                Defaults to ``os.cpu_count() - 1``.
            save_to (str | Path, optional): Directory where results will
                be saved. Defaults to ``<output_dir>/prolif_results``.

        Returns:
            pandas.DataFrame: Interaction fingerprint table with one row
            per ligand pose and columns for each residue–interaction pair.

        Raises:
            FileNotFoundError: If run() has not been called yet, the
                output directory does not exist, or the fixed PDB cannot
                be found.
        """
        if not hasattr(self, "output_dir") or not self.output_dir.exists():
            raise FileNotFoundError(
                "No output directory found. Call run() first.")

        # Derive the fixed PDB path from the original receptor PDBQT
        protein_pdb = (
            self._original_receptor.parent
            / "intermediate_proteins"
            / f"{self._original_receptor.stem}_prepared.pdb"
        )
        if not protein_pdb.is_file():
            raise FileNotFoundError(
                f"❌ Prepared protein PDB not found at expected location: {protein_pdb}. "
                "Make sure the receptor was prepared with Protein.prepare().")

        from .interaction_profiler import InteractionProfiler

        if save_to is None:
            save_to = self.output_dir / "prolif_results"

        self._profiler = InteractionProfiler(
            protein_pdb=protein_pdb,
            vina_output_pdbqt=self.output_pdbqt,
            _cpu=cpu,
        )
        return self._profiler.run(save_to=save_to)

    def view_interactions(self) -> None:
        """Display an interactive ProLIF interaction viewer in Jupyter.

        Renders a widget with Prev/Next buttons and a 2D/3D/Barcode toggle
        to browse protein–ligand interactions for each docked pose.

        Must be called after :meth:`interaction_profile`.

        Raises:
            AttributeError: If :meth:`interaction_profile` has not been
                called yet.
        """
        if not hasattr(self, "_profiler"):
            raise AttributeError(
                "No interaction profile found. "
                "Call interaction_profile() first."
            )
        self._profiler.view_interactions()