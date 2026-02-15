import subprocess
from pathlib import Path
import shutil
from typing import Union
import os
from .utils.viewer import view_docked_poses

AUTOGRID_EXE = (Path(__file__).parent / "bin" / "autodock" / "autogrid4.exe").resolve()
AUTODOCK_EXE = (Path(__file__).parent / "bin" / "autodock" / "autodock4.exe").resolve()


class AD4Docking:
    """
    AutoDock4 molecular docking interface.

    This module provides a Python wrapper for AutoDock4 and AutoGrid, implementing
    the classic genetic algorithm-based docking approach. AutoDock4 is widely used
    for protein-ligand docking with detailed energy calculations.

    The docking workflow:
        1. Directory setup: Creates output directory and copies input files.
        2. Grid Generation (AutoGrid):
            - Creates Grid Parameter File (GPF).
            - Runs AutoGrid to generate affinity maps (.map) and electrostatics.
        3. Docking (AutoDock):
            - Creates Docking Parameter File (DPF).
            - Runs AutoDock4 using the Lamarckian Genetic Algorithm.
        4. Result Processing:
            - Validates execution and output files.
            - Extracts docked conformations from the DLG log file.
            - Saves the best poses to a multi-model PDBQT file.

    AutoDock4 uses Lamarckian Genetic Algorithm (LGA) which combines genetic
    algorithm with local search for efficient conformational sampling.

    Note:
        Grid center coordinates should be determined from binding pocket
        prediction (e.g., using PocketFinder) or known binding site information.
    """

    def __init__(
        self,
        receptor: Union[str, Path],
        ligand: Union[str, Path],
        grid_center: tuple[float, float, float],
        grid_size: tuple[int, int, int] = (60, 60, 60),
        spacing: float = 0.375,
        dielectric: float = -0.1465,
        smooth: float = 0.5,
        # Genetic Algorithm Parameters
        ga_pop_size: int = 150,
        ga_num_evals: int = 2500000,
        ga_num_generations: int = 27000,
        ga_elitism: int = 1,
        ga_mutation_rate: float = 0.02,
        ga_crossover_rate: float = 0.8,
        ga_run: int = 10,
        rmstol: float = 2.0,
        seed: tuple[Union[int, str], Union[int, str]] = ("pid", "time")
    ):
        """Initializes an AutoDock4 docking run.

        Args:
            receptor (Union[str, Path]): Path to the receptor PDBQT file. Must be a prepared protein structure in PDBQT format.
            ligand (Union[str, Path]): Path to the ligand PDBQT file. Must be a prepared ligand structure in PDBQT format.
            grid_center (tuple[float, float, float]): Grid box center coordinates.
            grid_size (tuple[int, int, int], optional): Number of points in the grid box per axis. Defaults to (60, 60, 60).
            spacing (float, optional): Grid spacing in Å. Defaults to 0.375.
            dielectric (float, optional): Dielectric constant for electrostatics. Defaults to -0.1465.
            smooth (float, optional): Smoothing factor for potential maps. Defaults to 0.5.
            ga_pop_size (int, optional): Genetic algorithm population size. Defaults to 150.
            ga_num_evals (int, optional): Maximum number of energy evaluations in GA. Defaults to 2_500_000.
            ga_num_generations (int, optional): Maximum number of generations in GA. Defaults to 27_000.
            ga_elitism (int, optional): Number of top individuals preserved during GA. Defaults to 1.
            ga_mutation_rate (float, optional): Probability of mutation in GA. Defaults to 0.02.
            ga_crossover_rate (float, optional): Probability of crossover in GA. Defaults to 0.8.
            ga_run (int, optional): Number of independent GA runs. Defaults to 10.
            rmstol (float, optional): RMSD tolerance for clustering. Defaults to 2.0.
            seed (tuple[Union[int, str], Union[int, str]], optional): Seed for random number generation. Defaults to ("pid", "time").

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



        # Grid parameters
        self.grid_center = grid_center
        self.grid_size = grid_size
        self.spacing = spacing
        self.dielectric = dielectric
        self.smooth = smooth

        # Docking parameters
        self.ga_pop_size = ga_pop_size
        self.ga_num_evals = ga_num_evals
        self.ga_num_generations = ga_num_generations
        self.ga_elitism = ga_elitism
        self.ga_mutation_rate = ga_mutation_rate
        self.ga_crossover_rate = ga_crossover_rate
        self.ga_run = ga_run
        self.rmstol = rmstol
        self.seed = seed


    def _setup_environment(self):
        """Add AutoDock bin directory to system PATH.
        
        Ensures AutoDock executables can be found by adding the bin directory
        to the system PATH environment variable if not already present.
        """
        bin_dir = str((Path(__file__).parent / "bin" / "AutoDock").resolve())
        current_path = os.environ.get("PATH", "")
        if bin_dir not in current_path:
            os.environ["PATH"] = bin_dir + os.pathsep + current_path

    def _detect_atom_types(self, path):
        """Detect unique atom types from PDBQT file.
        
        Parses ATOM and HETATM records to extract atom type information
        from columns 77-78 of the PDBQT format.
        
        Args:
            path (Path): Path to PDBQT file.
            
        Returns:
            list: Sorted list of unique atom type strings.
        """
        atom_types = set()
        with path.open("r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    parts = line.split()
                    if len(parts) >= 3:
                        # atom_types.add(parts[-1])
                        atom_types.add(line[77:79].strip())
        return sorted(atom_types)


    def _create_gpf(self):
        """Generate AutoGrid grid parameter file (GPF).
        
        Creates a GPF file containing grid box dimensions, spacing, atom types,
        and map file specifications for AutoGrid to generate affinity maps.
        """
        maps_lines = "\n".join(
            f"map receptor.{t}.map" for t in self.ligand_types
        )
        content = f"""npts {self.grid_size[0]} {self.grid_size[1]} {self.grid_size[2]}
gridfld receptor.maps.fld
spacing {self.spacing}
receptor_types {' '.join(self.receptor_types)}
ligand_types {' '.join(self.ligand_types)}
receptor {self.receptor.name}
gridcenter {self.grid_center[0]} {self.grid_center[1]} {self.grid_center[2]}
smooth {self.smooth}
{maps_lines}
elecmap receptor.e.map
dsolvmap receptor.d.map
dielectric {self.dielectric}
"""
        self.gpf_file.write_text(content)

    def _create_dpf(self):
        """Generate AutoDock docking parameter file (DPF).
        
        Creates a DPF file containing genetic algorithm parameters, map file
        references, and search settings for AutoDock4 docking simulation.
        """
        maps_lines = "\n".join(
            f"map receptor.{t}.map" for t in self.ligand_types
        )
        seed_line = " ".join(str(s) for s in self.seed)
        content = f"""autodock_parameter_version 4.2
outlev 1
intelec
seed {seed_line}
ligand_types {' '.join(self.ligand_types)}
fld receptor.maps.fld
{maps_lines}
elecmap receptor.e.map
desolvmap receptor.d.map
move {self.ligand.name}

ga_pop_size {self.ga_pop_size}
ga_num_evals {self.ga_num_evals}
ga_num_generations {self.ga_num_generations}
ga_elitism {self.ga_elitism}
ga_mutation_rate {self.ga_mutation_rate}
ga_crossover_rate {self.ga_crossover_rate}
set_ga

sw_max_its 300
sw_max_succ 4 
sw_max_fail 4 
sw_rho 1.0
sw_lb_rho 0.01
ls_search_freq 0.06
set_psw1

ga_run {self.ga_run}
rmstol {self.rmstol}
analysis
"""
        self.dpf_file.write_text(content)

    def _extract_lowest_energy_conformations(self, dlg_file, output_pdbqt):
        """Extract docked conformations from DLG file to PDBQT format.
        
        Parses the AutoDock4 DLG (docking log) file and extracts all MODEL/ENDMDL
        blocks containing docked ligand conformations, writing them to a multi-model
        PDBQT file.
        
        Args:
            dlg_file (Path): Path to AutoDock4 DLG result file.
            output_pdbqt (Path): Path to output PDBQT file for docked poses.
        """
        with open(dlg_file, 'r') as f:
            lines = f.readlines()

        models = []
        capture = False
        current_model = []

        for line in lines:
            if line.startswith("MODEL"):
                capture = True
                current_model = [line]
            elif line.startswith("ENDMDL") and capture:
                current_model.append(line)
                models.append("".join(current_model))
                capture = False
            elif capture:
                current_model.append(line)

        if not models:
            return

        with open(output_pdbqt, 'w') as out:
            for model in models:
                out.write(model + "\n")



    def run(self, save_to: Union[str, Path] = None) -> Path:
        """Execute AutoDock4 docking simulation.

        Runs the complete AutoDock4 docking workflow, including grid map generation
        with AutoGrid and molecular docking with AutoDock4. It manages the creation
        of parameter files (GPF, DPF), executes the binaries, and processes the results.

        Args:
            save_to (Union[str, Path], optional): Directory path where docking
                results will be saved. If None, creates a directory named
                "ad4_docked_{receptor}_{ligand}_center_{x}_{y}_{z}" in the current directory.
                Defaults to None.

        Returns:
            Path: Absolute path to the output directory containing:
                - Receptor and ligand PDBQT files (copies of inputs).
                - receptor.gpf, receptor.glg: AutoGrid parameter and log files.
                - ligand.dpf, results.dlg: AutoDock parameter and log files.
                - receptor.*.map: Affinity maps generated by AutoGrid.
                - output.pdbqt: Extracted docked ligand poses (multi-model PDBQT).

        Raises:
            RuntimeError: If AutoGrid or AutoDock execution fails, or if expected
                output files (e.g., .fld, .dlg) are not created.
            subprocess.CalledProcessError: If the binary execution returns a non-zero exit code.
        """
        self._setup_environment()

        if save_to is None:
            center_str = "_".join(f"{c:.2f}" for c in self.grid_center)
            save_to = f"ad4_docked_{self.receptor.stem}_{self.ligand.stem}_center_{center_str}"
        self.output_dir = Path(save_to).resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)

        shutil.copy2(self.receptor, self.output_dir / self.receptor.name)
        shutil.copy2(self.ligand, self.output_dir / self.ligand.name)

        self.receptor = self.output_dir / self.receptor.name
        self.ligand = self.output_dir / self.ligand.name

        self.gpf_file = self.output_dir / "receptor.gpf"
        self.glg_file = self.output_dir / "receptor.glg"
        self.dpf_file = self.output_dir / "ligand.dpf"
        self.dlg_file = self.output_dir / "results.dlg"

        self.receptor_types = self._detect_atom_types(self.receptor)
        self.ligand_types = self._detect_atom_types(self.ligand)


        # Run AutoGrid
        self._create_gpf()
        autogrid_cmd = [str(AUTOGRID_EXE), "-p", str(self.gpf_file.name), "-l", str(self.glg_file.name)]
        result = subprocess.run(
            autogrid_cmd,
            cwd=str(self.output_dir),
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            if self.glg_file.exists():
                raise RuntimeError(f"❌ AutoGrid failed. Log file content:\n{self.glg_file.read_text()}")
            raise subprocess.CalledProcessError(result.returncode, autogrid_cmd, result.stdout, result.stderr)

        fld_file = self.output_dir / "receptor.maps.fld"
        if not fld_file.exists():
            raise RuntimeError("❌ AutoGrid did not create the .fld file")

        # Run AutoDock
        self._create_dpf()
        autodock_cmd = [str(AUTODOCK_EXE), "-p", str(self.dpf_file.name), "-l", str(self.dlg_file.name)]
        result = subprocess.run(
            autodock_cmd,
            cwd=str(self.output_dir),
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            if self.dlg_file.exists():
                raise RuntimeError(f"❌ AutoDock failed. Log file content:\n{self.dlg_file.read_text()}")
            raise subprocess.CalledProcessError(result.returncode, autodock_cmd, result.stdout, result.stderr)

        if not self.dlg_file.exists():
            raise RuntimeError("❌ Docking results are missing")

        self._extract_lowest_energy_conformations(self.dlg_file, Path(self.output_dir / "output.pdbqt"))

        print(f"✅ AutoDock4 completed. Results saved to: {self.output_dir}")
        return self.output_dir

    def parse_results(self, save_to: Union[str, Path, None] = None):
        """Parse the AutoDock4 DLG file and return a DataFrame of docking results.

        Delegates to :func:`docksuitex.utils.parser.parse_ad4_dlg_to_csv`
        to extract receptor/ligand names, grid parameters, GA settings,
        and cluster docking results from the DLG file generated by run().

        Args:
            save_to (str | Path, optional): Path to save the CSV file.
                Defaults to `<output_dir>/ad4_summary.csv`.

        Returns:
            pd.DataFrame: DataFrame with columns including
                Receptor, Ligand, Cluster_Rank, Binding_Energy, RMSD, etc.

        Raises:
            FileNotFoundError: If run() has not been called yet or the
                output directory does not exist.
        """
        if not hasattr(self, 'output_dir') or not self.output_dir.exists():
            raise FileNotFoundError("No output directory found. Call run() first.")

        from .utils.parser import parse_ad4_dlg_to_csv

        if save_to is None:
            save_to = self.output_dir / "ad4_summary.csv"

        return parse_ad4_dlg_to_csv(results_dir=self.output_dir, save_to=save_to)



    def view_docked_poses(self):
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
        view_docked_poses(protein_file=self.receptor, ligand_file=Path(self.output_dir / "output.pdbqt"))
