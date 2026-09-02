"""AutoDock4 docking module for DockSuiteX."""
import subprocess
from pathlib import Path
import shutil
import math
from typing import Union
import os
from .utils.viewer import view_docked_poses
from .platform_config import (
    AUTOGRID_EXE, AUTODOCK_EXE, MGLTOOLS_PATH, MGL_PYTHON_EXE,
    get_mgltools_env,
)

class AD4Docking:
    """AutoDock4 molecular docking interface.

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
        prediction (e.g., using GridCalculator) or known binding site information.
    """

    def __init__(
        self,
        receptor: Union[str, Path],
        ligand: Union[str, Path],
        grid_center: Union[tuple[float, float, float], str] = "blind_docking",
        grid_size: Union[tuple[int, int, int], str] = "blind_docking",
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
        seed: tuple[Union[int, str], Union[int, str]] = (27, 6)
    ):
        """Initialize an AutoDock4 docking run.

        Args:
            receptor (Union[str, Path]): Path to the receptor PDBQT file. Must be a prepared protein structure in PDBQT format.
            ligand (Union[str, Path]): Path to the ligand PDBQT file. Must be a prepared ligand structure in PDBQT format.
            grid_center (Union[tuple[float, float, float], str]): Grid box center
                coordinates, or ``"blind_docking"`` to auto-compute from receptor
                heavy atoms. Defaults to ``"blind_docking"``.
            grid_size (Union[tuple[float, float, float], str], optional): Grid box
                dimensions in Ångströms per axis, or ``"blind_docking"`` to
                auto-compute from receptor heavy atoms.  The conversion to
                AutoGrid ``npts`` (grid points) is done internally as
                ``npts = ceil(size_Å / spacing)``.  Defaults to ``"blind_docking"``.
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
            seed (tuple[Union[int, str], Union[int, str]], optional): Seed for random number generation. Defaults to (27, 6).

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
            import pprint
            print("Blind docking mode — auto-computed grid box:")
            pprint.pprint({"center": grid_center, "grid_size": grid_size}, sort_dicts=False)
        else:
            if not (isinstance(grid_center, tuple) and len(grid_center) == 3):
                raise ValueError("⚠️ 'grid_center' must be a 3-tuple of floats.")
            if not (isinstance(grid_size, tuple) and len(grid_size) == 3):
                raise ValueError("⚠️ 'grid_size' must be a 3-tuple of floats (Å).")
            if any(not isinstance(v, (float, int)) for v in grid_center + grid_size):
                raise TypeError(
                    "⚠️ Grid grid_center and grid_size values must be float or int.")

        # Grid parameters (grid_size is always in Ångströms)
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
        # Convert Å → grid points for AutoGrid (npts = ceil(size_Å / spacing))
        npts = tuple(int(math.ceil(s / self.spacing)) for s in self.grid_size)
        content = f"""npts {npts[0]} {npts[1]} {npts[2]}
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

        def format_atom_line(parts):
            atom_id = int(parts[1])
            atom_name = parts[2]
            res_name = parts[3]
            chain = parts[4]
            res_id = int(parts[5])
            x, y, z = map(float, parts[6:9])
            charge = float(parts[11])
            element = ''.join(filter(str.isalpha, atom_name))[0].upper()

            return (
                f"{'ATOM':<6}"
                f"{atom_id:>5} "
                f"{atom_name:>4} "
                f"{res_name:>3} "
                f"{chain:1}"
                f"{res_id:>4}    "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
                f"{0.00:>6.2f}{0.00:>6.2f}          "
                f"{element:>2}"
                f"{charge:>8.3f}\n"
            )

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
                current_model.append("ENDMDL\n")
                models.append(current_model)
                capture = False

            elif capture:

                # ✅ Convert USER → REMARK
                if line.startswith("USER"):
                    current_model.append("REMARK " + line[5:])

                # ✅ Fix ATOM/HETATM formatting
                elif line.startswith(("ATOM", "HETATM")):
                    parts = line.split()
                    try:
                        formatted = format_atom_line(parts)
                        current_model.append(formatted)
                    except Exception:
                        continue  # skip bad lines safely

                # ✅ Keep TER
                elif line.startswith("TER"):
                    current_model.append("TER\n")

        if not models:
            return

        with open(output_pdbqt, 'w') as out:
            for model in models:
                for line in model:
                    out.write(line)
                out.write("\n")



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

        print(f"✅ AutoDock4 completed. Results saved to: {self.output_dir}\n")
        return self.output_dir

    def parse_results(self, save_to: Union[str, Path, None] = None) -> "pd.DataFrame":
        """Parse the AutoDock4 DLG file and return a DataFrame of docking results.

        Delegates to :func:`docksuitex.utils.parser.parse_ad4_dlg`
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

        from .utils.parser import parse_ad4_dlg

        if save_to is None:
            save_to = self.output_dir / "ad4_summary.csv"

        return parse_ad4_dlg(dlg_file=self.dlg_file, save_to=save_to)



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
        view_docked_poses(protein_file=self.receptor, ligand_file=Path(self.output_dir / "output.pdbqt"))

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
                f"❌ Fixed protein PDB not found at expected location: {protein_pdb}. "
                "Make sure the receptor was prepared with Protein.prepare().")

        from .interaction_profiler import InteractionProfiler

        if save_to is None:
            save_to = self.output_dir / "prolif_results"

        self._profiler = InteractionProfiler(
            protein_pdb=protein_pdb,
            vina_output_pdbqt=self.output_dir / "output.pdbqt",
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
