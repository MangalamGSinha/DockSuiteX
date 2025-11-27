"""AutoDock4 molecular docking interface.

This module provides a Python wrapper for AutoDock4 and AutoGrid, implementing
the classic genetic algorithm-based docking approach. AutoDock4 is widely used
for protein-ligand docking with detailed energy calculations.

The docking workflow:
    1. Grid map generation with AutoGrid (receptor potential maps)
    2. Genetic algorithm-based conformational search with AutoDock4
    3. Clustering of docking results by RMSD
    4. Extraction of best binding poses

AutoDock4 uses Lamarckian Genetic Algorithm (LGA) which combines genetic
algorithm with local search for efficient conformational sampling.

Example:
    Basic AutoDock4 docking::

        from docksuitex import AD4Docking

        # Initialize docking
        docking = AD4Docking(
            receptor="protein.pdbqt",
            ligand="ligand.pdbqt",
            grid_center=(10.5, 15.2, 20.8),
            grid_size=(60, 60, 60),
            ga_run=50  # 50 independent runs
        )

        # Run docking
        results_dir = docking.run(save_to="ad4_results")

        # Visualize results
        docking.view_results()

    Advanced docking with custom GA parameters::

        docking = AD4Docking(
            receptor="receptor.pdbqt",
            ligand="ligand.pdbqt",
            grid_center=(25.0, 30.0, 15.0),
            ga_pop_size=300,
            ga_num_evals=25000000,
            ga_run=100
        )
        results = docking.run()

"""

import subprocess
from pathlib import Path
import shutil
from typing import Union
import uuid
import os
from .protein import Protein
from .ligand import Ligand
from .utils.viewer import view_results

# Paths to AutoDock4 executables bundled with DockSuiteX
AUTOGRID_EXE = (Path(__file__).parent / "bin" / "autodock" / "autogrid4.exe").resolve()
AUTODOCK_EXE = (Path(__file__).parent / "bin" / "autodock" / "autodock4.exe").resolve()



class AD4Docking:
    """
    A Python wrapper for AutoDock4 and AutoGrid to automate receptor–ligand docking.

    This class automates receptor–ligand docking using AutoDock4.
    It prepares grid parameter (GPF) and docking parameter (DPF) files,
    runs AutoGrid and AutoDock, and saves docking results.
    """

    def __init__(
        self,
        receptor: Union[str, Path, "Protein"],
        ligand: Union[str, Path, "Ligand"],
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
        """
        Initialize an AutoDock4 docking run.

        Parameters
        ----------
        receptor : str | Path
            Path to the receptor PDBQT file.
        ligand : str | Path
            Path to the ligand PDBQT file.
        grid_center : tuple[float, float, float], default=(0,0,0)
            Grid box center coordinates.
        grid_size : tuple[int, int, int], default=(60,60,60)
            Number of points in the grid box.
        spacing : float, default=0.375
            Grid spacing in Å.
        dielectric : float, default=-0.1465
            Dielectric constant for electrostatics.
        smooth : float, default=0.5
            Smoothing factor for potential maps.
        ga_pop_size : int, default=150
            Genetic algorithm population size.
        ga_num_evals : int, default=2_500_000
            Maximum number of energy evaluations in GA.
        ga_num_generations : int, default=27_000
            Maximum number of generations in GA.
        ga_elitism : int, default=1
            Number of top individuals preserved during GA.
        ga_mutation_rate : float, default=0.02
            Probability of mutation in GA.
        ga_crossover_rate : float, default=0.8
            Probability of crossover in GA.
        ga_run : int, default=10
            Number of independent GA runs.
        rmstol : float, default=2.0
            RMSD tolerance for clustering docking results.
        seed : tuple[int | str, int | str], default=("pid", "time")
            Each element can be an integer or the keywords "pid" or "time".
        """
        # normalize receptor

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
        """
        Runs AutoGrid and AutoDock for molecular docking.

        Raises:
        RuntimeError: If AutoGrid or AutoDock fails, or expected output
            files (.fld or .dlg) are missing.
        """
        self._setup_environment()

        if save_to is None:
            save_to = f"ad4_docked_{self.receptor.stem}_{self.ligand.stem}"
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



    def view_results(self):
        """
        Visualize docking results using NGLView.

        Opens the receptor and docked ligand in an interactive
        3D widget inside a Jupyter notebook.

        Returns:
            nglview.NGLWidget: Interactive visualization of receptor–ligand complex.
        """
        view_results(protein_file=self.receptor, ligand_file=Path(self.output_dir / "output.pdbqt"))
