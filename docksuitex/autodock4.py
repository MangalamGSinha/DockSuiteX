import subprocess
from pathlib import Path
import shutil
from typing import Union
import uuid
import os
from .protein import Protein
from .ligand import Ligand
from .utils.viewer import view_results

AUTOGRID_EXE = (Path(__file__).parent / "bin" / "autodock" / "autogrid4.exe").resolve()
AUTODOCK_EXE = (Path(__file__).parent / "bin" / "autodock" / "autodock4.exe").resolve()
TEMP_DIR = (Path(__file__).parent / "temp").resolve()
TEMP_DIR.mkdir(parents=True, exist_ok=True)


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

        # Temp directory
        self.temp_dir = TEMP_DIR / "ad4_results" / f"{self.receptor.stem}_{self.ligand.stem}_docked_ad4_{uuid.uuid4().hex[:8]}"
        self.temp_dir.mkdir(parents=True, exist_ok=True)

        shutil.copy2(self.receptor, self.temp_dir / self.receptor.name)
        shutil.copy2(self.ligand, self.temp_dir / self.ligand.name)

        self.receptor = self.temp_dir / self.receptor.name
        self.ligand = self.temp_dir / self.ligand.name

        self.gpf_file = self.temp_dir / "receptor.gpf"
        self.glg_file = self.temp_dir / "receptor.glg"
        self.dpf_file = self.temp_dir / "ligand.dpf"
        self.dlg_file = self.temp_dir / "results.dlg"

        self.receptor_types = self._detect_atom_types(self.receptor)
        self.ligand_types = self._detect_atom_types(self.ligand)


    def _setup_environment(self):
        bin_dir = str((Path(__file__).parent / "bin" / "AutoDock").resolve())
        current_path = os.environ.get("PATH", "")
        if bin_dir not in current_path:
            os.environ["PATH"] = bin_dir + os.pathsep + current_path

    def _detect_atom_types(self, path):
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
        maps_lines = "\n".join(
            f"map receptor.{t}.map" for t in self.ligand_types
        )
        content = f"""autodock_parameter_version 4.2
outlev 1
intelec
seed pid time
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





    def run(self):
        """
        Runs AutoGrid and AutoDock for molecular docking.

        Raises:
        RuntimeError: If AutoGrid or AutoDock fails, or expected output
            files (.fld or .dlg) are missing.
        """
        self._setup_environment()

        # Run AutoGrid
        self._create_gpf()
        autogrid_cmd = [str(AUTOGRID_EXE), "-p", str(self.gpf_file.name), "-l", str(self.glg_file.name)]
        result = subprocess.run(
            autogrid_cmd,
            cwd=str(self.temp_dir),
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            if self.glg_file.exists():
                raise RuntimeError(f"❌ AutoGrid failed. Log file content:\n{self.glg_file.read_text()}")
            raise subprocess.CalledProcessError(result.returncode, autogrid_cmd, result.stdout, result.stderr)

        fld_file = self.temp_dir / "receptor.maps.fld"
        if not fld_file.exists():
            raise RuntimeError("❌ AutoGrid did not create the .fld file")

        # Run AutoDock
        self._create_dpf()
        autodock_cmd = [str(AUTODOCK_EXE), "-p", str(self.dpf_file.name), "-l", str(self.dlg_file.name)]
        result = subprocess.run(
            autodock_cmd,
            cwd=str(self.temp_dir),
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            if self.dlg_file.exists():
                raise RuntimeError(f"❌ AutoDock failed. Log file content:\n{self.dlg_file.read_text()}")
            raise subprocess.CalledProcessError(result.returncode, autodock_cmd, result.stdout, result.stderr)

        self._extract_lowest_energy_conformations(self.dlg_file, Path(self.temp_dir / "output.pdbqt"))




    def view_results(self):
        """
        Visualize docking results using NGLView.

        Opens the receptor and docked ligand in an interactive
        3D widget inside a Jupyter notebook.

        Returns:
            nglview.NGLWidget: Interactive visualization of receptor–ligand complex.
        """
        view_results(protein_file=self.receptor, ligand_file=Path(self.temp_dir / "output.pdbqt"))




    def save_results(self, save_to: Union[str, Path] = Path("./ad4_results")):
        """Save docking results to a specified directory.

        Copies all temporary result files (GPF, DPF, DLG, PDBQT, logs)
        into a user-specified directory.

        Args:
            save_to (str | Path, optional): Destination directory. Defaults to "./ad4_results".

        Returns:
            Path: Resolved path to the saved results directory.

        Raises:
            RuntimeError: If docking has not been run or results are missing.
        """
        if not self.dlg_file.exists():
            raise RuntimeError(
                "❌ Docking results are missing. Please run docking before saving results.")

        save_to = Path(save_to).resolve()
        shutil.copytree(self.temp_dir, save_to, dirs_exist_ok=True)
        return save_to