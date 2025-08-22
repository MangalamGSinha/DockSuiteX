import subprocess
from pathlib import Path
import shutil
from typing import Union
import uuid
import os
import stat


AUTOGRID_EXE = (Path(__file__).parent / "bin" /
                "AutoDock" / "autogrid4.exe").resolve()
AUTODOCK_EXE = (Path(__file__).parent / "bin" /
                "AutoDock" / "autodock4.exe").resolve()
TEMP_DIR = (Path(__file__).parent / "temp").resolve()
TEMP_DIR.mkdir(parents=True, exist_ok=True)


class AD4Docking:
    """
    A Python wrapper for AutoDock4 and AutoGrid to automate receptor–ligand docking.

    This class prepares grid parameter files (GPF) and docking parameter files (DPF),
    runs AutoGrid and AutoDock, and organizes results in a temporary folder.
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
        self.receptor = Path(receptor).resolve()
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

        self.receptor_types = self._detect_receptor_types()
        self.ligand_types = self._detect_ligand_types()

    def _ensure_executable_permissions(self):
        for exe_path in [AUTOGRID_EXE, AUTODOCK_EXE]:
            if exe_path.exists():
                current_permissions = exe_path.stat().st_mode
                exe_path.chmod(current_permissions | stat.S_IEXEC |
                               stat.S_IXUSR | stat.S_IXGRP)

    def _setup_environment(self):
        bin_dir = str((Path(__file__).parent / "bin" / "AutoDock").resolve())
        current_path = os.environ.get("PATH", "")
        if bin_dir not in current_path:
            os.environ["PATH"] = bin_dir + os.pathsep + current_path

    def _detect_receptor_types(self):
        atom_types = set()
        with self.receptor.open("r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    parts = line.split()
                    if len(parts) >= 3:
                        # atom_types.add(parts[-1])
                        atom_types.add(line[77:79].strip())
        return sorted(atom_types)

    def _detect_ligand_types(self):
        atom_types = set()
        with self.ligand.open("r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    parts = line.split()
                    if len(parts) >= 3:
                        atom_types.add(parts[-1])
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

    def _run_with_error_handling(self, cmd, log_file):
        cmd_str = [str(x) for x in cmd]
        result = subprocess.run(
            cmd_str,
            cwd=str(self.temp_dir),
            capture_output=True,
            text=True,
            check=False
        )
        if result.returncode != 0:
            if log_file.exists():
                raise RuntimeError(
                    f"❌ Process failed. Log file content:\n{log_file.read_text()}")
            raise subprocess.CalledProcessError(
                result.returncode, cmd_str, result.stdout, result.stderr)
        return result

    def run(self):
        self._setup_environment()
        self._ensure_executable_permissions()

        if not AUTOGRID_EXE.exists():
            raise FileNotFoundError(
                f"❌ AutoGrid executable not found: {AUTOGRID_EXE}")
        if not AUTODOCK_EXE.exists():
            raise FileNotFoundError(
                f"❌ AutoDock executable not found: {AUTODOCK_EXE}")

        self._create_gpf()
        autogrid_cmd = [str(AUTOGRID_EXE), "-p",
                        str(self.gpf_file.name), "-l", str(self.glg_file.name)]
        self._run_with_error_handling(autogrid_cmd, self.glg_file)

        fld_file = self.temp_dir / "receptor.maps.fld"
        if not fld_file.exists():
            raise RuntimeError("❌ AutoGrid did not create the .fld file")

        self._create_dpf()
        autodock_cmd = [str(AUTODOCK_EXE), "-p",
                        str(self.dpf_file.name), "-l", str(self.dlg_file.name)]
        self._run_with_error_handling(autodock_cmd, self.dlg_file)

        self._extract_lowest_energy_conformations(
            self.dlg_file, Path(self.temp_dir / "output.pdbqt"))

    def save_results(self, save_dir: Union[str, Path] = Path("./ad4_results")):
        """
        Save docking results to a specified directory.

        Parameters
        ----------
        save_dir : str | Path, default="./ad4_results"
            Directory where results will be copied.

        Returns
        -------
        Path
            Path to the saved results.

        -------
        Raises:
            RuntimeError: If any result file is missing (docking not run or failed).
        """
        if not self.dlg_file.exists():
            raise RuntimeError(
                "❌ Docking results are missing. Please run docking before saving results.")

        save_dir = Path(save_dir).resolve()
        shutil.copytree(self.temp_dir, save_dir, dirs_exist_ok=True)
        return save_dir
