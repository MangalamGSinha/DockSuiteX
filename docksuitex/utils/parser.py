import re
import pandas as pd
from pathlib import Path
from typing import Union


def parse_vina_log_to_csv(
    results_dir: Union[str, Path],
    output_csv: Union[str, Path] = "vina_summary.csv"
) -> pd.DataFrame:
    """Parse AutoDock Vina log files into a structured CSV summary.

    This function scans a directory tree for Vina docking log files (`log.txt`),
    extracts receptor/ligand names, grid box parameters, and docking results
    (affinity and RMSD values), then writes the results to a CSV file.

    Args:
        results_dir (str | Path): Path to the directory containing docking result subfolders.
            Each folder should include a `log.txt` file from AutoDock Vina.
        output_csv (str | Path, optional): Path to save the generated CSV summary.
            Defaults to "vina_summary.csv" in the current directory.

    Returns:
        pd.DataFrame: DataFrame containing parsed docking results with columns:
            - Receptor, Ligand
            - Mode, Affinity (kcal/mol), RMSD LB, RMSD UB
            - Grid Center (X, Y, Z), Grid Size (X, Y, Z), Grid Spacing
            - Exhaustiveness

    Raises:
        FileNotFoundError: If no `log.txt` files are found under `results_dir`.
        RuntimeError: If parsing fails or no docking results are extracted.
    """
    results_dir = Path(results_dir).expanduser().resolve()
    log_files = list(results_dir.rglob("log.txt"))
    if not log_files:
        raise FileNotFoundError(f"No Vina log.txt files found in {results_dir}")

    results = []

    for log_file in log_files:
        with open(log_file, "r", encoding="utf-8") as f:
            text = f.read()

        # Extract receptor and ligand names
        receptor_match = re.search(r"Rigid receptor:\s*(.+\.pdbqt)", text)
        ligand_match = re.search(r"Ligand:\s*(.+\.pdbqt)", text)
        receptor_name = Path(receptor_match.group(1)).stem if receptor_match else "Unknown"
        ligand_name = Path(ligand_match.group(1)).stem if ligand_match else "Unknown"

        # Extract grid and parameters
        grid_center = re.search(r"Grid center:\s*X\s*([-\d.]+)\s*Y\s*([-\d.]+)\s*Z\s*([-\d.]+)", text)
        grid_size = re.search(r"Grid size\s*:\s*X\s*([-\d.]+)\s*Y\s*([-\d.]+)\s*Z\s*([-\d.]+)", text)
        grid_space = re.search(r"Grid space\s*:\s*([-\d.]+)", text)
        exhaustiveness = re.search(r"Exhaustiveness:\s*(\d+)", text)

        # Extract docking results table
        docking_results = re.findall(
            r"^\s*(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", text, re.MULTILINE
        )

        for mode, affinity, rmsd_lb, rmsd_ub in docking_results:
            results.append({
                "Receptor": receptor_name,
                "Ligand": ligand_name,
                "Grid Center X": float(grid_center.group(1)) if grid_center else None,
                "Grid Center Y": float(grid_center.group(2)) if grid_center else None,
                "Grid Center Z": float(grid_center.group(3)) if grid_center else None,
                "Grid Size X": float(grid_size.group(1)) if grid_size else None,
                "Grid Size Y": float(grid_size.group(2)) if grid_size else None,
                "Grid Size Z": float(grid_size.group(3)) if grid_size else None,
                "Grid Spacing": float(grid_space.group(1)) if grid_space else None,
                "Exhaustiveness": int(exhaustiveness.group(1)) if exhaustiveness else None,
                "Mode": int(mode),
                "Affinity (kcal/mol)": float(affinity),
                "RMSD LB": float(rmsd_lb),
                "RMSD UB": float(rmsd_ub),
            })

    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    return df


def parse_ad4_dlg_to_csv(
    results_dir: Union[str, Path],
    output_csv: Union[str, Path] = "ad4_summary.csv"
) -> pd.DataFrame:
    """Parse AutoDock4 DLG result files into a structured CSV summary.

    This function scans a directory tree for AutoDock4 docking log files
    (`results.dlg`), extracts receptor and ligand names, grid box parameters,
    genetic algorithm (GA) settings, and cluster docking results, then writes
    them to a CSV file.

    Args:
        results_dir (str | Path): Path to the directory containing docking result subfolders.
            Each folder should include a `results.dlg` file from AutoDock4.
        output_csv (str | Path, optional): Path to save the generated CSV summary.
            Defaults to "ad4_summary.csv" in the current directory.

    Returns:
        pd.DataFrame: DataFrame containing parsed docking results with columns:
            - Receptor, Ligand
            - Cluster_Rank, RMSD, Binding_Energy
            - Grid Center (X, Y, Z), Grid Size (X, Y, Z), Spacing
            - GA parameters (e.g., rmstol, ga_pop_size, ga_num_evals, etc.)

    Raises:
        FileNotFoundError: If no `results.dlg` files are found under `results_dir`.
        RuntimeError: If parsing fails or cluster information cannot be extracted.
    """
    results_dir = Path(results_dir).expanduser().resolve()
    dlg_files = list(results_dir.rglob("results.dlg"))
    if not dlg_files:
        raise FileNotFoundError(f"No results.dlg files found in {results_dir}")

    all_data = []

    for dlg_file in dlg_files:
        with open(dlg_file, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        receptor = ligand = None
        center = [None, None, None]
        size = [None, None, None]
        spacing = None
        ga_params = {
            "rmstol": None,
            "ga_pop_size": None,
            "ga_num_evals": None,
            "ga_num_generations": None,
            "ga_elitism": None,
            "ga_mutation_rate": None,
            "ga_crossover_rate": None,
            "ga_run": None,
        }

        in_cluster_section = False
        cluster_info = {}

        for i, line in enumerate(lines):
            line = line.strip()

            # Ligand and receptor
            if "Ligand PDBQT file" in line:
                match = re.search(r'"(.+?)"', line)
                ligand = Path(match.group(1)).stem if match else None
            if "Macromolecule file used to create Grid Maps" in line:
                receptor = Path(line.split("=")[-1].strip()).stem

            # Grid spacing
            if "Grid Point Spacing" in line:
                match = re.search(r"[\d.]+", line)
                spacing = float(match.group(0)) if match else None

            # Grid size (x/y/z points over next few lines)
            if "Even Number of User-specified Grid Points" in line:
                for j in range(i, i + 3):
                    s = lines[j]
                    if "x-points" in s:
                        size[0] = int(re.search(r"(\d+)", s).group(1))
                    if "y-points" in s:
                        size[1] = int(re.search(r"(\d+)", s).group(1))
                    if "z-points" in s:
                        size[2] = int(re.search(r"(\d+)", s).group(1))

            # Grid center
            if "Coordinates of Central Grid Point of Maps" in line:
                vals = re.findall(r"[-\d.]+", line)
                if len(vals) >= 3:
                    center = [float(v) for v in vals[:3]]

            # GA parameters
            for key in ga_params.keys():
                if line.startswith(f"DPF> {key}"):
                    match = re.search(r"[\d.]+", line)
                    if match:
                        ga_params[key] = float(match.group(0))

            # Cluster section
            if "LOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER" in line:
                in_cluster_section = True
                continue

            if in_cluster_section and line.startswith("MODEL"):
                cluster_info = {
                    "Receptor": receptor,
                    "Ligand": ligand,
                    "Center_X": center[0],
                    "Center_Y": center[1],
                    "Center_Z": center[2],
                    "Size_X": size[0],
                    "Size_Y": size[1],
                    "Size_Z": size[2],
                    "Spacing": spacing,
                    **ga_params,
                    "Cluster_Rank": None,
                    "RMSD": None,
                    "Binding_Energy": None,
                }

            if in_cluster_section and "Cluster Rank" in line:
                match = re.search(r"Cluster Rank\s*=\s*(\d+)", line)
                if match:
                    cluster_info["Cluster_Rank"] = int(match.group(1))

            if in_cluster_section and "RMSD from reference structure" in line:
                match = re.search(r"([\d.]+)", line)
                if match:
                    cluster_info["RMSD"] = float(match.group(1))

            if in_cluster_section and "Estimated Free Energy of Binding" in line:
                match = re.search(r"([-+]?\d*\.\d+|\d+)", line)
                if match:
                    cluster_info["Binding_Energy"] = float(match.group(1))

            if in_cluster_section and line.startswith("ENDMDL"):
                all_data.append(cluster_info)

    df = pd.DataFrame(all_data)
    df.to_csv(output_csv, index=False)
    return df
