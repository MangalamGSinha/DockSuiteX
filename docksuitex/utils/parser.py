import os
import re
import pandas as pd
from pathlib import Path


def parse_vina_log_to_csv(results_dir: str, output_csv: str = "vina_summary.csv") -> None:
    """
    Parse AutoDock Vina log files from multiple docking result folders and extract 
    docking results into a single CSV file.
    
    Parameters
    ----------
    results_dir : str
        Path to the parent directory containing all docking result folders.
        Each folder must contain a Vina log file (.txt).
    output_csv : str, optional
        Path to the output CSV file. Default is "vina_summary.csv".
    
    Output
    -------
    A CSV file containing receptor name, ligand name, docking scores, grid parameters, and
    search settings.
    """

    results = []

    for log_file in Path(results_dir).rglob("*.txt"):
        with open(log_file, "r") as f:
            text = f.read()

        # Extract receptor and ligand names from log
        receptor_match = re.search(r"Rigid receptor:\s*(.+\.pdbqt)", text)
        ligand_match = re.search(r"Ligand:\s*(.+\.pdbqt)", text)

        receptor_name = Path(receptor_match.group(1)).stem if receptor_match else "Unknown"
        ligand_name = Path(ligand_match.group(1)).stem if ligand_match else "Unknown"

        # Extract grid and parameters
        grid_center = re.search(r"Grid center:\s*X\s*([-\d.]+)\s*Y\s*([-\d.]+)\s*Z\s*([-\d.]+)", text)
        grid_size = re.search(r"Grid size\s*:\s*X\s*(\d+)\s*Y\s*(\d+)\s*Z\s*(\d+)", text)
        grid_space = re.search(r"Grid space\s*:\s*([-\d.]+)", text)
        exhaustiveness = re.search(r"Exhaustiveness:\s*(\d+)", text)

        # Extract docking results (mode, affinity, rmsd)
        docking_results = re.findall(r"^\s*(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", text, re.MULTILINE)

        for mode, affinity, rmsd_lb, rmsd_ub in docking_results:
            results.append({
                "Receptor": receptor_name,
                "Ligand": ligand_name,
                "Mode": int(mode),
                "Affinity (kcal/mol)": float(affinity),
                "RMSD LB": float(rmsd_lb),
                "RMSD UB": float(rmsd_ub),
                "Grid Center X": float(grid_center.group(1)) if grid_center else None,
                "Grid Center Y": float(grid_center.group(2)) if grid_center else None,
                "Grid Center Z": float(grid_center.group(3)) if grid_center else None,
                "Grid Size X": int(grid_size.group(1)) if grid_size else None,
                "Grid Size Y": int(grid_size.group(2)) if grid_size else None,
                "Grid Size Z": int(grid_size.group(3)) if grid_size else None,
                "Grid Spacing": float(grid_space.group(1)) if grid_space else None,
                "Exhaustiveness": int(exhaustiveness.group(1)) if exhaustiveness else None,
            })

    # Save to CSV
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    return df






def parse_ad4_dlg_to_csv(results_dir, output_csv="ad4_summary.csv"):
    """
    Parse AutoDock4 .dlg files from multiple docking result folders and extract 
    cluster-level docking results to a single CSV.

    Parameters
    ----------
    results_dir : str or Path
        Path to the parent directory containing all docking result folders.
        Each folder must contain a 'results.dlg' file.
    output_csv : str, optional
        Path to save the output CSV file. Default is 'ad4_summary.csv'.

    Output
    -------
    A CSV file containing cluster-level docking results including:
    Cluster Rank, receptor, ligand, binding energy, RMSD, grid center,
    grid size, spacing, and GA parameters.
    """
    
    results_dir = Path(results_dir)
    dlg_files = list(results_dir.glob("*/results.dlg"))
    if not dlg_files:
        raise FileNotFoundError(f"No results.dlg files found in subfolders of {results_dir}")

    all_data = []

    for dlg_file in dlg_files:
        with open(dlg_file, 'r') as f:
            lines = f.readlines()

        data = []

        # Initialize variables
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
            if 'Ligand PDBQT file' in line:
                ligand_full = re.search(r'"(.+?)"', line).group(1)
                ligand = Path(ligand_full).stem
            if 'Macromolecule file used to create Grid Maps' in line:
                receptor_full = line.split('=')[-1].strip()
                receptor = Path(receptor_full).stem

            # Grid spacing
            if 'Grid Point Spacing' in line:
                spacing = float(re.search(r'[\d.]+', line).group(0))

            # Grid size
            if 'Even Number of User-specified Grid Points' in line:
                for j in range(i, i+3):
                    s = lines[j]
                    if 'x-points' in s:
                        size[0] = int(re.search(r'(\d+)\s+x-points', s).group(1))
                    if 'y-points' in s:
                        size[1] = int(re.search(r'(\d+)\s+y-points', s).group(1))
                    if 'z-points' in s:
                        size[2] = int(re.search(r'(\d+)\s+z-points', s).group(1))
                    if all(x is not None for x in size):
                        break

            # Grid center
            if 'Coordinates of Central Grid Point of Maps' in line:
                center_vals = re.findall(r'[-\d.]+', line)
                if len(center_vals) >= 3:
                    center = [float(v) for v in center_vals[:3]]

            # GA parameters
            for key in ga_params.keys():
                if line.startswith(f'DPF> {key}'):
                    match = re.search(r'[\d.]+', line)
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
                    "Cluster_Rank": None,
                    "RMSD": None,
                    "Binding_Energy": None,
                    "Center_X": center[0],
                    "Center_Y": center[1],
                    "Center_Z": center[2],
                    "Size_X": size[0],
                    "Size_Y": size[1],
                    "Size_Z": size[2],
                    "Spacing": spacing,
                    **ga_params
                }

            # Extract Cluster Rank
            if in_cluster_section and 'Cluster Rank' in line:
                match = re.search(r'Cluster Rank\s*=\s*(\d+)', line)
                if match:
                    cluster_info["Cluster_Rank"] = int(match.group(1))

            # RMSD
            if in_cluster_section and 'RMSD from reference structure' in line:
                match = re.search(r'([\d.]+)', line)
                if match:
                    cluster_info["RMSD"] = float(match.group(1))

            # Binding energy
            if in_cluster_section and 'Estimated Free Energy of Binding' in line:
                match = re.search(r'([-+]?\d*\.\d+|\d+)', line)
                if match:
                    cluster_info["Binding_Energy"] = float(match.group(1))

            # End of model
            if in_cluster_section and line.startswith("ENDMDL"):
                data.append(cluster_info)

        all_data.extend(data)

    # Save to CSV
    df = pd.DataFrame(all_data)
    df.to_csv(output_csv, index=False)
    return df


