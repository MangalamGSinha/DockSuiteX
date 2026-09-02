"""Output parser utilities for DockSuiteX log files."""

import re
import pandas as pd
from pathlib import Path
from typing import Union, List



def parse_vina_log(
    log_file: Union[str, Path, List[Union[str, Path]]],
    save_to: Union[str, Path] = "vina_summary.csv"
) -> pd.DataFrame:
    """Parse AutoDock Vina log file(s) into a structured CSV summary.

    Extracts receptor/ligand names, grid box parameters, and docking results
    (affinity and RMSD values) from one or more Vina log files, then writes
    the results to a CSV file.

    Args:
        log_file (Union[str, Path, List[Union[str, Path]]]): Path to a single
            Vina log file, or a list of paths to multiple log files.
        save_to (Union[str, Path], optional): Path to save the generated CSV
            summary. Defaults to ``"vina_summary.csv"``.

    Returns:
        pd.DataFrame: DataFrame containing parsed docking results with columns:
            - Receptor, Ligand
            - Mode, Affinity (kcal/mol), RMSD LB, RMSD UB
            - Grid Center (X, Y, Z), Grid Size (X, Y, Z), Grid Spacing
            - Exhaustiveness

    Raises:
        FileNotFoundError: If any of the specified log files do not exist.
    """
    # Normalise to a list of Path objects
    if isinstance(log_file, (str, Path)):
        log_files = [Path(log_file).expanduser().resolve()]
    else:
        log_files = [Path(f).expanduser().resolve() for f in log_file]

    for lf in log_files:
        if not lf.is_file():
            raise FileNotFoundError(f"❌ Vina log file not found: {lf}")

    print(f"Starting Vina log parsing for {len(log_files)} file(s)...")
    results = []

    for lf in log_files:
        with open(lf, "r", encoding="utf-8") as f:
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

        print(f"✅ Parsed {lf} → {len(docking_results)} modes")

    df = pd.DataFrame(results)
    df.to_csv(save_to, index=False)
    print(f"✅ Parsing completed! Saved {len(results)} results to {save_to}\n")
    return df


def parse_ad4_dlg(
    dlg_file: Union[str, Path, List[Union[str, Path]]],
    save_to: Union[str, Path] = "ad4_summary.csv"
) -> pd.DataFrame:
    """Parse AutoDock4 DLG result file(s) into a structured CSV summary.

    Extracts receptor and ligand names, grid box parameters, genetic algorithm
    (GA) settings, and cluster docking results from one or more DLG files, then
    writes them to a CSV file.

    Args:
        dlg_file (Union[str, Path, List[Union[str, Path]]]): Path to a single
            AutoDock4 DLG file, or a list of paths to multiple DLG files.
        save_to (Union[str, Path], optional): Path to save the generated CSV
            summary. Defaults to ``"ad4_summary.csv"``.

    Returns:
        pd.DataFrame: DataFrame containing parsed docking results with columns:
            - Receptor, Ligand
            - Cluster_Rank, RMSD, Binding_Energy
            - Grid Center (X, Y, Z), Grid Size (X, Y, Z), Spacing
            - GA parameters (e.g., rmstol, ga_pop_size, ga_num_evals, etc.)

    Raises:
        FileNotFoundError: If any of the specified DLG files do not exist.
    """
    # Normalise to a list of Path objects
    if isinstance(dlg_file, (str, Path)):
        dlg_files = [Path(dlg_file).expanduser().resolve()]
    else:
        dlg_files = [Path(f).expanduser().resolve() for f in dlg_file]

    for df_path in dlg_files:
        if not df_path.is_file():
            raise FileNotFoundError(f"❌ DLG file not found: {df_path}")

    print(f"Starting AutoDock4 DLG parsing for {len(dlg_files)} file(s)...")
    all_data = []

    for df_path in dlg_files:
        with open(df_path, "r", encoding="utf-8", errors="ignore") as f:
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
                for j in range(i, min(i + 3, len(lines))):
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

        print(f"✅ Parsed {df_path} → {len([d for d in all_data if d.get('Ligand') == ligand])} clusters")

    df = pd.DataFrame(all_data)
    df.to_csv(save_to, index=False)
    print(f"✅ Parsing completed! Saved {len(all_data)} cluster results to {save_to}\n")
    return df
