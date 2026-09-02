"""Batch ProLIF interaction profiling module for DockSuiteX."""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Union, Sequence, Dict
import os
import pprint

import pandas as pd


def _profile_one(
    key: tuple[str, str, tuple[float, float, float]],
    result_path: Path,
    receptor_paths: list[Path],
    prolif_cpu: int,
    save_to: Union[str, Path, None],
) -> tuple[tuple[str, str, tuple[float, float, float]], Union[Path, str]]:
    """Profile a single docking result with ProLIF.

    Args:
        key (tuple): (receptor_name, ligand_name, center) identifier.
        result_path (Path): Directory containing the docking output.
        receptor_paths (list[Path]): Resolved receptor PDBQT paths.
        prolif_cpu (int): Number of CPUs for this ProLIF run.
        save_to (str | Path | None): Override base save directory.

    Returns:
        tuple: (key, prolif_dir_or_error_string)
    """
    from docksuitex.interaction_profiler import InteractionProfiler

    result_dir = Path(result_path)
    output_pdbqt = result_dir / "output.pdbqt"

    if not output_pdbqt.is_file():
        return key, "❌ No output.pdbqt found"

    # Derive the fixed PDB from the receptor name
    rec_name = key[0]

    original_receptor = None
    for receptor_path in receptor_paths:
        if receptor_path.name == rec_name:
            original_receptor = receptor_path
            break

    if original_receptor is None:
        return key, f"❌ Could not find original receptor for {rec_name}"

    protein_pdb = (
        original_receptor.parent
        / "intermediate_proteins"
        / f"{original_receptor.stem}_prepared.pdb"
    )

    if not protein_pdb.is_file():
        return key, (
            f"❌ Prepared protein PDB not found: {protein_pdb}. "
            "Make sure the receptor was prepared with Protein.prepare()."
        )

    try:
        prolif_dir = (Path(save_to) / result_dir.name / "prolif_results"
                      if save_to else result_dir / "prolif_results")

        profiler = InteractionProfiler(
            protein_pdb=protein_pdb,
            vina_output_pdbqt=output_pdbqt,
            _cpu=prolif_cpu,
        )
        profiler.run(save_to=prolif_dir)
        return key, prolif_dir
    except Exception as e:
        return key, f"❌ ProLIF failed: {e}"


def batch_interaction_profile(
    results: dict[tuple[str, str, tuple[float, float, float]], Union[Path, str]],
    receptors: Dict[Union[str, Path], Sequence[tuple[float, float, float]]],
    output_dir: Path,
    cpu: int = (os.cpu_count() or 2) - 1,
    save_to: Union[str, Path, None] = None,
) -> pd.DataFrame:
    """Compute ProLIF interaction fingerprints for all successful docking results.

    Iterates over each successful docking run and executes
    :class:`~docksuitex.interaction_profiler.InteractionProfiler` in parallel
    using a process pool.
    The protein PDB is automatically resolved from each receptor PDBQT
    using the ``intermediate_proteins/<stem>_prepared.pdb`` convention.

    After all runs complete, individual ``prolif_interactions.csv`` files
    are merged into a single ``prolif_interactions_combined.csv`` in the
    batch output directory, with Receptor, Ligand, and Center identifier
    columns prepended.

    Args:
        results (dict): Mapping from (receptor_name, ligand_name, center) to
            the docking result path (Path) or error message (str).
        receptors (dict): Dictionary mapping receptor PDBQT paths to their
            list of binding pocket centers (same as used in the batch class).
        output_dir (Path): Base output directory from the batch docking run.
        cpu (int, optional): Total number of CPU cores to use.
            Defaults to ``os.cpu_count() - 1``. CPUs are divided among workers,
            with each worker receiving multiple CPUs for ProLIF.
        save_to (str | Path, optional): Base directory where ProLIF
            results will be saved. Each result gets a ``prolif_results``
            subfolder inside its docking output directory. If provided,
            all results are saved under this directory instead.
            Defaults to None (saves inside each docking output directory).

    Returns:
        pd.DataFrame: Combined DataFrame of ProLIF interaction fingerprints
            from all successful docking results. Contains Receptor, Ligand,
            Center, and Pose identifier columns followed by interaction data.
            Returns an empty DataFrame if no results could be profiled.
    """
    profiling_results: dict[tuple[str, str, tuple[float, float, float]], Union[Path, str]] = {}

    # Separate successful docking results from failures
    tasks: list[tuple[tuple[str, str, tuple[float, float, float]], Path]] = []
    for key, result_path in results.items():
        if isinstance(result_path, str):
            profiling_results[key] = result_path
        else:
            tasks.append((key, Path(result_path)))

    if not tasks:
        print("⚠️ No successful docking results to profile.")
        return pd.DataFrame()

    # Pre-resolve receptor paths once for all workers
    receptor_paths = [Path(rp).expanduser().resolve() for rp in receptors]

    # Calculate workers and CPUs per worker
    total_tasks = len(tasks)
    max_workers = min(cpu, total_tasks)
    prolif_cpu = max(1, cpu // max_workers)

    print(f"Starting ProLIF interaction profiling for {total_tasks} results...")
    print(f"Using {max_workers} parallel workers, {prolif_cpu} CPUs per worker")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {}
        for key, result_path in tasks:
            future = executor.submit(
                _profile_one, key, result_path, receptor_paths, prolif_cpu, save_to,
            )
            futures[future] = key

        for future in as_completed(futures):
            key = futures[future]
            try:
                returned_key, result = future.result()
                profiling_results[returned_key] = result
            except Exception as e:
                profiling_results[key] = f"❌ ProLIF failed: {e}"

    print("✅ Batch interaction profiling completed!")
    pprint.pprint(profiling_results)

    # ------------------------------------------------------------------
    # Merge all individual prolif_interactions.csv into one combined CSV
    # ------------------------------------------------------------------
    frames: list[pd.DataFrame] = []

    for key, result in profiling_results.items():
        if not isinstance(result, Path):
            continue

        csv_path = Path(result) / "prolif_interactions.csv"
        if not csv_path.is_file():
            continue

        # Read the multi-level header (ligand, protein, interaction)
        df = pd.read_csv(csv_path, header=[0, 1, 2], index_col=0)

        # Add identifier columns
        rec_name, lig_name, center = key
        center_str = "_".join(f"{c:.2f}" for c in center)
        df.insert(0, ("Info", "Pose", ""), df.index)
        df.insert(0, ("Info", "Center", ""), center_str)
        df.insert(0, ("Info", "Ligand", ""), lig_name)
        df.insert(0, ("Info", "Receptor", ""), rec_name)

        frames.append(df)

    if frames:
        combined = pd.concat(frames, axis=0, ignore_index=True)
        combined = combined.fillna(False).infer_objects(copy=False)

        combined_path = output_dir / "prolif_interactions_combined.csv"
        combined.to_csv(combined_path, index=False)
        print(f"✅ Combined ProLIF CSV saved to: {combined_path}")
        return combined

    return pd.DataFrame()
