"""Grid box calculator for molecular docking.

Provides a unified :class:`GridCalculator` that supports two modes:

* ``"blind"``  — computes a bounding box encompassing all receptor heavy atoms.
* ``"p2rank"`` — runs P2Rank to predict binding pockets and computes a
  per-pocket bounding box from residue atom coordinates.
"""

import os
import subprocess
import math
import numpy as np
from pathlib import Path
from typing import List, Union, Dict, Optional
import sys
from .platform_config import P2RANK_PATH, IS_WINDOWS
import pandas as pd
from Bio import PDB


# P2RANK_PATH is imported from platform_config

# ── Geometry helpers (module-level, also used by docking engines) ─────────────

def compute_blind_box(
    pdb_file: Union[str, Path],
    padding: float = 10.0,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Compute a grid box enclosing all non-hydrogen atoms of a receptor.

    Reads ``ATOM`` / ``HETATM`` records directly from a PDB or PDBQT file
    (no BioPython needed — plain column parsing is sufficient and fast).

    Args:
        pdb_file: Path to a PDB or PDBQT file.
        padding: Extra space in Å added to each dimension. Defaults to 10.0.

    Returns:
        ``(center, size)`` — each a 3-tuple of floats in Å.

    Raises:
        ValueError: If no heavy-atom coordinates are found.
    """
    coords = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                element = line[76:78].strip()
                if element != "H":
                    try:
                        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    except ValueError:
                        continue

    if not coords:
        raise ValueError(
            f"No heavy-atom coordinates found in {pdb_file}. "
            "Cannot compute blind-docking grid box."
        )

    arr      = np.array(coords)
    min_vals = arr.min(axis=0)
    max_vals = arr.max(axis=0)
    center   = (min_vals + max_vals) / 2
    size     = (max_vals - min_vals) + padding

    return (
        tuple(round(float(v), 4) for v in center),
        tuple(round(float(v), 4) for v in size),
    )


def _get_all_atom_coords(pdb_path: Union[str, Path]) -> dict:
    """Parse a PDB file and return atom coordinates grouped by residue.

    Uses BioPython's PDBParser for robustness with insertion codes and
    non-standard residue numbering.

    Returns:
        dict mapping ``(chain_id, residue_number)`` → ``[[x, y, z], ...]``.
    """
    parser    = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))

    coords_map: dict = {}
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                res_num = residue.get_id()[1]
                key     = (chain_id, res_num)
                coords_map.setdefault(key, [])
                for atom in residue.get_atoms():
                    coords_map[key].append(atom.get_vector().get_array())
    return coords_map


def _compute_pocket_box(
    coords_map: dict,
    pocket_residues: list,
    padding: float,
) -> tuple[Optional[dict], list]:
    """Compute a bounding box for a specific set of pocket residues.

    Args:
        coords_map: Output of :func:`_get_all_atom_coords`.
        pocket_residues: List of ``(chain, residue_label)`` tuples.
        padding: Å of padding added to each side of the bounding box.

    Returns:
        ``(box_dict, missing)`` where *box_dict* has keys
        ``size_x``, ``size_y``, ``size_z`` (and ``min/max_x/y/z``),
        or ``None`` if no atoms were found.  *missing* lists residue IDs
        not found in *coords_map*.
    """
    all_coords: list = []
    missing:    list = []

    for chain, res_label in pocket_residues:
        try:
            key = (chain, int(res_label))
        except ValueError:
            continue
        atoms = coords_map.get(key)
        if atoms:
            all_coords.extend(atoms)
        else:
            missing.append(f"{chain}_{res_label}")

    if not all_coords:
        return None, missing

    arr        = np.array(all_coords)
    min_coords = arr.min(axis=0)
    max_coords = arr.max(axis=0)
    box_size   = (max_coords - min_coords) + padding

    return {
        "min_x":  round(float(min_coords[0]), 4),
        "min_y":  round(float(min_coords[1]), 4),
        "min_z":  round(float(min_coords[2]), 4),
        "max_x":  round(float(max_coords[0]), 4),
        "max_y":  round(float(max_coords[1]), 4),
        "max_z":  round(float(max_coords[2]), 4),
        "size_x": round(float(box_size[0]), 4),
        "size_y": round(float(box_size[1]), 4),
        "size_z": round(float(box_size[2]), 4),
    }, missing


# ── GridCalculator ────────────────────────────────────────────────────────────

class GridCalculator:
    """Unified grid box calculator for molecular docking.

    Supports two modes:

    ``"blind"``
        Computes a single box enclosing all receptor heavy atoms.

    ``"p2rank"``
        Runs P2Rank to predict binding pockets, then computes a
        per-pocket bounding box from residue atom coordinates.

    Both modes return a ``list[dict]``.  Each dict contains:
        - ``rank``        (int): Pocket rank (1 = best).
        - ``probability`` (float): Predicted binding probability (0–1).
          For blind mode this is always 1.0.
        - ``center``      (tuple[float, float, float]): Grid centre (x, y, z) Å.
        - ``grid_size``   (tuple[float, float, float] | None): Box dimensions Å.
    """

    DEFAULT_PADDING: float = 10.0

    def __init__(
        self,
        receptor: Union[str, Path],
        mode: str = "p2rank",
        _cpu: int = (os.cpu_count() or 2) - 1,
        padding: float = DEFAULT_PADDING,
    ):
        """Initialise GridCalculator.

        Args:
            receptor: Path to a receptor ``.pdb`` or ``.pdbqt`` file.
            mode: ``"blind"`` or ``"p2rank"``. Defaults to ``"p2rank"``.
            _cpu: CPU cores for P2Rank (ignored in blind mode).
                Defaults to ``cpu_count - 1``.
            padding: Padding in Å added to the bounding box. Defaults to 10.0.

        Raises:
            ValueError: If *mode* is not ``"blind"`` or ``"p2rank"``, or if
                the receptor file format is unsupported.
            FileNotFoundError: If the receptor file does not exist.
        """
        if mode not in ("blind", "p2rank"):
            raise ValueError(
                f"Invalid mode '{mode}'. Choose 'blind' or 'p2rank'."
            )

        self.receptor = Path(receptor).resolve()
        self.mode     = mode
        self.cpu      = _cpu
        self.padding  = padding

        if not self.receptor.is_file():
            raise FileNotFoundError(f"❌ Receptor file not found: {self.receptor}")

        if self.receptor.suffix.lower() not in (".pdb", ".pdbqt"):
            raise ValueError(
                "❌ Unsupported format. Only '.pdb' and '.pdbqt' are supported."
            )

    # ── Public API ────────────────────────────────────────────────────────────

    def run(self, save_to: Union[str, Path, None] = None) -> List[Dict]:
        """Compute the grid box(es).

        Args:
            save_to: Output directory for P2Rank files (P2Rank mode only).
                Defaults to ``./p2rank_results_<filename>`` in P2Rank mode.

        Returns:
            A ``list[dict]``, each with keys ``rank``, ``probability``,
            ``center``, ``grid_size``.  Blind mode returns a single-element
            list.
        """
        if self.mode == "blind":
            return self._run_blind()
        else:
            return self._run_p2rank(save_to)

    # ── Blind mode ────────────────────────────────────────────────────────────

    def _run_blind(self) -> List[Dict]:
        center, grid_size = compute_blind_box(self.receptor, self.padding)
        result = {
            "rank": 1,
            "probability": 1.0,
            "center": center,
            "grid_size": grid_size,
        }
        print("✅ Blind grid box computed:")
        import pprint
        pprint.pprint(result, sort_dicts=False)
        return [result]

    # ── P2Rank mode ───────────────────────────────────────────────────────────

    def _run_p2rank(self, save_to: Union[str, Path, None]) -> List[Dict]:
        if save_to is None:
            save_to = f"./p2rank_results_{self.receptor.name.replace('.', '_')}"

        output_dir = Path(save_to).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        result = subprocess.run(
            [
                str(P2RANK_PATH), "predict",
                "-f", str(self.receptor),
                "-o", str(output_dir),
                "-threads", str(self.cpu),
            ],
            shell=IS_WINDOWS,  # .bat files require shell=True on Windows
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"❌ P2Rank failed:\n{result.stderr}")

        pockets = self._parse_p2rank_output(output_dir)

        print(f"✅ P2Rank completed. Found {len(pockets)} pocket(s).")
        import pprint
        pprint.pprint(pockets, sort_dicts=False)
        return pockets

    def _resolve_pdb(self) -> Path:
        """Find the best PDB file for BioPython coordinate parsing.

        For ``.pdbqt`` inputs, tries the ``intermediate_proteins/<stem>_prepared.pdb``
        file produced by :meth:`Protein.prepare`. Falls back to the receptor file.
        """
        if self.receptor.suffix.lower() == ".pdbqt":
            fixed = (
                self.receptor.parent
                / "intermediate_proteins"
                / f"{self.receptor.stem}_prepared.pdb"
            )
            if fixed.is_file():
                return fixed
        return self.receptor

    def _parse_p2rank_output(self, output_dir: Path) -> List[Dict]:
        predictions_csv = output_dir / f"{self.receptor.name}_predictions.csv"
        residues_csv    = output_dir / f"{self.receptor.name}_residues.csv"

        if not predictions_csv.is_file():
            raise FileNotFoundError(
                f"❌ Prediction CSV not found: {predictions_csv}"
            )

        pred_df = pd.read_csv(predictions_csv, sep=r'\s*,\s*|\t', engine='python')
        pred_df.columns = pred_df.columns.str.strip()
        pred_df['name'] = pred_df['name'].astype(str).str.strip()

        res_df     = None
        coords_map = None

        if residues_csv.is_file():
            res_df = pd.read_csv(residues_csv, sep=r'\s*,\s*|\t', engine='python')
            res_df.columns          = res_df.columns.str.strip()
            res_df['chain']         = res_df['chain'].astype(str).str.strip()
            res_df['residue_label'] = res_df['residue_label'].astype(str).str.strip()
            res_df['pocket']        = res_df['pocket'].astype(str).str.strip()

            try:
                coords_map = _get_all_atom_coords(self._resolve_pdb())
            except Exception as e:
                print(f"  ⚠️ Could not parse PDB for grid-size computation: {e}")

        pockets: List[Dict] = []

        for idx, row in pred_df.iterrows():
            rank            = idx + 1
            pocket_name     = str(row.get('name', f'pocket{rank}')).strip()
            pocket_rank_str = pocket_name.replace('pocket', '')

            try:
                cx = float(row.get('center_x', 0))
                cy = float(row.get('center_y', 0))
                cz = float(row.get('center_z', 0))
            except (ValueError, TypeError) as e:
                raise ValueError(f"❌ Bad coordinates at row {rank}: {e}")

            probability = float(row.get('probability', 0))
            grid_size   = None

            if res_df is not None and coords_map is not None:
                pocket_rows     = res_df[res_df['pocket'] == pocket_rank_str]
                pocket_residues = list(zip(
                    pocket_rows['chain'],
                    pocket_rows['residue_label'],
                ))

                if pocket_residues:
                    box, missing = _compute_pocket_box(
                        coords_map, pocket_residues, self.padding
                    )
                    if missing:
                        print(f"  ⚠️ Pocket {rank}: {len(missing)} residues not found in PDB.")
                    if box is not None:
                        grid_size = (box['size_x'], box['size_y'], box['size_z'])

            pockets.append({
                "rank":        rank,
                "probability": probability,
                "center":      (cx, cy, cz),
                "grid_size":   grid_size,
            })

        if not pockets:
            raise ValueError(f"❌ No pockets found in: {predictions_csv}")

        return pockets
