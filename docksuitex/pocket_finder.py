import os
import subprocess
import csv
import shutil
from pathlib import Path
from typing import Optional, Tuple, List, Union, Dict
import uuid
from docksuitex.protein import Protein

# Paths
P2RANK_PATH = (Path(__file__).parent / "bin" / "p2rank" / "prank.bat").resolve()
TEMP_DIR = (Path(__file__).parent / "temp").resolve()
TEMP_DIR.mkdir(parents=True, exist_ok=True)


class PocketFinder:
    """
    Wrapper for running P2Rank to predict ligand-binding pockets in a protein.

    This class manages execution of the P2Rank tool on a given receptor
    structure (PDB/PDBQT format), parses the output and save results.
    """

    def __init__(self, receptor: Union[str, Path, "Protein"]):
        """
        Initialize the PocketFinder with a receptor structure.

        Args:
            receptor (Union[str, Path, Protein]):
                Either:
                - Path to a protein file in `.pdb` or `.pdbqt` format, or
                - A :class:`Protein` object prepared with ``Protein.prepare()``.

        Raises:
            FileNotFoundError: If the receptor file does not exist.
            ValueError: If the file is not a supported format.
            RuntimeError: If a Protein object is provided without being prepared.
        """
        if isinstance(receptor, Protein):
            if receptor.pdbqt_path is None or not Path(receptor.pdbqt_path).exists():
                raise RuntimeError("❌ Protein not prepared. Run Protein.prepare() first.")
            self.receptor = Path(receptor.pdbqt_path)
        else:
            self.receptor = Path(receptor).resolve()

        if not self.receptor.is_file():
            raise FileNotFoundError(
                f"❌ PDB file not found: {self.receptor}")

        if self.receptor.suffix.lower() not in [".pdb", ".pdbqt"]:
            raise ValueError(
                "❌ Unsupported file format. Only '.pdb' and 'pdbqt' is supported.")

        # Temp directories
        # Use receptor + uuid for uniqueness
        self.temp_dir = TEMP_DIR / "p2rank_results" / f"{self.receptor.stem}_pockets_{uuid.uuid4().hex[:8]}"
        self.temp_dir.mkdir(parents=True, exist_ok=True)


    def run(self) -> List[Dict[str, Union[int, Tuple[float, float, float]]]]:
        """
        Execute P2Rank to predict ligand-binding pockets.

        Returns:
            List[Dict[str, Union[int, Tuple[float, float, float]]]]:
                A list of pocket predictions, where each dictionary contains:
                - ``rank`` (int): Pocket ranking by confidence.
                - ``center`` (tuple[float, float, float]): (x, y, z) coordinates of the pocket center.

        Raises:
            RuntimeError: If P2Rank fails to run.
        """
        cmd = [
            str(P2RANK_PATH),
            "predict",
            "-f", str(self.receptor),
            "-o", str(self.temp_dir),
            "-threads", str(os.cpu_count() or 1)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"❌ Error running P2Rank:\n{result.stderr}")

        return self._parse_output()

    def _parse_output(self) -> List[Dict[str, Union[int, Tuple[float, float, float]]]]:
        """
        Parse P2Rank CSV output and extract predicted pocket centers.

        Returns:
            List[Dict[str, Union[int, Tuple[float, float, float]]]]:
                Pocket predictions with rank and coordinates.

        Raises:
            FileNotFoundError: If the P2Rank CSV output is missing.
            ValueError: If parsing fails or no pockets are found.
        """
        csv_path = self.temp_dir / \
            f"{self.receptor.name}_predictions.csv"
        if not csv_path.is_file():
            raise FileNotFoundError(f"❌ Prediction CSV not found: {csv_path}")

        pockets = []
        try:
            with open(csv_path, 'r') as f:
                reader = csv.DictReader(f)
                for idx, row in enumerate(reader, start=1):
                    try:
                        x = float(row.get("center_x", row.get(
                            "   center_x", "0")).strip())
                        y = float(row.get("center_y", row.get(
                            "   center_y", "0")).strip())
                        z = float(row.get("center_z", row.get(
                            "   center_z", "0")).strip())
                        pockets.append({"rank": idx, "center": (x, y, z)})
                    except Exception as e:
                        raise ValueError(
                            f"❌ Error parsing coordinates at row {idx}: {e}")
        except Exception as e:
            raise ValueError(f"❌ Failed to read prediction CSV: {e}")

        if not pockets:
            raise ValueError(f"❌ No pocket centers found in: {csv_path}")

        return pockets

    def save_report(self, save_to: Union[str, Path] = Path("./p2rank_results")) -> None:
        """
        Save P2Rank results to a user-specified directory.

        Args:
            save_to (Union[str, Path], optional):
                Destination folder to copy results into.
                Defaults to ``./p2rank_results``.

        Returns:
            Path: Path to the saved results directory.

        Raises:
            RuntimeError: If no results are available (``run()`` not called).
        """

        # Check if results exist before proceeding
        if not any(self.temp_dir.glob("*.csv")):
            raise RuntimeError(
                "❌ P2Rank results are missing. Please run PocketFinder before saving results."
            )
        save_to = Path(save_to).resolve()

        shutil.copytree(self.temp_dir, save_to, dirs_exist_ok=True)

        return save_to