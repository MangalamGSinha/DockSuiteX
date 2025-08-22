import os
import subprocess
import csv
import shutil
from pathlib import Path
from typing import Optional, Tuple, List, Union, Dict
import uuid

# Paths
P2RANK_PATH = (Path(__file__).parent / "bin" /
               "p2rank_2.5.1" / "p2rank_2.5.1" / "prank.bat").resolve()
TEMP_DIR = (Path(__file__).parent / "temp").resolve()
TEMP_DIR.mkdir(parents=True, exist_ok=True)


class PocketFinder:
    """
    A wrapper class for running P2Rank on a given protein structure (PDB format) to predict 
    potential ligand binding pockets. The class handles running the tool, parsing output, 
    and saving results.
    """

    def __init__(self, protein_path: Union[str, Path], threads: int = os.cpu_count() or 1):
        """
        Initialize the PocketFinder with a protein file path.

        Args:
            protein_path (Union[str, Path]): Path to the input protein file in PDB format.
            threads (int, optional): Number of threads to use for P2Rank. Defaults to number of CPU cores.

        Raises:
            FileNotFoundError: If the specified protein file does not exist.
            ValueError: If the input file is not a .pdb file.
        """
        self.protein_path = Path(protein_path).resolve()

        if not self.protein_path.is_file():
            raise FileNotFoundError(
                f"❌ PDB file not found: {self.protein_path}")

        if self.protein_path.suffix.lower() != ".pdb":
            raise ValueError(
                "❌ Unsupported file format. Only '.pdb' is supported.")

        # Temp directories
        # Use receptor + uuid for uniqueness
        self.temp_dir = TEMP_DIR / "p2rank_results" / f"{self.protein_path.stem}_pockets_{uuid.uuid4().hex[:8]}"
        self.temp_dir.mkdir(parents=True, exist_ok=True)

        self.threads = threads

    def run(self) -> List[Dict[str, Union[int, Tuple[float, float, float]]]]:
        """
        Run P2Rank to predict ligand-binding pockets in the protein structure.

        Returns:
            List[Dict[str, Union[int, Tuple[float, float, float]]]]: 
                A list of dictionaries, each containing:
                - 'rank': pocket ranking based on prediction confidence
                - 'center': a tuple of (x, y, z) coordinates of the pocket center

        Raises:
            RuntimeError: If P2Rank execution fails.
        """
        cmd = [
            str(P2RANK_PATH),
            "predict",
            "-f", str(self.protein_path),
            "-o", str(self.temp_dir),
            "-threads", str(self.threads)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"❌ Error running P2Rank:\n{result.stderr}")

        return self._parse_output()

    def _parse_output(self) -> List[Dict[str, Union[int, Tuple[float, float, float]]]]:
        """
        Parse the P2Rank-generated CSV output file and extract the pocket centers and ranks.

        Returns:
            List[Dict[str, Union[int, Tuple[float, float, float]]]]: 
                A list of predicted pockets with rank and center coordinates.

        Raises:
            FileNotFoundError: If the CSV output file is not found.
            ValueError: If coordinate parsing fails or no pockets are found.
        """
        csv_path = self.temp_dir / \
            f"{self.protein_path.name}_predictions.csv"
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

    def save_report(self, save_dir: Union[str, Path] = Path("./p2rank_results")) -> None:
        """
        Save the entire P2Rank output directory for the current protein to a user-specified location.

        Args:
            save_dir (Union[str, Path], optional): Directory to save the output folder.
                Defaults to "./p2rank_results/".

        Raises:
            RuntimeError: If the output directory has not been generated (i.e., `run()` was not called).
        """

        # Check if results exist before proceeding
        if not any(self.temp_dir.glob("*.csv")):
            raise RuntimeError(
                "❌ P2Rank results are missing. Please run PocketFinder before saving results."
            )
        save_dir = Path(save_dir).resolve()

        shutil.copytree(self.temp_dir, save_dir, dirs_exist_ok=True)

        return save_dir
