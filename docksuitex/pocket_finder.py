"""Binding pocket prediction using P2Rank.

This module provides a Python wrapper for P2Rank, a machine learning-based tool
for predicting ligand-binding pockets in protein structures. P2Rank uses a
template-free approach based on local chemical neighborhood features.

The pocket finding workflow:
    1. Run P2Rank prediction on protein structure (PDB or PDBQT)
    2. Parse CSV output to extract pocket centers and rankings
    3. Return pocket coordinates for use in docking grid setup

Example:
    Basic pocket finding::

        from docksuitex import PocketFinder

        # Find pockets in protein structure
        finder = PocketFinder("protein.pdb")
        pockets = finder.run(save_to="pocket_results")
        
        # Use top pocket for docking
        top_pocket_center = pockets[0]["center"]
        print(f"Top pocket at: {top_pocket_center}")

    Finding pockets for docking setup::

        finder = PocketFinder("prepared_receptor.pdbqt")
        pockets = finder.run()
        
        # Get centers for all predicted pockets
        centers = [pocket["center"] for pocket in pockets]

"""

import os
import subprocess
import csv
import shutil
from pathlib import Path
from typing import Optional, Tuple, List, Union, Dict
import uuid

# Path to P2Rank executable bundled with DockSuiteX
P2RANK_PATH = (Path(__file__).parent / "bin" / "p2rank" / "prank.bat").resolve()


class PocketFinder:
    """Ligand-binding pocket prediction using P2Rank.

    This class provides a wrapper for P2Rank, a machine learning-based tool
    for predicting ligand-binding sites in protein structures. It automates
    running P2Rank, parsing results, and extracting pocket center coordinates.

    P2Rank uses a template-free approach based on local chemical neighborhood
    features, making it fast and accurate for pocket prediction without
    requiring known ligand structures.

    Attributes:
        file_path (Path): Path to the input protein file (PDB or PDBQT).

    Example:
        Finding pockets for docking::

            from docksuitex import PocketFinder

            finder = PocketFinder("protein.pdb")
            pockets = finder.run(save_to="pocket_analysis")
            
            # Pockets are ranked by confidence
            for i, pocket in enumerate(pockets):
                rank = pocket["rank"]
                center = pocket["center"]
                print(f"Pocket {rank}: center at {center}")

        Using pocket centers for batch docking::

            finder = PocketFinder("receptor.pdbqt")
            pockets = finder.run()
            centers = [p["center"] for p in pockets[:3]]  # Top 3 pockets
            
            # Use centers for docking
            for center in centers:
                docking = VinaDocking(
                    receptor="receptor.pdbqt",
                    ligand="ligand.pdbqt",
                    grid_center=center
                )
                docking.run()

    Note:
        P2Rank generates multiple output files including visualizations and
        detailed predictions. The main CSV file is parsed to extract pocket
        center coordinates.
    """

    def __init__(self, input: Union[str, Path], _cpu: int = os.cpu_count() or 1,):
        """Initialize the PocketFinder with a receptor structure.

        Args:
            input (Union[str, Path]):
                Path to a protein file in `.pdb` or `.pdbqt` format.

        Raises:
            FileNotFoundError: If the receptor file does not exist.
            ValueError: If the file is not a supported format.
        """
        self.file_path = Path(input).resolve()
        self.cpu = _cpu

        if not self.file_path.is_file():
            raise FileNotFoundError(
                f"❌ PDB file not found: {self.file_path}")

        if self.file_path.suffix.lower() not in [".pdb", ".pdbqt"]:
            raise ValueError(
                "❌ Unsupported file format. Only '.pdb' and 'pdbqt' is supported.")

    def run(self, save_to: Union[str, Path] = None) -> List[Dict[str, Union[int, Tuple[float, float, float]]]]:
        """Execute P2Rank to predict ligand-binding pockets.
        
        Args:
            save_to (Union[str, Path], optional): Directory to save the P2Rank report. 
                Defaults to ".".

        Returns:
            List[Dict[str, Union[int, Tuple[float, float, float]]]]:
                A list of pocket predictions, where each dictionary contains:
                - ``rank`` (int): Pocket ranking by confidence.
                - ``center`` (tuple[float, float, float]): (x, y, z) coordinates of the pocket center.

        Raises:
            RuntimeError: If P2Rank fails to run.
        """

        if save_to is None:
            save_to = f"./p2rank_results_{self.file_path.name.replace('.', '_')}"

        output_dir = Path(save_to).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            str(P2RANK_PATH),
            "predict",
            "-f", str(self.file_path),
            "-o", str(output_dir),
            "-threads", str(self.cpu),  
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"❌ Error running P2Rank:\n{result.stderr}")

        pockets = self._parse_output(output_dir)
        print(f"✅ Pocket prediction completed. Found {len(pockets)} pockets.")
        for i, pocket in enumerate(pockets):
            print(f"Pocket {i+1} center: {pocket['center']}")
        return pockets

    def _parse_output(self, output_dir: Path) -> List[Dict[str, Union[int, Tuple[float, float, float]]]]:
        """Parse P2Rank CSV output and extract predicted pocket centers.
        
        Args:
            output_dir (Path): Directory where P2Rank results were saved.

        Returns:
            List[Dict[str, Union[int, Tuple[float, float, float]]]]:
                Pocket predictions with rank and coordinates.

        Raises:
            FileNotFoundError: If the P2Rank CSV output is missing.
            ValueError: If parsing fails or no pockets are found.
        """
        csv_path = output_dir / f"{self.file_path.name}_predictions.csv"
        
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