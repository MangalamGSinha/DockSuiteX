"""DockSuiteX: All-in-one Protein-Ligand Docking Package.

A comprehensive Python package for molecular docking that integrates multiple
tools including MGLTools, P2Rank, AutoDock Vina, and AutoDock4. DockSuiteX
provides a unified interface for protein and ligand preparation, binding pocket
prediction, and molecular docking simulations.

Key Features:
    - Automated protein and ligand preparation using PDBFixer and Open Babel
    - Binding pocket prediction with P2Rank
    - Molecular docking with AutoDock Vina and AutoDock4
    - Batch processing capabilities for high-throughput screening
    - Interactive 3D visualization with NGLView
    - Utility functions for fetching structures and parsing results

Note:
    This package is currently only supported on Windows platforms.

Example:
    Basic usage for single ligand docking::

        from docksuitex import Protein, Ligand, VinaDocking

        # Prepare protein
        protein = Protein("protein.pdb")
        protein.prepare(save_to="prepared_protein.pdbqt")

        # Prepare ligand
        ligand = Ligand("ligand.sdf")
        ligand.prepare(save_to="prepared_ligand.pdbqt")

        # Run docking
        docking = VinaDocking(
            receptor="prepared_protein.pdbqt",
            ligand="prepared_ligand.pdbqt",
            grid_center=(10.0, 15.0, 20.0)
        )
        results = docking.run()
"""

__version__ = "1.0.0"
__author__ = "DockSuiteX Team"

import platform

if platform.system() != "Windows":
    raise OSError(
        "DockSuiteX is only supported on Windows. "
        "Please install on a Windows environment."
    )

from .protein import Protein
from .ligand import Ligand
from .vina import VinaDocking
from .autodock4 import AD4Docking
from .pocket_finder import PocketFinder


__all__ = [
    "Protein",
    "Ligand",
    "VinaDocking",
    "AD4Docking",
    "PocketFinder",
]




import requests
import zipfile
import io
from pathlib import Path
from tqdm import tqdm
import shutil

# GitHub repository URL for binary dependencies
GITHUB_ZIP = "https://github.com/MangalamGSinha/DockSuiteX_Binaries/archive/refs/heads/main.zip"
BIN_DIR = Path(__file__).parent / "bin"


def download_binaries():
    """Download required binary executables from GitHub repository.

    This function automatically downloads and extracts the DockSuiteX binary
    dependencies (MGLTools, AutoDock Vina, AutoDock4, P2Rank, Open Babel)
    from the GitHub repository on first import. If binaries already exist,
    the download is skipped.

    The binaries are extracted to the `bin/` directory within the package
    installation directory.

    Raises:
        RuntimeError: If the operating system is not Windows.
        requests.HTTPError: If the download from GitHub fails.
        zipfile.BadZipFile: If the downloaded file is corrupted.

    Note:
        This function is automatically called when the package is imported.
        The download progress is displayed using tqdm.
    """
    if platform.system() != "Windows":
        raise RuntimeError("❌ DockSuiteX only supports Windows!")

    # If bin already exists, skip
    if BIN_DIR.exists() and any(BIN_DIR.iterdir()):
        # print(f"✅ Binaries already exist in {BIN_DIR}")
        return

    print("⬇️ Downloading DockSuiteX_Binaries ...")
    resp = requests.get(GITHUB_ZIP, stream=True)
    resp.raise_for_status()

    total_size = int(resp.headers.get('content-length', 0))
    zip_data = io.BytesIO()

    with tqdm(total=total_size, unit='B', unit_scale=True, desc='Downloading') as pbar:
        for chunk in resp.iter_content(chunk_size=1024*1024):
            if chunk:
                zip_data.write(chunk)
                pbar.update(len(chunk))
    zip_data.seek(0)

    with zipfile.ZipFile(zip_data) as zf:
        root = zf.namelist()[0].split("/")[0]  # DockSuiteX_Binaries-main
        for member in zf.namelist():
            if member.endswith("/"):  # skip directories for now
                continue
            if member.startswith(root):
                relative_path = member[len(root):].lstrip(
                    "/")  # remove root + leading slash
                target_path = BIN_DIR / relative_path

                target_path.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(member) as src, open(target_path, "wb") as dst:
                    dst.write(src.read())

    print(f"✅ All binaries saved in {BIN_DIR}")


download_binaries()