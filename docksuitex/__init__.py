"""
DockSuiteX: All-in-one Protein-Ligand Docking Package
Integrates MGLTools, P2Rank, and AutoDock Vina
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
from pathlib import Path


__all__ = [
    "Protein",
    "Ligand",
    "VinaDocking",
    "AD4Docking",
    "PocketFinder",
]



# Download Binaries with Download Progress and Extraction Notice
import platform
import requests
import zipfile
import io
from pathlib import Path
from tqdm import tqdm  # pip install tqdm

# Configuration
GITHUB_ZIP = "https://github.com/MangalamGSinha/DockSuiteX_Binaries/archive/refs/heads/main.zip"
BIN_DIR = Path(__file__).parent / "bin"  # target folder for all binaries
NEEDED_FOLDERS = ["mgltools", "p2rank_2.5.1", "OpenBabel-3.1.1", "AutoDock", "Vina"]

def download_binaries():
    if platform.system() != "Windows":
        raise RuntimeError("❌ DockSuiteX only supports Windows!")

    # Skip if already downloaded
    if all((BIN_DIR / f).exists() for f in NEEDED_FOLDERS):
        return

    print("⬇️ Downloading DockSuiteX_Binaries ...")
    resp = requests.get(GITHUB_ZIP, stream=True)
    resp.raise_for_status()

    total_size = int(resp.headers.get('content-length', 0))
    zip_data = io.BytesIO()
    
    # Download with progress bar
    with tqdm(total=total_size, unit='B', unit_scale=True, desc='Downloading') as pbar:
        for chunk in resp.iter_content(chunk_size=1024*1024):
            if chunk:
                zip_data.write(chunk)
                pbar.update(len(chunk))
    zip_data.seek(0)

    # Extraction
    print("📂 Extracting binaries ...")
    with zipfile.ZipFile(zip_data) as zf:
        root = zf.namelist()[0].split("/")[0]  # e.g., "DockSuiteX_Binaries-main/"
        for folder in NEEDED_FOLDERS:
            folder_prefix = f"{root}{folder}/"
            for member in zf.namelist():
                if member.startswith(folder_prefix):
                    target_path = BIN_DIR / member[len(root):]
                    if member.endswith("/"):
                        target_path.mkdir(parents=True, exist_ok=True)
                    else:
                        target_path.parent.mkdir(parents=True, exist_ok=True)
                        with zf.open(member) as src, open(target_path, "wb") as dst:
                            dst.write(src.read())

    print(f"✅ All binaries saved in {BIN_DIR}")

download_binaries()

