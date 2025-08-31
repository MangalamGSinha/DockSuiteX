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
from .batch_vina import BatchVinaDocking
from .autodock4 import AD4Docking
from .batch_autodock4 import BatchAD4Docking
from .pocket_finder import PocketFinder
from pathlib import Path


__all__ = [
    "Protein",
    "Ligand",
    "VinaDocking",
    "BatchVinaDocking",
    "AD4Docking",
    "BatchAD4Docking",
    "PocketFinder",
]







import platform
import requests
import zipfile
import io
from pathlib import Path
from tqdm import tqdm
import shutil

# GitHub repo zip
GITHUB_ZIP = "https://github.com/MangalamGSinha/DockSuiteX_Binaries/archive/refs/heads/main.zip"
BIN_DIR = Path(__file__).parent / "bin"


def download_binaries():
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