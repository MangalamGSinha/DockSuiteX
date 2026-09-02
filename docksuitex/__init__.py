"""DockSuiteX core package initialization."""

import platform
import warnings

# Suppress known deprecation warnings from third-party dependencies (MDAnalysis, ProLIF)
warnings.filterwarnings("ignore", message=".*topology.tables.*", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning, module="MDAnalysis")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="prolif")
warnings.filterwarnings("ignore", message=".*has been superseded by a new class.*", category=UserWarning)

try:
    import pandas as pd
    pd.set_option('future.no_silent_downcasting', True)
except (ImportError, AttributeError):
    pass

if platform.system() not in ("Windows", "Linux"):
    raise OSError(
        "DockSuiteX is only supported on Windows and Linux. "
        f"Detected platform: {platform.system()}"
    )

from .protein import Protein
from .ligand import Ligand
from .vina import VinaDocking
from .autodock4 import AD4Docking
from .grid_calculator import GridCalculator
from .interaction_profiler import InteractionProfiler


__all__ = [
    "Protein",
    "Ligand",
    "VinaDocking",
    "AD4Docking",
    "GridCalculator",
    "InteractionProfiler",
]


import requests
import zipfile
import io
from pathlib import Path
from tqdm import tqdm
import shutil
import stat

from .platform_config import GITHUB_ZIP, BIN_DIR, IS_LINUX


def download_binaries():
    """Download required binary executables from GitHub repository.

    This function automatically downloads and extracts the DockSuiteX binary
    dependencies (MGLTools, AutoDock Vina, AutoDock4, P2Rank)
    from the GitHub repository on first import. If binaries already exist,
    the download is skipped.

    On Linux, downloaded binaries are automatically marked as executable.

    The binaries are extracted to the `bin/` directory within the package
    installation directory.

    Raises:
        requests.HTTPError: If the download from GitHub fails.
        zipfile.BadZipFile: If the downloaded file is corrupted.

    Note:
        This function is automatically called when the package is imported.
        The download progress is displayed using tqdm.
    """
    # Ensure bin directory exists
    BIN_DIR.mkdir(exist_ok=True)

    # Check if binaries exist (ignore .keep file)
    existing_files = [f for f in BIN_DIR.iterdir() if f.name != '.keep']
    if existing_files:
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

    # On Linux, mark key binaries as executable
    if IS_LINUX:
        _set_linux_permissions()

    print(f"✅ All binaries saved in {BIN_DIR}")


def _set_linux_permissions():
    """Mark downloaded Linux binaries as executable (chmod +x)."""
    executable_patterns = [
        BIN_DIR / "vina" / "vina",
        BIN_DIR / "autodock" / "autodock4",
        BIN_DIR / "autodock" / "autogrid4",
        BIN_DIR / "p2rank" / "prank",
        BIN_DIR / "java" / "bin" / "java",
        BIN_DIR / "mgltools" / "bin" / "MGLpython2.7",
    ]
    for path in executable_patterns:
        if path.exists():
            path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    # Also make all .sh scripts and binaries inside mgltools executable
    mgltools_bin = BIN_DIR / "mgltools" / "bin"
    if mgltools_bin.is_dir():
        for f in mgltools_bin.iterdir():
            if f.is_file():
                f.chmod(f.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    # Make all .so shared libraries readable
    mgltools_lib = BIN_DIR / "mgltools" / "lib"
    if mgltools_lib.is_dir():
        for f in mgltools_lib.rglob("*"):
            if f.is_file() and f.suffix == ".so":
                f.chmod(f.stat().st_mode | stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)


download_binaries()