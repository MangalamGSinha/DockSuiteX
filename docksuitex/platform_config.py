"""Platform-specific path resolution for DockSuiteX.

This module centralises all binary executable paths and environment
configuration so that both Windows and Linux are supported from a
single codebase.  Every other module should import paths from here
instead of hard-coding ``.exe`` / ``.bat`` suffixes.

Supported platforms
-------------------
* **Windows** — binaries downloaded from ``DockSuiteX_Binaries_Windows``
* **Linux**   — binaries downloaded from ``DockSuiteX_Binaries_Linux``
"""

import os
import platform
from pathlib import Path

# ── Platform detection ────────────────────────────────────────────────────────

IS_WINDOWS = platform.system() == "Windows"
IS_LINUX = platform.system() == "Linux"

if not (IS_WINDOWS or IS_LINUX):
    raise OSError(
        "DockSuiteX is only supported on Windows and Linux. "
        f"Detected platform: {platform.system()}"
    )

# ── Root binary directory ─────────────────────────────────────────────────────

BIN_DIR = Path(__file__).parent / "bin"

# ── GitHub download URLs ──────────────────────────────────────────────────────

if IS_WINDOWS:
    GITHUB_ZIP = (
        "https://github.com/MangalamGSinha/DockSuiteX_Binaries_Windows"
        "/archive/refs/heads/main.zip"
    )
else:
    GITHUB_ZIP = (
        "https://github.com/MangalamGSinha/DockSuiteX_Binaries_Linux"
        "/archive/refs/heads/main.zip"
    )

# ── AutoDock Vina ─────────────────────────────────────────────────────────────

if IS_WINDOWS:
    VINA_PATH = (BIN_DIR / "vina" / "vina.exe").resolve()
else:
    VINA_PATH = (BIN_DIR / "vina" / "vina").resolve()

# ── AutoDock4 / AutoGrid4 ────────────────────────────────────────────────────

if IS_WINDOWS:
    AUTOGRID_EXE = (BIN_DIR / "autodock" / "autogrid4.exe").resolve()
    AUTODOCK_EXE = (BIN_DIR / "autodock" / "autodock4.exe").resolve()
else:
    AUTOGRID_EXE = (BIN_DIR / "autodock" / "autogrid4").resolve()
    AUTODOCK_EXE = (BIN_DIR / "autodock" / "autodock4").resolve()

# ── MGLTools ──────────────────────────────────────────────────────────────────

MGLTOOLS_PATH = (BIN_DIR / "mgltools").resolve()

if IS_WINDOWS:
    MGL_PYTHON_EXE = (MGLTOOLS_PATH / "python.exe").resolve()
    PREPARE_RECEPTOR_SCRIPT = (
        MGLTOOLS_PATH / "Lib" / "site-packages"
        / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py"
    ).resolve()
    PREPARE_LIGAND_SCRIPT = (
        MGLTOOLS_PATH / "Lib" / "site-packages"
        / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"
    ).resolve()
else:
    MGL_PYTHON_EXE = (MGLTOOLS_PATH / "bin" / "MGLpython2.7").resolve()
    PREPARE_RECEPTOR_SCRIPT = (
        MGLTOOLS_PATH / "MGLToolsPckgs"
        / "AutoDockTools" / "Utilities24" / "prepare_receptor4.py"
    ).resolve()
    PREPARE_LIGAND_SCRIPT = (
        MGLTOOLS_PATH / "MGLToolsPckgs"
        / "AutoDockTools" / "Utilities24" / "prepare_ligand4.py"
    ).resolve()

# ── P2Rank ────────────────────────────────────────────────────────────────────

if IS_WINDOWS:
    P2RANK_PATH = (BIN_DIR / "p2rank" / "prank.bat").resolve()
else:
    P2RANK_PATH = (BIN_DIR / "p2rank" / "prank").resolve()


# ── Environment helpers ──────────────────────────────────────────────────────

def get_mgltools_env() -> dict:
    """Return an environment dict with MGLTools paths configured.

    On Linux the bundled MGLTools Python 2.7 requires ``LD_LIBRARY_PATH``,
    ``PYTHONHOME``, and ``PYTHONPATH`` to be pointed at the extracted
    MGLTools tree.  On Windows the standard ``PATH`` is sufficient.

    Returns:
        dict: A copy of ``os.environ`` with the necessary overrides.
    """
    env = os.environ.copy()
    if IS_LINUX:
        base = str(MGLTOOLS_PATH)
        env["LD_LIBRARY_PATH"] = f"{base}/lib"
        env["PYTHONHOME"] = base
        env["PYTHONPATH"] = f"{base}/MGLToolsPckgs"
    return env
