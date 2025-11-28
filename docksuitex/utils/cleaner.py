import shutil
from pathlib import Path

BIN_DIR = (Path(__file__).parent.parent / "bin").resolve()

def delete_binaries():
    """Delete the DockSuiteX binaries directory and all its contents.

    This function removes the entire `bin/` directory containing all bundled
    executables (MGLTools, AutoDock Vina, AutoDock4, P2Rank, Open Babel).
    Use this to clean up disk space or force a fresh download of binaries.

    The binaries will be automatically re-downloaded on the next import of
    DockSuiteX.

    Raises:
        PermissionError: If the directory or files are in use or protected.

    Warning:
        This will delete all binary executables. They will need to be
        re-downloaded (~500MB) on next import.

    Example:
        ::
            from docksuitex.utils import delete_binaries

            # Remove all binaries
            delete_binaries()

    """
    if BIN_DIR.exists() and BIN_DIR.is_dir():
        shutil.rmtree(BIN_DIR)
        print(f"🗑️ Deleted {BIN_DIR}")
    else:
        print(f"⚠️ {BIN_DIR} does not exist.")
