import shutil
from pathlib import Path

# Go one level above the current script directory to access temp folder
TEMP_DIR = (Path(__file__).parent.parent / "temp").resolve()


# def clean_temp_folder():
#     """
#     Deletes and recreates the temporary folder used for intermediate files.

#     This function ensures the TEMP_DIR is clean by deleting all its contents
#     and recreating the directory. Useful for resetting the state before a new run.
#     """
#     if TEMP_DIR.exists() and TEMP_DIR.is_dir():
#         shutil.rmtree(TEMP_DIR)
#         TEMP_DIR.mkdir(exist_ok=True)
#         print(f"✅ Temp folder cleaned: {TEMP_DIR}")
#     else:
#         print("⚠️ Temp folder does not exist.")
#         TEMP_DIR.mkdir(exist_ok=True)



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

            # Binaries will be re-downloaded on next import
            import docksuitex  # Triggers download
    """
    if BIN_DIR.exists() and BIN_DIR.is_dir():
        shutil.rmtree(BIN_DIR)
        print(f"🗑️ Deleted {BIN_DIR}")
    else:
        print(f"⚠️ {BIN_DIR} does not exist.")
