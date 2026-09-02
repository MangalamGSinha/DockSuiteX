"""Command line interface for DockSuiteX."""

import warnings
import sys
import os
import stat
import shutil
from pathlib import Path

# Forcefully suppress all DeprecationWarnings at the environment level
# This ensures sub-processes (Streamlit) also inherit the suppression.
os.environ["PYTHONWARNINGS"] = "ignore::DeprecationWarning"

# Also apply local filters as a backup
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")

import streamlit.web.cli as stcli


def _remove_readonly(func, path, exc_info=None):
    """Clear the read-only / write-protection bit and retry removal."""
    try:
        os.chmod(path, stat.S_IWRITE | stat.S_IWUSR)
        func(path)
    except Exception:
        pass


def get_bin_dir() -> Path:
    """Get the binaries directory path without triggering package-level imports."""
    return Path(__file__).resolve().parent / "bin"


def clean_binaries() -> int:
    """Remove downloaded binaries."""
    bin_dir = get_bin_dir()
    
    if bin_dir.exists():
        try:
            if sys.version_info >= (3, 12):
                def _onexc(func, path, exc):
                    _remove_readonly(func, path)
                shutil.rmtree(bin_dir, onexc=_onexc)
            else:
                shutil.rmtree(bin_dir, onerror=_remove_readonly)
            print(f"✅ Binaries removed from {bin_dir}")
            print("Now run: pip uninstall docksuitex")
            return 0
        except Exception as e:
            print(f"❌ Failed to remove binaries: {e}")
            return 1
    else:
        print("⚠️ No binaries found")
        print("Now run: pip uninstall docksuitex")
        return 0


def launch_gui():
    """Launch the Streamlit GUI application."""
    from docksuitex import download_binaries
    
    # Ensure binaries exist before launching GUI
    download_binaries()

    base_path = Path(__file__).resolve().parent.parent
    app_path = base_path / "streamlit_app" / "app.py"

    sys.argv = [
        "streamlit",
        "run",
        str(app_path),
    ]

    stcli.main()


def main():
    """Run the DockSuiteX command-line interface."""
    if len(sys.argv) > 1 and sys.argv[1] in ("--clean", "-c"):
        sys.exit(clean_binaries())
    else:
        launch_gui()


if __name__ == "__main__":
    main()