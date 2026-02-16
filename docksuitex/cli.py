import streamlit.web.cli as stcli
import sys
import shutil
from pathlib import Path


def get_bin_dir():
    """Get the binaries directory path."""
    try:
        import docksuitex
        return Path(docksuitex.__file__).parent / "bin"
    except ImportError:
        return None


def clean_binaries():
    """Remove downloaded binaries."""
    bin_dir = get_bin_dir()
    
    if not bin_dir:
        print("❌ DockSuiteX is not installed")
        return 1
    
    if bin_dir.exists():
        try:
            shutil.rmtree(bin_dir)
            print(f"✅ Binaries removed from {bin_dir}")
            print("\nNow run: pip uninstall docksuitex")
            return 0
        except Exception as e:
            print(f"❌ Failed to remove binaries: {e}")
            return 1
    else:
        print("ℹ️  No binaries found")
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
    if len(sys.argv) > 1 and sys.argv[1] == "--clean":
        sys.exit(clean_binaries())
    else:
        launch_gui()


if __name__ == "__main__":
    main()