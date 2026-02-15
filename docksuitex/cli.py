import streamlit.web.cli as stcli
import sys
from pathlib import Path
from docksuitex import download_binaries


def main():
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
