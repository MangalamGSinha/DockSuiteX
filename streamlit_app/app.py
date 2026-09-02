"""
DockSuiteX — Streamlit App (Multi-Page)
Entry point: loads CSS, sets page config, and defines page navigation.
"""
import warnings

import os
os.environ["PYTHONWARNINGS"] = "ignore::DeprecationWarning"

# Suppress known third-party deprecation warnings before anything else is imported
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="prolif")

import sys
from pathlib import Path

# Ensure docksuitex is importable from the project root (one folder above streamlit_app/)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import streamlit as st


# ── Page config ────────────────────────────────────────────────
st.set_page_config(
    page_title="DockSuiteX",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Load CSS ───────────────────────────────────────────────────
css_path = Path(__file__).parent / "style.css"
if css_path.exists():
    st.markdown(f"<style>{css_path.read_text()}</style>", unsafe_allow_html=True)

# ── Page navigation ───────────────────────────────────────────
pages_dir = Path(__file__).parent / "pages"

home_page   = st.Page(str(pages_dir / "home.py"),             title="Home",             icon="🏠")
single_page = st.Page(str(pages_dir / "single_docking.py"),   title="Single Docking",   icon="🔬")
multi_page  = st.Page(str(pages_dir / "batch_docking.py"),    title="Batch Docking",    icon="⚡")
fetch_page  = st.Page(str(pages_dir / "fetch_molecules.py"),  title="Fetch Molecules",  icon="🌐")
conv_page   = st.Page(str(pages_dir / "format_converter.py"), title="Format Converter", icon="🧩")
molsynth_pg = st.Page(str(pages_dir / "molsynthai.py"),       title="MolSynthAI",       icon="🧪")
about_page  = st.Page(str(pages_dir / "about.py"),            title="About",            icon="📌")

nav = st.navigation([home_page, single_page, multi_page, fetch_page, conv_page, molsynth_pg, about_page])
nav.run()
