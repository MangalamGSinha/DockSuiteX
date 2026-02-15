"""
DockSuiteX — Streamlit App (Multi-Page)
Entry point: loads CSS, sets page config, and defines page navigation.
"""
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
about_page  = st.Page(str(pages_dir / "about.py"),            title="About",            icon="📌")

nav = st.navigation([home_page, single_page, multi_page, fetch_page, about_page])
nav.run()
