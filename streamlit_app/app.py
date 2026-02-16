"""
DockSuiteX — Streamlit App (Multi-Page)
Entry point: loads CSS, sets page config, and defines page navigation.
"""
import sys
from pathlib import Path

# Ensure docksuitex is importable from the project root (one folder above streamlit_app/)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import streamlit as st
import streamlit.components.v1 as components

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

# ── Theme toggle (auto-detects browser preference on first load) ──
if "dark_mode" not in st.session_state:
    # Check if browser preference was communicated via query param
    browser_pref = st.query_params.get("_theme")
    if browser_pref in ("dark", "light"):
        st.session_state["dark_mode"] = browser_pref == "dark"
        st.session_state["_theme_initialised"] = True
        # Clean up the query param so it doesn't linger in the URL
        st.query_params.pop("_theme", None)
    else:
        # Default to dark; JS below will detect & redirect once
        st.session_state["dark_mode"] = True
        st.session_state["_theme_initialised"] = False

with st.sidebar:
    dark_mode = st.toggle("🌙 Dark Mode", value=st.session_state["dark_mode"], key="dark_mode_toggle")
    st.session_state["dark_mode"] = dark_mode

theme = "dark" if st.session_state["dark_mode"] else "light"

# Inject JS: on first load detect browser preference & redirect once,
# afterwards just apply the chosen theme attribute.
if not st.session_state.get("_theme_initialised"):
    # First visit: detect prefers-color-scheme and reload with ?_theme=…
    components.html(
        """
        <script>
            (function() {
                const url = new URL(window.parent.location.href);
                if (!url.searchParams.has('_theme')) {
                    const prefersDark = window.parent.matchMedia('(prefers-color-scheme: dark)').matches;
                    url.searchParams.set('_theme', prefersDark ? 'dark' : 'light');
                    window.parent.location.replace(url.toString());
                }
            })();
        </script>
        """,
        height=0,
    )
    st.stop()  # halt rendering until redirect completes
else:
    # Normal run: apply the selected theme
    components.html(
        f"""
        <script>
            const root = window.parent.document.querySelector('.stApp');
            if (root) {{
                root.setAttribute('data-theme', '{theme}');
            }}
        </script>
        """,
        height=0,
    )

# ── Page navigation ───────────────────────────────────────────
pages_dir = Path(__file__).parent / "pages"

home_page   = st.Page(str(pages_dir / "home.py"),             title="Home",             icon="🏠")
single_page = st.Page(str(pages_dir / "single_docking.py"),   title="Single Docking",   icon="🔬")
multi_page  = st.Page(str(pages_dir / "batch_docking.py"),    title="Batch Docking",    icon="⚡")
fetch_page  = st.Page(str(pages_dir / "fetch_molecules.py"),  title="Fetch Molecules",  icon="🌐")
about_page  = st.Page(str(pages_dir / "about.py"),            title="About",            icon="📌")

nav = st.navigation([home_page, single_page, multi_page, fetch_page, about_page])
nav.run()
