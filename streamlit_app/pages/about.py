"""
DockSuiteX — About
Citations and Contributing information.
"""
import streamlit as st
from pathlib import Path


# ── Helper to read markdown ────────────────────────────────────
def read_md(filename):
    try:
        root = Path(__file__).resolve().parents[2]
        path = root / "docs" / "about" / filename
        if path.exists():
            return path.read_text(encoding="utf-8")
        else:
            return f"⚠️ Could not find `{filename}`."
    except Exception as e:
        return f"⚠️ Error reading `{filename}`: {e}"


# ── Header ─────────────────────────────────────────────────────
st.markdown(
    """
    <div class="hero-header animate-in">
        <h1>📌 About DockSuiteX</h1>
        <p>Citations, references, and how to contribute.</p>
    </div>
    <div class="styled-divider"></div>
    """,
    unsafe_allow_html=True,
)

# ── Tabs for Citations and Contributing ────────────────────────
tab_cite, tab_contrib = st.tabs(["📜 Citations", "🤝 Contributing"])

with tab_cite:
    st.markdown(read_md("cite.md"))

with tab_contrib:
    st.markdown(read_md("contributing.md"))

# ── Footer ─────────────────────────────────────────────────────
st.markdown("---")
st.caption("DockSuiteX © 2025")
