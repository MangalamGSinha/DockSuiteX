"""
DockSuiteX — Home / Landing Page
Premium landing page for the DockSuiteX Streamlit GUI.
"""
import streamlit as st

# ── Hero Section ───────────────────────────────────────────────
st.markdown(
    """
    <div class="hero-header animate-in" style="text-align:center; padding:3.5rem 2rem;">
        <div style="font-size:4rem; margin-bottom:0.5rem;">🧬</div>
        <h1 style="font-size:3.2rem; margin-bottom:0.3rem;">DockSuiteX</h1>
        <p class="hero-subtitle">All-in-one Protein-Ligand Docking Suite</p>
        <p style="max-width:640px; margin:1rem auto 0; font-size:1rem; opacity:0.8;">
            A unified graphical interface for protein &amp; ligand preparation,
            binding pocket prediction, and molecular docking simulations — powered
            by AutoDock Vina, AutoDock4, P2Rank, and more.
        </p>
    </div>
    """,
    unsafe_allow_html=True,
)

# ── Key Features ───────────────────────────────────────────────
st.markdown(
    '<div class="section-header"><h2>✨ Key Features</h2></div>',
    unsafe_allow_html=True,
)

FEATURES = [
    ("🔧", "Automated Preparation",
     "Streamlined protein and ligand preparation using PDBFixer, OBabel, and AutoDockTools."),
    ("🎯", "Pocket Prediction",
     "Predict binding pockets using P2Rank machine learning and select the best docking site."),
    ("⚡", "Docking Engines",
     "Run docking with AutoDock Vina or AutoDock4 with full parameter control from the UI."),
    ("📊", "Batch Processing",
     "Screen many ligands against many receptors simultaneously with parallel execution."),
    ("👁️", "Interactive 3D Viewer",
     "Explore proteins, ligands, and docked poses in an interactive NGL.js 3D molecular viewer."),
    ("🌐", "Molecule Fetcher",
     "Fetch protein (PDB) and ligand (SDF) structures directly from RCSB PDB and PubChem."),
]

# Render as a 3-column grid
rows = [FEATURES[i : i + 3] for i in range(0, len(FEATURES), 3)]
for row in rows:
    cols = st.columns(len(row), gap="medium")
    for col, (icon, title, desc) in zip(cols, row):
        with col:
            st.markdown(
                f"""
                <div class="glass-card" style="height:100%;">
                    <span class="card-icon">{icon}</span>
                    <h3>{title}</h3>
                    <p>{desc}</p>
                </div>
                """,
                unsafe_allow_html=True,
            )

st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)

# ── Docking Scenarios ──────────────────────────────────────────
st.markdown(
    '<div class="section-header"><h2>🚀 Flexible Docking Scenarios</h2></div>',
    unsafe_allow_html=True,
)
st.caption("DockSuiteX supports every possible docking combination for high-throughput screening.")

SCENARIOS = [
    ("1 → 1", "Single Docking",
     "One ligand against one receptor for detailed interaction analysis."),
    ("N → 1", "Virtual Screening",
     "Screen compound libraries against a single protein target."),
    ("1 → N", "Inverse Docking",
     "Test one compound against multiple protein variants or targets."),
    ("N → N", "Combinatorial",
     "All ligand-receptor combinations in a single parallel run."),
]

cols = st.columns(4, gap="medium")
for col, (ratio, title, desc) in zip(cols, SCENARIOS):
    with col:
        st.markdown(
            f"""
            <div class="scenario-card">
                <div class="scenario-ratio">{ratio}</div>
                <div class="scenario-title">{title}</div>
                <div class="scenario-desc">{desc}</div>
            </div>
            """,
            unsafe_allow_html=True,
        )

st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)

# ── How It Works ───────────────────────────────────────────────
st.markdown(
    '<div class="section-header"><h2>🛠️ How It Works</h2></div>',
    unsafe_allow_html=True,
)

STEPS = [
    ("1", "Prepare", "Upload your protein (PDB) and ligand files. The app handles format conversion, hydrogen addition, and charge assignment automatically."),
    ("2", "Predict", "Use P2Rank to predict binding pockets on your protein and select the best docking site."),
    ("3", "Dock", "Choose AutoDock Vina or AutoDock4, configure parameters, and launch the docking simulation."),
    ("4", "Analyze", "View ranked poses in a sortable table, inspect 3D docked conformations, and download results."),
]

cols = st.columns(4, gap="medium")
for col, (num, title, desc) in zip(cols, STEPS):
    with col:
        st.markdown(
            f"""
            <div class="glass-card" style="text-align:center; height:100%;">
                <div class="step-number">{num}</div>
                <h3>{title}</h3>
                <p>{desc}</p>
            </div>
            """,
            unsafe_allow_html=True,
        )

st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)

# ── Powered By ─────────────────────────────────────────────────
st.markdown(
    '<div class="section-header"><h2>🙏 Powered By</h2></div>',
    unsafe_allow_html=True,
)

TECH = [
    "AutoDock Vina", "AutoDock4", "MGLTools",
    "Open Babel", "PDBFixer", "P2Rank",
    "NGLView", "RCSB PDB", "PubChem",
]

badge_html = " ".join(
    f'<span class="tech-badge">{t}</span>' for t in TECH
)
st.markdown(
    f'<div style="display:flex; flex-wrap:wrap; gap:10px; justify-content:center; padding:1rem 0;">{badge_html}</div>',
    unsafe_allow_html=True,
)

st.markdown("")

# ── Footer ─────────────────────────────────────────────────────
st.caption("DockSuiteX © 2026")
