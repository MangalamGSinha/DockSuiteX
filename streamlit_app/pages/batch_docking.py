"""
DockSuiteX — Batch Docking Page
Batch protein + ligand workflow — no visualization.
"""
import streamlit as st
import io, os, traceback, contextlib, shutil, tempfile, zipfile
from pathlib import Path
import pandas as pd


# ── Helpers ────────────────────────────────────────────────────

class StreamlitSink:
    """Redirects writes to a Streamlit placeholder for real-time logging."""
    def __init__(self, key, placeholder):
        self.key = key
        self.placeholder = placeholder
        self.buffer = io.StringIO()

    def write(self, text):
        self.buffer.write(text)
        if self.placeholder:
            self.placeholder.code(self.buffer.getvalue(), language="text")

    def flush(self):
        self.buffer.flush()


@contextlib.contextmanager
def capture_log(key, placeholder=None):
    """Capture stdout to session state key and optionally stream to UI."""
    if placeholder:
        sink = StreamlitSink(key, placeholder)
        with contextlib.redirect_stdout(sink):
            yield
        st.session_state[key] = sink.buffer.getvalue()
    else:
        f = io.StringIO()
        with contextlib.redirect_stdout(f):
            yield
        st.session_state[key] = f.getvalue()


def parse_comma_list(raw: str):
    return [x.strip() for x in raw.split(",") if x.strip()] if raw else None


def copy_inputs_to_dir(items, input_dir: Path):
    """Copy uploaded/fetched input items into a target input directory."""
    for item in items:
        dst = input_dir / item["name"]
        if item["type"] == "fetch":
            shutil.copy(item["path"], dst)
        else:
            dst.write_bytes(item["data"])


def render_log_if_present(session_key: str, title: str = "Logs:"):
    if st.session_state.get(session_key):
        st.caption(title)
        st.code(st.session_state[session_key], language="text")


def render_input_section(title, key_prefix, storage_key, allowed_exts, uploader_label, help_text=None):
    """Renders a file upload section WITHOUT visualization."""
    st.markdown(f'<div class="section-header"><h2>{title}</h2></div>', unsafe_allow_html=True)
    uploaded = st.file_uploader(
        uploader_label,
        type=allowed_exts,
        accept_multiple_files=True,
        key=f"{key_prefix}_uploader",
        help=help_text,
    )
    if uploaded:
        st.session_state[storage_key] = [
            {"name": f.name, "type": "upload", "data": f.getvalue(), "ext": Path(f.name).suffix}
            for f in uploaded
        ]
    else:
        st.session_state[storage_key] = []


# ── Session state init ─────────────────────────────────────────
def init_session_state():
    try:
        for key in [
            "prepared_ligands", "prepared_proteins",
            "pocket_results", "docking_results",
            "ligand_inputs", "protein_inputs",
        ]:
            if key not in st.session_state:
                if "inputs" in key:
                    st.session_state[key] = []
                else:
                    st.session_state[key] = ""

        for key in ["log_ligand", "log_protein", "log_docking", "log_pockets"]:
            if key not in st.session_state:
                st.session_state[key] = ""

        if "dock_engine_used" not in st.session_state:
            st.session_state["dock_engine_used"] = None
        if "interaction_profile" not in st.session_state:
            st.session_state["interaction_profile"] = None
        if "pocket_mode" not in st.session_state:
            st.session_state["pocket_mode"] = "Grid Calculator"
    except Exception:
        pass


def _fmt_pocket(pocket: dict) -> str:
    """Format a pocket dict as (cx, cy, cz):(sx, sy, sz)."""
    c = pocket["center"]
    s = pocket["grid_size"]
    return f"({c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f}):({s[0]:.3f}, {s[1]:.3f}, {s[2]:.3f})"


def activate_manual_mode():
    """Callback: editing a pocket text area switches to Manual mode."""
    st.session_state["pocket_mode"] = "Manual Pockets"


def restore_gc_pockets():
    """Callback: switching radio to Grid Calculator restores predicted pockets."""
    if st.session_state["pocket_mode"] == "Grid Calculator":
        proteins = st.session_state.get("prepared_proteins", [])
        if proteins:
            for p in proteins:
                gc_key = f"gc_{p.stem}"
                mc_key = f"mc_{p.stem}"
                if gc_key in st.session_state:
                    st.session_state[mc_key] = st.session_state[gc_key]


init_session_state()

# ── Hero header ────────────────────────────────────────────────
st.markdown(
    """
    <div class="hero-header animate-in">
        <h1>⚡ Batch Docking</h1>
        <p>Batch molecular docking workflow &mdash; prepare multiple proteins &amp; 
        ligands, find pockets, dock, and profile all combinations.</p>
        <div style="margin-top: 15px; font-size: 0.9rem; opacity: 0.85; border-top: 1px solid rgba(255,255,255,0.1); padding-top: 10px;">
            <b> 1:N</b> (Virtual Screening) &nbsp; | &nbsp; 
            <b> N:1</b> (Reverse Docking) &nbsp; | &nbsp; 
            <b> N:N</b> (Cross Docking)
        </div>
    </div>
    """,
    unsafe_allow_html=True,
)

# ── Temp output directory ──────────────────────────────────────
if "m_temp_dir_obj" not in st.session_state:
    st.session_state["m_temp_dir_obj"] = tempfile.TemporaryDirectory(prefix="docksuitex_batch_")
    st.session_state["m_temp_dir"] = st.session_state["m_temp_dir_obj"].name
out_path = Path(st.session_state["m_temp_dir"])

# ── CPU Settings ─────────────────────────────────────────────
cpu_avail = os.cpu_count() or 2
global_max_cpu = st.slider("Max CPUs", 1, cpu_avail, max(cpu_avail - 1, 1), key="m_max_cpu")
st.caption("If you encounter an **Out of Memory (OOM)** error, try reducing the number of CPUs.")


# ── Files panel styling ────────────────────────────────────────
st.markdown("""
<style>
section[data-testid="stSidebar"] > div:first-child {
    overflow-y: auto;
    max-height: 100vh;
}
/* Light individual file download buttons (secondary only — ZIP is primary) */
section[data-testid="stSidebar"] div[data-testid="stDownloadButton"] button[kind="secondary"] {
    background: rgba(128, 128, 128, 0.08) !important;
    border: 1px solid rgba(128, 128, 128, 0.15) !important;
    box-shadow: none !important;
    color: inherit !important;
}
section[data-testid="stSidebar"] div[data-testid="stDownloadButton"] button[kind="secondary"]:hover {
    background: rgba(128, 128, 128, 0.18) !important;
}
</style>
""", unsafe_allow_html=True)


def _fmt_size(size_bytes: int) -> str:
    """Format file size to human-readable string."""
    if size_bytes < 1024:
        return f"{size_bytes} B"
    elif size_bytes < 1048576:
        return f"{size_bytes / 1024:.1f} KB"
    return f"{size_bytes / 1048576:.1f} MB"


def _render_dir_contents(current: Path, root: Path):
    """Recursively render directory contents as expanders and download buttons."""
    try:
        items = sorted(current.iterdir(), key=lambda x: (x.is_file(), x.name.lower()))
    except OSError:
        return
    for item in items:
        if item.is_dir():
            with st.expander(f"📁 {item.name}"):
                _render_dir_contents(item, root)
        else:
            sz = _fmt_size(item.stat().st_size)
            rel_key = str(item.relative_to(root)).replace(os.sep, "/")
            st.download_button(
                f"📄 {item.name}  ({sz})",
                data=item.read_bytes(),
                file_name=item.name,
                key=f"m_dl_{rel_key}",
                width='stretch',
            )


def render_files_panel(root_dir: Path):
    """Render a files panel in the sidebar showing actual folder structure."""
    with st.sidebar:
        st.markdown("### 📁 Generated Files")

        if not root_dir.exists():
            st.caption("No files generated yet.")
            return

        files_only = [f for f in root_dir.rglob("*") if f.is_file()]

        if not files_only:
            st.caption("No files generated yet.")
            return

        # Download All as ZIP
        zip_buf = io.BytesIO()
        with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
            for fp in files_only:
                zf.write(fp, fp.relative_to(root_dir))
        zip_buf.seek(0)

        st.download_button(
            "📦 Download All (ZIP)",
            data=zip_buf.getvalue(),
            file_name="docksuitex_results.zip",
            mime="application/zip",
            width='stretch',
            key="m_zip_dl_all",
            type="primary",
        )

        st.divider()

        # Render actual folder tree
        _render_dir_contents(root_dir, root_dir)


st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)

# ════════════════════════════════════════════════════════════════
# ROW 1: PROTEIN PREP (left)  |  LIGAND PREP (right)
# ════════════════════════════════════════════════════════════════
col_prot, col_lig = st.columns(2, gap="large")

# ── Protein Preparation ───────────────────────────────────────
with col_prot:
    render_input_section(
        "🧬 Protein Preparation", 
        "prot", 
        "protein_inputs", 
        ["pdb"],
        "Upload proteins",
        help_text="Please upload files containing only a single model. If you have multi-model files, use the Format Converter to split them first."
    )
    all_proteins = st.session_state["protein_inputs"]

    with st.expander("⚙️ Protein Options", expanded=True):
        p1, p2 = st.columns(2)
        with p1:
            prot_fix   = st.checkbox("Fix PDB(PDBFixer)", True, key="prot_fix")
            prot_het   = st.checkbox("Remove heterogens", True, key="prot_het")
            prot_water = st.checkbox("Remove water", True, key="prot_water")
            prot_chg = st.checkbox("Add Gasteiger charges", True, key="prot_chg")
            prot_h   = st.checkbox("Add hydrogens", True, key="prot_h")
        with p2:
            prot_ph = st.number_input("pH", value=7.4, min_value=0.0, max_value=14.0,
                                       step=0.1, format="%.1f", key="prot_ph",
                                       help="pH for protonation state assignment via PDBFixer")
            prot_preserve = st.text_input("Preserve charge types", placeholder="e.g. Zn, Fe", key="prot_pres")
    prot_preserve_list = parse_comma_list(prot_preserve)

    if st.button("🚀 Prepare Proteins", type="primary", width='stretch', key="btn_prot"):
        if not all_proteins:
            st.warning("Upload protein files first.")
        else:
            log_ph = st.empty()
            with st.spinner(f"Preparing {len(all_proteins)} protein(s)…"):
                with capture_log("log_protein", log_ph):
                    try:
                        prot_input_dir = out_path / "proteins"
                        prot_input_dir.mkdir(parents=True, exist_ok=True)
                        copy_inputs_to_dir(all_proteins, prot_input_dir)
                        from docksuitex.batch_docking.batch_protein import BatchProtein
                        prot_save_dir = out_path / "prepared_proteins"
                        batch = BatchProtein(inputs=prot_input_dir, fix_pdb=prot_fix, remove_heterogens=prot_het,
                                             remove_water=prot_water, add_hydrogens=prot_h, ph=prot_ph,
                                             add_charges=prot_chg,
                                             preserve_charge_types=prot_preserve_list)
                        results = batch.prepare_all(save_to=prot_save_dir, cpu=global_max_cpu)
                        success = [r for r in results if r["status"] == "success"]
                        if success:
                            st.session_state["prepared_proteins"] = [Path(r["pdbqt_path"]) for r in success]
                            st.session_state["protein_out_dir"] = prot_save_dir
                    except Exception as e:
                        st.error(str(e))
                        with st.expander("Traceback"): st.code(traceback.format_exc())
            log_ph.empty()
    render_log_if_present("log_protein", "Protein Preparation Log:")

# ── Ligand Preparation ────────────────────────────────────────
with col_lig:
    render_input_section(
        "💊 Ligand Preparation", 
        "lig", 
        "ligand_inputs", 
        ["sdf", "mol2"],
        "Upload ligands",
        help_text="Please upload files containing only a single molecule. If you have multi-molecule libraries, use the Format Converter to split them first."
    )
    all_ligands = st.session_state["ligand_inputs"]

    with st.expander("⚙️ Ligand Options", expanded=True):
        lig_minimize = st.selectbox("Energy minimization", [None,"mmff94","mmff94s","uff","gaff"],
                                    format_func=lambda x: "None (skip)" if x is None else x.upper(),
                                    help="Select a forcefield for energy minimization (2500 steps) for all ligands.",
                                    key="lig_min")
        l1, l2 = st.columns(2)
        with l1:
            lig_rw = st.checkbox("Remove water", True, key="lig_rw")
            lig_h  = st.checkbox("Add hydrogens", True, key="lig_h")
            lig_chg = st.checkbox("Add Gasteiger charges", True, key="lig_chg")
        with l2:
            lig_ph = st.number_input("pH", value=7.4, min_value=0.0, max_value=14.0,
                                      step=0.1, format="%.1f", key="lig_ph",
                                      help="pH for protonation state assignment via Open Babel")
            lig_preserve = st.text_input("Preserve charge types", placeholder="e.g. Zn, Fe", key="lig_pres")
    lig_preserve_list = parse_comma_list(lig_preserve)

    if st.button("🚀 Prepare Ligands", type="primary", width='stretch', key="btn_lig"):
        if not all_ligands:
            st.warning("Upload ligand files first.")
        else:
            log_ph = st.empty()
            with st.spinner(f"Preparing {len(all_ligands)} ligand(s)…"):
                with capture_log("log_ligand", log_ph):
                    try:
                        lig_input_dir = out_path / "ligands"
                        lig_input_dir.mkdir(parents=True, exist_ok=True)
                        copy_inputs_to_dir(all_ligands, lig_input_dir)
                        from docksuitex.batch_docking.batch_ligand import BatchLigand
                        lig_save_dir = out_path / "prepared_ligands"
                        batch = BatchLigand(inputs=lig_input_dir, minimize=lig_minimize, remove_water=lig_rw,
                                            add_hydrogens=lig_h, ph=lig_ph, add_charges=lig_chg,
                                            preserve_charge_types=lig_preserve_list)
                        results = batch.prepare_all(save_to=lig_save_dir, cpu=global_max_cpu)
                        success = [r for r in results if r["status"] == "success"]
                        if success:
                            st.session_state["prepared_ligands"] = [Path(r["pdbqt_path"]) for r in success]
                            st.session_state["ligand_out_dir"] = lig_save_dir
                    except Exception as e:
                        st.error(str(e))
                        with st.expander("Traceback"): st.code(traceback.format_exc())
            log_ph.empty()
    render_log_if_present("log_ligand", "Ligand Preparation Log:")

# ════════════════════════════════════════════════════════════════
# GRID CALCULATOR (full width, outside columns)
# ════════════════════════════════════════════════════════════════
st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)
st.markdown('<div class="section-header"><h2>📐 Grid Calculator</h2></div>', unsafe_allow_html=True)
st.caption("Compute grid box(es) for each protein — predict pockets with P2Rank or auto-compute blind boxes.")

gc_batch_mode = st.radio(
    "Mode",
    ["P2Rank (pocket prediction)", "Blind Docking (full receptor)"],
    horizontal=True,
    key="gc_batch_mode",
)
is_gc_batch_p2rank = gc_batch_mode == "P2Rank (pocket prediction)"

if is_gc_batch_p2rank:
    pf_max = st.number_input("Max pockets per protein (0 = all)", min_value=0, value=3, key="pf_max")
else:
    pf_max = 0

batch_btn_label = "🔍 Run P2Rank" if is_gc_batch_p2rank else "📦 Compute Blind Boxes"
if st.button(batch_btn_label, type="primary", width='stretch', key="btn_pf"):
    can_find = st.session_state["prepared_proteins"] not in (None, "", [])
    if not can_find:
        st.warning("Prepare proteins above first.")
    else:
        log_ph = st.empty()
        spinner_msg = "Running P2Rank…" if is_gc_batch_p2rank else "Computing blind docking box(es)…"
        with st.spinner(spinner_msg):
            with capture_log("log_pockets", log_ph):
                pf_out = out_path / "p2rank_outputs"
                try:
                    from docksuitex.batch_docking.batch_grid_calculator import BatchGridCalculator
                    batch = BatchGridCalculator(
                        inputs=[str(p) for p in st.session_state["prepared_proteins"]],
                        mode="p2rank" if is_gc_batch_p2rank else "blind",
                        max_pockets=pf_max if is_gc_batch_p2rank and pf_max > 0 else None,
                    )
                    results = batch.run_all(save_to=pf_out, cpu=global_max_cpu)
                    st.session_state["pocket_results"] = results

                    # Populate pocket text areas from Grid Calculator results
                    for receptor_path, pockets in results.items():
                        stem = Path(receptor_path).stem
                        text = "\n".join(_fmt_pocket(p) for p in pockets)
                        st.session_state[f"gc_{stem}"] = text   # original predictions
                        st.session_state[f"mc_{stem}"] = text   # editable copy
                except Exception as e:
                    st.error(str(e))
                    with st.expander("Traceback"): st.code(traceback.format_exc())
        log_ph.empty()
render_log_if_present("log_pockets", "Grid Calculator Log:")


# ════════════════════════════════════════════════════════════════
# BATCH DOCKING
# ════════════════════════════════════════════════════════════════
st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)
st.markdown('<div class="section-header"><h2>⚡ Batch Docking</h2></div>', unsafe_allow_html=True)

dock_engine = st.radio("Docking Engine", ["AutoDock Vina", "AutoDock4"], horizontal=True, key="dock_eng")

proteins_ready = st.session_state["prepared_proteins"] not in (None, "", [])
ligands_ready  = st.session_state["prepared_ligands"] not in (None, "", [])
pockets_ready  = st.session_state["pocket_results"] not in (None, "", {})

if not (proteins_ready and ligands_ready):
    st.info("Prepare both proteins and ligands above to enable docking.")
elif not pockets_ready:
    st.info("Run Pocket Finder above, or enter manual centers below.")

st.radio("Pocket Source", ["Grid Calculator", "Manual Pockets"],
         key="pocket_mode", horizontal=True, on_change=restore_gc_pockets)
use_manual = st.session_state["pocket_mode"] == "Manual Pockets"
if proteins_ready:
    st.caption("Enter/Adjust pockets for each receptor. Format: `(center_x, center_y, center_z):(size_x, size_y, size_z)` per line.")
    for p in st.session_state["prepared_proteins"]:
        st.text_area(f"Pockets for {p.stem}",
                     placeholder="e.g. (10.5, 20.3, 30.1):(20.0, 20.0, 20.0)",
                     key=f"mc_{p.stem}",
                     on_change=activate_manual_mode)


def build_rwc():
    """Build the receptors_with_pockets dict for batch docking engines.

    Converts GridCalculator pocket results (list of pocket dicts) or
    manual entries to the dictionary format expected by docking engines.
    """
    import re
    # Pattern: (cx, cy, cz):(sx, sy, sz)  — parens and colon-size are optional
    _pat = re.compile(
        r"\(?\s*([^,:\s\)]+)\s*,\s*([^,:\s\)]+)\s*,\s*([^,:\s\)]+)\s*\)?"
        r"(?:\s*:\s*\(?\s*([^,:\s\)]+)\s*,\s*([^,:\s\)]+)\s*,\s*([^,:\s\)]+)\s*\)?)?"
    )
    default_gs = (20.0, 20.0, 20.0) if dock_engine == "AutoDock Vina" else (22.5, 22.5, 22.5)

    if use_manual:
        m_rwc = {}
        for p in st.session_state["prepared_proteins"]:
            mv = st.session_state.get(f"mc_{p.stem}")
            if mv:
                try:
                    pockets = []
                    for line in re.split(r"[\n;]", mv):
                        line = line.strip()
                        if not line:
                            continue
                        m = _pat.search(line)
                        if m:
                            g = m.groups()
                            center = tuple(float(x) for x in g[:3])
                            size = tuple(float(x) for x in g[3:]) if g[3] is not None else default_gs
                            pockets.append({"center": center, "grid_size": size})
                    if pockets:
                        m_rwc[str(p)] = pockets
                except Exception:
                    continue
        return m_rwc
    elif pockets_ready:
        return st.session_state["pocket_results"]
    return {}


with st.expander("⚙️ Docking Parameters", expanded=True):
    if dock_engine == "AutoDock Vina":
        d1, d2, d3 = st.columns(3)
        with d1:
            v_exh = st.number_input("Exhaustiveness", value=8, min_value=1, key="v_exh", help="Search exhaustiveness (higher = more accurate but slower)")
        with d2:
            v_nm  = st.number_input("Num modes", value=9, min_value=1, key="v_nm", help="Max number of binding modes to generate")
        with d3:
            v_seed = st.number_input("Seed (0 = auto)", value=42, min_value=0, key="v_seed", help="Random seed for reproducibility")
    else:
        st.markdown("##### Grid Settings")
        g1, g2 = st.columns(2)
        with g1:
            a_sp  = st.number_input("Spacing (Å)", value=0.375, format="%.3f", key="a_sp", help="Spacing between grid points in Å")
            a_diel   = st.number_input("Dielectric", value=-0.1465, format="%.4f", key="a_diel",
                                        help="Dielectric constant for electrostatics")
        with g2:
            a_smooth = st.number_input("Smooth", value=0.5, format="%.2f", min_value=0.0, key="a_smooth",
                                        help="Smoothing factor for potential maps")

        st.markdown("##### Genetic Algorithm Settings")
        ga1, ga2, ga3 = st.columns(3)
        with ga1:
            a_run   = st.number_input("GA runs", value=10, min_value=1, key="a_run", help="Number of GA runs (independent dockings)")
            a_pop   = st.number_input("Population size", value=150, min_value=1, key="a_pop", help="Number of individuals in population")
        with ga2:
            a_evals = st.number_input("Num evals", value=2500000, min_value=1, key="a_evals", help="Max number of energy evaluations")
            a_gens  = st.number_input("Num generations", value=27000, min_value=1, key="a_gens",
                                      help="Maximum number of generations")
        with ga3:
            a_elit  = st.number_input("Elitism", value=1, min_value=0, key="a_elit",
                                       help="Number of top individuals preserved")
            a_mut   = st.number_input("Mutation rate", value=0.02, format="%.3f", min_value=0.0, max_value=1.0,
                                       key="a_mut", help="Probability of mutation")
            a_cross = st.number_input("Crossover rate", value=0.8, format="%.2f", min_value=0.0, max_value=1.0,
                                       key="a_cross", help="Probability of crossover")

        st.markdown("##### Other Settings")
        o1, o2 = st.columns(2)
        with o1:
            a_rmstol = st.number_input("RMSD tolerance", value=2.0, format="%.2f", min_value=0.0,
                                        key="a_rmstol", help="RMSD tolerance for clustering")
        with o2:
            a_seed_mode = st.selectbox("Seed", ["Custom", "Auto (pid, time)"], key="a_seed_mode")
            if a_seed_mode == "Custom":
                sc1, sc2 = st.columns(2)
                with sc1:
                    a_seed1 = st.number_input("Seed 1", value=27, min_value=0, key="a_seed1", help="First random seed")
                with sc2:
                    a_seed2 = st.number_input("Seed 2", value=6, min_value=0, key="a_seed2", help="Second random seed")

rwc = build_rwc()
can_dock = proteins_ready and ligands_ready and len(rwc) > 0

if st.button("🚀 Run Batch Docking", type="primary", width='stretch', key="btn_dock"):
    if not can_dock:
        if not proteins_ready or not ligands_ready:
            st.warning("Prepare both proteins and ligands first.")
        else:
            st.warning("Run Pocket Finder or enter manual pocket centers first.")
    else:
        total = sum(len(st.session_state["prepared_ligands"]) * len(c) for c in rwc.values())
        log_ph = st.empty()
        batch = None
        with st.spinner(f"Running {total} docking task(s) with {dock_engine}…"):
            with capture_log("log_docking", log_ph):
                dock_out = out_path / ("batch_vina_results" if dock_engine == "AutoDock Vina" else "batch_ad4_results")
                try:
                    lig_paths = [str(p) for p in st.session_state["prepared_ligands"]]
                    if dock_engine == "AutoDock Vina":
                        from docksuitex.batch_docking.batch_vina import BatchVinaDocking
                        batch = BatchVinaDocking(receptors_with_pockets=rwc, ligands=lig_paths,
                                                 exhaustiveness=v_exh, num_modes=v_nm,
                                                 seed=v_seed if v_seed > 0 else None)
                        results = batch.run_all(cpu=global_max_cpu, save_to=dock_out)
                    else:
                        from docksuitex.batch_docking.batch_autodock4 import BatchAD4Docking
                        ad4_seed = (27, 6)
                        if a_seed_mode == "Custom":
                            ad4_seed = (a_seed1, a_seed2)
                        elif a_seed_mode == "Auto (pid, time)":
                            ad4_seed = ("pid", "time")
                        batch = BatchAD4Docking(receptors_with_pockets=rwc, ligands=lig_paths,
                                                spacing=a_sp,
                                                dielectric=a_diel,
                                                smooth=a_smooth,
                                                ga_pop_size=a_pop,
                                                ga_num_evals=a_evals,
                                                ga_num_generations=a_gens,
                                                ga_elitism=a_elit,
                                                ga_mutation_rate=a_mut,
                                                ga_crossover_rate=a_cross,
                                                ga_run=a_run,
                                                rmstol=a_rmstol,
                                                seed=ad4_seed)
                        results = batch.run_all(cpu=global_max_cpu, save_to=dock_out)

                    # Parse results immediately after docking (like single_docking)
                    # csv_name = "vina_summary.csv" if dock_engine == "AutoDock Vina" else "ad4_summary.csv"
                    # csv_path = dock_out / csv_name
                    # df = batch.parse_results(save_to=csv_path)
                    df = batch.parse_results()

                    st.session_state["dock_engine_used"] = dock_engine
                    st.session_state["docking_results"] = dock_out

                    succeeded = {k: v for k, v in results.items() if isinstance(v, Path)}
                    failed    = {k: v for k, v in results.items() if isinstance(v, str)}
                    m1, m2, m3 = st.columns(3)
                    with m1: st.markdown(f'<div class="metric-card"><div class="metric-value">{total}</div><div class="metric-label">Total Tasks</div></div>', unsafe_allow_html=True)
                    with m2: st.markdown(f'<div class="metric-card"><div class="metric-value">{len(succeeded)}</div><div class="metric-label">Succeeded</div></div>', unsafe_allow_html=True)
                    with m3: st.markdown(f'<div class="metric-card"><div class="metric-value">{len(failed)}</div><div class="metric-label">Failed</div></div>', unsafe_allow_html=True)
                    if failed:
                        with st.expander("⚠️ Errors"):
                            for key, err in failed.items():
                                st.error(f"{key[0]} + {key[1]} @ {key[2]}: {err}")
                except Exception as e:
                    st.error(str(e))
                    with st.expander("Traceback"): st.code(traceback.format_exc())
        log_ph.empty()

        # Run interaction profiler automatically after docking
        if st.session_state.get("docking_results") and batch is not None:
            try:
                prolif_log_ph = st.empty()
                with st.spinner("Running interaction profiler…"):
                    with capture_log("log_prolif", prolif_log_ph):
                        ip_df = batch.interaction_profile(
                            cpu=global_max_cpu,
                            # save_to=dock_out / "prolif_results"
                        )
                        st.session_state["interaction_profile"] = ip_df
                prolif_log_ph.empty()
                # Append ProLIF log to Docking Log
                prolif_log = st.session_state.get("log_prolif", "")
                if prolif_log:
                    st.session_state["log_docking"] = st.session_state.get("log_docking", "") + "\n" + prolif_log
            except Exception as e:
                st.warning(f"Interaction profiler failed: {e}")
                with st.expander("Profiler Traceback"): st.code(traceback.format_exc())

render_log_if_present("log_docking", "Docking Log:")

# ════════════════════════════════════════════════════════════════
# RESULTS VIEWER  (outside button handler so it persists)
# ════════════════════════════════════════════════════════════════
if st.session_state.get("docking_results"):
    results_dir = st.session_state["docking_results"]
    engine_used = st.session_state.get("dock_engine_used", "AutoDock Vina")
    st.markdown("---")
    st.markdown("##### 📊 Docking Results")

    csv_name = "vina_summary.csv" if engine_used == "AutoDock Vina" else "ad4_summary.csv"
    csv_path = Path(results_dir) / csv_name

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        st.markdown(f'<span class="status-badge">✅ {len(df)} result(s) found</span>', unsafe_allow_html=True)

        ecol = "Affinity (kcal/mol)" if "Affinity (kcal/mol)" in df.columns else "Binding_Energy"

        if ecol in df.columns and not df[ecol].isna().all():
            mc1, mc2, mc3 = st.columns(3)
            with mc1: st.markdown(f'<div class="metric-card"><div class="metric-value">{len(df)}</div><div class="metric-label">Total Results</div></div>', unsafe_allow_html=True)
            with mc2: st.markdown(f'<div class="metric-card"><div class="metric-value">{df[ecol].min():.2f}</div><div class="metric-label">Best Energy</div></div>', unsafe_allow_html=True)
            with mc3: st.markdown(f'<div class="metric-card"><div class="metric-value">{df[ecol].mean():.2f}</div><div class="metric-label">Mean Energy</div></div>', unsafe_allow_html=True)

        st.dataframe(df, width='stretch', hide_index=True)

        # ── Interaction Profile ──
        if st.session_state.get("interaction_profile") is not None:
            st.markdown("##### 🧪 Interaction Profile")
            ip_display = st.session_state["interaction_profile"].reset_index()
            st.dataframe(ip_display, width='stretch', hide_index=True)


    else:
        st.info("Docking completed but no results CSV found yet.")



render_files_panel(out_path)

st.caption("DockSuiteX © 2026")
