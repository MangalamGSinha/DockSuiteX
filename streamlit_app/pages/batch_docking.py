"""
DockSuiteX — Batch Docking Page
Batch protein + ligand workflow — no visualization.
"""
import streamlit as st
import io, os, traceback, contextlib, shutil
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


def render_input_section(title, key_prefix, storage_key, allowed_exts):
    """Renders a file upload section WITHOUT visualization."""
    st.markdown(f'<div class="section-header"><h2>{title}</h2></div>', unsafe_allow_html=True)
    uploaded = st.file_uploader(
        f"Upload files ({' / '.join([x.upper() for x in allowed_exts])})",
        type=allowed_exts,
        accept_multiple_files=True,
        key=f"{key_prefix}_uploader",
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
    except Exception:
        # Ignore if running in a context without session state (e.g. initial import)
        pass

init_session_state()

# ── Hero header ────────────────────────────────────────────────
st.markdown(
    """
    <div class="hero-header animate-in">
        <h1>⚡ Batch Docking</h1>
        <p>Batch molecular docking workflow &mdash; prepare multiple ligands &amp;
        proteins, find pockets, dock all combinations, and analyze results.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

# ── Global Settings ────────────────────────────────────────────
if "global_out_dir" not in st.session_state:
    st.session_state["global_out_dir"] = str(Path.home() / "Desktop" / "docksuitex_results")


def select_folder():
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk(); root.withdraw(); root.wm_attributes('-topmost', 1)
        folder = filedialog.askdirectory(master=root); root.destroy()
        if folder:
            st.session_state["global_out_dir"] = folder
    except Exception as e:
        st.error(f"Could not open folder dialog: {e}")


st.markdown('<div class="section-header"><h2>🌍 Global Settings</h2></div>', unsafe_allow_html=True)
selected_out_dir = st.session_state["global_out_dir"]

with st.container():
    col_btn, col_txt, col_cpu = st.columns([0.15, 0.55, 0.3])
    with col_btn:
        st.write(""); st.write("")
        st.button("📂 Choose...", on_click=select_folder, help="Pick a folder", key="m_choose_dir")
    with col_txt:
        val = st.text_input("Output Directory", key="global_out_dir")
        if val: selected_out_dir = val
    with col_cpu:
        cpu_avail = os.cpu_count() or 1
        global_max_cpu = st.slider("Max CPUs", 1, cpu_avail, cpu_avail, key="m_max_cpu")

from datetime import datetime as _dt
if "m_run_dir" not in st.session_state:
    st.session_state["m_run_dir"] = _dt.now().strftime("%Y-%m-%d_%H-%M-%S")
out_path = Path(selected_out_dir).resolve() / st.session_state["m_run_dir"]
if not out_path.exists():
    try: out_path.mkdir(parents=True, exist_ok=True)
    except Exception as e: st.error(f"Could not create output directory: {e}")

st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)

# ════════════════════════════════════════════════════════════════
# ROW 1: PROTEIN PREP (left)  |  LIGAND PREP (right)
# ════════════════════════════════════════════════════════════════
col_prot, col_lig = st.columns(2, gap="large")

# ── Protein Preparation ───────────────────────────────────────
with col_prot:
    render_input_section("🧬 Protein Preparation", "prot", "protein_inputs", ["pdb"])
    all_proteins = st.session_state["protein_inputs"]

    with st.expander("⚙️ Protein Options", expanded=True):
        p1, p2 = st.columns(2)
        with p1:
            prot_fix   = st.checkbox("Fix PDB", True, key="prot_fix")
            prot_het   = st.checkbox("Remove heterogens", True, key="prot_het")
            prot_water = st.checkbox("Remove water", True, key="prot_water")
            prot_chg = st.checkbox("Add Gasteiger charges", True, key="prot_chg")
            prot_h   = st.checkbox("Add hydrogens", True, key="prot_h")
        with p2:
            prot_preserve = st.text_input("Preserve charge types", placeholder="e.g. Zn, Fe", key="prot_pres")
    prot_preserve_list = parse_comma_list(prot_preserve)

    if st.button("🚀 Prepare Proteins", type="primary", use_container_width=True, key="btn_prot"):
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
                                             remove_water=prot_water, add_hydrogens=prot_h, add_charges=prot_chg,
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
    render_input_section("💊 Ligand Preparation", "lig", "ligand_inputs", ["mol2","sdf","pdb","mol","smi"])
    all_ligands = st.session_state["ligand_inputs"]

    with st.expander("⚙️ Ligand Options", expanded=True):
        lig_minimize = st.selectbox("Energy minimization", [None,"mmff94","mmff94s","uff","gaff"],
                                    format_func=lambda x: "None (skip)" if x is None else x.upper(), key="lig_min")
        l1, l2 = st.columns(2)
        with l1:
            lig_rw = st.checkbox("Remove water", True, key="lig_rw")
            lig_h  = st.checkbox("Add hydrogens", True, key="lig_h")
            lig_chg = st.checkbox("Add Gasteiger charges", True, key="lig_chg")
        with l2:
            lig_preserve = st.text_input("Preserve charge types", placeholder="e.g. Zn, Fe", key="lig_pres")
    lig_preserve_list = parse_comma_list(lig_preserve)

    if st.button("🚀 Prepare Ligands", type="primary", use_container_width=True, key="btn_lig"):
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
                                            add_hydrogens=lig_h, add_charges=lig_chg,
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
# POCKET FINDER (full width, outside columns)
# ════════════════════════════════════════════════════════════════
st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)
st.markdown('<div class="section-header"><h2>🎯 Pocket Finder</h2></div>', unsafe_allow_html=True)
st.caption("Predict binding pockets on prepared proteins using P2Rank.")

pf_max = st.number_input("Max pockets per protein (0 = all)", min_value=0, value=3, key="pf_max")

if st.button("🔍 Find Pockets", type="primary", use_container_width=True, key="btn_pf"):
    can_find = st.session_state["prepared_proteins"] not in (None, "", [])
    if not can_find:
        st.warning("Prepare proteins above first.")
    else:
        log_ph = st.empty()
        with st.spinner("Running P2Rank…"):
            with capture_log("log_pockets", log_ph):
                pf_out = out_path / "pockets_output"
                try:
                    from docksuitex.batch_docking.batch_pocket_finder import BatchPocketFinder
                    batch = BatchPocketFinder(
                        inputs=[str(p) for p in st.session_state["prepared_proteins"]],
                        max_centres=pf_max if pf_max > 0 else None,
                    )
                    results = batch.run_all(save_to=pf_out, cpu=global_max_cpu)
                    st.session_state["pocket_results"] = results
                except Exception as e:
                    st.error(str(e))
                    with st.expander("Traceback"): st.code(traceback.format_exc())
        log_ph.empty()
render_log_if_present("log_pockets", "Pocket Finder Log:")


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

use_manual = st.checkbox("Enter pocket centers manually (skip Pocket Finder)", key="manual_centers")
manual_rwc = {}
if use_manual and proteins_ready:
    st.caption("Enter centers for each receptor. Format: `x, y, z` separated by `;` for multiple.")
    for p in st.session_state["prepared_proteins"]:
        mv = st.text_input(f"Centers for {p.stem}", placeholder="e.g. 10.5, 20.3, 30.1 ; 15.0, 25.0, 35.0", key=f"mc_{p.stem}")
        if mv:
            try:
                parsed = []
                for group in mv.split(";"):
                    parts = [float(x.strip()) for x in group.strip().split(",")]
                    if len(parts) == 3: parsed.append(tuple(parts))
                if parsed: manual_rwc[str(p)] = parsed
            except ValueError:
                st.warning(f"Could not parse centers for {p.stem}")


def build_rwc():
    if use_manual: return manual_rwc
    elif pockets_ready: return st.session_state["pocket_results"]
    return {}


with st.expander("⚙️ Docking Parameters", expanded=True):
    if dock_engine == "AutoDock Vina":
        d1, d2, d3 = st.columns(3)
        with d1:
            v_gsx = st.number_input("Grid Size X", value=20, min_value=1, key="v_gsx", help="Search space size in X dimension (Å)")
            v_gsy = st.number_input("Grid Size Y", value=20, min_value=1, key="v_gsy", help="Search space size in Y dimension (Å)")
            v_gsz = st.number_input("Grid Size Z", value=20, min_value=1, key="v_gsz", help="Search space size in Z dimension (Å)")
        with d2:
            v_exh = st.number_input("Exhaustiveness", value=8, min_value=1, key="v_exh", help="Search exhaustiveness (higher = more accurate but slower)")
            v_nm  = st.number_input("Num modes", value=9, min_value=1, key="v_nm", help="Max number of binding modes to generate")
        with d3:
            v_seed = st.number_input("Seed (0 = auto)", value=0, min_value=0, key="v_seed", help="Random seed for reproducibility")
    else:
        st.markdown("##### Grid Settings")
        g1, g2, g3 = st.columns(3)
        with g1:
            a_gsx = st.number_input("Grid Size X", value=60, min_value=1, key="a_gsx", help="Number of grid points in X dimension")
            a_gsy = st.number_input("Grid Size Y", value=60, min_value=1, key="a_gsy", help="Number of grid points in Y dimension")
        with g2:
            a_gsz = st.number_input("Grid Size Z", value=60, min_value=1, key="a_gsz", help="Number of grid points in Z dimension")
            a_sp  = st.number_input("Spacing (Å)", value=0.375, format="%.3f", key="a_sp", help="Spacing between grid points in Å")
        with g3:
            a_diel   = st.number_input("Dielectric", value=-0.1465, format="%.4f", key="a_diel",
                                        help="Dielectric constant for electrostatics")
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
            a_elit    = st.number_input("Elitism", value=1, min_value=0, key="a_elit",
                                        help="Number of top individuals preserved")
            a_mut     = st.number_input("Mutation rate", value=0.02, format="%.3f", min_value=0.0, max_value=1.0,
                                        key="a_mut", help="Probability of mutation")
            a_cross   = st.number_input("Crossover rate", value=0.8, format="%.2f", min_value=0.0, max_value=1.0,
                                        key="a_cross", help="Probability of crossover")

        st.markdown("##### Other Settings")
        o1, o2 = st.columns(2)
        with o1:
            a_rmstol = st.number_input("RMSD tolerance", value=2.0, format="%.2f", min_value=0.0,
                                        key="a_rmstol", help="RMSD tolerance for clustering")
        with o2:
            a_seed_mode = st.selectbox("Seed", ["Auto (pid, time)", "Custom"], key="a_seed_mode")
            if a_seed_mode == "Custom":
                sc1, sc2 = st.columns(2)
                with sc1:
                    a_seed1 = st.number_input("Seed 1", value=0, min_value=0, key="a_seed1", help="First random seed")
                with sc2:
                    a_seed2 = st.number_input("Seed 2", value=0, min_value=0, key="a_seed2", help="Second random seed")

rwc = build_rwc()
can_dock = proteins_ready and ligands_ready and len(rwc) > 0

if st.button("🚀 Run Batch Docking", type="primary", use_container_width=True, key="btn_dock"):
    if not can_dock:
        if not proteins_ready or not ligands_ready:
            st.warning("Prepare both proteins and ligands first.")
        else:
            st.warning("Run Pocket Finder or enter manual pocket centers first.")
    else:
        total = sum(len(st.session_state["prepared_ligands"]) * len(c) for c in rwc.values())
        log_ph = st.empty()
        with st.spinner(f"Running {total} docking task(s) with {dock_engine}…"):
            with capture_log("log_docking", log_ph):
                dock_out = out_path / ("batch_vina_results" if dock_engine == "AutoDock Vina" else "batch_ad4_results")
                try:
                    lig_paths = [str(p) for p in st.session_state["prepared_ligands"]]
                    if dock_engine == "AutoDock Vina":
                        from docksuitex.batch_docking.batch_vina import BatchVinaDocking
                        batch = BatchVinaDocking(receptors_with_centers=rwc, ligands=lig_paths,
                                                 grid_size=(v_gsx,v_gsy,v_gsz), exhaustiveness=v_exh,
                                                 num_modes=v_nm, seed=v_seed if v_seed>0 else None)
                        results = batch.run_all(cpu=global_max_cpu, save_to=dock_out)
                    else:
                        from docksuitex.batch_docking.batch_autodock4 import BatchAD4Docking
                        ad4_seed = ("pid", "time")
                        if a_seed_mode == "Custom":
                            ad4_seed = (a_seed1, a_seed2)
                        batch = BatchAD4Docking(receptors_with_centers=rwc, ligands=lig_paths,
                                                grid_size=(a_gsx,a_gsy,a_gsz),
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
                    csv_name = "vina_summary.csv" if dock_engine == "AutoDock Vina" else "ad4_summary.csv"
                    csv_path = dock_out / csv_name
                    df = batch.parse_results(save_to=csv_path)

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


    else:
        st.info("Docking completed but no results CSV found yet.")

    st.info(f"Docking complete! Results in `{results_dir}`")

st.caption("DockSuiteX © 2025")
