"""
DockSuiteX — Streamlit App (Single Page)
All-in-one batch molecular docking workflow.
Layout: Ligand Prep (left) | Protein Prep (right) | Pocket Finder | Docking | Results
"""
import streamlit as st
import tempfile, shutil, zipfile, io, os, traceback
from pathlib import Path
import pandas as pd

# ── Page config ────────────────────────────────────────────────
st.set_page_config(
    page_title="DockSuiteX",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ── Load CSS ───────────────────────────────────────────────────
css_path = Path(__file__).parent / "style.css"
if css_path.exists():
    st.markdown(f"<style>{css_path.read_text()}</style>", unsafe_allow_html=True)

# ── Persist state across reruns ────────────────────────────────
for key in [
    "prepared_ligands", "prepared_proteins",
    "pocket_results", "docking_results",
]:
    if key not in st.session_state:
        st.session_state[key] = None

if "protein_tmp_dir" not in st.session_state:
    st.session_state["protein_tmp_dir"] = None
if "ligand_tmp_dir" not in st.session_state:
    st.session_state["ligand_tmp_dir"] = None

# ── Hero ───────────────────────────────────────────────────────
st.markdown(
    """
    <div class="hero-header animate-in">
        <h1>🧬 DockSuiteX</h1>
        <p>All-in-one batch molecular docking workflow &mdash;
        prepare ligands &amp; proteins, find pockets, dock, and analyze results.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

# ════════════════════════════════════════════════════════════════
# ROW 1: LIGAND PREP (left)  |  PROTEIN PREP (right)
# ════════════════════════════════════════════════════════════════

col_lig, col_prot = st.columns(2, gap="large")

# ── LEFT: Ligand Preparation ──────────────────────────────────
with col_lig:
    st.markdown(
        '<div class="section-header"><h2>💊 Ligand Preparation</h2></div>',
        unsafe_allow_html=True,
    )

    lig_files = st.file_uploader(
        "Upload ligand files (MOL2 / SDF / PDB / MOL / SMI)",
        type=["mol2", "sdf", "pdb", "mol", "smi"],
        accept_multiple_files=True,
        key="lig_upload",
    )

    with st.expander("⚙️ Ligand Options", expanded=False):
        lig_minimize = st.selectbox(
            "Energy minimization",
            [None, "mmff94", "mmff94s", "uff", "gaff"],
            format_func=lambda x: "None (skip)" if x is None else x.upper(),
            key="lig_min",
        )
        l_c1, l_c2 = st.columns(2)
        with l_c1:
            lig_remove_water = st.checkbox("Remove water", value=True, key="lig_rw")
            lig_add_h = st.checkbox("Add hydrogens", value=True, key="lig_h")
        with l_c2:
            lig_add_charges = st.checkbox("Add Gasteiger charges", value=True, key="lig_chg")
            lig_preserve = st.text_input("Preserve charge types", placeholder="e.g. Zn, Fe", key="lig_pres")

    lig_preserve_list = (
        [x.strip() for x in lig_preserve.split(",") if x.strip()] if lig_preserve else None
    )

    if lig_files and st.button("🚀 Prepare Ligands", type="primary", use_container_width=True, key="btn_lig"):
        with st.spinner(f"Preparing {len(lig_files)} ligand(s)…"):
            tmp_dir = Path(tempfile.mkdtemp())
            input_dir = tmp_dir / "inputs"
            output_dir = tmp_dir / "outputs"
            input_dir.mkdir(); output_dir.mkdir()
            try:
                for f in lig_files:
                    (input_dir / f.name).write_bytes(f.getvalue())

                from docksuitex.batch_docking.batch_ligand import BatchLigand
                batch = BatchLigand(
                    inputs=input_dir,
                    minimize=lig_minimize,
                    remove_water=lig_remove_water,
                    add_hydrogens=lig_add_h,
                    add_charges=lig_add_charges,
                    preserve_charge_types=lig_preserve_list,
                )
                results = batch.prepare_all(save_to=output_dir, cpu=os.cpu_count() or 1)

                success = [r for r in results if r["status"] == "success"]
                errors = [r for r in results if r["status"] == "error"]

                st.markdown(
                    f'<span class="status-badge status-success">✅ {len(success)} prepared</span>'
                    + (f' <span class="status-badge status-error">❌ {len(errors)} failed</span>' if errors else ""),
                    unsafe_allow_html=True,
                )
                if errors:
                    for r in errors:
                        st.error(f"{Path(r['file']).name}: {r.get('error','')}")

                if success:
                    st.session_state["prepared_ligands"] = [Path(r["pdbqt_path"]) for r in success]
                    st.session_state["ligand_tmp_dir"] = tmp_dir
                    # Download ZIP
                    zbuf = io.BytesIO()
                    with zipfile.ZipFile(zbuf, "w", zipfile.ZIP_DEFLATED) as zf:
                        for r in success:
                            p = Path(r["pdbqt_path"])
                            if p.exists(): zf.write(p, p.name)
                    zbuf.seek(0)
                    st.download_button("⬇️ Download Ligands (ZIP)", zbuf, "prepared_ligands.zip", "application/zip")
            except Exception as e:
                st.error(str(e))
                with st.expander("Traceback"): st.code(traceback.format_exc())

    if st.session_state["prepared_ligands"]:
        st.caption(f"✅ {len(st.session_state['prepared_ligands'])} ligand(s) ready for docking")

# ── RIGHT: Protein Preparation ────────────────────────────────
with col_prot:
    st.markdown(
        '<div class="section-header"><h2>🧬 Protein Preparation</h2></div>',
        unsafe_allow_html=True,
    )

    prot_files = st.file_uploader(
        "Upload protein files (PDB)",
        type=["pdb"],
        accept_multiple_files=True,
        key="prot_upload",
    )

    with st.expander("⚙️ Protein Options", expanded=False):
        p_c1, p_c2 = st.columns(2)
        with p_c1:
            prot_fix = st.checkbox("Fix PDB", value=True, key="prot_fix")
            prot_het = st.checkbox("Remove heterogens", value=True, key="prot_het")
            prot_water = st.checkbox("Remove water", value=True, key="prot_water")
        with p_c2:
            prot_h = st.checkbox("Add hydrogens", value=True, key="prot_h")
            prot_chg = st.checkbox("Add Gasteiger charges", value=True, key="prot_chg")
            prot_preserve = st.text_input("Preserve charge types", placeholder="e.g. Zn, Fe", key="prot_pres")

    prot_preserve_list = (
        [x.strip() for x in prot_preserve.split(",") if x.strip()] if prot_preserve else None
    )

    if prot_files and st.button("🚀 Prepare Proteins", type="primary", use_container_width=True, key="btn_prot"):
        with st.spinner(f"Preparing {len(prot_files)} protein(s)…"):
            tmp_dir = Path(tempfile.mkdtemp())
            input_dir = tmp_dir / "inputs"
            output_dir = tmp_dir / "outputs"
            input_dir.mkdir(); output_dir.mkdir()
            try:
                for f in prot_files:
                    (input_dir / f.name).write_bytes(f.getvalue())

                from docksuitex.batch_docking.batch_protein import BatchProtein
                batch = BatchProtein(
                    inputs=input_dir,
                    fix_pdb=prot_fix,
                    remove_heterogens=prot_het,
                    remove_water=prot_water,
                    add_hydrogens=prot_h,
                    add_charges=prot_chg,
                    preserve_charge_types=prot_preserve_list,
                )
                results = batch.prepare_all(save_to=output_dir, cpu=os.cpu_count() or 1)

                success = [r for r in results if r["status"] == "success"]
                errors = [r for r in results if r["status"] == "error"]

                st.markdown(
                    f'<span class="status-badge status-success">✅ {len(success)} prepared</span>'
                    + (f' <span class="status-badge status-error">❌ {len(errors)} failed</span>' if errors else ""),
                    unsafe_allow_html=True,
                )
                if errors:
                    for r in errors:
                        st.error(f"{Path(r['file']).name}: {r.get('error','')}")

                if success:
                    st.session_state["prepared_proteins"] = [Path(r["pdbqt_path"]) for r in success]
                    st.session_state["protein_tmp_dir"] = tmp_dir
                    zbuf = io.BytesIO()
                    with zipfile.ZipFile(zbuf, "w", zipfile.ZIP_DEFLATED) as zf:
                        for r in success:
                            p = Path(r["pdbqt_path"])
                            if p.exists(): zf.write(p, p.name)
                    zbuf.seek(0)
                    st.download_button("⬇️ Download Proteins (ZIP)", zbuf, "prepared_proteins.zip", "application/zip")
            except Exception as e:
                st.error(str(e))
                with st.expander("Traceback"): st.code(traceback.format_exc())

    if st.session_state["prepared_proteins"]:
        st.caption(f"✅ {len(st.session_state['prepared_proteins'])} protein(s) ready")

    # ── POCKET FINDER (below proteins) ─────────────────────────
    st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="section-header"><h2>🎯 Pocket Finder</h2></div>',
        unsafe_allow_html=True,
    )
    st.caption("Predict binding pockets on prepared proteins using P2Rank.")

    pf_max = st.number_input("Max pockets per protein (0 = all)", min_value=0, value=3, key="pf_max")

    can_find = st.session_state["prepared_proteins"] is not None
    if not can_find:
        st.info("Prepare proteins above first.")

    if can_find and st.button("🔍 Find Pockets", type="primary", use_container_width=True, key="btn_pf"):
        with st.spinner("Running P2Rank…"):
            pf_tmp = Path(tempfile.mkdtemp())
            try:
                from docksuitex.batch_docking.batch_pocket_finder import BatchPocketFinder
                batch = BatchPocketFinder(
                    inputs=[str(p) for p in st.session_state["prepared_proteins"]],
                    max_centres=pf_max if pf_max > 0 else None,
                )
                results = batch.run_all(save_to=pf_tmp, cpu=os.cpu_count() or 1)

                st.session_state["pocket_results"] = results

                for fpath, centers in results.items():
                    name = Path(fpath).stem
                    st.markdown(f"**{name}** — {len(centers)} pocket(s)")
                    data = [{"Rank": i+1, "X": round(c[0],3), "Y": round(c[1],3), "Z": round(c[2],3)}
                            for i, c in enumerate(centers)]
                    if data:
                        st.dataframe(pd.DataFrame(data), use_container_width=True, hide_index=True)
            except Exception as e:
                st.error(str(e))
                with st.expander("Traceback"): st.code(traceback.format_exc())

    if st.session_state["pocket_results"]:
        st.caption(f"✅ Pockets found for {len(st.session_state['pocket_results'])} protein(s)")

# ════════════════════════════════════════════════════════════════
# ROW 2: DOCKING
# ════════════════════════════════════════════════════════════════

st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-header"><h2>⚡ Batch Docking</h2></div>',
    unsafe_allow_html=True,
)

dock_engine = st.radio("Docking Engine", ["AutoDock Vina", "AutoDock4"], horizontal=True, key="dock_eng")

# Check readiness
proteins_ready = st.session_state["prepared_proteins"] is not None
ligands_ready = st.session_state["prepared_ligands"] is not None
pockets_ready = st.session_state["pocket_results"] is not None

if not (proteins_ready and ligands_ready):
    st.info("Prepare both proteins and ligands above to enable docking.")
elif not pockets_ready:
    st.info("Run Pocket Finder above, or enter manual centers below.")

# Manual centers option
use_manual = st.checkbox("Enter pocket centers manually (skip Pocket Finder)", key="manual_centers")

manual_rwc = {}
if use_manual and proteins_ready:
    st.caption("Enter centers for each receptor. Format: `x, y, z` separated by `;` for multiple.")
    for p in st.session_state["prepared_proteins"]:
        val = st.text_input(
            f"Centers for {p.stem}",
            placeholder="e.g. 10.5, 20.3, 30.1 ; 15.0, 25.0, 35.0",
            key=f"mc_{p.stem}",
        )
        if val:
            try:
                parsed = []
                for group in val.split(";"):
                    parts = [float(x.strip()) for x in group.strip().split(",")]
                    if len(parts) == 3:
                        parsed.append(tuple(parts))
                if parsed:
                    manual_rwc[str(p)] = parsed
            except ValueError:
                st.warning(f"Could not parse centers for {p.stem}")

# Build receptors_with_centers
def build_rwc():
    if use_manual:
        return manual_rwc
    elif pockets_ready:
        return st.session_state["pocket_results"]
    return {}

# ── Docking parameters ─────────────────────────────────────────
with st.expander("⚙️ Docking Parameters", expanded=False):
    if dock_engine == "AutoDock Vina":
        d1, d2, d3 = st.columns(3)
        with d1:
            v_gsx = st.number_input("Grid Size X", value=20, min_value=1, key="v_gsx")
            v_gsy = st.number_input("Grid Size Y", value=20, min_value=1, key="v_gsy")
            v_gsz = st.number_input("Grid Size Z", value=20, min_value=1, key="v_gsz")
        with d2:
            v_exh = st.number_input("Exhaustiveness", value=8, min_value=1, key="v_exh")
            v_nm = st.number_input("Num modes", value=9, min_value=1, key="v_nm")
        with d3:
            v_seed = st.number_input("Seed (0 = auto)", value=0, min_value=0, key="v_seed")
    else:
        d1, d2, d3 = st.columns(3)
        with d1:
            a_gsx = st.number_input("Grid Size X", value=40, min_value=1, key="a_gsx")
            a_gsy = st.number_input("Grid Size Y", value=40, min_value=1, key="a_gsy")
            a_gsz = st.number_input("Grid Size Z", value=40, min_value=1, key="a_gsz")
        with d2:
            a_sp = st.number_input("Spacing (Å)", value=0.375, format="%.3f", key="a_sp")
            a_run = st.number_input("GA runs", value=10, min_value=1, key="a_run")
        with d3:
            a_pop = st.number_input("Pop size", value=150, min_value=1, key="a_pop")
            a_evals = st.number_input("Num evals", value=2500000, min_value=1, key="a_evals")

cpu_count = os.cpu_count() or 1
dock_cpus = st.slider("CPUs for docking", 1, cpu_count, cpu_count, key="dock_cpu")

rwc = build_rwc()
can_dock = proteins_ready and ligands_ready and len(rwc) > 0

if can_dock and st.button("🚀 Run Batch Docking", type="primary", use_container_width=True, key="btn_dock"):
    total = sum(len(st.session_state["prepared_ligands"]) * len(c) for c in rwc.values())
    with st.spinner(f"Running {total} docking task(s) with {dock_engine}…"):
        dock_tmp = Path(tempfile.mkdtemp())
        try:
            lig_paths = [str(p) for p in st.session_state["prepared_ligands"]]

            if dock_engine == "AutoDock Vina":
                from docksuitex.batch_docking.batch_vina import BatchVinaDocking
                batch = BatchVinaDocking(
                    receptors_with_centers=rwc,
                    ligands=lig_paths,
                    grid_size=(v_gsx, v_gsy, v_gsz),
                    exhaustiveness=v_exh,
                    num_modes=v_nm,
                    seed=v_seed if v_seed > 0 else None,
                )
                results = batch.run_all(cpu=dock_cpus, save_to=dock_tmp)
            else:
                from docksuitex.batch_docking.batch_autodock4 import BatchAD4Docking
                batch = BatchAD4Docking(
                    receptors_with_centers=rwc,
                    ligands=lig_paths,
                    grid_size=(a_gsx, a_gsy, a_gsz),
                    spacing=a_sp,
                    ga_run=a_run,
                    ga_pop_size=a_pop,
                    ga_num_evals=a_evals,
                )
                results = batch.run_all(cpu=dock_cpus, save_to=dock_tmp)

            st.session_state["docking_results"] = dock_tmp

            succeeded = {k: v for k, v in results.items() if isinstance(v, Path)}
            failed = {k: v for k, v in results.items() if isinstance(v, str)}

            m1, m2, m3 = st.columns(3)
            with m1:
                st.markdown(
                    f'<div class="metric-card"><div class="metric-value">{total}</div>'
                    f'<div class="metric-label">Total Tasks</div></div>',
                    unsafe_allow_html=True)
            with m2:
                st.markdown(
                    f'<div class="metric-card"><div class="metric-value">{len(succeeded)}</div>'
                    f'<div class="metric-label">Succeeded</div></div>',
                    unsafe_allow_html=True)
            with m3:
                st.markdown(
                    f'<div class="metric-card"><div class="metric-value">{len(failed)}</div>'
                    f'<div class="metric-label">Failed</div></div>',
                    unsafe_allow_html=True)

            if failed:
                with st.expander("⚠️ Errors"):
                    for key, err in failed.items():
                        st.error(f"{key[0]} + {key[1]} @ {key[2]}: {err}")

            if succeeded:
                zbuf = io.BytesIO()
                with zipfile.ZipFile(zbuf, "w", zipfile.ZIP_DEFLATED) as zf:
                    for f in dock_tmp.rglob("*"):
                        if f.is_file(): zf.write(f, f.relative_to(dock_tmp))
                zbuf.seek(0)
                st.download_button(
                    "⬇️ Download Docking Results (ZIP)", zbuf,
                    f"batch_{dock_engine.lower().replace(' ','_')}_results.zip",
                    "application/zip",
                )
        except Exception as e:
            st.error(str(e))
            with st.expander("Traceback"): st.code(traceback.format_exc())

# ════════════════════════════════════════════════════════════════
# ROW 3: RESULTS VIEWER
# ════════════════════════════════════════════════════════════════

st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)
st.markdown(
    '<div class="section-header"><h2>📊 View Results</h2></div>',
    unsafe_allow_html=True,
)

result_tab = st.radio("Parse results from", ["Docking output above", "Upload ZIP"], horizontal=True, key="res_src")

results_dir = None

if result_tab == "Docking output above":
    if st.session_state["docking_results"]:
        results_dir = st.session_state["docking_results"]
        st.caption(f"Scanning: `{results_dir}`")
    else:
        st.info("Run docking above first, or switch to 'Upload ZIP'.")
else:
    result_zip = st.file_uploader("Upload results ZIP", type=["zip"], key="res_zip")
    if result_zip:
        results_dir = Path(tempfile.mkdtemp()) / "extracted"
        results_dir.mkdir()
        with zipfile.ZipFile(io.BytesIO(result_zip.getvalue()), "r") as zf:
            zf.extractall(results_dir)

parser_type = st.radio("Parser", ["AutoDock Vina", "AutoDock4"], horizontal=True, key="parser_type")

if results_dir and st.button("🔍 Parse Results", type="primary", use_container_width=True, key="btn_parse"):
    with st.spinner("Parsing…"):
        try:
            csv_path = Path(tempfile.mkdtemp()) / "results.csv"

            if parser_type == "AutoDock Vina":
                from docksuitex.utils.parser import parse_vina_log_to_csv
                df = parse_vina_log_to_csv(results_dir=results_dir, output_csv=csv_path)
                energy_col = "Affinity (kcal/mol)"
            else:
                from docksuitex.utils.parser import parse_ad4_dlg_to_csv
                df = parse_ad4_dlg_to_csv(results_dir=results_dir, output_csv=csv_path)
                energy_col = "Binding_Energy"

            st.markdown(
                f'<span class="status-badge status-success">✅ {len(df)} result(s) parsed</span>',
                unsafe_allow_html=True,
            )

            # Metrics
            if energy_col in df.columns and not df[energy_col].isna().all():
                mc1, mc2, mc3 = st.columns(3)
                with mc1:
                    st.markdown(
                        f'<div class="metric-card"><div class="metric-value">{len(df)}</div>'
                        f'<div class="metric-label">Total Results</div></div>',
                        unsafe_allow_html=True)
                with mc2:
                    st.markdown(
                        f'<div class="metric-card"><div class="metric-value">{df[energy_col].min():.2f}</div>'
                        f'<div class="metric-label">Best Energy (kcal/mol)</div></div>',
                        unsafe_allow_html=True)
                with mc3:
                    st.markdown(
                        f'<div class="metric-card"><div class="metric-value">{df[energy_col].mean():.2f}</div>'
                        f'<div class="metric-label">Mean Energy (kcal/mol)</div></div>',
                        unsafe_allow_html=True)

            st.dataframe(df, use_container_width=True, hide_index=True)

            # Bar chart
            if energy_col in df.columns and "Ligand" in df.columns:
                chart_df = df.sort_values(energy_col).head(30).copy()
                if "Mode" in chart_df.columns:
                    chart_df["Label"] = chart_df["Ligand"] + " (m" + chart_df["Mode"].astype(str) + ")"
                elif "Cluster_Rank" in chart_df.columns:
                    chart_df["Label"] = chart_df["Ligand"] + " (r" + chart_df["Cluster_Rank"].astype(str) + ")"
                else:
                    chart_df["Label"] = chart_df["Ligand"]
                st.bar_chart(chart_df.set_index("Label")[energy_col])

            # Download CSV
            csv_data = csv_path.read_bytes() if csv_path.exists() else df.to_csv(index=False).encode()
            st.download_button("⬇️ Download CSV", csv_data, "docking_results.csv", "text/csv")

        except Exception as e:
            st.error(str(e))
            with st.expander("Traceback"): st.code(traceback.format_exc())

# ── Footer ─────────────────────────────────────────────────────
st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)
st.caption("Built with ❤️ using Streamlit · DockSuiteX © 2025 Mangalam Gautam Sinha")
