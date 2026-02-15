"""
DockSuiteX — Single Docking Page
One protein + one ligand workflow with NGL.js 3D visualization.
"""
import streamlit as st
import streamlit.components.v1 as components
import io, os, traceback, base64, contextlib
from pathlib import Path
import pandas as pd


# ── Helpers ────────────────────────────────────────────────────

def ngl_viewer(structures: list, height: int = 500, elem_id: str = None):
    """Render molecules in NGL.js embedded viewer.
    
    Args:
        structures: List of dicts, each with:
            - data: file content string
            - ext: file extension (e.g. "pdb")
            - type: "protein" or "ligand"
        height: Viewer height in pixels.
        elem_id: Optional fixed ID for the viewer container.
    """
    
    # Generate unique viewport ID
    if elem_id:
        viewport_id = elem_id
    else:
        viewport_id = "ngl-default-viewer"

    # Build JS for loading each structure
    load_scripts = []
    
    encoded_data_map = {}
    
    for i, struct in enumerate(structures):
        ext = struct['ext'].lower().replace('.', '')
        # Mapping common extensions
        if ext == "pdbqt": ext = "pdbqt"
        elif ext == "mol2": ext = "mol2" 
        elif ext == "sdf": ext = "sdf"
        else: ext = "pdb"
            
        b64_data = base64.b64encode(struct['data'].encode("utf-8", errors="replace")).decode()
        var_name = f"blob_{i}"
        encoded_data_map[var_name] = b64_data
        
        # Representation string
        if struct.get("type") == "ligand":
            repr_js = """c.addRepresentation("ball+stick", {multipleBond:true});"""
        else:
            repr_js = """
            c.addRepresentation("cartoon", {color:"chainid"});
            c.addRepresentation("ball+stick", {sele:"hetero and not water", multipleBond:true});
            """

        load_script = f"""
        var {var_name} = new Blob([atob("{b64_data}")], {{type:"text/plain"}});
        promises.push(stage.loadFile({var_name}, {{ext:"{ext}"}}).then(function(c){{
            {repr_js}
            return c;
        }}));
        """
        load_scripts.append(load_script)

    js_load_block = "\n".join(load_scripts)

    html_template = """
    <div id="__VIEWPORT_ID__" style="width:100%;height:__HEIGHT__px;border-radius:12px;background:#0e1117;border:1px solid rgba(255,255,255,0.08);"></div>

    <script>
    (function() {

        function initNGL() {
            var stage = new NGL.Stage("__VIEWPORT_ID__", {backgroundColor:"#0e1117"});
            var promises = [];

            __LOAD_BLOCK__

            Promise.all(promises).then(function() {
                stage.autoView();
                stage.handleResize();
            }).catch(function(e) {
                console.error(e);
                document.getElementById("__VIEWPORT_ID__").innerHTML =
                    "<p style='color:red;'>Error loading structure</p>";
            });

            window.addEventListener("resize", function() {
                stage.handleResize();
            });
        }

        if (typeof NGL === "undefined") {
            var script = document.createElement("script");
            script.src = "https://unpkg.com/ngl@2.3.1/dist/ngl.js";
            script.onload = initNGL;
            document.head.appendChild(script);
        } else {
            initNGL();
        }

    })();
    </script>
    """

    html = (
        html_template
        .replace("__VIEWPORT_ID__", viewport_id)
        .replace("__HEIGHT__", str(height))
        .replace("__LOAD_BLOCK__", js_load_block)
    )




    components.html(html, height=height + 10)



def extract_models(pdbqt_content: str) -> list:
    """Extract individual models from a multi-model PDBQT string."""
    models = []
    current_model = []
    in_model = False
    
    lines = pdbqt_content.splitlines()
    # If no MODEL/ENDMDL tags, update to treat whole file as one model
    has_models = any(l.startswith("MODEL") for l in lines)
    if not has_models:
        return [pdbqt_content]

    for line in lines:
        if line.startswith("MODEL"):
            in_model = True
            current_model = [line]
        elif line.startswith("ENDMDL"):
            if in_model:
                current_model.append(line)
                models.append("\n".join(current_model))
                in_model = False
        elif in_model:
            current_model.append(line)
            
    return models


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


def render_log_if_present(session_key: str, title: str = "Logs:"):
    if st.session_state.get(session_key):
        st.caption(title)
        st.code(st.session_state[session_key], language="text")


# ── Session state init ─────────────────────────────────────────
def init_session_state():
    try:
        for key in [
            "s_prepared_protein", "s_prepared_ligand",
            "s_pocket_results", "s_docking_results",
            "s_protein_input", "s_ligand_input",
        ]:
            if key not in st.session_state:
                st.session_state[key] = None

        for key in ["s_log_ligand", "s_log_protein", "s_log_docking", "s_log_pockets"]:
            if key not in st.session_state:
                st.session_state[key] = ""
    except Exception:
        # Ignore if running in a context without session state (e.g. initial import)
        pass

init_session_state()


# ── Hero header ────────────────────────────────────────────────
st.markdown(
    """
    <div class="hero-header animate-in">
        <h1>🔬 Single Docking</h1>
        <p>Prepare one protein &amp; one ligand, find pockets, dock, and visualize
        results with interactive 3D molecular views.</p>
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
    col_btn, col_txt = st.columns([0.15, 0.85])
    with col_btn:
        st.write(""); st.write("")
        st.button("📂 Choose...", on_click=select_folder, help="Pick a folder", key="s_choose_dir")
    with col_txt:
        val = st.text_input("Output Directory", key="global_out_dir")
        if val: selected_out_dir = val

from datetime import datetime as _dt
if "s_run_dir" not in st.session_state:
    st.session_state["s_run_dir"] = _dt.now().strftime("%Y-%m-%d_%H-%M-%S")
out_path = Path(selected_out_dir).resolve() / st.session_state["s_run_dir"]
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
    st.markdown('<div class="section-header"><h2>🧬 Protein Preparation</h2></div>', unsafe_allow_html=True)
    prot_file = st.file_uploader("Upload protein (PDB)", type=["pdb"], key="s_prot_uploader")

    if prot_file:
        st.session_state["s_protein_input"] = {"name": prot_file.name, "data": prot_file.getvalue(), "ext": Path(prot_file.name).suffix}
        try:
            ngl_viewer(
                structures=[{"data": prot_file.getvalue().decode("utf-8", errors="replace"), "ext": Path(prot_file.name).suffix, "type": "protein"}],
                height=300,
                elem_id="ngl-protein-input"
            )
        except Exception as e:
            st.error("Viewer Error:")
            st.code(traceback.format_exc())

    else:
        st.session_state["s_protein_input"] = None

    with st.expander("⚙️ Protein Options", expanded=True):
        p1, p2 = st.columns(2)
        with p1:
            prot_fix   = st.checkbox("Fix PDB", True, key="s_prot_fix")
            prot_het   = st.checkbox("Remove heterogens", True, key="s_prot_het")
            prot_water = st.checkbox("Remove water", True, key="s_prot_water")
            prot_chg = st.checkbox("Add Gasteiger charges", True, key="s_prot_chg")
            prot_h   = st.checkbox("Add hydrogens", True, key="s_prot_h")
        with p2:
            prot_preserve = st.text_input("Preserve charge types", placeholder="e.g. Zn, Fe", key="s_prot_pres")
    prot_preserve_list = parse_comma_list(prot_preserve)

    if st.button("🚀 Prepare Protein", type="primary", use_container_width=True, key="s_btn_prot"):
        inp = st.session_state["s_protein_input"]
        if not inp:
            st.warning("Upload a protein file first.")
        else:
            log_ph = st.empty()
            with st.spinner("Preparing protein…"):
                with capture_log("s_log_protein", log_ph):
                    ipath = out_path / inp["name"]; ipath.parent.mkdir(parents=True, exist_ok=True); ipath.write_bytes(inp["data"])
                    try:
                        from docksuitex import Protein
                        prot = Protein(input=str(ipath), fix_pdb=prot_fix, remove_heterogens=prot_het,
                                       remove_water=prot_water, add_hydrogens=prot_h, add_charges=prot_chg,
                                       preserve_charge_types=prot_preserve_list)
                        st.session_state["s_prepared_protein"] = Path(prot.prepare(save_to=out_path))
                    except Exception as e:
                        st.error(str(e))
                        with st.expander("Traceback"): st.code(traceback.format_exc())
            log_ph.empty()
    render_log_if_present("s_log_protein", "Protein Preparation Log:")

    if st.session_state.get("s_prepared_protein"):
        prepared_prot_path = st.session_state["s_prepared_protein"]
        if prepared_prot_path and prepared_prot_path.exists():
            try:
                ngl_viewer(
                    structures=[{"data": prepared_prot_path.read_text(), "ext": prepared_prot_path.suffix, "type": "protein"}],
                    height=300,
                    elem_id="ngl-protein-prepared"
                )
            except Exception as e:
                st.error("Viewer Error:")
                st.code(traceback.format_exc())


# ── Ligand Preparation ────────────────────────────────────────
with col_lig:
    st.markdown('<div class="section-header"><h2>💊 Ligand Preparation</h2></div>', unsafe_allow_html=True)
    lig_file = st.file_uploader("Upload ligand (MOL2/SDF/PDB/MOL/SMI)", type=["mol2","sdf","pdb","mol","smi"], key="s_lig_uploader")

    if lig_file:
        st.session_state["s_ligand_input"] = {"name": lig_file.name, "data": lig_file.getvalue(), "ext": Path(lig_file.name).suffix}
        try:
            ngl_viewer(
                structures=[{"data": lig_file.getvalue().decode("utf-8", errors="replace"), "ext": Path(lig_file.name).suffix, "type": "ligand"}],
                height=300,
                elem_id="ngl-ligand-input"
            )
        except Exception as e:
            st.error("Viewer Error:")
            st.code(traceback.format_exc())

    else:
        st.session_state["s_ligand_input"] = None

    with st.expander("⚙️ Ligand Options", expanded=True):
        lig_minimize = st.selectbox("Energy minimization", [None,"mmff94","mmff94s","uff","gaff"],
                                    format_func=lambda x: "None (skip)" if x is None else x.upper(), key="s_lig_min")
        l1, l2 = st.columns(2)
        with l1:
            lig_rw  = st.checkbox("Remove water", True, key="s_lig_rw")
            lig_h   = st.checkbox("Add hydrogens", True, key="s_lig_h")
            lig_chg = st.checkbox("Add Gasteiger charges", True, key="s_lig_chg")
        with l2:
            lig_preserve = st.text_input("Preserve charge types", placeholder="e.g. Zn, Fe", key="s_lig_pres")
    lig_preserve_list = parse_comma_list(lig_preserve)

    if st.button("🚀 Prepare Ligand", type="primary", use_container_width=True, key="s_btn_lig"):
        inp = st.session_state["s_ligand_input"]
        if not inp:
            st.warning("Upload a ligand file first.")
        else:
            log_ph = st.empty()
            with st.spinner("Preparing ligand…"):
                with capture_log("s_log_ligand", log_ph):
                    ipath = out_path / inp["name"]; ipath.parent.mkdir(parents=True, exist_ok=True); ipath.write_bytes(inp["data"])
                    try:
                        from docksuitex import Ligand
                        lig = Ligand(input=str(ipath), minimize=lig_minimize, remove_water=lig_rw,
                                     add_hydrogens=lig_h, add_charges=lig_chg,
                                     preserve_charge_types=lig_preserve_list)
                        st.session_state["s_prepared_ligand"] = Path(lig.prepare(save_to=out_path))
                    except Exception as e:
                        st.error(str(e))
                        with st.expander("Traceback"): st.code(traceback.format_exc())
            log_ph.empty()
    render_log_if_present("s_log_ligand", "Ligand Preparation Log:")

    if st.session_state.get("s_prepared_ligand"):
        lp = st.session_state["s_prepared_ligand"]
        if lp.exists():
            try:
                ngl_viewer(
                    structures=[{"data": lp.read_text(errors="replace"), "ext": lp.suffix, "type": "ligand"}],
                    height=300,
                    elem_id="ngl-ligand-prepared"
                )
            except Exception as e:
                st.error("Viewer Error:")
                st.code(traceback.format_exc())


# ════════════════════════════════════════════════════════════════
# POCKET FINDER
# ════════════════════════════════════════════════════════════════
st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)
st.markdown('<div class="section-header"><h2>🎯 Pocket Finder</h2></div>', unsafe_allow_html=True)
st.caption("Predict binding pockets on the prepared protein using P2Rank.")

pf_max = st.number_input("Max pockets (0 = all)", min_value=0, value=3, key="s_pf_max")
can_find = st.session_state.get("s_prepared_protein") is not None
if not can_find:
    st.info("Prepare a protein above first.")

if st.button("🔍 Find Pockets", type="primary", use_container_width=True, key="s_btn_pf"):
    if not can_find:
        st.warning("Prepare a protein above first.")
    else:
        log_ph = st.empty()
        with st.spinner("Running P2Rank…"):
            with capture_log("s_log_pockets", log_ph):
                try:
                    from docksuitex import PocketFinder
                    pf = PocketFinder(input=str(st.session_state["s_prepared_protein"]))
                    pockets = pf.run(save_to=out_path / f"p2rank_results_{st.session_state['s_prepared_protein'].name.replace('.', '_')}")
                    if pf_max > 0:
                        pockets = pockets[:pf_max]
                    st.session_state["s_pocket_results"] = pockets
                except Exception as e:
                    st.error(str(e))
                    with st.expander("Traceback"): st.code(traceback.format_exc())
        log_ph.empty()

render_log_if_present("s_log_pockets", "Pocket Finder Log:")




# ════════════════════════════════════════════════════════════════
# DOCKING
# ════════════════════════════════════════════════════════════════
st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)
st.markdown('<div class="section-header"><h2>⚡ Docking</h2></div>', unsafe_allow_html=True)

dock_engine = st.radio("Docking Engine", ["AutoDock Vina", "AutoDock4"], horizontal=True, key="s_dock_eng")
protein_ready = st.session_state.get("s_prepared_protein") is not None
ligand_ready  = st.session_state.get("s_prepared_ligand") is not None
pockets_ready = st.session_state.get("s_pocket_results") is not None

if not (protein_ready and ligand_ready):
    st.info("Prepare both protein and ligand above to enable docking.")
elif not pockets_ready:
    st.info("Run Pocket Finder above, or enter a manual center below.")

use_manual = st.checkbox("Enter pocket center manually", key="s_manual_center")
grid_center = None
if use_manual:
    mc_val = st.text_input("Pocket center (x, y, z)", placeholder="e.g. 32.14, 47.70, 13.75", key="s_mc_input")
    if mc_val:
        try:
            parts = [float(x.strip()) for x in mc_val.split(",")]
            if len(parts) == 3: grid_center = tuple(parts)
        except ValueError:
            st.warning("Could not parse center coordinates.")
elif pockets_ready:
    pockets = st.session_state["s_pocket_results"]
    labels = [f"Pocket {i+1}  ({p['center'][0]:.2f}, {p['center'][1]:.2f}, {p['center'][2]:.2f})" for i, p in enumerate(pockets)]
    sel = st.selectbox("Select pocket", range(len(labels)), format_func=lambda i: labels[i], key="s_pocket_sel")
    grid_center = tuple(pockets[sel]["center"])

with st.expander("⚙️ Docking Parameters", expanded=True):
    if dock_engine == "AutoDock Vina":
        d1, d2, d3 = st.columns(3)
        with d1:
            v_gsx = st.number_input("Grid Size X", value=20, min_value=1, key="s_v_gsx", help="Search space size in X dimension (Å)")
            v_gsy = st.number_input("Grid Size Y", value=20, min_value=1, key="s_v_gsy", help="Search space size in Y dimension (Å)")
            v_gsz = st.number_input("Grid Size Z", value=20, min_value=1, key="s_v_gsz", help="Search space size in Z dimension (Å)")
        with d2:
            v_exh = st.number_input("Exhaustiveness", value=8, min_value=1, key="s_v_exh", help="Search exhaustiveness (higher = more accurate but slower)")
            v_nm  = st.number_input("Num modes", value=9, min_value=1, key="s_v_nm", help="Max number of binding modes to generate")
        with d3:
            v_seed = st.number_input("Seed (0 = auto)", value=0, min_value=0, key="s_v_seed", help="Random seed for reproducibility")
    else:
        st.markdown("##### Grid Settings")
        g1, g2, g3 = st.columns(3)
        with g1:
            a_gsx = st.number_input("Grid Size X", value=60, min_value=1, key="s_a_gsx", help="Number of grid points in X dimension")
            a_gsy = st.number_input("Grid Size Y", value=60, min_value=1, key="s_a_gsy", help="Number of grid points in Y dimension")
        with g2:
            a_gsz = st.number_input("Grid Size Z", value=60, min_value=1, key="s_a_gsz", help="Number of grid points in Z dimension")
            a_sp  = st.number_input("Spacing (Å)", value=0.375, format="%.3f", key="s_a_sp", help="Spacing between grid points in Å")
        with g3:
            a_diel   = st.number_input("Dielectric", value=-0.1465, format="%.4f", key="s_a_diel",
                                        help="Dielectric constant for electrostatics")
            a_smooth = st.number_input("Smooth", value=0.5, format="%.2f", min_value=0.0, key="s_a_smooth",
                                        help="Smoothing factor for potential maps")

        st.markdown("##### Genetic Algorithm Settings")
        ga1, ga2, ga3 = st.columns(3)
        with ga1:
            a_run   = st.number_input("GA runs", value=10, min_value=1, key="s_a_run", help="Number of GA runs (independent dockings)")
            a_pop   = st.number_input("Population size", value=150, min_value=1, key="s_a_pop", help="Number of individuals in population")
        with ga2:
            a_evals = st.number_input("Num evals", value=2500000, min_value=1, key="s_a_evals", help="Max number of energy evaluations")
            a_gens  = st.number_input("Num generations", value=27000, min_value=1, key="s_a_gens",
                                      help="Maximum number of generations")
        with ga3:
            a_elit    = st.number_input("Elitism", value=1, min_value=0, key="s_a_elit",
                                        help="Number of top individuals preserved")
            a_mut     = st.number_input("Mutation rate", value=0.02, format="%.3f", min_value=0.0, max_value=1.0,
                                        key="s_a_mut", help="Probability of mutation")
            a_cross   = st.number_input("Crossover rate", value=0.8, format="%.2f", min_value=0.0, max_value=1.0,
                                        key="s_a_cross", help="Probability of crossover")

        st.markdown("##### Other Settings")
        o1, o2 = st.columns(2)
        with o1:
            a_rmstol = st.number_input("RMSD tolerance", value=2.0, format="%.2f", min_value=0.0,
                                        key="s_a_rmstol", help="RMSD tolerance for clustering")
        with o2:
            a_seed_mode = st.selectbox("Seed", ["Auto (pid, time)", "Custom"], key="s_a_seed_mode")
            if a_seed_mode == "Custom":
                sc1, sc2 = st.columns(2)
                with sc1:
                    a_seed1 = st.number_input("Seed 1", value=0, min_value=0, key="s_a_seed1", help="First random seed")
                with sc2:
                    a_seed2 = st.number_input("Seed 2", value=0, min_value=0, key="s_a_seed2", help="Second random seed")

can_dock = protein_ready and ligand_ready and grid_center is not None

if st.button("🚀 Run Docking", type="primary", use_container_width=True, key="s_btn_dock"):
    if not can_dock:
        st.warning("Prepare both molecules and select/enter a pocket center first.")
    else:
        log_ph = st.empty()
        with st.spinner(f"Running {dock_engine}…"):
            with capture_log("s_log_docking", log_ph):
                dock_out = out_path
                try:
                    rec = str(st.session_state["s_prepared_protein"])
                    lig = str(st.session_state["s_prepared_ligand"])
                    d = None
                    if dock_engine == "AutoDock Vina":
                        from docksuitex import VinaDocking
                        d = VinaDocking(
                            receptor=rec, ligand=lig,
                            grid_center=grid_center,
                            grid_size=(v_gsx, v_gsy, v_gsz),
                            exhaustiveness=v_exh, num_modes=v_nm,
                            seed=v_seed if v_seed > 0 else None
                        )
                    else:
                        from docksuitex import AD4Docking
                        ad4_seed = ("pid", "time")
                        if a_seed_mode == "Custom":
                            ad4_seed = (a_seed1, a_seed2)
                        d = AD4Docking(
                            receptor=rec, ligand=lig,
                            grid_center=grid_center,
                            grid_size=(a_gsx, a_gsy, a_gsz),
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
                            seed=ad4_seed
                        )
                    
                    # run() returns the output directory Path
                    result_dir = d.run(save_to=str(dock_out))

                    # Parse results using the class method
                    df = d.parse_results()
                    csv_name = "vina_summary.csv" if dock_engine == "AutoDock Vina" else "ad4_summary.csv"
                    df.to_csv(dock_out / csv_name, index=False)

                    st.session_state["s_docking_results"] = dock_out
                except Exception as e:
                    st.error(str(e))
                    with st.expander("Traceback"): st.code(traceback.format_exc())
        log_ph.empty()

render_log_if_present("s_log_docking", "Docking Log:")

# ════════════════════════════════════════════════════════════════
# RESULTS VIEWER  (outside button handler so it persists)
# ════════════════════════════════════════════════════════════════
if st.session_state.get("s_docking_results"):
    results_dir = st.session_state["s_docking_results"]
    st.markdown("---")
    st.markdown("##### 📊 Docking Results")

    csv_path = results_dir / ("vina_summary.csv" if dock_engine == "AutoDock Vina" else "ad4_summary.csv")

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        st.markdown(f'<span class="status-badge">✅ {len(df)} pose(s) found</span>', unsafe_allow_html=True)

        ecol = "Affinity (kcal/mol)" if "Affinity (kcal/mol)" in df.columns else "Binding_Energy"

        if ecol in df.columns and not df[ecol].isna().all():
            mc1, mc2, mc3 = st.columns(3)
            with mc1: st.markdown(f'<div class="metric-card"><div class="metric-value">{len(df)}</div><div class="metric-label">Total Poses</div></div>', unsafe_allow_html=True)
            with mc2: st.markdown(f'<div class="metric-card"><div class="metric-value">{df[ecol].min():.2f}</div><div class="metric-label">Best Energy</div></div>', unsafe_allow_html=True)
            with mc3: st.markdown(f'<div class="metric-card"><div class="metric-value">{df[ecol].mean():.2f}</div><div class="metric-label">Mean Energy</div></div>', unsafe_allow_html=True)

        st.dataframe(df, width='stretch', hide_index=True)

    else:
        st.info("Docking completed but no results CSV found yet.")

    # ── Interactive 3D Viewer ──
    st.markdown("##### 🧬 Interactive 3D View")
    
    try:
        # 1. Get Receptor content — find by prepared protein's filename
        prep_prot = st.session_state.get("s_prepared_protein")
        rec_path = results_dir / prep_prot.name if prep_prot else None
        
        rec_content = ""
        rec_ext = "pdbqt"
        if rec_path and rec_path.exists():
            rec_content = rec_path.read_text()
            rec_ext = rec_path.suffix
        elif prep_prot and prep_prot.exists():
            rec_content = prep_prot.read_text()
            rec_ext = prep_prot.suffix
        
        # 2. Get Ligand Output content
        lig_out_path = results_dir / "output.pdbqt"
        # Some ad4 outputs might differ, checking generic first
        if not lig_out_path.exists():
            # Try finding any pdbqt that isn't the receptor
            rec_name = prep_prot.name if prep_prot else ""
            candidates = list(results_dir.glob("*.pdbqt"))
            for c in candidates:
                if c.name != rec_name and c.name != "output.pdbqt":
                    lig_out_path = c
                    break
        
        if rec_content and lig_out_path.exists():
            lig_content = lig_out_path.read_text()
            models = extract_models(lig_content)
            
            if models:
                num_models = len(models)

                # ── Pose navigation controls ──
                ctrl_cols = st.columns([0.08, 0.08, 0.84])
                with ctrl_cols[0]:
                    if st.button("◀", key="s_pose_prev", help="Previous pose"):
                        cur = st.session_state.get("s_pose_slider", 1)
                        st.session_state["s_pose_slider"] = max(1, cur - 1)
                with ctrl_cols[1]:
                    if st.button("▶", key="s_pose_next", help="Next pose"):
                        cur = st.session_state.get("s_pose_slider", 1)
                        st.session_state["s_pose_slider"] = min(num_models, cur + 1)
                with ctrl_cols[2]:
                    model_idx = st.slider("Select Pose", 1, num_models,
                                          st.session_state.get("s_pose_slider", 1),
                                          key="s_pose_slider") - 1

                selected_model = models[model_idx]
                
                # Build structures list
                view_structs = [
                    {"data": rec_content, "ext": rec_ext, "type": "protein"},
                    {"data": selected_model, "ext": "pdbqt", "type": "ligand"}
                ]
                
                st.caption(f"Showing Pose {model_idx+1} of {num_models}")
                viewer_id = f"ngl-docking-viewer-{model_idx}"

                ngl_viewer(
                    structures=view_structs,
                    height=500,
                    elem_id=viewer_id
                )


            else:
                st.warning("No models found in output PDBQT.")
        else:
            st.warning("Could not find receptor or ligand output files for visualization.")

    except Exception as e:
        st.error(f"Visualization error: {e}")
        with st.expander("Details"): st.code(traceback.format_exc())

    st.info(f"Docking complete! Results in `{results_dir}`")

st.caption("DockSuiteX © 2025")
