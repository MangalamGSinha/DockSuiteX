"""
DockSuiteX — Single Docking Page
One protein + one ligand workflow with NGL.js 3D visualization.
"""
import streamlit as st

import io, os, traceback, base64, contextlib, tempfile, zipfile, time
from pathlib import Path
import pandas as pd


# ── Helpers ────────────────────────────────────────────────────

def ngl_viewer(structures: list, height: int = 500, elem_id: str = None, grid_box: dict = None):
    """Render molecules in NGL.js embedded viewer.
    
    Args:
        structures: List of dicts, each with:
            - data: file content string
            - ext: file extension (e.g. "pdb")
            - type: "protein" or "ligand"
        height: Viewer height in pixels.
        elem_id: Optional fixed ID for the viewer container.
        grid_box: Optional dict with grid box params to draw:
            - center: [x, y, z] grid center coordinates
            - size: [sx, sy, sz] box dimensions in Å
            - color: [r, g, b] box edge color (0-1), default red
            - center_color: [r, g, b] center cross color (0-1), default green
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
        elif ext in ["cif", "mmcif"]: ext = "mmcif"
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

    # ── Grid box overlay (NGL.js Shape API — ADT-style RGB) ──
    grid_box_js = ""
    if grid_box and grid_box.get("center") and grid_box.get("size"):
        cx, cy, cz = grid_box["center"]
        gsx, gsy, gsz = grid_box["size"]

        # Half‑sizes
        hx, hy, hz = gsx / 2, gsy / 2, gsz / 2

        grid_box_js = f"""
            // ── ADT-style grid box ──
            var shape = new NGL.Shape("gridbox");
            var cx={cx}, cy={cy}, cz={cz};
            var hx={hx}, hy={hy}, hz={hz};
            var RED=[1,0,0], GREEN=[0,1,0], BLUE=[0,0,1];

            // 8 corners
            var corners = [];
            var dxs = [-hx, hx], dys = [-hy, hy], dzs = [-hz, hz];
            for (var di=0; di<2; di++)
              for (var dj=0; dj<2; dj++)
                for (var dk=0; dk<2; dk++)
                  corners.push([cx+dxs[di], cy+dys[dj], cz+dzs[dk]]);

            // Edges coloured by parallel axis
            // X-parallel → RED
            var xEdges=[[0,4],[1,5],[2,6],[3,7]];
            for(var e=0;e<xEdges.length;e++) shape.addCylinder(corners[xEdges[e][0]],corners[xEdges[e][1]],RED,0.2);
            // Y-parallel → GREEN
            var yEdges=[[0,2],[1,3],[4,6],[5,7]];
            for(var e=0;e<yEdges.length;e++) shape.addCylinder(corners[yEdges[e][0]],corners[yEdges[e][1]],GREEN,0.2);
            // Z-parallel → BLUE
            var zEdges=[[0,1],[2,3],[4,5],[6,7]];
            for(var e=0;e<zEdges.length;e++) shape.addCylinder(corners[zEdges[e][0]],corners[zEdges[e][1]],BLUE,0.2);

            // Center cross
            var cs = 2.0;
            shape.addCylinder([cx-cs,cy,cz],[cx+cs,cy,cz], RED, 0.3);
            shape.addCylinder([cx,cy-cs,cz],[cx,cy+cs,cz], GREEN, 0.3);
            shape.addCylinder([cx,cy,cz-cs],[cx,cy,cz+cs], BLUE, 0.3);

            stage.addComponentFromObject(shape).addRepresentation("buffer");

            // ── Translucent panels (ADT RGB per axis) ──
            var panelShape = new NGL.Shape("gridbox_panels");
            // 6 faces: [corner quad indices, colour]
            var faceDefs = [
              [[0,1,3,2], RED],  [[4,5,7,6], RED],     // X faces
              [[0,1,5,4], GREEN],[[2,3,7,6], GREEN],   // Y faces
              [[0,2,6,4], BLUE], [[1,3,7,5], BLUE]     // Z faces
            ];
            var pos=[], col=[], idx=[], vo=0;
            for(var f=0; f<faceDefs.length; f++) {{
              var q=faceDefs[f][0], c=faceDefs[f][1];
              for(var v=0; v<4; v++) {{
                pos.push(corners[q[v]][0], corners[q[v]][1], corners[q[v]][2]);
                col.push(c[0], c[1], c[2]);
              }}
              idx.push(vo,vo+1,vo+2, vo,vo+2,vo+3);
              vo+=4;
            }}
            panelShape.addMesh(
              new Float32Array(pos), new Float32Array(col), new Uint32Array(idx)
            );
            stage.addComponentFromObject(panelShape).addRepresentation("buffer", {{opacity: 0.15, side: "double"}});
        """

    html_template = """
    <style>
        body { margin: 0; padding: 0; overflow: hidden; background-color: transparent; }
    </style>
    <div id="__VIEWPORT_ID__" style="width:100%;height:__HEIGHT__px;border-radius:12px;border:1px solid rgba(128,128,128,0.15);"></div>

    <script>
    (function() {
        // Detect browser light/dark mode from parent window
        function getThemeColors() {
            var isDark = window.parent.matchMedia('(prefers-color-scheme: dark)').matches;
            return {
                bg: isDark ? '#0e1117' : '#f8fafc',
                border: isDark ? 'rgba(255,255,255,0.08)' : 'rgba(0,0,0,0.1)'
            };
        }

        function initNGL() {
            var colors = getThemeColors();
            var container = document.getElementById("__VIEWPORT_ID__");
            container.style.background = colors.bg;
            container.style.borderColor = colors.border;

            var stage = new NGL.Stage("__VIEWPORT_ID__", {backgroundColor: colors.bg});
            var promises = [];

            __LOAD_BLOCK__

            Promise.all(promises).then(function() {
                __GRID_BOX_BLOCK__
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

            // Listen for theme changes and update background
            window.parent.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', function() {
                var newColors = getThemeColors();
                container.style.background = newColors.bg;
                container.style.borderColor = newColors.border;
                stage.setParameters({backgroundColor: newColors.bg});
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
        .replace("__GRID_BOX_BLOCK__", grid_box_js)
    )




    st.iframe(html, height=height + 10)




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
            "s_interaction_profile",
        ]:
            if key not in st.session_state:
                st.session_state[key] = None

        for key in ["s_log_ligand", "s_log_protein", "s_log_docking", "s_log_pockets", "s_log_prolif"]:
            if key not in st.session_state:
                st.session_state[key] = ""

        if "s_grid_box_params" not in st.session_state:
            st.session_state["s_grid_box_params"] = None
        
        if "s_pose_slider" not in st.session_state:
            st.session_state["s_pose_slider"] = 1
            
    except Exception:
        # Ignore if running in a context without session state (e.g. initial import)
        pass

init_session_state()


# ── Hero header ────────────────────────────────────────────────
st.markdown(
    """
    <div class="hero-header animate-in">
        <h1>🔬 Single Docking</h1>
        <p>Prepare protein &amp; ligand, find pockets, dock, and profile
        results with interactive molecular views.</p>
    </div>
    """,
    unsafe_allow_html=True,
)


# ── Temp output directory ──────────────────────────────────────
if "s_temp_dir_obj" not in st.session_state:
    st.session_state["s_temp_dir_obj"] = tempfile.TemporaryDirectory(prefix="docksuitex_single_")
    st.session_state["s_temp_dir"] = st.session_state["s_temp_dir_obj"].name
out_path = Path(st.session_state["s_temp_dir"])


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
                key=f"s_dl_{rel_key}",
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
            key="s_zip_dl_all",
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
    st.markdown('<div class="section-header"><h2>🧬 Protein Preparation</h2></div>', unsafe_allow_html=True)
    prot_file = st.file_uploader(
        "Upload protein", 
        type=["pdb"], 
        key="s_prot_uploader",
        help="Please upload a file containing only a single model. If you have a multi-model file, use the Format Converter to split it first."
    )

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
            prot_fix   = st.checkbox("Fix PDB(PDBFixer)", True, key="s_prot_fix")
            prot_het   = st.checkbox("Remove heterogens", True, key="s_prot_het")
            prot_water = st.checkbox("Remove water", True, key="s_prot_water")
            prot_chg = st.checkbox("Add Gasteiger charges", True, key="s_prot_chg")
            prot_h   = st.checkbox("Add hydrogens", True, key="s_prot_h")
        with p2:
            prot_ph = st.number_input("pH", value=7.4, min_value=0.0, max_value=14.0,
                                       step=0.1, format="%.1f", key="s_prot_ph",
                                       help="pH for protonation state assignment via PDBFixer")
            prot_preserve = st.text_input("Preserve charge types", placeholder="e.g. Zn, Fe", key="s_prot_pres")
    prot_preserve_list = parse_comma_list(prot_preserve)

    if st.button("🚀 Prepare Protein", type="primary", width='stretch', key="s_btn_prot"):
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
                                       remove_water=prot_water, add_hydrogens=prot_h, ph=prot_ph,
                                       add_charges=prot_chg,
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
                    structures=[{"data": prepared_prot_path.read_text(errors="replace"), "ext": prepared_prot_path.suffix, "type": "protein"}],
                    height=300,
                    elem_id="ngl-protein-prepared"
                )
            except Exception as e:
                st.error("Viewer Error:")
                st.code(traceback.format_exc())


# ── Ligand Preparation ────────────────────────────────────────
with col_lig:
    st.markdown('<div class="section-header"><h2>💊 Ligand Preparation</h2></div>', unsafe_allow_html=True)
    lig_file = st.file_uploader(
        "Upload ligand", 
        type=["sdf", "mol2"], 
        key="s_lig_uploader",
        help="Please upload a file containing only a single molecule. If you have a multi-molecule library, use the Format Converter to split it first."
    )

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
                                    format_func=lambda x: "None (skip)" if x is None else x.upper(), 
                                    help="Select a forcefield for energy minimization (2500 steps) via Open Babel.",
                                    key="s_lig_min")
        l1, l2 = st.columns(2)
        with l1:
            lig_rw  = st.checkbox("Remove water", True, key="s_lig_rw")
            lig_h   = st.checkbox("Add hydrogens", True, key="s_lig_h")
            lig_chg = st.checkbox("Add Gasteiger charges", True, key="s_lig_chg")
        with l2:
            lig_ph = st.number_input("pH", value=7.4, min_value=0.0, max_value=14.0,
                                      step=0.1, format="%.1f", key="s_lig_ph",
                                      help="pH for protonation state assignment via Open Babel")
            lig_preserve = st.text_input("Preserve charge types", placeholder="e.g. Zn, Fe", key="s_lig_pres")
    lig_preserve_list = parse_comma_list(lig_preserve)

    if st.button("🚀 Prepare Ligand", type="primary", width='stretch', key="s_btn_lig"):
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
                                     add_hydrogens=lig_h, ph=lig_ph, add_charges=lig_chg,
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
# GRID CALCULATOR
# ════════════════════════════════════════════════════════════════
st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)
st.markdown('<div class="section-header"><h2>📐 Grid Calculator</h2></div>', unsafe_allow_html=True)
st.caption("Compute grid box(es) for docking — predict pockets with P2Rank or auto-compute a blind box from the receptor.")

gc_mode = st.radio(
    "Mode",
    ["P2Rank (pocket prediction)", "Blind Docking (full receptor)"],
    horizontal=True,
    key="s_gc_mode",
)
is_gc_p2rank = gc_mode == "P2Rank (pocket prediction)"

can_find = st.session_state.get("s_prepared_protein") is not None
if not can_find:
    st.info("Prepare a protein above first.")


btn_label = "🔍 Run P2Rank" if is_gc_p2rank else "📦 Compute Blind Box"
if st.button(btn_label, type="primary", width='stretch', key="s_btn_pf"):
    if not can_find:
        st.warning("Prepare a protein above first.")
    else:
        log_ph = st.empty()
        spinner_msg = "Running P2Rank…" if is_gc_p2rank else "Computing blind docking box…"
        with st.spinner(spinner_msg):
            with capture_log("s_log_pockets", log_ph):
                try:
                    from docksuitex.grid_calculator import GridCalculator
                    prep_prot = st.session_state["s_prepared_protein"]

                    if is_gc_p2rank:
                        gc = GridCalculator(receptor=str(prep_prot), mode="p2rank")
                        pockets = gc.run(
                            save_to=out_path / f"p2rank_results_{prep_prot.name.replace('.', '_')}"
                        )
                        st.session_state["s_pocket_results"] = pockets
                    else:
                        gc = GridCalculator(receptor=str(prep_prot), mode="blind")
                        pockets = gc.run()
                        st.session_state["s_pocket_results"] = pockets
                except Exception as e:
                    st.error(str(e))
                    with st.expander("Traceback"): st.code(traceback.format_exc())
        log_ph.empty()

render_log_if_present("s_log_pockets", "Grid Calculator Log:")




# ════════════════════════════════════════════════════════════════
# DOCKING
# ════════════════════════════════════════════════════════════════
st.markdown('<div class="styled-divider"></div>', unsafe_allow_html=True)
st.markdown('<div class="section-header"><h2>⚡ Docking</h2></div>', unsafe_allow_html=True)

dock_engine = st.radio("Docking Engine", ["AutoDock Vina", "AutoDock4"], horizontal=True, key="s_dock_eng")
protein_ready  = st.session_state.get("s_prepared_protein") is not None
ligand_ready   = st.session_state.get("s_prepared_ligand") is not None
pockets_ready  = st.session_state.get("s_pocket_results") is not None

if not (protein_ready and ligand_ready):
    st.info("Prepare both protein and ligand above to enable docking.")
elif not pockets_ready:
    st.info("Run the Grid Calculator above, or choose Manual Entry below.")

# ── Session State Initialization ──────────────────────────────
# We initialize these so widgets with 'key' don't conflict with 'value'
st.session_state.setdefault("s_cx", 0.0)
st.session_state.setdefault("s_cy", 0.0)
st.session_state.setdefault("s_cz", 0.0)
st.session_state.setdefault("s_v_gsx", 20.0)
st.session_state.setdefault("s_v_gsy", 20.0)
st.session_state.setdefault("s_v_gsz", 20.0)
st.session_state.setdefault("s_a_gsx", 22.5)
st.session_state.setdefault("s_a_gsy", 22.5)
st.session_state.setdefault("s_a_gsz", 22.5)

# ── helper: switch to Manual Entry when user edits center inputs ──
def _switch_to_manual():
    st.session_state["s_center_mode"] = "Manual Entry"

# ── Build the unified grid options list ──
# Each entry: {"label": str, "center": (x,y,z), "grid_size": (sx,sy,sz)|None}
_grid_options = []

if pockets_ready:
    for i, p in enumerate(st.session_state["s_pocket_results"]):
        cx, cy, cz = p["center"]
        gs = p.get("grid_size")
        gs_str = f"({gs[0]:.3f}, {gs[1]:.3f}, {gs[2]:.3f})" if gs else "N/A"
        prob = p.get('probability', 0)
        _grid_options.append({
            "label": f"Rank: {p.get('rank', i+1)}, Probability: {prob:.3f}, Center: ({cx:.4f}, {cy:.4f}, {cz:.4f}), Grid Size: {gs_str}",
            "center": p["center"],
            "grid_size": gs,
        })

# ── Grid Mode selector ──
center_mode = st.radio(
    "Grid Mode",
    ["Grid Selection", "Manual Entry"],
    horizontal=True,
    key="s_center_mode",
    help=(
        "**Grid Selection** – choose from P2Rank pockets or a pre-computed blind box.  "
        "**Manual Entry** – type in custom X / Y / Z coordinates directly."
    ),
)

# ── Resolve selected grid option and update state ──
if center_mode == "Grid Selection":
    if _grid_options:
        _opt_idx = st.selectbox(
            "Select grid",
            range(len(_grid_options)),
            format_func=lambda i: _grid_options[i]["label"],
            key="s_grid_sel",
        )
        _selected_opt = _grid_options[_opt_idx]
        
        # Update center
        _new_c = _selected_opt["center"]
        if (round(st.session_state["s_cx"], 4) != round(float(_new_c[0]), 4) or
            round(st.session_state["s_cy"], 4) != round(float(_new_c[1]), 4) or
            round(st.session_state["s_cz"], 4) != round(float(_new_c[2]), 4)):
            st.session_state["s_cx"] = float(_new_c[0])
            st.session_state["s_cy"] = float(_new_c[1])
            st.session_state["s_cz"] = float(_new_c[2])
            
        # Update sizes if available
        _gs = _selected_opt.get("grid_size")
        if _gs:
            if dock_engine == "AutoDock Vina":
                st.session_state["s_v_gsx"] = float(_gs[0])
                st.session_state["s_v_gsy"] = float(_gs[1])
                st.session_state["s_v_gsz"] = float(_gs[2])
            else:
                st.session_state["s_a_gsx"] = float(_gs[0])
                st.session_state["s_a_gsy"] = float(_gs[1])
                st.session_state["s_a_gsz"] = float(_gs[2])
    else:
        st.info("Run the Grid Calculator above to populate this list.")

# ── Grid Parameters ──
with st.expander("📐 Grid Parameters", expanded=True):

    # ── Grid Center (always shown as 3 inputs) ──
    # Note: 'value' is omitted because we use 'key' and set session_state manually above
    st.markdown("###### Grid Center (Å)")
    cc1, cc2, cc3 = st.columns(3)
    with cc1:
        cx = st.number_input("Center X", format="%.4f", step=0.1,
                             key="s_cx", on_change=_switch_to_manual,
                             help="Grid box center X coordinate (Å)")
    with cc2:
        cy = st.number_input("Center Y", format="%.4f", step=0.1,
                             key="s_cy", on_change=_switch_to_manual,
                             help="Grid box center Y coordinate (Å)")
    with cc3:
        cz = st.number_input("Center Z", format="%.4f", step=0.1,
                             key="s_cz", on_change=_switch_to_manual,
                             help="Grid box center Z coordinate (Å)")

    grid_center = (cx, cy, cz)

    st.markdown("")  # spacer

    # ── Grid Size ──
    # Note: 'value' is omitted because we use 'key' and initialize session_state
    if dock_engine == "AutoDock Vina":
        gc1, gc2, gc3 = st.columns(3)
        with gc1:
            v_gsx = st.number_input("Grid Size X (Å)", min_value=1.0, step=1.0,
                                    format="%.3f", key="s_v_gsx",
                                    help="Search space size in X dimension (Å)")
        with gc2:
            v_gsy = st.number_input("Grid Size Y (Å)", min_value=1.0, step=1.0,
                                    format="%.3f", key="s_v_gsy",
                                    help="Search space size in Y dimension (Å)")
        with gc3:
            v_gsz = st.number_input("Grid Size Z (Å)", min_value=1.0, step=1.0,
                                    format="%.3f", key="s_v_gsz",
                                    help="Search space size in Z dimension (Å)")
    else:
        gc1, gc2, gc3 = st.columns(3)
        with gc1:
            a_gsx = st.number_input("Grid Size X (Å)", min_value=1.0, step=1.0,
                                    format="%.3f", key="s_a_gsx",
                                    help="Search space size in X (Å) — npts = size ÷ spacing")
        with gc2:
            a_gsy = st.number_input("Grid Size Y (Å)", min_value=1.0, step=1.0,
                                    format="%.3f", key="s_a_gsy",
                                    help="Search space size in Y dimension (Å)")
        with gc3:
            a_gsz = st.number_input("Grid Size Z (Å)", min_value=1.0, step=1.0,
                                    format="%.3f", key="s_a_gsz",
                                    help="Search space size in Z dimension (Å)")

        gc4, gc5, gc6 = st.columns(3)
        with gc4:
            a_sp = st.number_input("Spacing (Å)", value=0.375, format="%.3f", key="s_a_sp",
                                   help="Spacing between grid points in Å")
        with gc5:
            a_diel = st.number_input("Dielectric", value=-0.1465, format="%.4f", key="s_a_diel",
                                     help="Dielectric constant for electrostatics")
        with gc6:
            a_smooth = st.number_input("Smooth", value=0.5, format="%.2f", min_value=0.0,
                                       key="s_a_smooth", help="Smoothing factor for potential maps")

can_dock = protein_ready and ligand_ready and grid_center is not None

# ── Grid Box Preview ──
if protein_ready and grid_center is not None:
    prep_prot_path = st.session_state["s_prepared_protein"]
    preview_center    = list(grid_center)
    if dock_engine == "AutoDock Vina":
        preview_grid_size = [v_gsx, v_gsy, v_gsz]
    else:
        preview_grid_size = [a_gsx, a_gsy, a_gsz]

    if preview_center and preview_grid_size:
        st.markdown("##### 📦 Grid Box Preview")
        st.caption(
            f"Center: ({preview_center[0]:.2f}, {preview_center[1]:.2f}, {preview_center[2]:.2f})  |  "
            f"Size: {preview_grid_size[0]:.1f} × {preview_grid_size[1]:.1f} × {preview_grid_size[2]:.1f} Å"
        )
        try:
            if prep_prot_path and prep_prot_path.exists():
                ngl_viewer(
                    structures=[{"data": prep_prot_path.read_text(errors="replace"), "ext": prep_prot_path.suffix, "type": "protein"}],
                    height=400,
                    elem_id="ngl-gridbox-preview",
                    grid_box={"center": preview_center, "size": preview_grid_size},
                )
        except Exception as e:
            st.error(f"Grid box preview error: {e}")

# ── Docking Parameters ──
with st.expander("⚙️ Docking Parameters", expanded=True):
    if dock_engine == "AutoDock Vina":
        dp1, dp2, dp3 = st.columns(3)
        with dp1:
            v_exh = st.number_input("Exhaustiveness", value=8, min_value=1, key="s_v_exh", help="Search exhaustiveness (higher = more accurate but slower)")
        with dp2:
            v_nm  = st.number_input("Num modes", value=9, min_value=1, key="s_v_nm", help="Max number of binding modes to generate")
        with dp3:
            v_seed = st.number_input("Seed (0 = auto)", value=42, min_value=0, key="s_v_seed", help="Random seed for reproducibility")
    else:
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
            a_elit  = st.number_input("Elitism", value=1, min_value=0, key="s_a_elit",
                                       help="Number of top individuals preserved")
            a_mut   = st.number_input("Mutation rate", value=0.02, format="%.3f", min_value=0.0, max_value=1.0,
                                       key="s_a_mut", help="Probability of mutation")
            a_cross = st.number_input("Crossover rate", value=0.8, format="%.2f", min_value=0.0, max_value=1.0,
                                       key="s_a_cross", help="Probability of crossover")

        st.markdown("##### Other Settings")
        o1, o2 = st.columns(2)
        with o1:
            a_rmstol = st.number_input("RMSD tolerance", value=2.0, format="%.2f", min_value=0.0,
                                        key="s_a_rmstol", help="RMSD tolerance for clustering")
        with o2:
            a_seed_mode = st.selectbox("Seed", ["Custom", "Auto (pid, time)"], key="s_a_seed_mode")
            if a_seed_mode == "Custom":
                sc1, sc2 = st.columns(2)
                with sc1:
                    a_seed1 = st.number_input("Seed 1", value=27, min_value=0, key="s_a_seed1", help="First random seed")
                with sc2:
                    a_seed2 = st.number_input("Seed 2", value=6, min_value=0, key="s_a_seed2", help="Second random seed")

if st.button("🚀 Run Docking", type="primary", width='stretch', key="s_btn_dock"):
    if not can_dock:
        st.warning("Prepare both molecules and select/enter a pocket center first.")
    else:
        log_ph = st.empty()
        with st.spinner(f"Running {dock_engine} docking…"):
            with capture_log("s_log_docking", log_ph):
                dock_out = out_path
                try:
                    rec = str(st.session_state["s_prepared_protein"])
                    lig = str(st.session_state["s_prepared_ligand"])
                    d = None
                    if dock_engine == "AutoDock Vina":
                        from docksuitex.vina import VinaDocking
                        d = VinaDocking(
                            receptor=rec, ligand=lig,
                            grid_center=grid_center,
                            grid_size=(v_gsx, v_gsy, v_gsz),
                            exhaustiveness=v_exh, num_modes=v_nm,
                            seed=v_seed if v_seed > 0 else None
                        )
                        # run() returns the output directory Path
                        d.run(save_to=dock_out / "vina_results")
                        # Parse results using the class method
                        df = d.parse_results()
                        st.session_state["s_docking_results"] = dock_out / "vina_results"

                    else:
                        from docksuitex.autodock4 import AD4Docking
                        ad4_seed = (27, 6)
                        if a_seed_mode == "Custom":
                            ad4_seed = (a_seed1, a_seed2)
                        elif a_seed_mode == "Auto (pid, time)":
                            ad4_seed = ("pid", "time")
                        d = AD4Docking(
                            receptor=rec, ligand=lig,
                            grid_center=grid_center,
                            grid_size=(a_gsx, a_gsy, a_gsz),  # Å — npts computed internally
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
                        d.run(save_to=dock_out / "ad4_results")
                        # Parse results using the class method
                        df = d.parse_results()
                        st.session_state["s_docking_results"] = dock_out / "ad4_results"
                    



                    st.session_state["s_pose_slider"] = 1  # Reset pose slider

                    # Store grid box params for results viewer
                    # grid_size is always in Å for both engines
                    gb_center = list(d.grid_center)
                    gb_size = list(d.grid_size)
                    st.session_state["s_grid_box_params"] = {"center": gb_center, "size": gb_size}
                except Exception as e:
                    st.error(str(e))
                    with st.expander("Traceback"): st.code(traceback.format_exc())
        log_ph.empty()

        # Run interaction profiler automatically after docking
        if st.session_state.get("s_docking_results") and d is not None:
            try:
                prolif_log_ph = st.empty()
                with st.spinner("Running interaction profiler…"):
                    with capture_log("s_log_prolif", prolif_log_ph):
                        ip_df = d.interaction_profile()
                        st.session_state["s_interaction_profile"] = ip_df
                prolif_log_ph.empty()
                # Append ProLIF log to Docking Log
                prolif_log = st.session_state.get("s_log_prolif", "")
                if prolif_log:
                    st.session_state["s_log_docking"] = st.session_state.get("s_log_docking", "") + "\n" + prolif_log
            except Exception as e:
                st.warning(f"Interaction profiler failed: {e}")
                with st.expander("Profiler Traceback"): st.code(traceback.format_exc())

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
        st.markdown(f'<span class="status-badge">✅ {len(df)} result(s) found</span>', unsafe_allow_html=True)

        ecol = "Affinity (kcal/mol)" if "Affinity (kcal/mol)" in df.columns else "Binding_Energy"

        if ecol in df.columns and not df[ecol].isna().all():
            mc1, mc2, mc3 = st.columns(3)
            with mc1: st.markdown(f'<div class="metric-card"><div class="metric-value">{len(df)}</div><div class="metric-label">Total Results</div></div>', unsafe_allow_html=True)
            with mc2: st.markdown(f'<div class="metric-card"><div class="metric-value">{df[ecol].min():.2f}</div><div class="metric-label">Best Energy</div></div>', unsafe_allow_html=True)
            with mc3: st.markdown(f'<div class="metric-card"><div class="metric-value">{df[ecol].mean():.2f}</div><div class="metric-label">Mean Energy</div></div>', unsafe_allow_html=True)

        st.dataframe(df, width='stretch', hide_index=True)

        # ── Interactive 3D Viewer ──
        st.markdown("##### 🧬 Interactive 3D View")

        try:
            # 1. Get Receptor content — find by prepared protein's filename
            prep_prot = st.session_state.get("s_prepared_protein")
            rec_path = results_dir / prep_prot.name if prep_prot else None

            rec_content = ""
            rec_ext = "pdbqt"
            if rec_path and rec_path.exists():
                rec_content = rec_path.read_text(errors="replace")
                rec_ext = rec_path.suffix
            elif prep_prot and prep_prot.exists():
                rec_content = prep_prot.read_text(errors="replace")
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
                lig_content = lig_out_path.read_text(errors="replace")
                models = extract_models(lig_content)

                if models:
                    num_models = len(models)

                    # ── Pose mode: Show All vs Single ──
                    st.session_state.setdefault("s_show_all_poses", False)
                    st.session_state.setdefault("s_play_poses", False)
                    st.session_state.setdefault("s_play_speed", 1.0)

                    # Ensure pose index is within bounds
                    if st.session_state["s_pose_slider"] > num_models:
                        st.session_state["s_pose_slider"] = 1

                    # ── Controls row ──
                    ctrl_cols = st.columns([0.07, 0.07, 0.07, 0.79])
                    with ctrl_cols[0]:
                        if st.button("◀", key="s_pose_prev", help="Previous pose"):
                            st.session_state["s_play_poses"] = False
                            st.session_state["s_show_all_poses"] = False
                            st.session_state["s_pose_slider"] = max(1, st.session_state["s_pose_slider"] - 1)
                    with ctrl_cols[1]:
                        if st.button("⏯️", key="s_play_btn", help="Play / Pause animation"):
                            st.session_state["s_play_poses"] = not st.session_state.get("s_play_poses", False)
                            if st.session_state["s_play_poses"]:
                                st.session_state["s_show_all_poses"] = False
                    with ctrl_cols[2]:
                        if st.button("▶", key="s_pose_next", help="Next pose"):
                            st.session_state["s_play_poses"] = False
                            st.session_state["s_show_all_poses"] = False
                            st.session_state["s_pose_slider"] = min(num_models, st.session_state["s_pose_slider"] + 1)
                    with ctrl_cols[3]:
                        st.session_state["s_play_speed"] = st.slider(
                            "Speed", 0.5, 5.0, st.session_state.get("s_play_speed", 1.0),
                            step=0.5, key="s_speed_slider", help="Animation speed (poses/sec)"
                        )

                    show_all = st.toggle("Show All Poses", key="s_show_all_poses", help="Show all poses simultaneously")
                    if show_all:
                        st.session_state["s_play_poses"] = False

                    # ── Auto-play logic ──
                    if st.session_state.get("s_play_poses") and not show_all:
                        next_pose = st.session_state["s_pose_slider"] + 1
                        if next_pose > num_models:
                            next_pose = 1
                        st.session_state["s_pose_slider"] = next_pose

                    model_idx = st.session_state["s_pose_slider"] - 1

                    # ── Build structures list ──
                    view_structs = [
                        {"data": rec_content, "ext": rec_ext, "type": "protein"},
                    ]

                    if show_all:
                        # Add all poses as separate ligand entries
                        for i, m in enumerate(models):
                            view_structs.append({"data": m, "ext": "pdbqt", "type": "ligand"})
                        st.caption(f"Showing all {num_models} pose(s)")
                        viewer_id = "ngl-docking-viewer-all"
                    else:
                        selected_model = models[model_idx]
                        view_structs.append({"data": selected_model, "ext": "pdbqt", "type": "ligand"})
                        st.caption(f"Showing Pose {model_idx+1} of {num_models}")
                        viewer_id = f"ngl-docking-viewer-{model_idx}"

                    # ── Grid box toggle ──
                    show_box = st.checkbox("📦 Show Grid Box", value=False, key="s_show_gridbox",
                                           help="Toggle the docking grid box overlay")

                    active_grid_box = st.session_state.get("s_grid_box_params") if show_box else None

                    ngl_viewer(
                        structures=view_structs,
                        height=500,
                        elem_id=viewer_id,
                        grid_box=active_grid_box,
                    )

                    # ── Trigger rerun for animation after viewer is rendered ──
                    if st.session_state.get("s_play_poses") and not show_all:
                        speed = st.session_state.get("s_play_speed", 1.0)
                        time.sleep(1.0 / speed)
                        st.rerun()


                else:
                    st.warning("No models found in output PDBQT.")
            else:
                st.warning("Could not find receptor or ligand output files for visualization.")

        except Exception as e:
            st.error(f"Visualization error: {e}")
            with st.expander("Details"): st.code(traceback.format_exc())

        # ── Interaction Profile ──
        if st.session_state.get("s_interaction_profile") is not None:
            st.markdown("##### 🧪 Interaction Profile")
            ip_display = st.session_state["s_interaction_profile"].reset_index()
            st.dataframe(ip_display, width='stretch', hide_index=True)

            # ── ProLIF Visualizations (2D / 3D / Barcode) ──
            prolif_dir = results_dir / "prolif_results"
            if prolif_dir.exists():
                st.markdown("##### 🔬 Interaction Visualizations")

                viz_tab = st.radio(
                    "View", ["2D Network", "3D Complex", "Barcode"],
                    horizontal=True, key="s_prolif_viz_mode",
                    label_visibility="collapsed",
                )

                if viz_tab == "Barcode":
                    barcode_path = prolif_dir / "prolif_barcode.png"
                    if barcode_path.exists():
                        st.image(str(barcode_path), caption="Interaction Barcode (all poses)", width='stretch')
                    else:
                        st.info("Barcode image not generated.")

                else:
                    # Count available poses
                    viz_suffix = "2D" if viz_tab == "2D Network" else "3D"
                    viz_files = sorted(prolif_dir.glob(f"prolif_pose_*_interactions_{viz_suffix}.html"))
                    n_viz = len(viz_files)

                    if n_viz == 0:
                        st.info(f"No {viz_tab} visualizations found.")
                    else:
                        # Pose selector
                        st.session_state.setdefault("s_prolif_pose", 1)
                        if st.session_state["s_prolif_pose"] > n_viz:
                            st.session_state["s_prolif_pose"] = 1

                        pc1, pc2, pc3, pc4, pc5 = st.columns([3, 1, 1.5, 1, 3])
                        with pc2:
                            if st.button("◀ Prev", key=f"s_prolif_prev_{viz_suffix}", width='stretch'):
                                st.session_state["s_prolif_pose"] = max(1, st.session_state["s_prolif_pose"] - 1)
                        with pc3:
                            st.markdown(f"<div style='display: flex; align-items: center; justify-content: center; height: 38px;'><span style='font-size: 0.9rem; opacity: 0.8;'>Pose {st.session_state['s_prolif_pose']} of {n_viz}</span></div>", unsafe_allow_html=True)
                        with pc4:
                            if st.button("Next ▶", key=f"s_prolif_next_{viz_suffix}", width='stretch'):
                                st.session_state["s_prolif_pose"] = min(n_viz, st.session_state["s_prolif_pose"] + 1)

                        pose_idx = st.session_state["s_prolif_pose"]
                        html_file = prolif_dir / f"prolif_pose_{pose_idx}_interactions_{viz_suffix}.html"

                        if html_file.exists():
                            html_content = html_file.read_text(encoding="utf-8", errors="replace")
                            iframe_height = 650 if viz_suffix == "3D" else 600
                            _spacer_l, _viz_center, _spacer_r = st.columns([1, 6, 1])
                            with _viz_center:
                                st.iframe(html_content, height=iframe_height)
                        else:
                            st.warning(f"File not found: {html_file.name}")

    else:
        st.info("Docking completed but no results CSV found yet.")




render_files_panel(out_path)

st.caption("DockSuiteX © 2026")
