"""
DockSuiteX — Fetch Molecules Page
Download PDB files from RCSB and SDF files from PubChem.
"""
import streamlit as st
import io, traceback, contextlib
from pathlib import Path


# ── Helpers ────────────────────────────────────────────────────

class StreamlitSink:
    """Redirects writes to a Streamlit placeholder for real-time logging."""
    def __init__(self, placeholder):
        self.placeholder = placeholder
        self.buffer = io.StringIO()

    def write(self, text):
        self.buffer.write(text)
        if self.placeholder:
            self.placeholder.code(self.buffer.getvalue(), language="text")

    def flush(self):
        self.buffer.flush()


@contextlib.contextmanager
def capture_log(placeholder=None):
    import sys
    sink = StreamlitSink(placeholder)
    old_stdout = sys.stdout
    sys.stdout = sink
    try:
        yield sink
    finally:
        sys.stdout = old_stdout


# ── Session state init ─────────────────────────────────────────
for key in ["fetch_pdb_log", "fetch_sdf_log", "fetch_pdb_results", "fetch_sdf_results"]:
    if key not in st.session_state:
        st.session_state[key] = ""

# ── Hero header ────────────────────────────────────────────────
st.markdown(
    """
    <div class="hero-header animate-in">
        <h1>🌐 Fetch Molecules</h1>
        <p>Download protein structures (PDB) from RCSB and
        ligand structures (SDF) from PubChem.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

# ── Output directory ───────────────────────────────────────────
if "fetch_out_dir" not in st.session_state:
    st.session_state["fetch_out_dir"] = str(Path.home() / "Desktop" / "docksuitex_fetched")

st.text_input(
    "📂 Save directory",
    key="fetch_out_dir",
    help="Directory where downloaded files will be saved.",
)

st.markdown("---")

# ── Layout: Two columns ───────────────────────────────────────
col_pdb, col_sdf = st.columns(2)

# ════════════════════════════════════════════════════════════════
# FETCH PDB (from RCSB)
# ════════════════════════════════════════════════════════════════
with col_pdb:
    st.markdown('<div class="section-header"><h2>🧬 Fetch PDB</h2></div>', unsafe_allow_html=True)
    st.caption("Download protein structures from the **RCSB Protein Data Bank**.")

    pdb_input = st.text_input(
        "PDB IDs (comma-separated)",
        placeholder="e.g. 1HVR, 4LGS, 2HBA",
        key="pdb_id_input",
        help="Enter one or more 4-character PDB IDs separated by commas.",
    )

    if st.button("🚀 Fetch PDB", type="primary", key="btn_fetch_pdb"):
        if not pdb_input.strip():
            st.warning("Please enter at least one PDB ID.")
        else:
            ids = [x.strip() for x in pdb_input.split(",") if x.strip()]
            save_dir = Path(st.session_state["fetch_out_dir"])
            save_dir.mkdir(parents=True, exist_ok=True)

            log_ph = st.empty()
            try:
                with capture_log(log_ph) as sink:
                    from docksuitex.utils.fetcher import _download_pdb
                    results = []
                    for pid in ids:
                        results.append(_download_pdb(pid, save_to=str(save_dir)))

                st.session_state["fetch_pdb_log"] = sink.buffer.getvalue()
                st.session_state["fetch_pdb_results"] = [str(r) for r in results]
                st.rerun()
            except Exception as e:
                st.error(f"❌ {e}")
                with st.expander("Traceback"):
                    st.code(traceback.format_exc())

    # Persistent log & results
    if st.session_state.get("fetch_pdb_log"):
        st.caption("Download Log:")
        st.code(st.session_state["fetch_pdb_log"], language="text")

    if st.session_state.get("fetch_pdb_results"):
        paths = st.session_state["fetch_pdb_results"]


# ════════════════════════════════════════════════════════════════
# FETCH SDF (from PubChem)
# ════════════════════════════════════════════════════════════════
with col_sdf:
    st.markdown('<div class="section-header"><h2>💊 Fetch SDF</h2></div>', unsafe_allow_html=True)
    st.caption("Download 3D ligand structures from **PubChem**.")

    sdf_input = st.text_input(
        "PubChem CIDs (comma-separated)",
        placeholder="e.g. 2244, 5988, 3672",
        key="sdf_cid_input",
        help="Enter one or more numeric PubChem Compound IDs separated by commas.",
    )

    if st.button("🚀 Fetch SDF", type="primary", key="btn_fetch_sdf"):
        if not sdf_input.strip():
            st.warning("Please enter at least one CID.")
        else:
            cids = [x.strip() for x in sdf_input.split(",") if x.strip()]
            save_dir = Path(st.session_state["fetch_out_dir"])
            save_dir.mkdir(parents=True, exist_ok=True)

            log_ph = st.empty()
            try:
                with capture_log(log_ph) as sink:
                    from docksuitex.utils.fetcher import _download_sdf
                    results = []
                    for cid in cids:
                        results.append(_download_sdf(cid, save_to=str(save_dir)))

                st.session_state["fetch_sdf_log"] = sink.buffer.getvalue()
                st.session_state["fetch_sdf_results"] = [str(r) for r in results]
                st.rerun()
            except Exception as e:
                st.error(f"❌ {e}")
                with st.expander("Traceback"):
                    st.code(traceback.format_exc())

    # Persistent log & results
    if st.session_state.get("fetch_sdf_log"):
        st.caption("Download Log:")
        st.code(st.session_state["fetch_sdf_log"], language="text")

    if st.session_state.get("fetch_sdf_results"):
        paths = st.session_state["fetch_sdf_results"]


# ── Footer ─────────────────────────────────────────────────────
st.caption("DockSuiteX © 2026")
