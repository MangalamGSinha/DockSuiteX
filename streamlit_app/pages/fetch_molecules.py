"""
DockSuiteX — Fetch Molecules Page
Download protein structures (PDB) from RCSB and ligand structures (SDF)
from PubChem or ChEMBL.
Files are fetched in-memory and delivered directly to the browser via download buttons.
"""
import streamlit as st
import requests
import traceback


# ── Session state init ─────────────────────────────────────────
for key in ["pdb_downloads", "sdf_downloads"]:
    if key not in st.session_state:
        st.session_state[key] = []  # list of (filename, bytes) tuples


# ── In-memory fetch helpers ────────────────────────────────────

def _fetch_pdb_bytes(pid: str) -> tuple[str, bytes]:
    """Fetch a PDB file from RCSB and return (filename, content_bytes)."""
    pid = pid.upper().strip()
    if len(pid) != 4 or not pid.isalnum():
        raise ValueError(f"Invalid PDB ID '{pid}'. Must be 4-character alphanumeric.")

    url = f"https://files.rcsb.org/download/{pid}.pdb"
    resp = requests.get(url)
    if resp.status_code != 200:
        raise RuntimeError(f"Failed to download PDB from: {url}")

    return f"{pid}.pdb", resp.content


def _fetch_sdf_bytes(molecule_id: str) -> tuple[str, bytes]:
    """
    Fetch an SDF file from PubChem or ChEMBL and return (filename, content_bytes).

    Auto-detects the source:
    - Numeric IDs → PubChem (3D conformer)
    - IDs starting with 'CHEMBL' → ChEMBL
    """
    molecule_id = str(molecule_id).strip()

    # ── ChEMBL ─────────────────────────────────────────────
    if molecule_id.upper().startswith("CHEMBL"):
        chembl_id = molecule_id.upper()
        if not chembl_id[6:].isdigit():
            raise ValueError(
                f"Invalid ChEMBL ID '{chembl_id}'. "
                "Must be in the format CHEMBLnnn (e.g. CHEMBL25)."
            )
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.sdf"
        filename = f"{chembl_id}.sdf"
        source = "ChEMBL"

    # ── PubChem ────────────────────────────────────────────
    elif molecule_id.isdigit():
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{molecule_id}/SDF?record_type=3d"
        filename = f"{molecule_id}.sdf"
        source = "PubChem"

    else:
        raise ValueError(
            f"Unrecognised ID '{molecule_id}'. "
            "Provide a numeric PubChem CID or a ChEMBL ID (e.g. CHEMBL25)."
        )

    resp = requests.get(url)
    if resp.status_code != 200 or not resp.text.strip():
        raise RuntimeError(f"Failed to download SDF from {source}: {url}")

    return filename, resp.content


# ── Hero header ────────────────────────────────────────────────
st.markdown(
    """
    <div class="hero-header animate-in">
        <h1>🌐 Fetch Molecules</h1>
        <p>Download protein structures (PDB) from RCSB and
        ligand structures (SDF) from PubChem or ChEMBL.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

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
            downloads = []
            errors = []

            with st.spinner(f"Fetching {len(ids)} PDB file(s)..."):
                for pid in ids:
                    try:
                        fname, data = _fetch_pdb_bytes(pid)
                        downloads.append((fname, data))
                    except (ValueError, RuntimeError) as e:
                        errors.append(str(e))

            if errors:
                for err in errors:
                    st.error(f"❌ {err}")

            st.session_state["pdb_downloads"] = downloads
            if downloads:
                st.rerun()

    # ── Download buttons ────────────────────────────────────────
    if st.session_state["pdb_downloads"]:
        st.success(f"✅ {len(st.session_state['pdb_downloads'])} file(s) ready to download.")
        for fname, data in st.session_state["pdb_downloads"]:
            st.download_button(
                label=f"⬇️ {fname}",
                data=data,
                file_name=fname,
                mime="chemical/x-pdb",
                key=f"dl_pdb_{fname}",
            )


# ════════════════════════════════════════════════════════════════
# FETCH SDF (from PubChem or ChEMBL)
# ════════════════════════════════════════════════════════════════
with col_sdf:
    st.markdown('<div class="section-header"><h2>💊 Fetch SDF</h2></div>', unsafe_allow_html=True)
    st.caption("Download ligand structures from **PubChem** or **ChEMBL**.")

    sdf_input = st.text_input(
        "PubChem CIDs / ChEMBL IDs (comma-separated)",
        placeholder="e.g. 2244, CHEMBL25, 5988, CHEMBL1201559",
        key="sdf_id_input",
        help="Enter PubChem CIDs (numeric) and/or ChEMBL IDs — you can mix both.",
    )

    if st.button("🚀 Fetch SDF", type="primary", key="btn_fetch_sdf"):
        if not sdf_input.strip():
            st.warning("Please enter at least one CID or ChEMBL ID.")
        else:
            ids = [x.strip() for x in sdf_input.split(",") if x.strip()]
            downloads = []
            errors = []

            with st.spinner(f"Fetching {len(ids)} SDF file(s)..."):
                for mid in ids:
                    try:
                        fname, data = _fetch_sdf_bytes(mid)
                        downloads.append((fname, data))
                    except (ValueError, RuntimeError) as e:
                        errors.append(str(e))

            if errors:
                for err in errors:
                    st.error(f"❌ {err}")

            st.session_state["sdf_downloads"] = downloads
            if downloads:
                st.rerun()

    # ── Download buttons ────────────────────────────────────────
    if st.session_state["sdf_downloads"]:
        st.success(f"✅ {len(st.session_state['sdf_downloads'])} file(s) ready to download.")
        for fname, data in st.session_state["sdf_downloads"]:
            st.download_button(
                label=f"⬇️ {fname}",
                data=data,
                file_name=fname,
                mime="chemical/x-mdl-sdfile",
                key=f"dl_sdf_{fname}",
            )


# ── Footer ─────────────────────────────────────────────────────
st.caption("DockSuiteX © 2026")

