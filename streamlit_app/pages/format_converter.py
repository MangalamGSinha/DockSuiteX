"""
DockSuiteX — Format Converter Page
Simple molecular file format conversion using the bundled Open Babel.
"""
import streamlit as st
import traceback, tempfile, zipfile, io, sys
from pathlib import Path

# ── Format presets ─────────────────────────────────────────────
FORMAT_PRESETS = {
    "MOL2   — Tripos Mol2":              "mol2",
    "SDF    — MDL SD File":              "sdf",
    "PDB    — Protein Data Bank":        "pdb",
    "PDBQT  — AutoDock":                 "pdbqt",
    "MOL    — MDL MOL":                  "mol",
    "XYZ    — Cartesian coordinates":    "xyz",
    "SMI    — SMILES":                   "smi",
    "CAN    — Canonical SMILES":         "can",
    "INCHI  — IUPAC InChI":              "inchi",
    "CIF    — Crystallographic":         "cif",
    "MMCIF  — Macromolecular CIF":       "mmcif",
    "GRO    — GROMACS format":           "gro",
    "FASTA  — Sequence format":          "fasta",
    "RXN    — MDL Reaction":             "rxn",
    "CDX    — ChemDraw":                 "cdx",
}


# ── Session state init ─────────────────────────────────────────
for key in ["conv_result_name", "conv_output_bytes"]:
    if key not in st.session_state:
        st.session_state[key] = None

if "conv_tmp_dir" not in st.session_state:
    st.session_state["conv_tmp_dir"] = tempfile.mkdtemp(prefix="docksuitex_conv_")


# ── Hero header ────────────────────────────────────────────────
st.markdown(
    """
    <div class="hero-header animate-in">
        <h1>🧩 Format Converter</h1>
        <p>Convert molecular files between formats using Open Babel.</p>
    </div>
    """,
    unsafe_allow_html=True,
)


# ── Upload & Options ───────────────────────────────────────────
uploaded = st.file_uploader(
    "Upload a molecular file",
    type=list(FORMAT_PRESETS.values()),
    key="conv_upload",
    help="Supports MOL2, SDF, PDB, MOL, XYZ, PDBQT, SMILES, CIF, InChI.",
)

out_format_label = st.selectbox(
    "🎯 Target format",
    options=list(FORMAT_PRESETS.keys()),
    index=0,
    key="conv_target",
    help="Select the output file format.",
)
out_format = FORMAT_PRESETS[out_format_label]

gen3d = st.checkbox(
    "Generate 3D coordinates",
    value=False,
    key="conv_gen3d",
    help="Enable when converting from SMILES or other 2D/text formats to 3D formats like PDB, MOL2, SDF.",
)

split = st.checkbox(
    "Split multi-model files",
    value=False,
    key="conv_split",
    help="Split multi-molecule/model files (e.g. multi-molecule SDF) into individual files.",
)


# ── Convert Button ─────────────────────────────────────────────
st.markdown("---")

if st.button("🚀 Convert", type="primary", key="btn_convert", disabled=uploaded is None):
    tmp_dir = Path(st.session_state["conv_tmp_dir"])
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Save uploaded file to disk
    input_path = tmp_dir / uploaded.name
    input_path.write_bytes(uploaded.getvalue())

    try:
        with st.spinner("Converting..."):
            from docksuitex.utils.converter import convert
            result = convert(
                input=input_path,
                output_format=out_format,
                gen3d=gen3d,
                split=split,
            )

        if isinstance(result, list):
            if not result:
                st.error("❌ Open Babel produced no output files despite reporting success.")
                st.stop()

            # Create a ZIP of all converted files
            zip_buf = io.BytesIO()
            with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
                for fp in result:
                    if fp.exists():
                        zf.write(fp, fp.name)
            
            st.session_state["conv_result_name"] = f"{input_path.stem}_converted.zip"
            st.session_state["conv_output_bytes"] = zip_buf.getvalue()
        else:
            st.session_state["conv_result_name"] = result.name
            st.session_state["conv_output_bytes"] = result.read_bytes()
        
        st.rerun()
    except Exception as e:
        st.error(f"❌ {e}")
        with st.expander("Traceback"):
            st.code(traceback.format_exc())


# ── Download result ────────────────────────────────────────────
if st.session_state.get("conv_output_bytes") and st.session_state.get("conv_result_name"):
    st.success(f"✅ Conversion complete!")
    st.download_button(
        label=f"⬇️ Download {st.session_state['conv_result_name']}",
        data=st.session_state["conv_output_bytes"],
        file_name=st.session_state["conv_result_name"],
        mime="application/octet-stream",
        key="dl_conv",
    )


# ── Footer ─────────────────────────────────────────────────────
st.caption("DockSuiteX © 2026")
