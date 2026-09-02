"""ProLIF interaction profiling module for DockSuiteX."""

import subprocess
import shutil
import os
from pathlib import Path
from typing import Union, Optional
import warnings
import logging

# Aggressive warning suppression across all modules and sub-processes
warnings.simplefilter("ignore", category=DeprecationWarning)
warnings.simplefilter("ignore", category=UserWarning)
warnings.filterwarnings("ignore", module="MDAnalysis")
os.environ['PYTHONWARNINGS'] = 'ignore'

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
logging.getLogger("MDAnalysis").setLevel(logging.ERROR)

import MDAnalysis as mda
import prolif as plf
import pandas as pd

# Locate obabel from PATH (installed via: pip install openbabel-wheel)
OBABEL_EXE = shutil.which("obabel")


class InteractionProfiler:
    """Protein–ligand interaction fingerprinting using ProLIF.

    This module wraps ProLIF to analyse the non-covalent interactions between
    a protein and docked ligand poses produced by AutoDock Vina or AutoDock4.

    The analysis workflow:
        1. Converting the multi-model PDBQT output → SDF using Open Babel.
        2. Loading the protein with MDAnalysis / ProLIF.
        3. Loading all ligand poses from the SDF via ``plf.sdf_supplier()``.
        4. Computing interaction fingerprints with ProLIF.
        5. Saving the resulting DataFrame, 2D/3D HTML, and barcode PNG.

    Supported Input:
        - Protein: PDB (.pdb) — typically the ``_prepared.pdb`` from Protein.prepare()
        - Ligand poses: multi-model PDBQT (.pdbqt) from Vina / AD4

    Note:
        The intermediate SDF file is saved in ``prolif_intermediates/``
        inside the output directory.
    """

    def __init__(
        self,
        protein_pdb: Union[str, Path],
        vina_output_pdbqt: Union[str, Path],
        _cpu: int = (os.cpu_count() or 2) - 1,
    ):
        """Initialize an InteractionProfiler with input paths and settings.

        Args:
            protein_pdb (str | Path): Path to the protein PDB file.
                This should be the original or prepared PDB (not PDBQT).
            vina_output_pdbqt (str | Path): Path to the multi-model
                ``output.pdbqt`` produced by AutoDock Vina or AutoDock4.
            _cpu (int, optional): Number of CPUs for ProLIF
                fingerprint calculation. Defaults to ``os.cpu_count() - 1``.

        Raises:
            FileNotFoundError: If either input file does not exist.
            ValueError: If file extensions are incorrect.
        """
        self.protein_pdb = Path(protein_pdb).resolve()
        self.vina_output_pdbqt = Path(vina_output_pdbqt).resolve()
        self.cpu = _cpu

        if not self.protein_pdb.is_file():
            raise FileNotFoundError(
                f"❌ Protein PDB file not found: {self.protein_pdb}")

        if not self.vina_output_pdbqt.is_file():
            raise FileNotFoundError(
                f"❌ Vina output PDBQT file not found: {self.vina_output_pdbqt}")

        if self.protein_pdb.suffix.lower() != ".pdb":
            raise ValueError(
                f"❌ Protein file must be .pdb, got '{self.protein_pdb.suffix}'")

        if self.vina_output_pdbqt.suffix.lower() != ".pdbqt":
            raise ValueError(
                f"❌ Ligand poses file must be .pdbqt, got '{self.vina_output_pdbqt.suffix}'")

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _convert_pdbqt_to_sdf(input_pdbqt: Path, output_sdf: Path) -> None:
        """Convert a multi-model PDBQT file to SDF using Open Babel.

        Args:
            input_pdbqt (Path): Input PDBQT file (may contain multiple models).
            output_sdf (Path): Output SDF file.

        Raises:
            RuntimeError: If Open Babel fails or produces an empty file.
        """
        cmd = [
            OBABEL_EXE,
            str(input_pdbqt),
            "-O", str(output_sdf),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise RuntimeError(
                f"❌ Open Babel PDBQT → SDF conversion failed:\n{result.stderr}")

        if not output_sdf.exists() or output_sdf.stat().st_size == 0:
            raise RuntimeError(
                f"❌ Open Babel produced an empty SDF file: {output_sdf}")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(
        self,
        save_to: Union[str, Path, None] = None,
    ) -> pd.DataFrame:
        """Run the ProLIF interaction fingerprint analysis.

        Performs the full workflow: convert PDBQT → SDF → load structures
        → compute fingerprints → save CSV, 2D/3D HTML, and barcode.

        Folder structure created inside ``save_to``::

            save_to/
            ├── prolif_intermediates/
            │   ├── <protein>.pdb                     # copy of input protein
            │   ├── <vina_output>.pdbqt               # copy of input poses
            │   └── ligands.sdf                       # converted ligand poses
            ├── prolif_interactions.csv
            ├── prolif_barcode.png
            ├── prolif_pose_1_interactions_2D.html
            ├── prolif_pose_1_interactions_3D.html
            └── ...

        Args:
            save_to (str | Path, optional): Directory where outputs will
                be saved. Defaults to ``prolif_results``.

        Returns:
            pandas.DataFrame: Interaction fingerprint table with one row
            per ligand pose and columns for each residue–interaction pair.

        Raises:
            RuntimeError: If Open Babel conversion or ProLIF analysis fails.
        """
        # Enforce global warning suppression specifically for Jupyter environments
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        warnings.filterwarnings("ignore", category=UserWarning)

        save_to = Path(save_to or "prolif_results").expanduser().resolve()
        self.save_dir = save_to
        save_to.mkdir(parents=True, exist_ok=True)

        intermediates_dir = save_to / "prolif_intermediates"
        intermediates_dir.mkdir(parents=True, exist_ok=True)

        # Copy input files into intermediates for reference
        shutil.copy2(self.protein_pdb, intermediates_dir / self.protein_pdb.name)
        shutil.copy2(self.vina_output_pdbqt, intermediates_dir / self.vina_output_pdbqt.name)

        # 1. Convert multi-model PDBQT → SDF
        sdf_path = intermediates_dir / "ligands.sdf"
        self._convert_pdbqt_to_sdf(self.vina_output_pdbqt, sdf_path)

        # 2. Load protein
        u_prot = mda.Universe(str(self.protein_pdb))
        protein_mol = plf.Molecule.from_mda(u_prot, inferrer=None)

        # 3. Load ligand poses from SDF
        from rdkit import Chem

        rdkit_mols = [
            mol for mol in Chem.SDMolSupplier(str(sdf_path), removeHs=False)
            if mol is not None
        ]

        # Remove residue metadata (fix HOH issue)
        for mol in rdkit_mols:
            for atom in mol.GetAtoms():
                atom.SetMonomerInfo(None)

        lig_mols = [plf.Molecule.from_rdkit(mol) for mol in rdkit_mols]

        # Run ProLIF safely
        fp = plf.Fingerprint()
        fp.run_from_iterable(
            lig_mols,
            protein_mol,
            residues=list(protein_mol.residues),
            n_jobs=self.cpu
        )

        # Store for visualization methods
        self.fp = fp
        self.lig_mols = lig_mols
        self.protein_mol = protein_mol

        # 5. Save CSV
        df = fp.to_dataframe(index_col="Pose")
        df.index = df.index + 1  # 1-indexed poses

        csv_path = save_to / "prolif_interactions.csv"
        df.drop("parent_indices", axis=1, level=1, errors="ignore").to_csv(csv_path)

        # 6. Save 2D and 3D Diagrams for all poses
        for pose_idx in range(1, len(lig_mols) + 1):
            try:
                self.save_2d_view(save_to / f"prolif_pose_{pose_idx}_interactions_2D.html", pose_idx=pose_idx)
                self.save_3d_view(save_to / f"prolif_pose_{pose_idx}_interactions_3D.html", pose_idx=pose_idx)
            except Exception as e:
                print(f"  ⚠️ Visualization failed for pose {pose_idx}: {e}")

        # 7. Save Barcode plot
        try:
            import matplotlib.pyplot as plt
            ax = self.fp.plot_barcode(xlabel="Pose")

            # Matplotlib plots 0-based trajectories by default; offset labels
            ticks = ax.get_xticks()
            ax.set_xticks(ticks)
            ax.set_xticklabels([str(int(t) + 1) for t in ticks])

            ax.figure.savefig(save_to / "prolif_barcode.png", bbox_inches="tight", dpi=300)
            plt.close(ax.figure)
        except Exception as e:
            print(f"  ⚠️ Barcode plot failed: {e}")

        print(f"✅ ProLIF analysis complete. Results: {save_to}")
        return df

    # ------------------------------------------------------------------
    # Visualization helpers
    # ------------------------------------------------------------------

    def save_2d_view(self, save_path: Union[str, Path], pose_idx: int = 1) -> None:
        """Save a 2D network diagram of the protein-ligand interactions.
        
        Args:
            save_path (Union[str, Path]): Path to save the 2D HTML file.
            pose_idx (int, optional): Pose index (1-based). Defaults to 1.
        """
        from prolif.plotting.network import LigNetwork
        frame = pose_idx - 1
        net = LigNetwork.from_fingerprint(self.fp, self.lig_mols[frame], kind="frame", frame=frame)
        net.save(str(save_path))

    def save_3d_view(self, save_path: Union[str, Path], pose_idx: int = 1) -> None:
        """Save an interactive 3D view of the protein-ligand complex interactions.

        Args:
            save_path (Union[str, Path]): Path to save the 3D HTML file.
            pose_idx (int, optional): Pose index (1-based). Defaults to 1.
        """
        frame = pose_idx - 1
        view = self.fp.plot_3d(self.lig_mols[frame], self.protein_mol, frame=frame, display_all=False)
        with open(save_path, "w", encoding="utf-8") as f:
            f.write(str(view._make_html()))

    def view_interactions(self):
        """Render an interactive viewer with Prev/Next buttons and a 2D/3D/Barcode toggle.

        This method is designed for use in Jupyter Notebooks. It displays
        a widget with controls to switch between different ligand poses
        and different visualization modes (2D network, 3D complex, or
        the interaction barcode).
        """
        import ipywidgets as widgets
        from IPython.display import display, clear_output, HTML
        import base64

        save_dir = getattr(self, "save_dir", Path("prolif_results").resolve())
        max_poses = len(self.lig_mols)
        state = {"pose": 1, "mode": "3D"}

        # --- Widgets ---
        btn_prev  = widgets.Button(description="◄ Prev", button_style="info", layout=widgets.Layout(width="80px"))
        btn_next  = widgets.Button(description="Next ►", button_style="info", layout=widgets.Layout(width="80px"))
        lbl_pose  = widgets.Label(value=f"Pose 1 of {max_poses}", layout=widgets.Layout(width="70px"))
        toggle    = widgets.ToggleButtons(
            options=["2D", "3D", "Barcode"], value="3D", button_style="",
            style={"button_width": "80px", "font_weight": "bold"},
            layout=widgets.Layout(margin="0px")
        )
        out       = widgets.Output()

        def render():
            with out:
                clear_output(wait=True)

                if state["mode"] == "Barcode":
                    # Hide pose navigation for barcode (shows all poses)
                    btn_prev.layout.visibility = "hidden"
                    btn_next.layout.visibility = "hidden"
                    lbl_pose.value = "All Poses"

                    barcode_path = save_dir / "prolif_barcode.png"
                    if barcode_path.exists():
                        b64 = base64.b64encode(barcode_path.read_bytes()).decode("utf-8")
                        display(HTML(f'<img src="data:image/png;base64,{b64}" style="width:100%;"/>'))
                    else:
                        print("❌ prolif_barcode.png not found. Run profiler.run() first.")
                else:
                    # Show pose navigation for 2D/3D
                    btn_prev.layout.visibility = "visible"
                    btn_next.layout.visibility = "visible"
                    lbl_pose.value = f"Pose {state['pose']} of {max_poses}"

                    html_path = save_dir / f"prolif_pose_{state['pose']}_interactions_{state['mode']}.html"
                    if html_path.exists():
                        b64 = base64.b64encode(html_path.read_bytes()).decode("utf-8")
                        display(HTML(f'<iframe src="data:text/html;base64,{b64}" width="100%" height="620px" style="border:none;"></iframe>'))
                    else:
                        print(f"❌ {html_path.name} not found. Run profiler.run() first.")

        def on_prev(b):
            state["pose"] = (state["pose"] - 2) % max_poses + 1
            render()

        def on_next(b):
            state["pose"] = state["pose"] % max_poses + 1
            render()

        def on_toggle(change):
            state["mode"] = change["new"]
            render()

        btn_prev.on_click(on_prev)
        btn_next.on_click(on_next)
        toggle.observe(on_toggle, names="value")

        spacer = widgets.HBox(layout=widgets.Layout(flex="1 1 auto"))
        view_label = widgets.Label("View:", layout=widgets.Layout(width="auto", margin="0 4px 0 0"))
        nav = widgets.HBox(
            [btn_prev, lbl_pose, btn_next, spacer, view_label, toggle],
            layout=widgets.Layout(width="100%", align_items="center")
        )
        display(widgets.VBox([nav, out]))
        render()