import nglview as nv
import ipywidgets as widgets
import tempfile
import threading
import time
from IPython.display import display
from pathlib import Path


def view_molecule(file_path: str | Path) -> nv.NGLWidget:
    """
    Renders a molecular structure in Jupyter Notebook using NGLView.

    Args:
        file_path (str): Path to the molecular file (.pdb, .pdbqt, .mol2, or .sdf).

    Returns:
        nv.NGLWidget: An NGLView widget.

    Raises:
        FileNotFoundError: If the specified file does not exist.
    """
    file_path = Path(file_path).resolve()
    if not file_path.exists():
        raise FileNotFoundError(f"❌ File not found: {file_path}")

    view = nv.show_file(str(file_path))   # replace with your protein file
    return view



def view_results(protein_file: str | Path, ligand_file: str | Path) -> None:
    """
    Visualize docking results (multiple poses) of a ligand with a protein
    using NGLView and interactive Jupyter widgets.

    Features:
        - Step through individual docking poses.
        - Toggle between showing one pose at a time or all poses simultaneously.
        - Play/Pause automatic animation of poses.
        - Adjust animation speed with a slider.

    Args:
        protein_file (str | Path): Path to the receptor protein file (e.g., .pdb).
        ligand_file (str | Path): Path to the ligand docking results file (.pdbqt).
            The file should contain multiple docking poses in MODEL/ENDMDL blocks.

    Returns:
        None: Displays the visualization and interactive controls directly
        in the Jupyter Notebook.

    Raises:
        FileNotFoundError: If the protein or ligand file does not exist.
        ValueError: If the ligand file does not contain valid MODEL/ENDMDL blocks.
    """
    protein_file = str(Path(protein_file).resolve())
    ligand_file = str(Path(ligand_file).resolve())

    # Extract ligand poses into temp files
    poses, current = [], 0
    with open(ligand_file) as f:
        pose = []
        for line in f:
            if line.startswith("MODEL"):
                pose = [line]
            elif line.startswith("ENDMDL"):
                pose.append(line)
                tmp = tempfile.NamedTemporaryFile(
                    delete=False, suffix=".pdbqt", mode="w")
                tmp.write("".join(pose))
                tmp.close()
                poses.append(tmp.name)
            else:
                pose.append(line)

    playing, play_speed = [False], [1.0]

    # NGL Viewer
    view = nv.NGLWidget()
    protein = view.add_component(protein_file)
    protein.add_representation(
        "cartoon", selection="protein")

    # Widgets
    pose_label = widgets.Label()
    show_all = widgets.ToggleButton(description="Show All Poses")
    play_btn = widgets.ToggleButton(description="Play", icon="play")
    prev_btn = widgets.Button(description="◀️ Prev")
    next_btn = widgets.Button(description="Next ▶️")
    speed_slider = widgets.FloatSlider(
        value=1.0, min=0.2, max=5, step=0.1, description="Speed:")

    # Track ligand components
    ligand_components = []

    def update(_=None):
        # Remove old ligands
        for lig in ligand_components:
            try:
                view.remove_component(lig)
            except Exception:
                pass
        ligand_components.clear()

        if show_all.value:
            for lig in poses:
                comp = view.add_component(lig)
                comp.add_representation("ball+stick")
                ligand_components.append(comp)
            pose_label.value = f"All poses ({len(poses)})"
        else:
            comp = view.add_component(poses[current])
            comp.add_representation("ball+stick")
            ligand_components.append(comp)
            pose_label.value = f"Pose: {current+1}/{len(poses)}"

    def step(d):
        nonlocal current
        if not show_all.value:
            current = (current + d) % len(poses)
            update()

    def toggle(change):
        playing[0] = change["new"]
        play_btn.description, play_btn.icon = (
            "Pause", "pause") if playing[0] else ("Play", "play")
        if playing[0]:
            threading.Thread(target=loop, daemon=True).start()

    def loop():
        while playing[0]:
            time.sleep(1 / play_speed[0])
            if not show_all.value:
                step(1)

    # Widget callbacks
    show_all.observe(update, "value")
    play_btn.observe(toggle, "value")
    prev_btn.on_click(lambda _: step(-1))
    next_btn.on_click(lambda _: step(1))
    speed_slider.observe(
        lambda c: play_speed.__setitem__(0, c["new"]), "value")

    # Initial update
    update()

    # Display
    controls = widgets.VBox([
        widgets.HBox([prev_btn, pose_label, next_btn]),
        widgets.HBox([play_btn, speed_slider]),
        show_all
    ])
    display(controls, view)
