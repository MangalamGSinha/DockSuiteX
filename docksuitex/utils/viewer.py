"""Molecular visualization utilities for DockSuiteX."""

import nglview as nv
import ipywidgets as widgets
import tempfile
import threading
import time
from IPython.display import display
from pathlib import Path



def view_molecule(file_path: str | Path) -> nv.NGLWidget:
    """Render a molecular structure in Jupyter Notebook using NGLView.

    Args:
        file_path (str | Path): Path to the molecular file (.pdb, .pdbqt, .mol2, or .sdf).

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


def view_grid_box(
    protein_file: str | Path,
    grid_center: tuple[float, float, float],
    grid_size: tuple[float, float, float] = (20, 20, 20),
    color: list[float] | None = None,
    center_color: list[float] | None = None,
    radius: float = 0.2,
    cross_size: float = 2.0,
) -> nv.NGLWidget:
    """Visualize a protein structure with a docking grid box overlay.

    Renders the protein and draws a wireframe box (12 edges as cylinders)
    around the specified grid center, plus a green 3-axis cross marking
    the pocket center.  Useful for inspecting the search space before or
    after a docking run.

    Args:
        protein_file (str | Path): Path to the protein structure file
            (.pdb, .pdbqt, .mol2, or .sdf).
        grid_center (tuple[float, float, float]): (x, y, z) coordinates
            of the grid box center in Ångströms.
        grid_size (tuple[float, float, float], optional): (sx, sy, sz)
            dimensions of the grid box in Ångströms. Defaults to
            (20, 20, 20).
        color (list[float] | None, optional): RGB color for the box
            edges, each component in [0, 1].
            Defaults to red ``[1, 0, 0]``.
        center_color (list[float] | None, optional): RGB color for the
            center cross. Defaults to green ``[0, 1, 0]``.
        radius (float, optional): Cylinder radius for the box edges.
            Defaults to 0.2.
        cross_size (float, optional): Half-length of each axis line
            in the center cross, in Ångströms. Defaults to 2.0.

    Returns:
        nv.NGLWidget: An interactive NGLView widget (display in Jupyter).

    Raises:
        FileNotFoundError: If *protein_file* does not exist.

    Example::

        from docksuitex.utils import view_grid_box
        view_grid_box("receptor.pdbqt", grid_center=(10, 15, 20))
    """
    if color is None:
        color = [1, 0, 0]
    if center_color is None:
        center_color = [0, 1, 0]

    protein_file = Path(protein_file).resolve()
    if not protein_file.exists():
        raise FileNotFoundError(f"❌ File not found: {protein_file}")

    view = nv.show_file(str(protein_file))

    cx, cy, cz = grid_center
    sx, sy, sz = [s / 2 for s in grid_size]

    # Eight corners of the box
    # Index: 0(-,-,-)  1(-,-,+)  2(-,+,-)  3(-,+,+)
    #        4(+,-,-)  5(+,-,+)  6(+,+,-)  7(+,+,+)
    corners = [
        (cx + dx, cy + dy, cz + dz)
        for dx in (-sx, sx)
        for dy in (-sy, sy)
        for dz in (-sz, sz)
    ]

    RED   = [1, 0, 0]
    GREEN = [0, 1, 0]
    BLUE  = [0, 0, 1]

    # Twelve edges, coloured by their parallel axis (ADT style)
    # X-parallel edges → RED
    for i, j in [(0, 4), (1, 5), (2, 6), (3, 7)]:
        view.shape.add_cylinder(list(corners[i]), list(corners[j]), RED, radius)
    # Y-parallel edges → GREEN
    for i, j in [(0, 2), (1, 3), (4, 6), (5, 7)]:
        view.shape.add_cylinder(list(corners[i]), list(corners[j]), GREEN, radius)
    # Z-parallel edges → BLUE
    for i, j in [(0, 1), (2, 3), (4, 5), (6, 7)]:
        view.shape.add_cylinder(list(corners[i]), list(corners[j]), BLUE, radius)

    # Center cross — axis indicator at the pocket center
    cross_radius = radius * 1.5
    view.shape.add_cylinder(
        [cx - cross_size, cy, cz], [cx + cross_size, cy, cz],
        RED, cross_radius,
    )
    view.shape.add_cylinder(
        [cx, cy - cross_size, cz], [cx, cy + cross_size, cz],
        GREEN, cross_radius,
    )
    view.shape.add_cylinder(
        [cx, cy, cz - cross_size], [cx, cy, cz + cross_size],
        BLUE, cross_radius,
    )

    # ── Translucent panels via NGL.js injection ──
    # nglview's Python wrapper doesn't apply opacity to buffer meshes
    # correctly, so we inject the same JS used in the Streamlit viewer.
    corners_js = "[" + ",".join(
        f"[{c[0]},{c[1]},{c[2]}]" for c in corners
    ) + "]"

    panel_js = f"""
    var corners = {corners_js};
    var RED=[1,0,0], GREEN=[0,1,0], BLUE=[0,0,1];
    var panelShape = new NGL.Shape("gridbox_panels");
    var faceDefs = [
      [[0,1,3,2], RED],  [[4,5,7,6], RED],
      [[0,1,5,4], GREEN],[[2,3,7,6], GREEN],
      [[0,2,6,4], BLUE], [[1,3,7,5], BLUE]
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
    this.stage.addComponentFromObject(panelShape).addRepresentation(
      "buffer", {{opacity: 0.15, side: "double"}}
    );
    """
    view._js(panel_js)

    return view


def view_docked_poses(protein_file: str | Path, ligand_file: str | Path) -> None:
    """Visualize docked poses (multiple poses) of a ligand with a protein.

    This uses NGLView and interactive Jupyter widgets.

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
