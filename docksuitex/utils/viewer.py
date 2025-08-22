from pathlib import Path
from typing import Literal
import py3Dmol

def view_molecule(
    file_path: str,
    style: Literal["stick", "line", "sphere", "cartoon", "surface"] = "stick",
    background: str = "white",
    color: str = "spectrum",
    width: int = 500,
    height: int = 500,
):
    """
    Renders a molecular structure in a Jupyter Notebook using py3Dmol.

    Args:
        file_path (str): Path to the molecular file (.pdb, .pdbqt, .mol2, or .sdf).
        style (Literal): Visualization style ('stick', 'line', 'sphere', 'cartoon', 'surface').
        background (str): Background color of the 3Dmol viewer (e.g., 'white', 'black').
        color (str): Coloring method for atoms or residues (e.g., 'spectrum', 'chain').
        width (int): Width of the rendered viewer window in pixels.
        height (int): Height of the rendered viewer window in pixels.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file format is unsupported.

    Returns:
        py3Dmol.view: A 3Dmol viewer object rendered in a Jupyter notebook.
    """
    path = Path(file_path)

    if not path.exists():
        raise FileNotFoundError(f"❌ File not found: {file_path}")

    ext = path.suffix.lower().lstrip(".")
    supported_formats = {"pdb", "pdbqt", "mol2", "sdf"}

    if ext not in supported_formats:
        raise ValueError(f"❌ Unsupported format '.{ext}'. Supported formats: {supported_formats}")

    with open(path, "r") as file:
        molecule_data = file.read()

    # .pdbqt format behaves like .pdb in py3Dmol
    mol_format = "pdb" if ext in {"pdb", "pdbqt"} else ext

    viewer = py3Dmol.view(width=width, height=height)
    viewer.setBackgroundColor(background)
    viewer.addModel(molecule_data, mol_format)
    viewer.setStyle({ }, { style: { "color": color } })
    viewer.zoomTo()

    return viewer.show()
