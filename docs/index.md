# DockSuiteX

**All-in-one Protein-Ligand Docking Package**

DockSuiteX is a comprehensive Python package for molecular docking that integrates multiple tools including MGLTools, P2Rank, AutoDock Vina, and AutoDock4. It provides a unified, user-friendly interface for protein and ligand preparation, binding pocket prediction, and molecular docking simulations.

## Key Features

- **Automated Preparation**: Streamlined protein and ligand preparation using PDBFixer and Open Babel
- **Pocket Prediction**: Binding pocket prediction with P2Rank machine learning
- **Multiple Docking Engines**: Support for both AutoDock Vina and AutoDock4
- **Batch Processing**: High-throughput screening with parallel execution
- **Interactive Visualization**: 3D molecular visualization with NGLView
- **Utility Functions**: Built-in tools for fetching structures and parsing results

## Platform Support

!!! warning "Windows Only"
    DockSuiteX currently supports **Windows platforms only**. The package includes pre-compiled binaries for Windows.

## Quick Example

```python
from docksuitex import Protein, Ligand, VinaDocking

# Prepare protein
protein = Protein("protein.pdb")
protein.prepare(save_to="prepared_protein.pdbqt")

# Prepare ligand
ligand = Ligand("ligand.sdf")
ligand.prepare(save_to="prepared_ligand.pdbqt")

# Run docking
docking = VinaDocking(
    receptor="prepared_protein.pdbqt",
    ligand="prepared_ligand.pdbqt",
    grid_center=(10.0, 15.0, 20.0)
)
results = docking.run()
```

## Integrated Tools

DockSuiteX bundles and integrates the following tools:

- **MGLTools**: Protein and ligand preparation (prepare_receptor4.py, prepare_ligand4.py)
- **AutoDock Vina**: Fast gradient-based docking
- **AutoDock4**: Classic genetic algorithm-based docking
- **P2Rank**: Machine learning pocket prediction
- **Open Babel**: Format conversion and energy minimization
- **PDBFixer**: Protein structure repair

## Getting Started

Check out the [Installation Guide](getting-started/installation.md) to get started, or jump straight to the [Quick Start Tutorial](getting-started/quickstart.md) to see DockSuiteX in action!

---

## Acknowledgments

This package builds upon and automates workflows using:

- [AutoDock4 & AutoDock Vina](http://autodock.scripps.edu/)
- [MGLTools](http://mgltools.scripps.edu/)
- [Open Babel](https://openbabel.org/)
- [PDBFixer](http://openmm.org/)
- [P2Rank](https://github.com/rdk/p2rank)
- [RCSB PDB](https://www.rcsb.org/) & [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- [NGLView](https://pypi.org/project/nglview/)

For detailed citations, see the [Citations](about/cite.md) page.
