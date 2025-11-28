# DockSuiteX

**All-in-one Protein-Ligand Docking Package**

DockSuiteX is a comprehensive Python package for molecular docking that integrates multiple tools including MGLTools, P2Rank, AutoDock Vina, and AutoDock4. It provides a unified, user-friendly interface for protein and ligand preparation, binding pocket prediction, and molecular docking simulations.

## ✨ Key Features

- 🔧 **Automated Preparation**: Streamlined protein and ligand preparation using PDBFixer, Open Babel, and AutoDockTools
- 🎯 **Pocket Prediction**: Binding pocket prediction with P2Rank machine learning
- ⚡ **Multiple Docking Engines**: Support for both AutoDock Vina and AutoDock4
- 📊 **Batch Processing**: High-throughput screening with parallel execution
- 👁️ **Interactive Visualization**: 3D molecular visualization with NGLView
- 🛠️ **Utility Functions**: Built-in tools for fetching structures and parsing results

## 🚀 Flexible Batch Docking

!!! success "Complete Flexibility in Docking Scenarios"
    **DockSuiteX is uniquely designed to handle ALL possible docking combinations**, making it the most versatile docking package for high-throughput screening. Whether you're working with a single molecule or thousands, DockSuiteX adapts to your needs:

    ✅ **1 Ligand → 1 Receptor**  
    Standard molecular docking for detailed interaction analysis
    
    ✅ **Many Ligands → 1 Receptor**  
    Virtual screening - screen compound libraries against a single protein target
    
    ✅ **1 Ligand → Many Receptors**  
    Inverse docking - test one compound against multiple protein variants or different targets
    
    ✅ **Many Ligands → Many Receptors**  
    Complete combinatorial screening - explore all possible ligand-receptor interactions in one run

    All batch operations leverage **parallel execution** to maximize computational efficiency, enabling high-throughput virtual screening on your local machine without requiring cluster computing resources.

## 💻 Platform Support

!!! warning "Windows Only"
    DockSuiteX currently supports **Windows platforms only**. The package includes pre-compiled binaries for Windows.

## 📝 Quick Example

```python
from docksuitex import Protein, Ligand, PocketFinder, VinaDocking
from docksuitex.utils import fetch_pdb, fetch_sdf

# 1. Fetch and prepare receptor
protein_path = fetch_pdb("1HVR")
protein = Protein(protein_path)
protein.prepare()

# 2. Fetch and prepare ligand (Aspirin)
ligand_path = fetch_sdf("2244")
ligand = Ligand(ligand_path)
ligand.prepare(minimize="mmff94")

# 3. Find binding pockets
pocket_finder = PocketFinder("1HVR.pdbqt")
pockets = pocket_finder.run()
best_pocket = pockets[0]['center']

# 4. Run Vina docking
vina = VinaDocking(
    receptor="1HVR.pdbqt",
    ligand="2244.pdbqt",
    grid_center=best_pocket,
    grid_size=(20, 20, 20)
)
vina.run()
```

## 🚦 Getting Started

Check out the [Installation Guide](getting-started/installation.md) and Examples [Single Docking](examples/single_docking.ipynb), [Multiple Docking](examples/multiple_docking.ipynb) to get started.


## 🙏 Acknowledgments

This package builds upon and automates workflows using:

- [AutoDock4 & AutoDock Vina](http://autodock.scripps.edu/)
- [MGLTools](http://mgltools.scripps.edu/)
- [Open Babel](https://openbabel.org/)
- [PDBFixer](http://openmm.org/)
- [P2Rank](https://github.com/rdk/p2rank)
- [RCSB PDB](https://www.rcsb.org/) & [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- [NGLView](https://pypi.org/project/nglview/)

For detailed citations, see the [Citations](about/cite.md) page.
