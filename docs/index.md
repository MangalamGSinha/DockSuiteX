# DockSuiteX

**All-in-one Protein-Ligand Docking Package**

DockSuiteX is a comprehensive Python package and GUI application for molecular docking that integrates multiple tools including MGLTools, P2Rank, AutoDock Vina, and AutoDock4. It provides a unified, user-friendly interface for protein and ligand preparation, binding pocket prediction, and molecular docking simulations.

## ✨ Key Features

- 🔧 **Automated Preparation**: Streamlined protein and ligand preparation using PDBFixer, Open Babel, and AutoDockTools
- 🎯 **Pocket Prediction**: Binding pocket prediction with P2Rank machine learning
- ⚡ **Multiple Docking Engines**: Support for both AutoDock Vina and AutoDock4
- 📊 **Batch Processing**: High-throughput screening with parallel execution
- 👁️ **Interactive Visualization**: 3D molecular visualization with NGLView
- 🛠️ **Utility Functions**: Built-in tools for fetching structures and parsing results
- 🖥️ **GUI Application**: Complete desktop application — no coding required

## 🚀 Batch Docking

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

    All batch operations leverage parallel execution to maximize computational efficiency, with the option to set the number of cpu threads to use. 

## 💻 Platform Support

!!! warning "Windows Only"
    DockSuiteX currently supports **Windows platforms only**. The package includes pre-compiled binaries for Windows.


## 🚦 Getting Started

=== "GUI Application"

    DockSuiteX includes a full-featured **GUI application** built with Streamlit.
    
    To use the GUI:
    ```bash
    pip install docksuitex
    docksuitex
    ```
    
    See the [GUI Application Guide](getting-started/gui.md).

=== "Python Package"

    Check out the [Installation Guide](getting-started/python_package.md) and Examples: [Single Docking](api/examples/single_docking.ipynb) | [Batch Docking](api/examples/batch_docking.ipynb) to get started with the Python interface.

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
