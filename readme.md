# 🧬 DockSuiteX

[![PyPI version](https://badge.fury.io/py/docksuitex.svg)](https://badge.fury.io/py/docksuitex)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**All-in-one Protein-Ligand Docking Package**

DockSuiteX is a comprehensive Python package for molecular docking that integrates multiple tools including MGLTools, P2Rank, AutoDock Vina, and AutoDock4. It provides a unified, user-friendly interface for protein and ligand preparation, binding pocket prediction, and molecular docking simulations.

---

## ✨ Key Features

- **Automated Preparation**: Streamlined protein and ligand preparation using PDBFixer and Open Babel
- **Pocket Prediction**: Binding pocket prediction with P2Rank machine learning
- **Multiple Docking Engines**: Support for both AutoDock Vina and AutoDock4
- **Batch Processing**: High-throughput screening with parallel execution
- **Interactive Visualization**: 3D molecular visualization with NGLView
- **Utility Functions**: Built-in tools for fetching structures and parsing results

---

## 📦 Installation

### Requirements

- **Python**: 3.8 or higher
- **Operating System**: Windows (currently Windows-only)
- **Java**: Required for P2Rank ([Download here](https://adoptium.net/temurin/releases))

### Install from PyPI

```bash
pip install docksuitex
```

### Install from Source

```bash
git clone https://github.com/MangalamGSinha/DockSuiteX.git
cd DockSuiteX
pip install .
```

---

## 🚀 Quick Start

```python
from docksuitex import Protein, Ligand, PocketFinder, VinaDocking
from docksuitex.utils import fetch_pdb, fetch_sdf

# Fetch & prepare protein
protein_file = fetch_pdb("1HVR")
prot = Protein(protein_file)
prot.prepare()

# Fetch & prepare ligand
ligand_file = fetch_sdf(2244)  # Aspirin
lig = Ligand(ligand_file)
lig.prepare(minimize="mmff94")

# Predict binding pockets
finder = PocketFinder(prot)
pockets = finder.run()
center = pockets[0]['center']

# Run docking
vina = VinaDocking(
    receptor=prot,
    ligand=lig,
    grid_center=center,
    grid_size=(20, 20, 20)
)
vina.run()
vina.save_results("vina_results")
```

---

## 📂 Project Structure

```
docksuitex/
├── protein.py          # Protein preparation (PDBFixer, Open Babel, AutoDockTools)
├── ligand.py           # Ligand preparation (Open Babel, minimization, AutoDockTools)
├── pocket_finder.py    # Pocket detection with P2Rank
├── autodock4.py        # AutoDock4 docking wrapper
├── vina.py             # AutoDock Vina docking wrapper
├── batch_autodock4.py  # Parallel batch docking with AutoDock4
├── batch_vina.py       # Parallel batch docking with AutoDock Vina
├── utils/
│   ├── fetcher.py      # Fetch PDB (RCSB) & SDF (PubChem)
│   ├── viewer.py       # NGLView visualization
│   ├── parser.py       # Parse logs to CSV summaries
│   └── cleaner.py      # Reset temp folders and delete bin folder
│
└── bin/                # Auto-downloaded on first run
    ├── mgltools/       # MGLTools binaries and scripts
    ├── obabel/         # Open Babel executables
    ├── vina/           # AutoDock Vina executable (vina.exe)
    ├── p2rank/         # P2Rank executable and scripts
    └── autodock/       # AutoDock4 & AutoGrid executables
```

---

## 📚 Documentation

**Full documentation is available at: [https://yourusername.github.io/DockSuiteX](https://yourusername.github.io/DockSuiteX)**

- [Installation Guide](https://yourusername.github.io/DockSuiteX/getting-started/installation/)
- [Quick Start Tutorial](https://yourusername.github.io/DockSuiteX/getting-started/quickstart/)
- [Comprehensive Examples](https://yourusername.github.io/DockSuiteX/getting-started/examples/)
- [API Reference](https://yourusername.github.io/DockSuiteX/api/protein/)

---

## 📂 Examples

Runnable examples are available in the [examples](https://github.com/MangalamGSinha/DockSuiteX/tree/main/examples) folder.

---

## 🙏 Acknowledgments

This package builds upon:

- [AutoDock4 & AutoDock Vina](http://autodock.scripps.edu/)
- [MGLTools](http://mgltools.scripps.edu/)
- [Open Babel](https://openbabel.org/)
- [PDBFixer](http://openmm.org/)
- [P2Rank](https://github.com/rdk/p2rank)
- [RCSB PDB](https://www.rcsb.org/) & [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- [NGLView](https://pypi.org/project/nglview/)

---

## 📜 License

This project is licensed under the GNU GPL v3 License - see the [LICENSE](LICENSE) file for details.
