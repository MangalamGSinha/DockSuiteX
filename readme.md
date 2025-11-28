# DockSuiteX

[![PyPI version](https://badge.fury.io/py/docksuitex.svg)](https://badge.fury.io/py/docksuitex)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Documentation](https://img.shields.io/badge/docs-mkdocs-blue.svg)](https://mangalamgsinha.github.io/DockSuiteX/)

**All-in-one Protein-Ligand Docking Package**

DockSuiteX is a comprehensive Python package for molecular docking that integrates multiple tools including MGLTools, P2Rank, AutoDock Vina, and AutoDock4. It provides a unified, user-friendly interface for protein and ligand preparation, binding pocket prediction, and molecular docking simulations.

## ✨ Key Features

- 🔧 **Automated Preparation**: Streamlined protein and ligand preparation using PDBFixer and Open Babel
- 🎯 **Pocket Prediction**: Binding pocket prediction with P2Rank machine learning
- ⚡ **Multiple Docking Engines**: Support for both AutoDock Vina and AutoDock4
- 📊 **Batch Processing**: High-throughput screening with parallel execution
- 👁️ **Interactive Visualization**: 3D molecular visualization with NGLView
- 🛠️ **Utility Functions**: Built-in tools for fetching structures and parsing results

## 🚀 Flexible Batch Docking - A Unique Feature

**DockSuiteX is uniquely designed to handle ALL possible docking combinations**, making it the most versatile docking package for high-throughput screening. Whether you're working with a single molecule or thousands, DockSuiteX adapts to your needs:

**✅ 1 Ligand → 1 Receptor**
Standard molecular docking for detailed interaction analysis

**✅ Many Ligands → 1 Receptor**
Virtual screening - screen compound libraries against a single protein target

✅ **1 Ligand → Many Receptors**
Inverse docking - test one compound against multiple protein variants or different targets

✅ **Many Ligands → Many Receptors**
Complete combinatorial screening - explore all possible ligand-receptor interactions in one run

All batch operations leverage **parallel execution** to maximize computational efficiency, enabling high-throughput virtual screening on your local machine without requiring cluster computing resources.

## 📚 Documentation

**📖 For complete documentation, API reference, and examples, visit:**

### **[https://mangalamgsinha.github.io/DockSuiteX/](https://mangalamgsinha.github.io/DockSuiteX/)**

---

## 📦 Quick Install

```bash
pip install docksuitex
```

**Requirements:** Python 3.8+, Windows OS, Java (for P2Rank)

---

## 📜 License

GNU GPL v3 - see [LICENSE](LICENSE) for details.

---

## 🙏 Acknowledgments

This package builds upon and automates workflows using:

- [AutoDock4 &amp; AutoDock Vina](http://autodock.scripps.edu/)
- [MGLTools](http://mgltools.scripps.edu/)
- [Open Babel](https://openbabel.org/)
- [PDBFixer](http://openmm.org/)
- [P2Rank](https://github.com/rdk/p2rank)
- [RCSB PDB](https://www.rcsb.org/) & [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- [NGLView](https://pypi.org/project/nglview/)
