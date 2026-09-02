# API Overview

The API reference documents the Python modules that power both the **Python package** and the **GUI application**. Whether you're writing scripts or using the GUI, the same core engine runs under the hood.

!!! info "Using the GUI?"
    If you're using the GUI application, you don't need to interact with the API directly. The GUI wraps all these modules in a graphical interface. This reference is for Python developers who want to use DockSuiteX programmatically.

## Project Structure

The DockSuiteX package is organized into core modules, batch processing modules, and utility functions:

```
docksuitex/
├── protein.py        # Protein preparation (PDBFixer, Open Babel, AutoDockTools)
├── ligand.py         # Ligand preparation (Open Babel, AutoDockTools)
├── grid_calculator.py# Grid box calculation and pocket prediction using P2Rank
├── autodock4.py              # AutoDock4 docking wrapper
├── vina.py                   # AutoDock Vina docking wrapper
├── interaction_profiler.py   # ProLIF interaction fingerprinting
├── cli.py                    # Command Line Interface (CLI) entry point
│ 
├── batch_docking/             # Batch processing modules
│   ├── batch_protein.py       # Parallel batch protein preparation
│   ├── batch_ligand.py        # Parallel batch ligand preparation
│   ├── batch_grid_calculator.py # Parallel batch grid calculation
│   ├── batch_autodock4.py     # Parallel batch docking with AutoDock4
│   ├── batch_vina.py                    # Parallel batch docking with AutoDock Vina
│   └── batch_interaction_profiler.py    # Parallel batch interaction profiling
│ 
├── utils/
│   ├── fetcher.py      # Fetch PDB (RCSB) & SDF (PubChem)
│   ├── viewer.py       # NGLView visualization
│   ├── parser.py       # Parse logs to CSV summaries
│   └── converter.py    # Format conversion using Open Babel
│ 
├── bin/                # Auto-downloaded on first run
│    ├── mgltools/       # MGLTools(AutoDockTools) binaries and scripts
│    ├── obabel/         # Open Babel executables
│    ├── vina/           # AutoDock Vina executable
│    ├── p2rank/         # P2Rank executable and scripts
│    ├── java/           # Java Runtime Environment (JRE) for P2Rank
│    └── autodock/       # AutoDock4 & AutoGrid executables
│
└── streamlit_app/          # Streamlit GUI Application
    ├── app.py              # Main GUI entry point
    ├── style.css           # Custom CSS for the GUI
    └── pages/              # Individual GUI pages (Docking, Fetch, etc.)
```

## Module Categories

### Core Modules

Fundamental classes for molecular docking workflows:

- **[Protein](protein.md)** - Protein structure preparation and conversion to PDBQT 
- **[Ligand](ligand.md)** - Ligand preparation with energy minimization options and conversion to PDBQT
- **[GridCalculator](grid-calculator.md)** - Grid box calculation and binding pocket prediction using P2Rank
- **[VinaDocking](vina.md)** - AutoDock Vina docking interface
- **[AD4Docking](autodock4.md)** - AutoDock4 docking interface
- **[InteractionProfiler](interaction-profiler.md)** - Post-docking interaction fingerprinting with ProLIF

### Batch Processing

High-throughput docking with parallel execution:

- **[BatchProtein](batch-protein.md)** - Batch protein preparation for multiple receptors
- **[BatchLigand](batch-ligand.md)** - Batch ligand preparation for multiple ligands
- **[BatchGridCalculator](batch-grid-calculator.md)** - Batch grid box calculation for multiple receptors
- **[BatchVinaDocking](batch-vina.md)** - Parallel Vina docking for multiple receptors, ligands and pockets
- **[BatchAD4Docking](batch-autodock4.md)** - Parallel AutoDock4 docking for multiple receptors, ligands and pockets
- **[BatchInteractionProfiler](batch-interaction-profiler.md)** - Parallel interaction profiling for batch docking results



### Utilities

Helper functions for common tasks:

- **[Viewer](viewer.md)** - 3D molecular visualization with NGLView
- **[Parser](parser.md)** - Parse docking outputs to CSV format
- **[Fetcher](fetcher.md)** - Download structures from RCSB PDB and PubChem
- **[Converter](converter.md)** - Format conversion using Open Babel

## Binary Dependencies

DockSuiteX automatically downloads and manages the following third-party tools:

| Tool | Purpose | Location |
|------|---------|----------|
| MGLTools | Protein/ligand preparation scripts | `bin/mgltools/` |
| AutoDock Vina | Fast gradient-based docking | `bin/vina/` |
| AutoDock4 | Classic genetic algorithm docking | `bin/autodock/` |
| P2Rank | Machine learning pocket prediction | `bin/p2rank/` |
| Java JRE | Runtime for P2Rank | `bin/java/` |
| Open Babel | Format conversion and energy minimization | `bin/obabel/` |

All binaries are downloaded from the [DockSuiteX_Binaries](https://github.com/MangalamGSinha/DockSuiteX_Binaries) repository on first import.
