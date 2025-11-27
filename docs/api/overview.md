# API Overview

## Project Structure

The DockSuiteX package is organized into core modules, batch processing modules, and utility functions:

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

## Module Categories

### Core Modules

Fundamental classes for molecular docking workflows:

- **[Protein](protein.md)** - Protein structure preparation and conversion to PDBQT format
- **[Ligand](ligand.md)** - Ligand preparation with energy minimization options
- **[VinaDocking](vina.md)** - AutoDock Vina docking interface
- **[AD4Docking](autodock4.md)** - AutoDock4 docking interface
- **[PocketFinder](pocket-finder.md)** - Binding pocket prediction using P2Rank

### Batch Processing

High-throughput docking with parallel execution:

- **[BatchVinaDocking](batch-vina.md)** - Parallel Vina docking for multiple ligands and pockets
- **[BatchAD4Docking](batch-autodock4.md)** - Parallel AutoDock4 docking
- **[BatchProtein](batch-protein.md)** - Batch protein preparation
- **[BatchLigand](batch-ligand.md)** - Batch ligand preparation
- **[BatchPocketFinder](batch-pocket-finder.md)** - Batch pocket prediction

### Utilities

Helper functions for common tasks:

- **[Viewer](viewer.md)** - 3D molecular visualization with NGLView
- **[Parser](parser.md)** - Parse docking logs to CSV format
- **[Fetcher](fetcher.md)** - Download structures from RCSB PDB and PubChem
- **[Cleaner](cleaner.md)** - Clean temporary files and manage binaries

## Binary Dependencies

DockSuiteX automatically downloads and manages the following third-party tools:

| Tool | Purpose | Location |
|------|---------|----------|
| MGLTools | Protein/ligand preparation scripts | `bin/mgltools/` |
| AutoDock Vina | Fast gradient-based docking | `bin/vina/` |
| AutoDock4 | Classic genetic algorithm docking | `bin/autodock/` |
| P2Rank | Machine learning pocket prediction | `bin/p2rank/` |
| Open Babel | Format conversion and minimization | `bin/obabel/` |

All binaries are downloaded from the [DockSuiteX_Binaries](https://github.com/MangalamGSinha/DockSuiteX_Binaries) repository on first import.
