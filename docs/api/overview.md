# API Overview

## Project Structure

The DockSuiteX package is organized into core modules, batch processing modules, and utility functions:

```
docksuitex/
├── protein.py        # Protein preparation (PDBFixer, Open Babel, AutoDockTools)
├── ligand.py         # Ligand preparation (Open Babel, AutoDockTools)
├── pocket_finder.py  # Pocket detection with P2Rank
├── autodock4.py      # AutoDock4 docking wrapper
├── vina.py           # AutoDock Vina docking wrapper
│ 
├── batch_docking/             # Batch processing modules
│   ├── batch_protein.py       # Parallel batch protein preparation
│   ├── batch_ligand.py        # Parallel batch ligand preparation
│   ├── batch_pocket_finder.py # Parallel batch pocket detection
│   ├── batch_autodock4.py     # Parallel batch docking with AutoDock4
│   └── batch_vina.py          # Parallel batch docking with AutoDock Vina
│ 
├── utils/
│   ├── fetcher.py      # Fetch PDB (RCSB) & SDF (PubChem)
│   ├── viewer.py       # NGLView visualization
│   ├── parser.py       # Parse logs to CSV summaries
│   └── cleaner.py      # Delete bin folder
│ 
└── bin/                # Auto-downloaded on first run
    ├── mgltools/       # MGLTools(AutoDockTools) binaries and scripts
    ├── obabel/         # Open Babel executables
    ├── vina/           # AutoDock Vina executable
    ├── p2rank/         # P2Rank executable and scripts
    └── autodock/       # AutoDock4 & AutoGrid executables
```

## Module Categories

### Core Modules

Fundamental classes for molecular docking workflows:

- **[Protein](protein.md)** - Protein structure preparation and conversion to PDBQT 
- **[Ligand](ligand.md)** - Ligand preparation with energy minimization options and conversion to PDBQT
- **[PocketFinder](pocket-finder.md)** - Binding pocket prediction using P2Rank
- **[VinaDocking](vina.md)** - AutoDock Vina docking interface
- **[AD4Docking](autodock4.md)** - AutoDock4 docking interface

### Batch Processing

High-throughput docking with parallel execution:

- **[BatchProtein](batch-protein.md)** - Batch protein preparation for multiple receptors
- **[BatchLigand](batch-ligand.md)** - Batch ligand preparation for multiple ligands
- **[BatchPocketFinder](batch-pocket-finder.md)** - Batch pocket prediction for multiple receptors
- **[BatchVinaDocking](batch-vina.md)** - Parallel Vina docking for multiple receptors, ligands and pockets
- **[BatchAD4Docking](batch-autodock4.md)** - Parallel AutoDock4 docking for multiple receptors, ligands and pockets



### Utilities

Helper functions for common tasks:

- **[Viewer](viewer.md)** - 3D molecular visualization with NGLView
- **[Parser](parser.md)** - Parse docking outputs to CSV format
- **[Fetcher](fetcher.md)** - Download structures from RCSB PDB and PubChem
- **[Cleaner](cleaner.md)** - Delete binaries, useful when you want to re-download fresh binaries.

## Binary Dependencies

DockSuiteX automatically downloads and manages the following third-party tools:

| Tool | Purpose | Location |
|------|---------|----------|
| MGLTools | Protein/ligand preparation scripts | `bin/mgltools/` |
| AutoDock Vina | Fast gradient-based docking | `bin/vina/` |
| AutoDock4 | Classic genetic algorithm docking | `bin/autodock/` |
| P2Rank | Machine learning pocket prediction | `bin/p2rank/` |
| Open Babel | Format conversion and energy minimization | `bin/obabel/` |

All binaries are downloaded from the [DockSuiteX_Binaries](https://github.com/MangalamGSinha/DockSuiteX_Binaries) repository on first import.
