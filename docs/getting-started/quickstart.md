# Quick Start Guide

Get started with DockSuiteX in minutes! This guide walks you through a complete docking workflow from fetching molecules to analyzing results.

## Complete Workflow Example

```python
from docksuitex import Protein, Ligand, PocketFinder, VinaDocking
from docksuitex.utils import clean_temp_folder, fetch_pdb, fetch_sdf, parse_vina_log_to_csv

# Clean temporary folder from previous runs
clean_temp_folder()

# 1. Fetch & prepare protein
protein_file = fetch_pdb("1HVR")  # HIV-1 Protease
prot = Protein(protein_file)
prot.prepare()

# 2. Fetch & prepare ligand
ligand_file = fetch_sdf(2244)  # Aspirin from PubChem
lig = Ligand(ligand_file)
lig.prepare(minimize="mmff94")

# 3. Predict binding pockets
finder = PocketFinder(prot)
pockets = finder.run()
center = pockets[0]['center']  # Use the top-ranked pocket

# 4. Run docking (using Vina)
vina = VinaDocking(
    receptor=prot,
    ligand=lig,
    grid_center=center,
    grid_size=(20, 20, 20),
    exhaustiveness=16
)
vina.run()
vina.save_results("vina_results")

# 5. Parse and combine results from multiple runs
parse_vina_log_to_csv("vina_results", "vina_results/vina_summary.csv")
```

## What This Does

1. **Fetches molecules**: Downloads HIV-1 Protease structure from RCSB PDB and Aspirin from PubChem
2. **Prepares protein**: Fixes missing residues, adds hydrogens, assigns charges, converts to PDBQT
3. **Prepares ligand**: Generates 3D coordinates, minimizes energy with MMFF94, converts to PDBQT
4. **Predicts pockets**: Uses P2Rank machine learning to find potential binding sites
5. **Runs docking**: Performs AutoDock Vina docking with the top-ranked pocket
6. **Parses results**: Extracts binding affinities and poses into a CSV file

## Next Steps

- Check out the [Examples](examples.md) page for more use cases
- Explore the [API Reference](../api/protein.md) for detailed documentation
- View results in Jupyter Notebook using the visualization utilities

## Visualization (Jupyter Only)

If you're working in Jupyter Notebook, you can visualize results:

```python
# Visualize the protein
prot.view()

# Visualize the ligand
lig.view()

# Visualize docking results with interactive controls
vina.view_results()
```
