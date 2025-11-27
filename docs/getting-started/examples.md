# Examples

Comprehensive examples demonstrating all major features of DockSuiteX.

## Protein Preparation

Prepare a protein structure for docking:

```python
from docksuitex import Protein

# Load protein
prot = Protein("protein.pdb")

# Prepare protein for docking
prot.prepare(
    fix_pdb=True,
    remove_heterogens=True,
    remove_water=True,
    add_hydrogens=True,
    add_charges=True
)

# Save the final PDBQT file
prot.save_pdbqt("protein_prepared.pdbqt")

# Visualize protein (in Jupyter Notebook)
prot.view()
```

### Preserving Specific Charge Types

If your protein contains metal ions or other atoms whose charges should be preserved:

```python
prot.prepare(
    add_charges=True,
    preserve_charge_types=["Zn", "Fe", "Mg"]  # Preserve charges for these atoms
)
```

---

## Ligand Preparation

Prepare a ligand with energy minimization:

```python
from docksuitex import Ligand

# Load ligand
lig = Ligand("ligand.sdf")

# Prepare ligand with MMFF94 minimization
lig.prepare(
    minimize="mmff94",  # Options: "mmff94", "mmff94s", "uff", "gaff", or None
    remove_water=True,
    add_hydrogens=True,
    add_charges=True
)

# Save the final PDBQT file
lig.save_pdbqt("ligand_prepared.pdbqt")

# Visualize ligand (in Jupyter Notebook)
lig.view()
```

---

## Pocket Detection

Predict binding sites using P2Rank:

```python
from docksuitex import PocketFinder

# Initialize pocket finder
pf = PocketFinder("protein.pdb")

# Run P2Rank to predict pockets
pockets = pf.run()

# Display predicted pockets
for pocket in pockets:
    print(f"Rank {pocket['rank']}: Center at {pocket['center']}")

# Save full P2Rank output folder
pf.save_report("p2rank_report")
```

---

## AutoDock Vina Docking

### Basic Vina Docking

```python
from docksuitex import VinaDocking

# Initialize Vina docking
vina = VinaDocking(
    receptor="protein_prepared.pdbqt",
    ligand="ligand_prepared.pdbqt",
    grid_center=(10.0, 12.5, 8.0),
    grid_size=(20, 20, 20),
    exhaustiveness=8,
    num_modes=9
)

# Run docking
vina.run()

# Visualize results (in Jupyter Notebook)
vina.view_results()

# Save results
vina.save_results("vina_docking")
```

### Using Protein/Ligand Objects

You can pass `Protein` and `Ligand` objects directly:

```python
from docksuitex import Protein, Ligand, VinaDocking

prot = Protein("protein.pdb")
prot.prepare()

lig = Ligand("ligand.sdf")
lig.prepare(minimize="mmff94")

vina = VinaDocking(
    receptor=prot,  # Pass Protein object
    ligand=lig,     # Pass Ligand object
    grid_center=(10.0, 12.5, 8.0)
)
vina.run()
```

---

## AutoDock4 Docking

```python
from docksuitex import AD4Docking

# Initialize AutoDock4 docking
ad4 = AD4Docking(
    receptor="protein_prepared.pdbqt",
    ligand="ligand_prepared.pdbqt",
    grid_center=(10.0, 12.5, 8.0),
    grid_size=(60, 60, 60),
    spacing=0.375,
    ga_run=10,
    ga_num_evals=2500000
)

# Run AutoGrid + AutoDock
ad4.run()

# Visualize results (in Jupyter Notebook)
ad4.view_results()

# Save results
ad4.save_results("ad4_docking")
```

---

## Batch Docking with Vina

Run multiple ligands against multiple binding pockets in parallel:

```python
from docksuitex import BatchVinaDocking

# Input
receptor = "protein_prepared.pdbqt"
ligands = ["lig1_prepared.pdbqt", "lig2_prepared.pdbqt", "lig3_prepared.pdbqt"]
centers = [
    (10.0, 12.5, 8.0),
    (-8.2, 14.6, 25.3),
    (-12.2, -10.1, 8.3)
]

# Initialize batch docking
batch = BatchVinaDocking(
    receptor=receptor,
    ligand_list=ligands,
    center_list=centers,
    grid_size=(20, 20, 20),
    exhaustiveness=8,
    num_modes=9,
    seed=42
)

# Run all jobs in parallel
results = batch.run_all(save_to="batch_vina")

# Check results
for (lig, center), res in results.items():
    print(f"Ligand {lig} at center {center} → {res}")
```

---

## Batch Docking with AutoDock4

```python
from docksuitex import BatchAD4Docking

# Input
receptor = "protein_prepared.pdbqt"
ligands = ["lig1_prepared.pdbqt", "lig2_prepared.pdbqt"]
centers = [(10.0, 12.5, 8.0), (-8.2, 14.6, 25.3)]

# Initialize batch docking
batch = BatchAD4Docking(
    receptor=receptor,
    ligand_list=ligands,
    center_list=centers,
    grid_size=(60, 60, 60),
    ga_run=10
)

# Run all jobs in parallel
results = batch.run_all(save_to="batch_ad4")

for (lig, center), res in results.items():
    print(f"Ligand {lig} at center {center} → {res}")
```

---

## Parsing Docking Results

### Parse Vina Logs to CSV

```python
from docksuitex.utils import parse_vina_log_to_csv

# Parse all Vina log files in a directory
df = parse_vina_log_to_csv(
    results_dir="vina_docking",
    output_csv="vina_summary.csv"
)

print(df.head())
```

### Parse AutoDock4 DLG Files to CSV

```python
from docksuitex.utils import parse_ad4_dlg_to_csv

# Parse all AutoDock4 DLG files in a directory
df = parse_ad4_dlg_to_csv(
    results_dir="ad4_docking",
    output_csv="ad4_summary.csv"
)

print(df.head())
```

---

## Fetching Molecules

### Fetch Protein from RCSB PDB

```python
from docksuitex.utils import fetch_pdb

# Download protein structure by PDB ID
pdb_file = fetch_pdb("1UBQ", save_to="pdbs")
print(f"Downloaded: {pdb_file}")
```

### Fetch Ligand from PubChem

```python
from docksuitex.utils import fetch_sdf

# Download ligand by PubChem CID
sdf_file = fetch_sdf(2244, save_to="ligands")  # Aspirin
print(f"Downloaded: {sdf_file}")
```

---

## Visualization (Jupyter Notebook)

### View Single Molecule

```python
from docksuitex.utils import view_molecule

# View a protein or ligand
view_molecule("protein.pdbqt")
```

### View Docking Results

Interactive viewer with pose controls:

```python
from docksuitex.utils import view_results

# Visualize protein with docking poses
# Features: step through poses, play/pause animation, show all poses
view_results(
    protein_file="protein.pdb",
    ligand_file="ligand_poses.pdbqt"
)
```

---

## Utility Functions

### Clean Temporary Folder

```python
from docksuitex.utils import clean_temp_folder

# Clean up temporary files from previous runs
clean_temp_folder()
```

### Delete Binaries

```python
from docksuitex.utils import delete_binaries

# Delete the bin directory to re-download fresh binaries
delete_binaries()
```

---

## Complete End-to-End Workflow

Putting it all together:

```python
from docksuitex import Protein, Ligand, PocketFinder, VinaDocking
from docksuitex.utils import (
    clean_temp_folder,
    fetch_pdb,
    fetch_sdf,
    parse_vina_log_to_csv
)

# Clean up
clean_temp_folder()

# 1. Fetch molecules
protein_file = fetch_pdb("1HVR")  # HIV-1 Protease
ligand_file = fetch_sdf(2244)     # Aspirin

# 2. Prepare protein
prot = Protein(protein_file)
prot.prepare()

# 3. Prepare ligand
lig = Ligand(ligand_file)
lig.prepare(minimize="mmff94")

# 4. Find binding pockets
finder = PocketFinder(prot)
pockets = finder.run()
center = pockets[0]['center']

# 5. Run docking
vina = VinaDocking(
    receptor=prot,
    ligand=lig,
    grid_center=center,
    grid_size=(20, 20, 20),
    exhaustiveness=16
)
vina.run()
vina.save_results("vina_results")

# 6. Parse results
df = parse_vina_log_to_csv("vina_results", "vina_results/summary.csv")
print(df)
```

---

## More Examples

For additional runnable examples, check out the [examples folder on GitHub](https://github.com/MangalamGSinha/DockSuiteX/tree/main/examples).
