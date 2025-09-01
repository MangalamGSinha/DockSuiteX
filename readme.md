# 🧬 DockSuiteX – Automated Protein–Ligand Docking with AutoDock4 \& Vina

DockSuiteX is a Python package that automates the **end-to-end workflow of molecular docking** using [AutoDock4](http://autodock.scripps.edu/) and [AutoDock Vina](http://vina.scripps.edu/).
It integrates **protein and ligand preparation**, **binding site prediction**, **docking execution**, **results parsing**, and **3D visualization** into a seamless pipeline.

---

## ✨ Features

* **Protein Preparation**

  * Input formats: `.pdb`, `.mol2`, `.sdf`, `.pdbqt`, `.cif`, `.ent`,  `.xyz`
  * Fix missing residues/atoms with **PDBFixer**
  * Remove water, remove heterogens, add charges and add polar hydrogens
  * Convert to `.pdbqt` using **AutoDockTools**
* **Ligand Preparation**

  * Input formats: `.mol2`, `.sdf`, `.pdb`, `.mol`, `.smi`
  * Automatic 3D generation and optional energy minimization (**MMFF94**, **MMFF94s**, **UFF**, **GAFF**) with **obabel**
  * Remove water, add charges and add polar hydrogens
  * Convert to `.pdbqt` using **AutoDockTools**
* **Pocket Detection**

  * Predict binding sites using **P2Rank**
  * Get ranked pockets and 3D coordinates of centers
* **Docking**

  * **AutoDock4**

    * Genetic algorithm with customizable docking parameters
  * **AutoDock Vina**

    * Fast and efficient search with adjustable exhaustiveness and CPU usage
  * **Parallel Batch Docking**

    * Run multiple ligands against multiple predicted binding pocket centers in parallel
    * Supports both **AutoDock4** and **AutoDock Vina**
    * Fully utilizes available CPU cores for large-scale virtual screening
* **Utilities**

  * Fetch proteins from **RCSB PDB** (`.pdb`) and Fetch ligands from **PubChem** (`.sdf`)
  * Parse multiple docking logs and combine results into a single CSV file.
  * Visualize molecules in Jupyter with **NGLView**
  * Visualize docking results (multiple poses, step-through, play/pause, show all).

---

## 📦 Installation

DockSuiteX currently works **only on Windows**.

#### Prerequisites

  **Java (JDK/JRE)** is required for **P2Rank**.
  You can download the latest LTS release from [Adoptium Temurin](https://adoptium.net/temurin/releases).
  After installation, make sure Java is added to your **System PATH** so it can be accessed from the command line.
  You can verify by running:

```bash
java -version
```

#### Install from PyPI

```bash
pip install docksuitex
```

#### Install from GitHub

```bash
git clone https://github.com/MangalamGSinha/DockSuiteX.git
cd DockSuiteX
pip install .
```

---

## 🚀 Quickstart Example

```python
from docksuitex import Protein, Ligand, PocketFinder, VinaDocking
from docksuitex.utils import clean_temp_folder, fetch_pdb, fetch_sdf, parse_vina_log_to_csv
clean_temp_folder()

# 1. Fetch & prepare protein
protein_file = fetch_pdb("1HVR")
prot = Protein(protein_file)
prot.prepare()

# 1. Fetch & prepare ligand
ligand_file = fetch_sdf(2244)  # Aspirin
lig = Ligand(ligand_file)
lig.prepare(minimize="mmff94")

# 3. Predict binding pockets
finder = PocketFinder(prot)
pockets = finder.run()
grid_center = pockets[0]['center'] #First Pocket

# 4. Run docking (using Vina)
vina = VinaDocking(
    receptor=prot,
    ligand=lig,
    grid_center=grid_center,
    grid_size=(20,20,20),
    exhaustiveness=16
)
vina.run()
vina.save_results(f"vina_results")

# 5. Parse and combine results from multiple runs
parse_vina_log_to_csv("vina_results", "vina_results/vina_summary.csv")

```

Detailed, runnable examples are available in the [examples](https://github.com/MangalamGSinha/DockSuiteX/tree/main/examples) folder.

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
│   └── cleaner.py      # Reset temp folders
│
└── bin/                # Auto-downloaded on first run
    ├── mgltools/       # MGLTools binaries and scripts
    ├── obabel/         # Open Babel executables
    ├── vina/           # AutoDock Vina executable (vina.exe)
    ├── p2rank/         # P2Rank executable and scripts
    └── autodock/       # AutoDock4 & AutoGrid executables (autodock4.exe, autogrid4.exe)

```

---

## 🧩 Module Documentation

### 1. `protein`

#### Class: `Protein`

| Parameter     | Type     | Default | Description                                                                                           |
| ------------- | -------- | ------- | ----------------------------------------------------------------------------------------------------- |
| `file_path` | str/Path | —      | Path to input protein file (`.pdb`, `.mol2`, `.sdf`, `.pdbqt`, `.cif`, `.ent`, `.xyz`). |

#### Method: `prepare()`

| Parameter                 | Type                 | Default | Description                                                                                                                     |
| ------------------------- | -------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `fix_pdb`               | bool                 | True    | Run PDBFixer cleanup (fix missing residues/atoms, replace non-standard residues).                                               |
| `remove_heterogens`     | bool                 | True    | Remove heterogens (non-protein residues, ligands, etc.) in PDBFixer step.                                                       |
| `remove_water`          | bool                 | True    | Remove water molecules.                                                                                                         |
| `add_hydrogens`         | bool                 | True    | Add polar hydrogens.                                                                                                            |
| `add_charges`           | bool                 | True    | Assign Gasteiger charges; if False, input charges are preserved.                                                                |
| `preserve_charge_types` | list\[str], optional | None    | Atom types (e.g.,`["Zn", "Fe"]`) whose charges are preserved; others get Gasteiger charges; ignored if `add_charges=False`. |

#### Method: `view()`

| Returns                                                                          |
| -------------------------------------------------------------------------------- |
| `nglview.NGLWidget` – Interactive 3D view of the protein in Jupyter Notebook. |

#### Method: `save_pdbqt()`

| Parameter   | Type     | Default | Description                                                                  |
| ----------- | -------- | ------- | ---------------------------------------------------------------------------- |
| `save_to` | str/Path | "."     | Destination file or directory where the prepared PDBQT file should be saved. |

**Returns:**
Absolute path (`Path`) to the saved PDBQT file.

Example:

```python
from docksuitex import Protein

# Load protein
prot = Protein("protein.pdb")

# Prepare protein for docking
prot.prepare(fix_pdb=True, add_hydrogens=True)

# Visualize protein (in Jupyter)
prot.view()

# Save the final PDBQT file
prot.save_pdbqt("protein_prepared.pdbqt")

```

---

### 2. `ligand`

#### Class: `Ligand`

| Parameter     | Type     | Default | Description                                                                    |
| ------------- | -------- | ------- | ------------------------------------------------------------------------------ |
| `file_path` | str/Path | —      | Path to input ligand file (`.mol2`, `.sdf`, `.pdb`, `.mol`, `.smi`). |

#### Method: `prepare()`

| Parameter                 | Type                 | Default | Description                                                                                                                     |
| ------------------------- | -------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `minimize`              | str, optional        | None    | Forcefield for energy minimization (`"mmff94"`, `"mmff94s"`, `"uff"`, `"gaff"`). If None, no minimization is performed. |
| `remove_water`          | bool                 | True    | Remove water molecules.                                                                                                         |
| `add_hydrogens`         | bool                 | True    | Add polar hydrogens.                                                                                                            |
| `add_charges`           | bool                 | True    | Assign Gasteiger charges; if False, input charges are preserved.                                                                |
| `preserve_charge_types` | list\[str], optional | None    | Atom types (e.g.,`["Zn", "Fe"]`) whose charges are preserved; others get Gasteiger charges; ignored if `add_charges=False`. |

#### Method: `view()`

| Returns                                                                         |
| ------------------------------------------------------------------------------- |
| `nglview.NGLWidget` – Interactive 3D view of the ligand in Jupyter Notebook. |

#### Method: `save_pdbqt()`

| Parameter   | Type     | Default | Description                                                                  |
| ----------- | -------- | ------- | ---------------------------------------------------------------------------- |
| `save_to` | str/Path | "."     | Destination file or directory where the prepared PDBQT file should be saved. |

**Returns:**
Absolute path (`Path`) to the saved PDBQT file.

Example:

```python
from docksuitex import Ligand

# Load ligand
lig = Ligand("ligand.sdf")

# Prepare ligand for docking
lig.prepare(minimize="mmff94", remove_water=True, add_hydrogens=True)

# Visualize ligand (in Jupyter)
lig.view()

# Save the final PDBQT file
lig.save_pdbqt("ligand_prepared.pdbqt")
```

---

### 3. `pocket_finder`

#### Class: `PocketFinder`

| Parameter    | Type             | Default | Description                                                             |
| ------------ | ---------------- | ------- | ----------------------------------------------------------------------- |
| `receptor` | str/Path/Protein | —      | Input receptor. Can be a PDB file, PDBQT file, or a `Protein` object. |

#### Method: `run()`

| Parameter | Type | Default | Description                                                                                                               |
| --------- | ---- | ------- | ------------------------------------------------------------------------------------------------------------------------- |
| —        | —   | —      | Runs P2Rank to predict ligand-binding pockets in the protein. Returns a list of pockets with rank and center coordinates. |

#### Method: `save_report()`

| Parameter   | Type     | Default          | Description                                              |
| ----------- | -------- | ---------------- | -------------------------------------------------------- |
| `save_to` | str/Path | ./p2rank_results | Directory where the P2Rank output folder will be copied. |

**Returns:**
Absolute path (`Path`) to the directory containing the copied P2Rank report.

Example:

```python
from docksuitex import PocketFinder

# Initialize pocket finder
pf = PocketFinder("protein.pdb")

# Run P2Rank to predict pockets
pockets = pf.run()

for pocket in pockets:
    print(f"Rank {pocket['rank']}: Center at {pocket['center']}")

# Save full P2Rank output folder
pf.save_report("p2rank_report")
```

---

### 4. `vina`

#### Class: `VinaDocking`

| Parameter          | Type                        | Default         | Description                                                                                             |
| ------------------ | --------------------------- | --------------- | ------------------------------------------------------------------------------------------------------- |
| `receptor`       | str/Path/Protein            | —              | Input receptor. Can be a PDBQT file or a `Protein` object.                                            |
| `ligand`         | str/Path/Protein            | —              | Input ligand. Can be a PDBQT file or a `Ligand` object.                                               |
| `grid_center`    | tuple\[float, float, float] | —              | Center coordinates (x, y, z) of the docking grid in Ångström.                                         |
| `grid_size`      | tuple\[int, int, int]       | (20, 20, 20)    | Physical length of the grid box along (x, y, z) in Ångström. Spacing is fixed internally at 0.375 Å. |
| `exhaustiveness` | int                         | 8               | How exhaustively Vina searches conformational space.                                                    |
| `num_modes`      | int                         | 9               | Maximum number of binding modes to output.                                                              |
| `cpu`            | int                         | os.cpu\_count() | Number of CPU cores to use.                                                                             |
| `verbosity`      | int                         | 1               | Level of console output (0=quiet, 1=normal, 2=verbose).                                                 |
| `seed`           | int, optional               | None            | Random seed for reproducibility.                                                                        |

#### Method: `run()`

| Parameter | Type | Default | Description                                                                                  |
| --------- | ---- | ------- | -------------------------------------------------------------------------------------------- |
| —        | —   | —      | Runs AutoDock Vina with the given parameters. Outputs docking results to a temporary folder. |

#### Method: `view_results()`

| Parameter | Type | Default | Description                                                                                                             |
| --------- | ---- | ------- | ----------------------------------------------------------------------------------------------------------------------- |
| —        | —   | —      | Opens an**interactive 3D visualization** of receptor + docked ligand using **NGLView** in Jupyter Notebook. |

#### Method: `save_results()`

| Parameter   | Type     | Default        | Description                                                    |
| ----------- | -------- | -------------- | -------------------------------------------------------------- |
| `save_to` | str/Path | ./vina_results | Directory where docking results (.pdbqt, .log) will be copied. |

**Returns:**
Absolute path (`Path`) to the directory containing all saved Vina results.

Example:

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

# Visualize ligand (in Jupyter)
vina.view_results()

# Save results
vina.save_results("vina_docking")
```

---

### 5. `autodock4`

#### Class: `AD4Docking`

| Parameter              | Type                        | Default      | Description                                                                                |
| ---------------------- | --------------------------- | ------------ | ------------------------------------------------------------------------------------------ |
| `receptor`           | str/Path/Protein            | —           | Input receptor. Can be a PDBQT file or a `Protein` object.                               |
| `ligand`             | str/Path/Protein            | —           | Input ligand. Can be a PDBQT file or a `Ligand` object.                                  |
| `grid_center`        | tuple\[float, float, float] | —           | Center coordinates (x, y, z) of the docking grid in Ångström.                            |
| `grid_size`          | tuple\[int, int, int]       | (60, 60, 60) | Number of grid points along (x, y, z) axes. Effective box size =`grid_size × spacing`. |
| `spacing`            | float                       | 0.375        | Distance between adjacent grid points (Å). Controls the resolution of the grid.           |
| `dielectric`         | float                       | -0.1465      | Dielectric constant for electrostatics.                                                    |
| `smooth`             | float                       | 0.5          | Smoothing factor for potential maps.                                                       |
| `ga_pop_size`        | int                         | 150          | Genetic algorithm population size.                                                         |
| `ga_num_evals`       | int                         | 2,500,000    | Maximum number of energy evaluations in GA.                                                |
| `ga_num_generations` | int                         | 27,000       | Maximum number of generations in GA.                                                       |
| `ga_elitism`         | int                         | 1            | Number of top individuals preserved during GA.                                             |
| `ga_mutation_rate`   | float                       | 0.02         | Probability of mutation in GA.                                                             |
| `ga_crossover_rate`  | float                       | 0.8          | Probability of crossover in GA.                                                            |
| `ga_run`             | int                         | 10           | Number of independent GA runs.                                                             |
| `rmstol`             | float                       | 2.0          | RMSD tolerance for clustering docking results.                                             |

#### Method: `run()`

| Parameter | Type | Default | Description                                                                                                    |
| --------- | ---- | ------- | -------------------------------------------------------------------------------------------------------------- |
| —        | —   | —      | Runs AutoGrid and AutoDock docking with the given parameters. Generates output files in a temporary directory. |

#### Method: `view_results()`

| Parameter | Type | Default | Description                                                                                                             |
| --------- | ---- | ------- | ----------------------------------------------------------------------------------------------------------------------- |
| —        | —   | —      | Opens an**interactive 3D visualization** of receptor + docked ligand using **NGLView** in Jupyter Notebook. |

#### Method: `save_results()`

| Parameter   | Type     | Default       | Description                                                                      |
| ----------- | -------- | ------------- | -------------------------------------------------------------------------------- |
| `save_to` | str/Path | ./ad4_results | Directory where the docking results (DLG, GPF, DPF, PDBQT files) will be copied. |

**Returns:**
Absolute path (`Path`) to the directory containing all saved AutoDock4 results.

Example:

```python
from docksuitex import AD4Docking

# Initialize docking
ad4 = AD4Docking(
    receptor="protein_prepared.pdbqt",
    ligand="ligand_prepared.pdbqt",
    grid_center=(10.0, 12.5, 8.0),
    grid_size=(60, 60, 60),
    ga_run=10
)

# Run AutoGrid + AutoDock
ad4.run()

# Visualize results (in Jupyter)
ad4.view_results()

# Save results
ad4.save_results("ad4_docking")
```

---

### 6. `batch_vina`

#### Class: `BatchVinaDocking`

| Parameter          | Type                                   | Default      | Description                                                                             |
| ------------------ | -------------------------------------- | ------------ | --------------------------------------------------------------------------------------- |
| `receptor`       | str/Path                               | —           | Path to receptor PDBQT file.                                                            |
| `ligand_list`    | Sequence\[str \| Path]                 | —           | List of ligand PDBQT files.                                                             |
| `center_list`    | Sequence\[tuple\[float, float, float]] | —           | List of docking box centers `(x, y, z)` in Ångström.                                |
| `grid_size`      | tuple\[int, int, int]                  | (20, 20, 20) | Dimensions of the docking search box (Å). Spacing fixed internally at 0.375 Å.        |
| `exhaustiveness` | int                                    | 8            | How exhaustively Vina searches conformational space. Higher = slower but more accurate. |
| `num_modes`      | int                                    | 9            | Maximum number of binding modes per ligand.                                             |
| `verbosity`      | int                                    | 1            | Console output level (0=quiet, 1=normal, 2=verbose).                                    |
| `seed`           | int, optional                          | None         | Random seed for reproducibility. If `None`, Vina chooses automatically.               |

#### Method: `run_all()`

| Parameter   | Type     | Default              | Description                                                       |
| ----------- | -------- | -------------------- | ----------------------------------------------------------------- |
| `cpu`     | int      | os.cpu_count()       | Number of CPU cores to use. Distributes CPUs across jobs.         |
| `save_to` | str/Path | ./batch_vina_results | Directory where all docking results (.pdbqt, .log) will be saved. |

**Returns:**
Dictionary mapping `(ligand_name, center)` → absolute path (`Path`) of result file if docking succeeded, or error message (`str`) if failed.

Example:

```python
from docksuitex import BatchVinaDocking

# Input
receptor = "protein_prepared.pdbqt" #Receptor path
ligands = ["lig1_prepared.pdbqt", "lig2_prepared.pdbqt"] #Ligand paths list
centers = [(10.0, 12.5, 8.0), (-8.2, 14.6, 25.3), (-12.2, -10.1, 8.3)] #Pocket center list

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

for (lig, center), res in results.items():
    print(f"Ligand {lig} at center {center} → results in {res}")
```

### 7. `batch_autodock4`

#### Class: `BatchAD4Docking`

| Parameter              | Type                                   | Default      | Description                                                                               |
| ---------------------- | -------------------------------------- | ------------ | ----------------------------------------------------------------------------------------- |
| `receptor`           | str/Path                               | —           | Path to receptor PDBQT file.                                                              |
| `ligand_list`        | Sequence\[str/Path]                    | —           | List of ligand PDBQT files.                                                               |
| `center_list`        | Sequence\[tuple\[float, float, float]] | —           | List of docking box centers (grid centers in Å).                                         |
| `grid_size`          | tuple\[int, int, int]                  | (60, 60, 60) | Number of grid points along (x, y, z) axes. Effective box size =`grid_size × spacing`. |
| `spacing`            | float                                  | 0.375        | Distance between adjacent grid points (Å). Controls grid resolution.                     |
| `dielectric`         | float                                  | -0.1465      | Dielectric constant for electrostatics.                                                   |
| `smooth`             | float                                  | 0.5          | Smoothing factor for potential maps.                                                      |
| `ga_pop_size`        | int                                    | 150          | Genetic algorithm population size.                                                        |
| `ga_num_evals`       | int                                    | 2,500,000    | Maximum number of energy evaluations in GA.                                               |
| `ga_num_generations` | int                                    | 27,000       | Maximum number of generations in GA.                                                      |
| `ga_elitism`         | int                                    | 1            | Number of top individuals preserved during GA.                                            |
| `ga_mutation_rate`   | float                                  | 0.02         | Probability of mutation in GA.                                                            |
| `ga_crossover_rate`  | float                                  | 0.8          | Probability of crossover in GA.                                                           |
| `ga_run`             | int                                    | 10           | Number of independent GA runs.                                                            |
| `rmstol`             | float                                  | 2.0          | RMSD tolerance for clustering docking results.                                            |

#### Method: `run_all()`

| Parameter   | Type     | Default             | Description                                            |
| ----------- | -------- | ------------------- | ------------------------------------------------------ |
| `cpu`     | int      | os.cpu_count()      | Number of CPU cores to use. Defaults to all available. |
| `save_to` | str/Path | ./batch_ad4_results | Directory where all docking results will be stored.    |

**Returns:**
Dictionary mapping `(ligand_name, center)` → absolute path (`Path`) of result file if docking succeeded, or error message (`str`) if failed.

---

Example:

```python
from docksuitex import BatchAD4Docking

# Input
receptor = "protein_prepared.pdbqt" #Receptor path
ligands = ["lig1_prepared.pdbqt", "lig2_prepared.pdbqt"] #Ligand paths list
centers = [(10.0, 12.5, 8.0), (-8.2, 14.6, 25.3), (-12.2, -10.1, 8.3)] #Pocket center list

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
    print(f"Ligand {lig} at center {center} → results in {res}")
```

---

### 8. `utils/parser`

#### Method: `parse_vina_log_to_csv()`

| Parameter       | Type | Default            | Description                                                                        |
| --------------- | ---- | ------------------ | ---------------------------------------------------------------------------------- |
| `results_dir` | str  | —                 | Parent directory containing docking result folders with Vina log files (`.txt`). |
| `output_csv`  | str  | "vina_summary.csv" | Path to save the output CSV file.                                                  |

**Returns:**
`pandas.DataFrame` containing parsed docking scores and poses.

Example:

```python
from docksuitex.utils import parse_vina_log_to_csv

df = parse_vina_log_to_csv(results_dir="vina_docking", output_csv="vina_summary.csv")
print(df.head())
```

#### Method: `parse_ad4_dlg_to_csv()`

| Parameter       | Type       | Default         | Description                                                              |
| --------------- | ---------- | --------------- | ------------------------------------------------------------------------ |
| `results_dir` | str / Path | —              | Parent directory containing docking result folders with `results.dlg`. |
| `output_csv`  | str        | ad4_summary.csv | Path to save the output CSV file.                                        |

**Returns:**
`pandas.DataFrame` containing parsed docking scores and poses.

Example:

```python
from docksuitex.utils import parse_ad4_dlg_to_csv

df = parse_ad4_dlg_to_csv(results_dir="ad4_docking", output_csv="ad4_summary.csv")
print(df.head())
```

---

### 9. `utils/fetcher`

#### **Method:** `fetch_pdb()`

| Parameter   | Type       | Default | Description                                    |
| ----------- | ---------- | ------- | ---------------------------------------------- |
| `pdbid`   | str        | —      | 4-character alphanumeric PDB ID (e.g., '1UBQ') |
| `save_to` | str / Path | "."     | Directory to save the `.pdb` file            |

**Returns:**
Absolute path (`Path`) to the downloaded `.pdb` file.

Example:

```python
from docksuitex.utils import fetch_pdb

pdb_file = fetch_pdb("1UBQ", save_to="pdbs")
```

#### Method: `fetch_sdf()`

| Parameter   | Type       | Default | Description                                       |
| ----------- | ---------- | ------- | ------------------------------------------------- |
| `cid`     | str / int  | —      | PubChem Compound ID (CID), e.g., 2244 for Aspirin |
| `save_to` | str / Path | "."     | Directory to save the `.sdf` file               |

**Returns:**
Absolute path (`Path`) to the downloaded `.sdf` file.

Example:

```python
from docksuitex.utils import fetch_sdf

sdf_file = fetch_sdf(2244, save_to="sdfs")
```

---

### 10. `utils/viewer`

#### Method: `view_molecule()`

| Parameter     | Type     | Description                                                                |
| ------------- | -------- | -------------------------------------------------------------------------- |
| `file_path` | str/Path | Path to the molecular file (`.pdb`, `.pdbqt`, `.mol2`, or `.sdf`). |

**Returns:**
An interactive NGLView widget.

Example:

```python
from docksuitex.utils import view_molecule

# View a protein or ligand
view_molecule(file_path="protein.pdbqt")
```

#### Method: `view_results()`

| Parameter        | Type     | Description                                                           |
| ---------------- | -------- | --------------------------------------------------------------------- |
| `protein_file` | str/Path | Path to the receptor protein file (e.g.,`.pdb`).                    |
| `ligand_file`  | str/Path | Path to ligand docking results file (`.pdbqt` with multiple poses). |

**Returns:**
Displays visualization and interactive controls (step-through, play/pause, show all) directly in Jupyter Notebook.

Example:

```python
from docksuitex.utils.viewer import view_results

# Visualize protein with docking poses
view_results("protein.pdb", "ligand_poses.pdbqt")
```

---

### 11. `utils/cleaner`

#### Method: `clean_temp_folder()`

```python
from docksuitex.utils import clean_temp_folder

clean_temp_folder()
```

#### Method: `delete_binaries()`

```python
from docksuitex.utils import delete_binaries

# Deletes the `bin` directory.
# Useful when you want to re-download fresh binaries.
delete_binaries()
```

---

## 🙏 Acknowledgments

This package builds upon and automates workflows using:

* [AutoDock4 \& AutoDock Vina](http://autodock.scripps.edu/)
* [MGLTools](http://mgltools.scripps.edu/)
* [Open Babel](https://openbabel.org/)
* [PDBFixer](http://openmm.org/)
* [P2Rank](https://github.com/rdk/p2rank)
* [RCSB PDB](https://www.rcsb.org/) \& [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
* [NGLView](https://pypi.org/project/nglview/)

---

## 📜 License

This project is licensed under the GNU GPL v3 License - see the [LICENSE](https://github.com/MangalamGSinha/DockSuiteX/blob/main/LICENSE) file for details.
