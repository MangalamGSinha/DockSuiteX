# 🧬 DockSuiteX – Automated Protein–Ligand Docking with AutoDock4 \& Vina

DockSuiteX is a Python package that automates the **end-to-end workflow of molecular docking** using [AutoDock4](http://autodock.scripps.edu/) and [AutoDock Vina](http://vina.scripps.edu/).
It integrates **protein and ligand preparation**, **binding site prediction**, **docking execution**, **results parsing**, and **3D visualization** into a seamless pipeline.

---

## ✨ Features

* **Protein Preparation**

  * Input formats: `.pdb`, `.mol2`, `.sdf`, `.pdbqt`, `.cif`, `.ent`,  `.xyz`
  * Fix missing residues/atoms with **PDBFixer**
  * Remove water, add charges and add polar hydrogens
  * Convert to `.pdbqt` using **AutoDockTools**
* **Ligand Preparation**

  * Input formats: `.mol2`, `.sdf`, `.pdb`, `.mol`, `.smi`
  * Automatic 3D generation and optional energy minimization (**MMFF94**, **UFF**, **GAFF**) with **obabel**
  * Remove water, add charges and add polar hydrogens
  * Convert to `.pdbqt` using **AutoDockTools**
* **Pocket Detection**

  * Predict binding sites using **P2Rank**
  * Get ranked pockets and 3D coordinates of centers
* **Docking**

  * **AutoDock4** (genetic algorithm with customizable parameters)
  * **AutoDock Vina** (fast \& efficient search with flexible exhaustiveness/CPU usage)
* **Utilities**

  * Fetch proteins from **RCSB PDB** (`.pdb`) and Fetch ligands from **PubChem** (`.sdf`)
  * Parse docking logs → single CSV summary (`vina_summary.csv`, `ad4_summary.csv`)
  * Visualize molecules in Jupyter with **py3Dmol**
  * Clean temp folders between runs

---

## 📦 Installation

DockSuiteX currently works **only on Windows**.

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

# 1. Fetch protein & ligand
protein_file = fetch_pdb("1UBQ")
ligand_file = fetch_sdf(2244)  # Aspirin

# 2. Prepare protein & ligand
prot = Protein(protein_file)
prot.prepare()
prot.save_pdbqt()

lig = Ligand(ligand_file)
lig.prepare(minimize="mmff94")
lig.save_pdbqt()

# 3. Predict binding pockets
finder = PocketFinder(prot.pdb_path)
pockets = finder.run()
grid_center = pockets[0]['center'] #First Pocket

# 4. Run docking (using Vina)
vina = VinaDocking(
    receptor=prot.pdbqt_path,
    ligand=lig.pdbqt_path,
    grid_center=grid_center,
    grid_size=(20,20,20),
    exhaustiveness=16
)
vina.run()
vina.save_results(f"vina_results/{prot.pdbqt_path.stem}_{lig.pdbqt_path.stem}_pocket_1_docking")

# 5. Parse and combine results from multiple runs
parse_vina_log_to_csv("vina_results", "vina_results/vina_summary.csv")

```

Detailed, runnable examples are available in the [`examples/`](./examples) folder.

---

## 📂 Project Structure

```
docksuitex/
├── protein.py        # Protein preparation (PDBFixer, Open Babel, MGLTools)
├── ligand.py         # Ligand preparation (Open Babel, minimization, MGLTools)
├── pocket_finder.py   # Pocket detection with P2Rank
├── autodock4.py      # AutoDock4 docking wrapper
├── vina.py           # AutoDock Vina docking wrapper
├── utils/
│   ├── fetcher.py    # Fetch PDB (RCSB) & SDF (PubChem)
│   ├── viewer.py     # Py3Dmol visualization
│   ├── parser.py     # Parse logs to CSV summaries
│   └── cleaner.py    # Reset temp folders
|
└── bin/              # Will be auto-downloaded on first run
    ├── mgltools/     # MGLTools binaries and scripts
    ├── obabel/       # Open Babel executables
    ├── vina/         # AutoDock Vina executable (vina.exe)
    ├── p2rank_2.5/   # P2Rank executable and scripts
    └── AutoDock/     # AutoDock4 & AutoGrid executables (autodock4.exe, autogrid4.exe)

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

#### Method: `save_pdbqt()`

| Parameter     | Type     | Default | Description                                                                  |
| ------------- | -------- | ------- | ---------------------------------------------------------------------------- |
| `save_path` | str/Path | "."     | Destination file or directory where the prepared PDBQT file should be saved. |

Example:

```python
from docksuitex import Protein

# Load protein
prot = Protein("protein.pdb")

# Prepare protein for docking
prot.prepare(fix_pdb=True, add_hydrogens=True)

# Save the final PDBQT file
prot.save_pdbqt("outputs/protein_prepared.pdbqt")
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

#### Method: `save_pdbqt()`

| Parameter     | Type     | Default | Description                                                                  |
| ------------- | -------- | ------- | ---------------------------------------------------------------------------- |
| `save_path` | str/Path | "."     | Destination file or directory where the prepared PDBQT file should be saved. |

Example:

```python
from docksuitex import Ligand

# Load ligand
lig = Ligand("ligand.pdb")

# Prepare ligand for docking
lig.prepare(minimize="mmff94", remove_water=True, add_hydrogens=True)

# Save the final PDBQT file
lig.save_pdbqt("outputs/ligand_prepared.pdbqt")
```

---

### 3. `pocket_finder`

#### Class: `PocketFinder`

| Parameter        | Type     | Default         | Description                               |
| ---------------- | -------- | --------------- | ----------------------------------------- |
| `protein_path` | str/Path | —              | Path to input protein file in PDB format. |
| `threads`      | int      | os.cpu\_count() | Number of threads to use for P2Rank.      |

#### Method: `run()`

| Parameter | Type | Default | Description                                                                                                               |
| --------- | ---- | ------- | ------------------------------------------------------------------------------------------------------------------------- |
| —        | —   | —      | Runs P2Rank to predict ligand-binding pockets in the protein. Returns a list of pockets with rank and center coordinates. |

#### Method: `save_report()`

| Parameter    | Type     | Default              | Description                                              |
| ------------ | -------- | -------------------- | -------------------------------------------------------- |
| `save_dir` | str/Path | `./p2rank_results` | Directory where the P2Rank output folder will be copied. |

Example:

```python
from docksuitex import PocketFinder

# Initialize pocket finder
pf = PocketFinder("protein.pdb", threads=4)

# Run P2Rank to predict pockets
pockets = pf.run()
for pocket in pockets:
    print(f"Rank {pocket['rank']}: Center at {pocket['center']}")

# Save full P2Rank output folder
pf.save_report("outputs/pockets")
```

---

### 4. `autodock4`

#### Class: `AD4Docking`

| Parameter              | Type                        | Default      | Description                                                                                |
| ---------------------- | --------------------------- | ------------ | ------------------------------------------------------------------------------------------ |
| `receptor`           | str/Path                    | —           | Path to receptor PDBQT file.                                                               |
| `ligand`             | str/Path                    | —           | Path to ligand PDBQT file.                                                                 |
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

#### Method: `save_results()`

| Parameter    | Type     | Default           | Description                                                                      |
| ------------ | -------- | ----------------- | -------------------------------------------------------------------------------- |
| `save_dir` | str/Path | `./ad4_results` | Directory where the docking results (DLG, GPF, DPF, PDBQT files) will be copied. |

Example:

```python
from docksuitex import AD4Docking

# Initialize docking
docking = AD4Docking(
    receptor="outputs/protein_prepared.pdbqt",
    ligand="outputs/ligand_prepared.pdbqt",
    grid_center=(10.0, 12.5, 8.0),
    grid_size=(60, 60, 60),
    ga_run=10
)

# Run AutoGrid + AutoDock
docking.run()

# Save results
docking.save_results("outputs/ad4_docking")
```

---

### 5. `vina`

#### Class: `VinaDocking`

| Parameter          | Type                        | Default         | Description                                                                                                        |
| ------------------ | --------------------------- | --------------- | ------------------------------------------------------------------------------------------------------------------ |
| `receptor`       | str/Path                    | —              | Path to receptor PDBQT file.                                                                                       |
| `ligand`         | str/Path                    | —              | Path to ligand PDBQT file.                                                                                         |
| `grid_center`    | tuple\[float, float, float] | —              | Center coordinates (x, y, z) of the docking grid in Ångström.                                                    |
| `grid_size`      | tuple\[int, int, int]       | (20, 20, 20)    | Physical length of the grid box along (x, y, z) in **Ångström**. Spacing is fixed internally at 0.375 Å. |
| `exhaustiveness` | int                         | 8               | How exhaustively Vina searches conformational space.                                                               |
| `num_modes`      | int                         | 9               | Maximum number of binding modes to output.                                                                         |
| `cpu`            | int                         | os.cpu\_count() | Number of CPU cores to use.                                                                                        |
| `verbosity`      | int                         | 1               | Level of console output (0=quiet, 1=normal, 2=verbose).                                                            |
| `seed`           | int, optional               | None            | Random seed for reproducibility.                                                                                   |

#### Method: `run()`

| Parameter | Type | Default | Description                                                                                  |
| --------- | ---- | ------- | -------------------------------------------------------------------------------------------- |
| —        | —   | —      | Runs AutoDock Vina with the given parameters. Outputs docking results to a temporary folder. |

#### Method: `save_results()`

| Parameter    | Type     | Default            | Description                                                    |
| ------------ | -------- | ------------------ | -------------------------------------------------------------- |
| `save_dir` | str/Path | `./vina_results` | Directory where docking results (.pdbqt, .log) will be copied. |

Example:

```python
from docksuitex import VinaDocking

# Initialize Vina docking
vina = VinaDocking(
    receptor="outputs/protein_prepared.pdbqt",
    ligand="outputs/ligand_prepared.pdbqt",
    grid_center=(10.0, 12.5, 8.0),
    exhaustiveness=8,
    num_modes=9
)

# Run docking
vina.run()

# Save results
vina.save_results("outputs/vina_docking")
```

---

### 6. `utils/parser`

#### Method: `parse_vina_log_to_csv()`

| Parameter       | Type | Default                | Description                                                                        |
| --------------- | ---- | ---------------------- | ---------------------------------------------------------------------------------- |
| `results_dir` | str  | —                     | Parent directory containing docking result folders with Vina log files (`.txt`). |
| `output_csv`  | str  | `"vina_summary.csv"` | Path to save the output CSV file.                                                  |

Example:

```python
from docksuitex.utils import parse_vina_log_to_csv

df = parse_vina_log_to_csv(results_dir="outputs", output_csv="outputs/vina_summary.csv")
print(df.head())
```

#### Method: `parse_ad4_dlg_to_csv()`

| Parameter       | Type       | Default               | Description                                                              |
| --------------- | ---------- | --------------------- | ------------------------------------------------------------------------ |
| `results_dir` | str / Path | —                    | Parent directory containing docking result folders with `results.dlg`. |
| `output_csv`  | str        | `"ad4_summary.csv"` | Path to save the output CSV file.                                        |

Example:

```python
from docksuitex.utils import parse_ad4_dlg_to_csv

df = parse_ad4_dlg_to_csv(results_dir="outputs", output_csv="outputs/ad4_summary.csv")
print(df.head())
```

---

### 7. `utils/fetcher`

#### **Method:** `fetch_pdb()`

| Parameter    | Type       | Default | Description                                    |
| ------------ | ---------- | ------- | ---------------------------------------------- |
| `pdbid`    | str        | —      | 4-character alphanumeric PDB ID (e.g., '1UBQ') |
| `save_dir` | str / Path | `.`   | Directory to save the `.pdb` file            |

Example:

```python
from docksuitex.utils import fetch_pdb

pdb_file = fetch_pdb("1UBQ", save_dir="data/pdbs")
```

#### Method: `fetch_sdf()`

| Parameter    | Type       | Default | Description                                       |
| ------------ | ---------- | ------- | ------------------------------------------------- |
| `cid`      | str / int  | —      | PubChem Compound ID (CID), e.g., 2244 for Aspirin |
| `save_dir` | str / Path | `.`   | Directory to save the `.sdf` file               |

Example:

```python
from docksuitex.utils import fetch_sdf

sdf_file = fetch_sdf(2244, save_dir="data/sdfs")
```

---

### 8. `utils/viewer`

#### Method: `view_molecule()`

| Parameter      | Type                                                     | Default    | Description                                                              |
| -------------- | -------------------------------------------------------- | ---------- | ------------------------------------------------------------------------ |
| `file_path`  | str                                                      | —         | Path to the molecular file (`.pdb`,`.pdbqt`,`.mol2`, or `.sdf`). |
| `style`      | Literal["stick", "line", "sphere", "cartoon", "surface"] | "stick"    | Visualization style of the molecule.                                     |
| `background` | str                                                      | "white"    | Background color of the 3Dmol viewer (e.g.,`"white"`,`"black"`).     |
| `color`      | str                                                      | "spectrum" | Coloring method for atoms or residues (e.g.,`"spectrum"`,`"chain"`). |
| `width`      | int                                                      | 500        | Width of the rendered viewer window in pixels.                           |
| `height`     | int                                                      | 500        | Height of the rendered viewer window in pixels.                          |

Example:

```python
from docksuitex.utils import view_molecule

view_molecule(file_path="protein.pdbqt", style="cartoon")
```

---

### 9. `utils/cleaner`

#### Method: `clean_temp_folder()`

```python
from docksuitex.utils import clean_temp_folder

clean_temp_folder()
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
* [py3Dmol](https://pypi.org/project/py3Dmol/)

---

## 📜 License

This project is licensed under the GNU GPL v3 License - see the [LICENSE](LICENSE) file for details.
