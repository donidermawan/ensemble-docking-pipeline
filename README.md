# Ensemble-Docking-Pipeline

A reproducible computational workflow for **Receptor structure-based virtual screening, ensemble molecular docking, redocking validation, RMSD analysis, consensus ranking, and publication-ready visualization**.

This repository contains Bash and R scripts designed to automate a multi-stage molecular docking workflow using multiple receptor conformations. The pipeline supports ligand preparation, receptor preparation, docking validation against experimentally determined ligand poses, ensemble docking across five receptor structures, consensus-based ligand prioritization, and generation of publication-quality figures.

---

## Overview

Protein flexibility can influence molecular docking performance because a single receptor conformation may not adequately represent the conformational landscape of a protein target.

This repository implements an **ensemble docking strategy** using five receptor structures:

* **1VYW**
* **7RWF**
* **5A14**
* **5IF1**
* **1W98**

The workflow combines multiple receptor conformations to evaluate ligand binding across different structural states.

The overall workflow is:

```text
Input Ligands
     │
     ▼
Ligand Preparation
prepare_ligands.sh
     │
     ▼
PDBQT Ligands
     │
     ├──────────────────────────────┐
     │                              │
     ▼                              ▼
Receptor Preparation          Redocking Validation
prepare_receptors.sh          redocking_rmsd.sh
     │                              │
     ▼                              ▼
PDBQT Receptors               RMSD + Affinity Validation
     │
     ▼
Vina Configuration
generate_vina_configs.sh
     │
     ▼
Ensemble Docking
ensemble_docking.sh
     │
     ▼
Individual Receptor Docking Results
     │
     ▼
Ensemble Statistics
Mean / Median / Best / Worst / SD
     │
     ▼
Consensus Ranking
     │
     ▼
ensemble_docking_ranked.csv
     │
     ▼
Publication-Ready Visualization
receptor_ensemble_docking_visualization.R
     │
     ▼
High-Resolution Figures
PNG + PDF
```

---

## Main Features

* Automated ligand conversion from **3D SDF to PDBQT**
* Parallel ligand preparation using **GNU Parallel**
* Automated preparation of multiple receptor structures
* Redocking validation against experimental crystal ligand poses
* Heavy-atom RMSD calculation using **Open Babel `obrms`**
* Mean RMSD calculation across generated docking poses
* Automatic receptor–ligand complex generation
* Automated generation of AutoDock Vina configuration files
* Ensemble docking against five receptor structures
* Parallelized docking using multiple simultaneous Vina processes
* Initial ligand filtering based on docking affinity
* Individual receptor ranking
* Ensemble mean, median, best, worst, and standard deviation calculations
* Consensus classification of ligands
* Consensus ranking of ligands successfully docked against all five receptor states
* Export of complete and long-format docking datasets
* Publication-ready visualization using **R/ggplot2**
* High-resolution figure export in PNG and PDF formats

---

# Repository Contents

```text
Receptor-Ensemble-Docking-Pipeline/
│
├── README.md
│
├── prepare_ligands.sh
├── prepare_receptors.sh
├── redocking_rmsd.sh
├── generate_vina_configs.sh
├── ensemble_docking.sh
│
├── receptor_ensemble_docking_visualization.R
│
├── ligands/
│   └── *.sdf
│
├── 1VYW/
│   ├── 1VYW.pdb
│   ├── 1VYW.pdbqt
│   └── ...
│
├── 7RWF/
│   ├── 7RWF.pdb
│   ├── 7RWF.pdbqt
│   └── ...
│
├── 5A14/
│   ├── 5A14.pdb
│   ├── 5A14.pdbqt
│   └── ...
│
├── 5IF1/
│   ├── 5IF1.pdb
│   └── 5IF1.pdbqt
│
├── 1W98/
│   ├── 1W98.pdb
│   └── 1W98.pdbqt
│
├── pdbqt/
│
├── logs/
│
├── failed/
│
├── results/
│
├── ensemble_results/
│
└── ensemble_figures/
```

The exact directory contents will depend on the dataset and docking run.

---

# Scripts

## 1. `prepare_ligands.sh`

### Purpose

Prepares ligand structures for AutoDock Vina by converting existing **3D SDF files into PDBQT format**.

The script:

1. Creates output directories.
2. Detects the available CPU cores.
3. Searches the `ligands/` directory for `.sdf` files.
4. Processes ligands in parallel using GNU Parallel.
5. Converts 3D SDF structures to PDBQT using Open Babel.
6. Adds hydrogens.
7. Applies protonation conditions at approximately pH 7.4.
8. Skips ligands that have already been processed.
9. Records failed ligand preparation attempts.
10. Reports the number of input, successfully converted, and failed ligands.

The script intentionally converts the existing 3D structures rather than generating new 3D coordinates.

### Input

```text
ligands/*.sdf
```

### Output

```text
pdbqt/*.pdbqt
logs/*.log
failed/failed_ligands.txt
```

### Run

```bash
chmod +x prepare_ligands.sh
./prepare_ligands.sh
```

---

# 2. `prepare_receptors.sh`

### Purpose

Prepares the five receptor structures for docking.

The receptors included in the workflow are:

```text
1VYW
7RWF
5A14
5IF1
1W98
```

The script uses:

```text
mk_prepare_receptor.py
```

and processes each receptor structure individually.

The receptor preparation command includes:

```text
--read_pdb
--allow_bad_res
-p
```

### Input

The script expects receptor PDB files named:

```text
1VYW.pdb
7RWF.pdb
5A14.pdb
5IF1.pdb
1W98.pdb
```

### Output

For each receptor:

```text
1VYW.pdbqt
7RWF.pdbqt
5A14.pdbqt
5IF1.pdbqt
1W98.pdbqt
```

Preparation logs are also generated.

### Run

```bash
chmod +x prepare_receptors.sh
./prepare_receptors.sh
```

---

# 3. `redocking_rmsd.sh`

## Purpose

Performs **reference-ligand redocking and RMSD validation** for selected receptor–ligand systems.

The validation systems included in the script are:

| Receptor | Reference ligand |
| -------- | ---------------- |
| 7RWF     | TW8672           |
| 1VYW     | PNU-292137       |
| 5A14     | K03861           |

The script uses:

* AutoDock Vina
* Open Babel
* `obrms`

The redocking procedure evaluates the generated Vina poses against experimentally determined crystal ligand poses.

### Validation workflow

```text
Crystal ligand
      │
      ├───────────────┐
      │               │
      ▼               ▼
Experimental      Vina docking
pose                  │
                      ▼
                 Docked poses
                      │
                      ▼
                RMSD calculation
                      │
                      ▼
              Mean RMSD calculation
```

The script evaluates all reported Vina poses and calculates the mean RMSD.

It also extracts the best Vina pose (`MODEL 1`) and combines it with the receptor to generate a receptor–ligand complex PDB file.

### Software check

The script verifies the availability of:

```bash
vina
obabel
obrms
```

and reports their installed versions.

### Validation systems

#### TW8672 – 7RWF

```text
Receptor:
7RWF/7RWF.pdbqt

Ligand:
7RWF/TW8672.pdbqt

Experimental ligand:
7RWF/TW8672_crystal.pdb
```

#### PNU-292137 – 1VYW

```text
Receptor:
1VYW/1VYW.pdbqt

Ligand:
1VYW/PNU-292137.pdbqt

Experimental ligand:
1VYW/PNU-292137_crystal.pdb
```

#### K03861 – 5A14

```text
Receptor:
5A14/5A14.pdbqt

Ligand:
5A14/K03861.pdbqt

Experimental ligand:
5A14/K03861_crystal.pdb
```

### Output

```text
results/
├── 7RWF_TW8672/
├── 1VYW_PNU292137/
├── 5A14_K03861/
└── redocking_summary.csv
```

The summary contains:

```text
System
Receptor
Ligand
Vina_Affinity_kcal_mol
Mean_RMSD_A
```

The workflow also produces receptor–ligand complex PDB files and RMSD output files.

### Run

```bash
chmod +x redocking_rmsd.sh
./redocking_rmsd.sh
```

---

# 4. `generate_vina_configs.sh`

## Purpose

Generates AutoDock Vina configuration files for the five receptor structures used in ensemble docking.

The configured receptor structures are:

```text
1VYW
7RWF
5A14
5IF1
1W98
```

The default docking parameters implemented in the script are:

```text
Grid size:
24 × 24 × 24 Å

Exhaustiveness:
8

Number of modes:
3

Energy range:
3 kcal/mol
```

### Configuration files

The script generates:

```text
vina_1VYW.txt
vina_7RWF.txt
vina_5A14.txt
vina_5IF1.txt
vina_1W98.txt
```

### Docking centers

The receptor-specific centers are defined in the script.

#### 1VYW

```text
center_x = -8.929
center_y = 205.749
center_z = 111.606
```

#### 7RWF

```text
center_x = -22.195
center_y = 2.569
center_z = 17.819
```

#### 5A14

```text
center_x = 22.287
center_y = 8.979
center_z = 10.665
```

#### 5IF1

```text
center_x = -73.670
center_y = -55.376
center_z = 19.600
```

#### 1W98

```text
center_x = 24.006
center_y = 20.821
center_z = -14.838
```

### Run

```bash
chmod +x generate_vina_configs.sh
./generate_vina_configs.sh
```

---

# 5. `ensemble_docking.sh`

## Purpose

Performs the main **receptor ensemble docking and consensus ranking workflow**.

The script docks selected ligands against five receptor conformations:

```text
1VYW
7RWF
5A14
5IF1
1W98
```

The initial ligand selection threshold is:

```text
Affinity ≤ -10.00 kcal/mol
```

The script is designed to perform parallel docking using:

```text
12 simultaneous Vina processes
1 CPU per Vina process
```

This parameter can be modified through:

```bash
JOBS=12
```

---

## Ensemble Docking Workflow

```text
Initial screening
       │
       ▼
Affinity ≤ -10 kcal/mol
       │
       ▼
Selected ligand list
       │
       ▼
 ┌─────────────┬─────────────┬─────────────┬─────────────┬
 │             │             │             │             │             
 ▼             ▼             ▼             ▼             ▼
1VYW         7RWF          5A14          5IF1          1W98
 │             │             │             │             │
 ▼             ▼             ▼             ▼             ▼
Vina          Vina          Vina          Vina          Vina
 │             │             │             │             │
 └─────────────┴─────────────┴─────────────┴─────────────┴
                              │
                              ▼
                   Ensemble statistics
                              │
              ┌───────────────┼────────────────┐
              ▼               ▼                ▼
             Mean           Median              SD
              │               │                │
              └───────────────┼────────────────┘
                              ▼
                     Consensus ranking
```

---

## Initial ligand selection

The script reads:

```text
vina_results_1VYW_1-7629.csv
```

and selects ligands with docking affinity:

```text
≤ -10.00 kcal/mol
```

The selected ligand information is stored as:

```text
ensemble_results/selected_ligands/
```

including:

```text
selected_ligands.csv
selected_ligands.txt
```

---

# Ensemble Receptor Docking

For every selected ligand, the script performs docking against each of the five receptor structures.

Each receptor generates an individual result file:

```text
vina_results_1VYW.csv
vina_results_7RWF.csv
vina_results_5A14.csv
vina_results_5IF1.csv
vina_results_1W98.csv
```

Ranked versions are also generated.

Docking poses are stored as PDBQT files.

Receptor–ligand complexes are generated as PDB files using the first Vina pose.

---

# Ensemble Statistics

For each ligand, the workflow calculates:

### Mean affinity

The arithmetic mean of the available docking affinities across receptor states.

### Median affinity

The median docking affinity across available receptor states.

### Best affinity

The most negative affinity value.

### Worst affinity

The least negative affinity value.

### Standard deviation

The standard deviation of affinity values across receptor states.

### Successful receptors

The number of receptor states for which docking produced a numerical affinity.

### Strong hits

The number of receptor states in which:

```text
Affinity ≤ -10.0 kcal/mol
```

---

# Consensus Classification

The pipeline classifies ligands according to the number of receptor states meeting the `≤ -10 kcal/mol` threshold.

The classification scheme implemented in the script is:

| Condition                        | Classification        |
| -------------------------------- | --------------------- |
| 0 successful receptors           | `FAILED`              |
| 5/5 successful + 5/5 strong hits | `Strong_consensus`    |
| ≥4 strong hits                   | `High_consensus`      |
| ≥3 strong hits                   | `Moderate_consensus`  |
| ≥2 strong hits                   | `Weak_consensus`      |
| 1 strong hit                     | `Single_receptor_hit` |
| 0 strong hits                    | `No_consensus`        |

More negative docking affinity is treated as a more favorable predicted binding score.

---

# Consensus Ranking

The primary consensus ranking is restricted to ligands successfully docked against all five receptor states.

Ranking is performed using:

1. **Mean affinity**
2. **Median affinity**
3. **Standard deviation**

More negative mean affinity receives higher priority.

This provides a simple consensus framework for prioritizing ligands that maintain favorable predicted binding across multiple receptor conformations.

---

# Ensemble Outputs

The main output directory is:

```text
ensemble_results/
```

with the following structure:

```text
ensemble_results/
│
├── selected_ligands/
│   ├── selected_ligands.csv
│   └── selected_ligands.txt
│
├── docking/
│   ├── 1VYW/
│   ├── 7RWF/
│   ├── 5A14/
│   ├── 5IF1/
│   └── 1W98/
│
├── complexes/
│   ├── 1VYW/
│   ├── 7RWF/
│   ├── 5A14/
│   ├── 5IF1/
│   └── 1W98/
│
├── receptor_pdb/
│
├── logs/
│   ├── 1VYW/
│   ├── 7RWF/
│   ├── 5A14/
│   ├── 5IF1/
│   └── 1W98/
│
└── summaries/
    ├── vina_results_1VYW.csv
    ├── vina_results_7RWF.csv
    ├── vina_results_5A14.csv
    ├── vina_results_5IF1.csv
    ├── vina_results_1W98.csv
    ├── vina_results_1VYW_ranked.csv
    ├── vina_results_7RWF_ranked.csv
    ├── vina_results_5A14_ranked.csv
    ├── vina_results_5IF1_ranked.csv
    ├── vina_results_1W98_ranked.csv
    ├── ensemble_docking_summary.csv
    ├── ensemble_docking_ranked.csv
    └── ensemble_docking_all_ligands_ranked.csv
```

---

# 6. `receptor_ensemble_docking_visualization.R`

## Purpose

Generates **publication-ready visualizations** from:

```text
ensemble_docking_ranked.csv
```

The visualization script uses:

* `readr`
* `dplyr`
* `tidyr`
* `ggplot2`
* `scales`

The script automatically checks for missing R packages and installs them when required.

---

# Input Requirements

The main input file is:

```text
ensemble_docking_ranked.csv
```

The script expects the following columns:

```text
Rank
Ligand
Initial_1VYW_Affinity
Mean_Affinity
Median_Affinity
Best_Affinity
Worst_Affinity
SD_Affinity
Successful_Receptors
1VYW
7RWF
5A14
5IF1
1W98
```

Only ligands with complete docking results across all five receptor states are retained for the primary visualization analysis.

---

# Generated Figures

The R script generates seven main figures.

## Figure 1 — Top-20 Consensus Ligands

Displays the highest-ranked ligands according to the ensemble consensus ranking.

## Figure 2 — Top-20 Ligand × Receptor Affinity Heatmap

Visualizes docking affinities of the top-ranked ligands across the five receptor conformations.

## Figure 3 — Receptor-Wise Affinity Distributions

Shows the distribution of docking affinities for each receptor state.

## Figure 4 — Ensemble Mean Affinity vs Variability

Compares mean ensemble affinity against affinity variability.

## Figure 5 — Consensus Rank vs Mean Affinity

Visualizes the relationship between consensus rank and ensemble mean affinity.

## Figure 6 — Top-10 Ligand Affinity Profiles

Shows affinity profiles of the top 10 ligands across the five receptor conformations.

## Figure 7 — Distribution of Ensemble Mean Affinities

Shows the overall distribution of ensemble mean docking affinities and highlights the `-10 kcal/mol` threshold and dataset mean.

All figures are exported at high resolution.

---

# Running the Visualization

Make sure:

```text
ensemble_docking_ranked.csv
```

is located in the working directory, or modify:

```r
INPUT_FILE <- "ensemble_docking_ranked.csv"
```

The output directory is:

```r
OUTPUT_DIR <- "ensemble_figures"
```

Run:

```bash
Rscript receptor_ensemble_docking_visualization.R
```

or:

```bash
chmod +x receptor_ensemble_docking_visualization.R
./receptor_ensemble_docking_visualization.R
```

---

# Visualization Outputs

The generated figures are stored in:

```text
ensemble_figures/
```

The script exports publication-quality PNG and PDF figures.

Additional exported datasets include:

```text
top20_consensus_ligands.csv
ensemble_docking_long_format.csv
ensemble_summary_statistics.csv
```

These files can be used for downstream statistical analysis, figure refinement, or supplementary data preparation.

---

# Installation

## 1. AutoDock Vina

Install AutoDock Vina and ensure that:

```bash
vina
```

is available from the command line.

Verify:

```bash
vina --version
```

The redocking validation script explicitly checks for AutoDock Vina and reports its version.

---

## 2. Open Babel

Open Babel is required for:

* SDF → PDBQT conversion
* PDBQT → PDB conversion
* PDBQT → SDF conversion
* RMSD preparation

Verify:

```bash
obabel -V
```

---

## 3. Open Babel `obrms`

The redocking validation workflow requires:

```bash
obrms
```

Verify:

```bash
obrms --version
```

---

## 4. GNU Parallel

GNU Parallel is required for parallel ligand preparation and ensemble docking.

Verify:

```bash
parallel --version
```

---

## 5. Receptor Preparation Tool

The receptor preparation workflow requires:

```bash
mk_prepare_receptor.py
```

Make sure the command is available in your environment.

Verify:

```bash
which mk_prepare_receptor.py
```

---

## 6. R

R is required for publication-ready visualization.

Required packages:

```r
readr
dplyr
tidyr
ggplot2
scales
```

The visualization script automatically attempts to install missing packages.

---

# Recommended Execution Order

For a new docking project, execute the scripts in the following order:

### Step 1 — Prepare ligands

```bash
./prepare_ligands.sh
```

### Step 2 — Prepare receptors

```bash
./prepare_receptors.sh
```

### Step 3 — Validate the docking protocol

```bash
./redocking_rmsd.sh
```

### Step 4 — Generate Vina configurations

```bash
./generate_vina_configs.sh
```

### Step 5 — Perform ensemble docking

```bash
./ensemble_docking.sh
```

### Step 6 — Generate publication figures

```bash
Rscript receptor_ensemble_docking_visualization.R
```

---

# Reproducibility

The pipeline is designed to keep the major computational parameters explicitly defined in the scripts.

Important parameters include:

### Ensemble receptors

```text
1VYW
7RWF
5A14
5IF1
1W98
```

### Initial affinity threshold

```text
≤ -10.00 kcal/mol
```

### Vina search parameters

```text
Grid size: 24 × 24 × 24 Å
Exhaustiveness: 8
Number of modes: 3
Energy range: 3 kcal/mol
```

### Parallelization

```text
12 simultaneous Vina processes
1 CPU per Vina process
```

These settings can be modified directly in the corresponding scripts.

---

# Important Notes

## 1. Docking affinity is a computational score

The docking affinity values generated by AutoDock Vina represent predicted binding energies and should not be interpreted as experimentally measured binding affinities.

More negative values indicate more favorable predicted docking scores within the scoring framework.

---

## 2. Ensemble consensus does not establish biological activity

A ligand receiving a strong ensemble consensus classification indicates consistent predicted docking performance across the receptor conformations used in this workflow.

It does not by itself establish:

* biochemical activity
* cellular activity
* pharmacological efficacy
* selectivity
* experimental binding affinity
* therapeutic potential

Experimental validation remains necessary.

---

## 3. Receptor structures are project-specific

The docking coordinates and receptor structures in the current scripts are specifically configured for the five receptor structures:

```text
1VYW
7RWF
5A14
5IF1
1W98
```

If different proteins or receptor structures are used, the receptor list, docking centers, input files, and configuration files must be modified accordingly.

---

## 4. Input structures must be prepared appropriately

The ligand preparation script assumes that the input SDF structures already contain appropriate 3D coordinates.

The workflow converts existing 3D SDF structures rather than generating new 3D conformations during ligand preparation.

---

## 5. Hardware considerations

Ensemble docking can be computationally intensive, particularly when thousands of ligands are screened against multiple receptor conformations.

The parameter:

```bash
JOBS=12
```

controls the number of simultaneous docking processes.

This should be adjusted according to the available CPU resources.

---

# Example Minimal Usage

After cloning the repository:

```bash
git clone https://github.com/donidermawan/ensemble-docking-pipeline.git
cd ensemble-docking-pipeline
```

Make the scripts executable:

```bash
chmod +x scripts/*.sh
chmod +x scripts/*.R
```

Prepare ligands:

```bash
./scripts/prepare_ligands.sh
```

Prepare receptors:

```bash
./scripts/prepare_receptors.sh
```

Perform redocking validation:

```bash
./scripts/redocking_rmsd.sh
```

Generate Vina configurations:

```bash
./scripts/generate_vina_configs.sh
```

Run ensemble docking:

```bash
./scripts/ensemble_docking.sh
```

Generate figures:

```bash
Rscript scripts/receptor_ensemble_docking_visualization.R
```

---

# Output Interpretation

The principal outputs of the pipeline are:

| Output                                    | Purpose                                                             |
| ----------------------------------------- | ------------------------------------------------------------------- |
| `selected_ligands.csv`                    | Ligands passing the initial affinity threshold                      |
| `vina_results_<receptor>.csv`             | Individual receptor docking results                                 |
| `vina_results_<receptor>_ranked.csv`      | Ranked receptor-specific results                                    |
| `ensemble_docking_summary.csv`            | Complete ensemble statistics                                        |
| `ensemble_docking_ranked.csv`             | Consensus ranking for complete 5/5 docking                          |
| `ensemble_docking_all_ligands_ranked.csv` | Statistical ranking including ligands with partial receptor results |
| `redocking_summary.csv`                   | Redocking affinity and mean RMSD summary                            |
| `*_RMSD.txt`                              | Individual redocking RMSD values                                    |
| `*_complex.pdb`                           | Receptor–ligand complexes                                           |
| `top20_consensus_ligands.csv`             | Top consensus-ranked ligands                                        |
| `ensemble_docking_long_format.csv`        | Long-format receptor/ligand affinity dataset                        |
| `ensemble_summary_statistics.csv`         | Overall ensemble summary statistics                                 |
| `Figure_*.png`                            | High-resolution visualization figures                               |
| `Figure_*.pdf`                            | Publication-ready PDF figures                                       |

---

# Scientific Scope

This repository is intended for computational studies involving:

* Structure-based drug discovery
* Molecular docking
* Virtual screening
* Receptor inhibitor discovery
* Ensemble docking
* Consensus ligand prioritization
* Docking validation
* RMSD-based redocking assessment
* Computational medicinal chemistry
* Structure-based bioinformatics
* Publication-ready docking visualization

---

# Receptor Ensemble

The current implementation uses five receptor structures:

```text
1VYW
7RWF
5A14
5IF1
1W98
```

The ensemble strategy is intended to reduce dependence on a single receptor conformation by evaluating ligand docking performance across multiple receptor structural states.

---

# Citation

If this repository is used in a scientific publication, please cite the associated research article and the software tools used in the workflow.

Relevant computational software includes:

* AutoDock Vina
* Open Babel
* GNU Parallel
* R
* ggplot2

Please also cite the original structural databases and molecular structures used in the analysis where appropriate.

---

# Disclaimer

This repository contains computational scripts for research and methodological purposes.

The docking scores, RMSD values, consensus classifications, and rankings generated by these scripts are computational predictions and should not be considered experimental evidence of molecular binding or biological efficacy.

Users are responsible for validating the input structures, docking parameters, receptor preparation, and resulting predictions for their specific research applications.

---

# Author

**Doni Dermawan**

Computational Biology / Bioinformatics / Structure-Based Drug Discovery

GitHub: https://github.com/donidermawan
---
