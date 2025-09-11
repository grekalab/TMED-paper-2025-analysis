# Code for: A therapeutic node for diverse proteinopathies centers on a modular cargo receptor system that controls secretory pathway traffic (2025, under review)

This repository contains the R scripts and data necessary to reproduce the phylogenetic tree and mass spectrometry analysis figures for the manuscript "A therapeutic node for diverse proteinopathies centers on a modular cargo receptor system that controls secretory pathway traffic".

## Repository Structure

- `data`: Contains the raw input data for the analyses.
- `scripts`: Contains the R scripts used for data processing, analysis, and figure generation.
- `output`: The destination for all generated figures and data tables.

---

## System Requirements

To run this analysis, you will need **Conda** installed on your system. The project files can be downloaded directly as a ZIP file or cloned using Git.

* **Conda:** Required for managing the R environment and all dependencies. We recommend installing **Miniconda**, a lightweight version of Conda.
    * **Installation Guide:** [https://www.anaconda.com/docs/getting-started/miniconda/install](https://www.anaconda.com/docs/getting-started/miniconda/install)

* **Git (Optional):** If you prefer to clone the repository using version control, you will need Git.
    * **Installation Guide:** [https://github.com/git-guides/install-git](https://github.com/git-guides/install-git)

The specific software versions, including R and all required packages for this project, are defined in the `environment.yml` file.

This analysis was performed using **R version 4.4.3** on Linux. The following R packages are required. They can be installed automatically by following the installation instructions below.

- `r-base>=4.3`
- `r-httr`
- `r-stringr`
- `r-phangorn`
- `r-ggplot2`
- `r-readr`
- `r-dplyr`
- `r-missforest`
- `r-reshape2`
- `r-tibble`
- `r-ggrepel`
- `bioconductor-msa`
- `bioconductor-biostrings`
- `bioconductor-ggtree`
- `bioconductor-limma`

---

## Installation

The installation process involves two main steps: getting the project files and then building the Conda environment.

### Step 1: Get the Project Files

You can get the files by either downloading a ZIP archive or cloning the repository with Git.

#### Option A: Download ZIP File (Easiest)

Download the ZIP file from here: [TMED-paper-2025-analysis-main.zip](https://github.com/grekalab/TMED-paper-2025-analysis/archive/refs/heads/main.zip).

#### Option B: Clone with Git (Recommended for developers)

Open your terminal and run the following command to clone the repository:
```bash
git clone [https://github.com/grekalab/TMED-paper-2025-analysis.git](https://github.com/grekalab/TMED-paper-2025-analysis.git)
```

### Step 2: Create and Activate the Conda Environment

1.  **Navigate into the project directory** using your terminal.
    ```bash
    # If you downloaded the ZIP, the folder name might end in "-main"
    cd path/to/TMED-paper-2025-analysis-main
    ```

2.  **Create the Conda Environment:** Use the `environment.yml` file to create a self-contained environment with all the necessary dependencies. This command will download and install the correct version of R and all required packages. This step may take several minutes.
    ```bash
    conda env create -n tmed-analysis -f environment.yml
    ```

3.  **Activate the Environment:** Before running any scripts, you must activate the newly created environment. This makes the specific R version and packages available in your terminal session.
    ```bash
    conda activate tmed-analysis
    ```
    Your terminal prompt should now be prefixed with `(tmed-analysis)`. You must do this every time you open a new terminal to work on this project.


---

## Usage

The scripts are located in the `scripts` directory and are numbered in the order they should be run.

### 1. TMED Family Phylogenetic Tree (`1_tmed_phylogenetic_tree.R`)

This script downloads human TMED protein sequences from UniProt, performs a multiple sequence alignment, and generates a maximum likelihood phylogenetic tree with bootstrap support values.

**To run from the terminal:**

```bash
# Make sure you have run 'conda activate tmed-analysis' first!
# Run from the project root directory.
Rscript scripts/2_tmed_phylogenetic_tree.R
```

**Description:**
The script performs the following steps:
1.  Downloads reviewed, canonical human TMED protein sequences from UniProt.
2.  Performs multiple sequence alignment using the `msa` package.
3.  Constructs a phylogenetic tree using a maximum likelihood model (`phangorn` package).
4.  Calculates bootstrap support with 1,000 replicates.
5.  Plots the final tree using `ggtree` and saves it as `output/tmed_tree_with_bootstrap.pdf`.

### 2. Co-IP Mass Spectrometry Analysis (`2_limma_volcano_plot.R`)

This script performs differential abundance analysis between TMED5 and TMED7 co-immunoprecipitation samples.

**To run from the terminal:**

```bash
# Make sure you have run 'conda activate tmed-analysis' first!
# Run from the project root directory.
Rscript scripts/1_limma_volcano_plot.R
```
**Description:**
The script performs the following steps:
1.  Reads the raw intensity data from `data/20240826_GrekaLab_MKG_original.csv`.
2.  Performs log2 transformation and median normalization.
3.  Imputes missing values using the `missForest` package.
4.  Conducts differential expression analysis using the `limma` package.
5.  Generates a volcano plot and saves it to `output/volcano_plot.png`.
6.  Saves the full results table to `output/fold_change_data.csv`.

---

## Citation

If you use this code or data, please cite our paper:

```<TBA>```

## Contact

For questions about the code, please contact John Carlos Ignacio at jignacio@broadinstitute.org.