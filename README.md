# Code for: A therapeutic node for diverse proteinopathies centers on a modular cargo receptor system that controls secretory pathway traffic (2025, under review)

This repository contains the R scripts and data necessary to reproduce the phylogenetic tree and mass spectrometry analysis figures for the manuscript "A therapeutic node for diverse proteinopathies centers on a modular cargo receptor system that controls secretory pathway traffic".

## Repository Structure

- `/data`: Contains the raw input data for the analyses.
- `/scripts`: Contains the R scripts used for data processing, analysis, and figure generation.
- `/output`: The destination for all generated figures and data tables.

---

## System Requirements

This analysis was performed using **R version 4.0.0** on Linux. The following R packages are required. They can be installed automatically by following the installation instructions below.

- `tidyverse`
- `limma`
- `ggrepel`
- `reshape2`
- `missForest`
- `httr`
- `Biostrings`
- `msa`
- `phangorn`
- `ape`
- `ggtree`
- `renv` (for environment management)

---

## Installation

This project uses the `renv` package to ensure reproducibility. To set up the correct environment:

1.  Clone this repository:
    ```bash
    git clone [https://github.com/grekalab/TMED-paper-2025-analysis.git](https://github.com/grekalab/TMED-paper-2025-analysis.git)
    cd TMED-paper-2025-analysis
    ```

2.  Open R or RStudio in this directory. The `renv` package will bootstrap itself and prompt you to restore the environment. Type `y` and press Enter.
    ```r
    renv::restore()
    ```
    This will install all the exact package versions used in this analysis into a project-specific library.

---

## Usage

The scripts are located in the `/scripts` directory and are numbered in the order they should be run.

### 1. TMED Family Phylogenetic Tree (`1_tmed_phylogenetic_tree.R`)

This script downloads human TMED protein sequences from UniProt, performs a multiple sequence alignment, and generates a maximum likelihood phylogenetic tree with bootstrap support values.

**To run:**

```r
# From the project root directory
source("scripts/1_tmed_phylogenetic_tree.R")
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

**To run:**

```r
# From the project root directory
source("scripts/2_limma_volcano_plot.R")
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

```<To be filled later>```

## Contact

For questions about the code, please contact John Carlos Ignacio at [jignacio@broadinstitute.org].