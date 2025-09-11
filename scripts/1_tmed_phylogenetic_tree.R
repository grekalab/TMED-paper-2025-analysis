# ---------------------------------------------------------------------------
# Script: 1_tmed_phylogenetic_tree.R
#
# Description:
# Downloads human TMED protein sequences from UniProt, performs a multiple
# sequence alignment, and generates a maximum likelihood phylogenetic tree
# with 1,000 bootstrap replicates.
#
# Usage:
# This script is intended to be run from the root of the project directory.
# Example: source("scripts/1_tmed_phylogenetic_tree.R")
# ---------------------------------------------------------------------------


# -----------------------------
# 0. Load Libraries
# -----------------------------
library(httr)
library(Biostrings)
library(stringr)
library(msa)
library(phangorn)
library(ggplot2)
library(ggtree)


# -----------------------------
# 1. Download and Prepare Data
# -----------------------------

# Ensure the output directory exists
dir.create("output", showWarnings = FALSE)

# Define UniProt API URL for canonical human TMED protein isoforms
url <- "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=(gene:TMED*)%20AND%20(organism_id:9606)%20AND%20(reviewed:true)"

# Download and save FASTA to the output directory
fasta_file <- "output/tmeds_uniprot.fasta"
GET(url, write_disk(fasta_file, overwrite = TRUE))

# Read protein sequences
seqs <- readAAStringSet(fasta_file)

# Remove TMED8 as it is a pseudogene in humans
seqs <- seqs[!grepl("TMED8", names(seqs))]

# Simplify sequence names (e.g., from "sp|Q13445|TMED1_HUMAN..." to "TMED1")
old_labels <- names(seqs)
names(seqs) <- str_match(old_labels, "TMED[0-9]+")[, 1]

# View basic info
print(seqs)


# -----------------------------
# 2. Multiple Sequence Alignment (MSA)
# -----------------------------
alignment <- msa(seqs, method = "Muscle")

# Convert to phangorn-compatible format
set.seed(20250701) # For reproducibility
alignment_phangorn <- msaConvert(alignment, type = "seqinr::alignment")
phy_data <- as.phyDat(alignment_phangorn, type = "AA")

# -----------------------------
# 3. Phylogenetic Tree Construction
# -----------------------------

# Create a distance matrix and a neighbor-joining starting tree
dist_matrix <- dist.ml(phy_data, model = "WAG")
start_tree <- nj(dist_matrix)

# Fit a maximum likelihood model using the neighbor-joining tree as a start point
# The 'WAG' model is a common substitution model for globular proteins
fit <- pml(start_tree, data = phy_data)
fit <- optim.pml(fit, model = "WAG", rearrangement = "stochastic")


# -----------------------------
# 4. Bootstrap Analysis
# -----------------------------
# Perform bootstrap analysis with 1,000 replicates for statistical support
set.seed(2025)  # For reproducibility
bs <- bootstrap.pml(fit, bs = 1000, optNni = TRUE)

# Get the final maximum likelihood tree
ml_tree <- fit$tree
# Ladderize the tree for cleaner visualization
ml_tree <- ladderize(ml_tree, right = FALSE)
# Associate the bootstrap values with the tree nodes
tree_with_bs <- plotBS(ml_tree, bs, type = "p")


# -----------------------------
# 5. Plot and Save Tree
# -----------------------------
# Plot tree
ggtree(tree_with_bs) +
  geom_tiplab(hjust = -0.1, linetype = "dotted") +
  geom_text2(aes(subset = !isTip, label = label), hjust = -0.3) +
  scale_x_continuous(limits = c(0, 1.8)) +
  geom_treescale(x=1.1, y=2, offset=-0.7, width=0.5)

# Save as PDF
ggsave("output/tmed_phylogenetic_tree.pdf", height = 3, width = 6)

print("Phylogenetic tree construction complete. Check the output/ directory for the PDF.")