# ---------------------------------------------------------------------------
# Script: 2_limma_volcano_plot.R
#
# Description:
# Performs differential abundance analysis for Co-IP mass spectrometry data.
# This script reads raw intensity data, performs normalization and imputation,
# runs a limma contrast analysis, and generates a volcano plot.
#
# Usage:
# This script is intended to be run from the root of the project directory.
# Example: source("scripts/2_limma_volcano_plot.R")
# ---------------------------------------------------------------------------


# -----------------------------
# 0. Load Libraries
# -----------------------------
# Clear workspace
rm(list = ls())

# Load libraries
library(readr)
library(dplyr)
library(missForest)
library(limma)
library(reshape2)
library(tibble)
library(ggrepel)

# -----------------------------
# 1. Read and Prepare Data
# -----------------------------

# Ensure the output directory exists
dir.create("output", showWarnings = FALSE)

# Read input data from the /data directory
input_file <- "data/20240826_GrekaLab_MKG_original.csv"
data_raw <- read_csv(input_file) %>%
  distinct(Genes, .keep_all = TRUE)  # Remove duplicate gene rows

# Save gene names (as uppercase)
gene_names <- toupper(data_raw$Genes)

# Remove metadata columns, keeping only numeric intensities
data_numeric <- data_raw %>%
  select(-Genes, -contains("Protein")) %>%
  as.matrix()

# Replace any missing values (NA) with 0
data_numeric[is.na(data_numeric)] <- 0

# -----------------------------
# 2. Log2 Transformation
# -----------------------------

data_log2 <- data_numeric
# Apply log2 transformation only to non-zero values to avoid -Inf
data_log2[data_numeric != 0] <- log2(data_numeric[data_numeric != 0])

# -----------------------------
# 3. Median Normalization (Optional, based on experimental design)
# -----------------------------
# Note: For Co-IP, normalization strategies can vary. This is one common approach.
# We are assuming columns 5-12 are the samples of interest for global median calculation.
# Adjust this range if your data format differs.
global_median <- median(data_log2[,5:12][data_log2[,5:12] != 0])

data_norm <- apply(data_log2, 2, function(sample) {
  sample_median <- median(sample[sample != 0])
  # Avoid division by zero for empty samples
  if (sample_median != 0) {
    return(sample / sample_median * global_median)
  } else {
    return(sample)
  }
})

# Add gene names back as rownames
rownames(data_norm) <- gene_names

# -----------------------------
# 4. Missing Value Imputation
# -----------------------------
# Impute zero values using missForest. Zeros are converted to NA for imputation.
set.seed(2025) # For reproducibility

# Define experimental groups and column indices
group_names <- c("EV", "TMED5", "TMED7")
reps <- 4
groups <- rep(group_names, each = reps)
treatment_cols_list <- list(EV = 1:4, TMED5 = 5:8, TMED7 = 9:12)

# Convert 0s to NAs for imputation
data_norm_na <- data_norm
data_norm_na[data_norm == 0] <- NA

# Function to impute missing data within a subset of columns
impute_within_subset <- function(data_subset) {
  if (any(is.na(data_subset))) {
    # missForest requires a data frame
    imputed_result <- missForest(as.data.frame(data_subset), verbose = FALSE)
    return(imputed_result$ximp)
  } else {
    return(data_subset)
  }
}

# Apply imputation to each treatment group separately
imputed_treatments_list <- lapply(treatment_cols_list, function(cols) {
  subset_data <- data_norm_na[, cols, drop = FALSE]
  cat("Imputing missing data for", names(treatment_cols_list), "...")
  imputed_data <- impute_within_subset(subset_data)
  cat("done!\n")
  names(imputed_data) <- NULL
  return(imputed_data)
})

# Combine the imputed data frames back into a single matrix
imputed_data <- as.matrix(do.call(cbind, imputed_treatments_list))
write_csv(as.data.frame(imputed_data) %>% rownames_to_column("Gene"), "output/imputed_data.csv")
colnames(imputed_data) <- groups

# -----------------------------
# 5. Plot Data Distribution
# -----------------------------
df_long <- melt(imputed_data, varnames = c("Gene", "Sample"), value.name = "Intensity")

dist_plot <- ggplot(df_long, aes(x = Intensity, color = Sample)) +
  geom_density() +
  theme_minimal() +
  labs(
    title = "Post-Normalization and Imputation Density",
    x = "Log2-Transformed Intensity",
    y = "Density"
  )
ggsave("output/density_plot.png", plot = dist_plot, width = 8, height = 5)


# -----------------------------
# 6. Differential Expression with limma
# -----------------------------

# Define group labels based on the column order of `imputed_data`
group <- factor(rep(c("EV", "TMED5", "TMED7"), each = 4))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Fit linear model
fit <- lmFit(imputed_data, design)

# Define contrast between TMED7 and TMED5
contrast.matrix <- makeContrasts(TMED7 - TMED5, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract differential expression results
results <- topTable(fit2, number = Inf, adjust.method = "BH")
results$Protein <- rownames(results)

# Label significant hits for plotting
target.genes <- c("TMED5", "TMED7", "GORASP2", "BLZF1")
rename.map <- c("GORASP2" = "GRASP55", "BLZF1" = "Golgin-45")
results <- results %>%
  mutate(
    Significant = case_when(
      Protein %in% target.genes     ~ "Gene of interest",
      adj.P.Val < 0.05 & logFC > 1   ~ "Upregulated in TMED7",
      adj.P.Val < 0.05 & logFC < -1  ~ "Upregulated in TMED5",
      TRUE                          ~ "Not significant"
    ),
    Protein.renamed = recode(Protein, !!!rename.map),
    Significant = factor(Significant, levels = c("Gene of interest", "Upregulated in TMED7", "Upregulated in TMED5", "Not significant"))
  )

write_csv(results, "output/fold_change_data_TMED7_vs_TMED5.csv")

# -----------------------------
# 7. Volcano Plot
# -----------------------------

volcano <- ggplot(results, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.5) +
  geom_point(data = ~filter(., Protein %in% target.genes), size = 2.5) +
  geom_text_repel(
    data = ~filter(., Protein %in% target.genes),
    aes(label = Protein.renamed),
    size = 4,
    fontface = "bold",
    color = "black",
    segment.color = "black",
    box.padding = 0.5
  ) +
  scale_color_manual(values = c("Gene of interest" = "black",
                                "Upregulated in TMED7" = "#6666FF",
                                "Upregulated in TMED5" = "#FF6666",
                                "Not significant" = "grey70")) +
  labs(
    title = "TMED7 vs TMED5 Co-IP",
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("p-value"))
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.title = element_blank()
  )

ggsave("output/volcano_plot_TMED7_vs_TMED5.png", plot = volcano, width = 8, height = 6)

print("Analysis complete. Check the output/ directory for results.")