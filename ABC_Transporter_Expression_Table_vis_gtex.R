# ABC Transporter Expression Visualization
# Inputs:  ABC_Transporter_Expression_Table.xlsx
# Author: Atonu Chakrabortty
# ABC Transporter Baseline Expression Analysis

## Overview
This script visualizes the baseline expression of ABC transporter genes across 40 human tissues
using normalized TPM values from the Human Protein Atlas. The goal is to characterize
tissue-specific expression patterns before any drug treatment or perturbation.

---

## Input Data
| File | Description |
|------|-------------|
| `ABC_Transporter_Expression_Table.xlsx` | ABC transporter gene expression table. Rows = genes, Columns = tissues. Values are TPM (Transcripts Per Million). |

### Column structure
- `Ensemble` — Ensembl gene ID
- `Gene` — HGNC gene symbol
- `MAX Tissue` — Tissue with highest expression
- `MAX TPM` — Highest TPM value across tissues
- `AVERAGE TPM` — Mean TPM across all tissues
- Remaining 40 columns — TPM per tissue (e.g. `liver`, `kidney`, `lung`, etc.)

---

## Script
| File | Description |
|------|-------------|
| `ABC_Transporter_Visualization.R` | Main R script. Loads data, transforms values, and generates 3 plots. |

### What the script does
1. Loads the Excel table using `readxl`
2. Applies **log2(TPM + 1)** transformation to normalize skewed expression values
3. Removes genes with zero variance (uninformative across tissues)
4. Produces three visualizations (printed to RStudio plot pane)

---

## Outputs (printed to plot pane)
| Plot | Description |
|------|-------------|
| Heatmap | Genes × 40 tissues, log2(TPM+1), hierarchical clustering (Ward.D2, Euclidean distance) |

---

## Dependencies
```r
# Install all with:
install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, tibble, ggplot2, pheatmap, RColorBrewer, scales, stringr)
```

---

## Notes
- TPM values are from the **Human Protein Atlas** (normal tissue expression)
- This is a **baseline / reference** analysis — no differential expression or treatment involved
- log2(TPM + 1) was chosen over z-score to preserve absolute expression differences
  while compressing the dynamic range of highly skewed TPM distributions

# Install / load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, ggplot2, pheatmap, RColorBrewer, scales, stringr)
install.packages("tibble")
library("tibble")
# Load data 
file_path <- "~/Library/CloudStorage/OneDrive-UniversityofIowa/ABC_Transporter_Expression_Table.xlsx"

df_raw <- read_excel(file_path)

# Peek at structure
cat("Dimensions:", nrow(df_raw), "x", ncol(df_raw), "\n")
cat("Columns:   ", paste(colnames(df_raw), collapse = ", "), "\n\n")

# Define tissue columns 
# All columns after the metadata columns are tissue TPM values
meta_cols <- c("Ensemble", "Gene", "MAX Tissue", "MAX TPM", "AVERAGE TPM")
tissue_cols <- setdiff(colnames(df_raw), meta_cols)

cat("Tissue columns found:", length(tissue_cols), "\n")
cat(paste(tissue_cols, collapse = ", "), "\n\n")

# Build matrix for heatmap 
mat <- df_raw %>%
  select(Gene, all_of(tissue_cols)) %>%
  column_to_rownames("Gene") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# Log2(TPM + 1) transformation — compresses skew, preserves zeros
mat_log <- log2(mat + 1)

# Remove genes with zero variance (uninformative rows)
row_var <- apply(mat_log, 1, var, na.rm = TRUE)
mat_log <- mat_log[row_var > 0, ]

cat("Genes retained after variance filter:", nrow(mat_log), "\n")

# Tidy tissue column names for display 
colnames(mat_log) <- str_replace_all(colnames(mat_log), "_", " ") %>%
  str_to_title()

# Clustered Heatmap 
cat("\nGenerating heatmap...\n")

# Colour palette: white → steel blue → dark navy
heatmap_colours <- colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100)

pheatmap(
  mat_log,
  color            = heatmap_colours,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D2",
  scale            = "none",          # already log2-transformed
  fontsize_row     = 8,
  fontsize_col     = 8,
  angle_col        = 45,
  border_color     = NA,
  main             = "ABC Transporter Expression Across Tissues\nlog2(TPM + 1)"
)


