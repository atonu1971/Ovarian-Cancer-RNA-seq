# ABC Transporter Expression Visualization
# Inputs:  ABC_Transporter_Expression_Table.xlsx
# Author: Atonu Chakrabortty

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


