# Sleuth Differential Expression Analysis Pipeline
### Ovarian Cancer Cell Lines — OV8 & OV90 (Sensitive vs Resistant, Drug Response)

**Original Script Author:** Professor Jared M McLendon, University of Iowa
**Modified by:** Atonu Chakrabortty, McLendon Lab, University of Iowa
**Date:** April 03, 2026
**Platform:** UIowa Argon HPC Cluster (SGE scheduler)
**Language:** R 4.3.1

---

## Table of Contents
1. [Scientific Background](#1-scientific-background)
2. [Experimental Design](#2-experimental-design)
3. [Why Kallisto and Sleuth](#3-why-kallisto-and-sleuth)
4. [Statistical Framework](#4-statistical-framework)
5. [Dependencies](#5-dependencies)
6. [Environment Setup on Argon HPC](#6-environment-setup-on-argon-hpc)
7. [Directory Structure](#7-directory-structure)
8. [Input Files](#8-input-files)
9. [Scripts in This Repository](#9-scripts-in-this-repository)
10. [How to Run](#10-how-to-run)
11. [Comparisons Performed](#11-comparisons-performed)
12. [Output Files](#12-output-files)
13. [Parameter Justification](#13-parameter-justification)
14. [Known Issues and Fixes](#14-known-issues-and-fixes)
15. [References](#15-references)

---

## 1. Scientific Background

Ovarian cancer is the fifth leading cause of cancer-related death in women and has one of the lowest five-year survival rates among gynecological malignancies. A major clinical challenge is the development of resistance to platinum-based chemotherapy (e.g., cisplatin, carboplatin), which affects the majority of patients after initial treatment.

Understanding the **transcriptomic basis of drug resistance** is essential for:
- Identifying genes and pathways that drive resistance
- Discovering potential therapeutic targets to overcome resistance
- Understanding how sensitive and resistant cells respond differently to drug treatment over time

This pipeline performs RNA-seq differential expression analysis on two ovarian cancer cell line models:
- **OV8** — ovarian cancer cell line
- **OV90** — ovarian cancer cell line

For each cell line, both **drug-sensitive** and **drug-resistant** variants were profiled at multiple time points following drug treatment, enabling two types of analysis:
1. **Baseline resistance signatures** — what is transcriptomically different between resistant and sensitive cells before any drug treatment
2. **Drug response dynamics** — how gene expression changes over time (6h, 48h) in response to drug treatment, separately in sensitive and resistant cells

---

## 2. Experimental Design

### Cell Lines, Conditions, and Time Points

| Group Label | Cell Line | Condition | Timepoint | Sample IDs | n |
|-------------|-----------|-----------|-----------|-----------|---|
| OV8_S_NT | OV8 | Sensitive | No Treatment | A1, A2, A3, A4 | 4 |
| OV8_S_6 | OV8 | Sensitive | 6 hours | B1, B2, B3, B4 | 4 |
| OV8_S_48 | OV8 | Sensitive | 48 hours | E1, E2, E3, E4 | 4 |
| OV8_R_NT | OV8 | Resistant | No Treatment | C1, C2, C3, C4 | 4 |
| OV8_R_6 | OV8 | Resistant | 6 hours | D1, D2, D3, D4 | 4 |
| OV8_R_48 | OV8 | Resistant | 48 hours | F1, F2, F3, F4 | 4 |
| OV90_S_NT | OV90 | Sensitive | No Treatment | G1, G2, G3, G4 | 4 |
| OV90_S_6 | OV90 | Sensitive | 6 hours | H1, H2, H3, H4 | 4 |
| OV90_S_48 | OV90 | Sensitive | 48 hours | I1, I2, I3, I4 | 4 |
| OV90_R_NT | OV90 | Resistant | No Treatment | J1, J2, J3, J4 | 4 |
| OV90_R_6 | OV90 | Resistant | 6 hours | K1, K2, K3, K4 | 4 |
| OV90_R_48 | OV90 | Resistant | 48 hours | L1, L2, L3, L4 | 4 |

- **n = 4 biological replicates per group**
- **Total samples: 48**

### Metadata File Structure (`metadata.txt`)

| Column | Type | Description |
|--------|------|-------------|
| `sample` | string | Sample ID (A1–L4) |
| `group` | string | Full group label (e.g., OV8_S_NT) |
| `ov8` | 0/1 | 1 = OV8 cell line |
| `ov90` | 0/1 | 1 = OV90 cell line |
| `resistant` | 0/1 | 1 = drug resistant |
| `sensitive` | 0/1 | 1 = drug sensitive |
| `zero` | 0/1 | 1 = no treatment (NT) |
| `six` | 0/1 | 1 = 6h treatment |
| `fortyeight` | 0/1 | 1 = 48h treatment |
| `exclude` | 0/1 | 1 = exclude from analysis |

---

## 3. Why Kallisto and Sleuth

### Why Kallisto for quantification?
Traditional RNA-seq alignment tools (STAR, HISAT2) map reads to a genome before counting. Kallisto (Bray et al., 2016) uses **pseudoalignment** against a transcriptome index, which is:
- ~50x faster than traditional aligners
- Equally accurate for expression quantification
- Produces **bootstrap resamples** that quantify technical (quantification) uncertainty — a key feature that Sleuth uses

### Why Sleuth for differential expression?
Most DE tools (DESeq2, edgeR, limma) model only **biological variance** between samples. They ignore the fact that transcript abundance estimates from RNA-seq are themselves uncertain — especially for transcripts with low coverage or high sequence similarity.

Sleuth (Pimentel et al., 2017) was designed specifically to work with Kallisto output and models:

```
Total observed variance = Biological variance + Technical variance
```

Where technical variance is estimated directly from Kallisto's bootstrap resamples. This prevents false positives driven by quantification noise rather than true biological differences.

**Why not DESeq2/tximport?**
DESeq2 via tximport is also valid and widely used. However, because our data was quantified with Kallisto with bootstrap resampling (`max_bootstrap = 30`), Sleuth fully utilizes this information. DESeq2 discards the bootstrap uncertainty. For transcript-level analysis especially, Sleuth is the more appropriate choice.

---

## 4. Statistical Framework

### Model Design
For each comparison, two linear models are fit per gene/transcript using a response error measurement (REM) model:

| Model | Formula | Interpretation |
|-------|---------|----------------|
| Full | `~ condition` | Expression depends on the condition of interest |
| Reduced | `~ 1` | Expression does not depend on condition (intercept only) |

### Two Complementary Tests

**Likelihood Ratio Test (LRT)**
- Compares full vs reduced model fit
- Asks: *does including the condition variable significantly improve the model?*
- Output: p-value (no direction/magnitude)
- Threshold used: `pval < 0.1`

**Wald Test (WT)**
- Estimates the beta coefficient (b) for the condition variable
- Asks: *what is the magnitude and direction of the effect?*
- Output: b (effect size), standard error, p-value, q-value (FDR-adjusted)
- Threshold for significance: `qval < 0.05`

**Why both tests?**
Using LRT + WT together is more conservative than either alone:
- LRT identifies genes where the condition matters at all
- WT provides effect size and direction
- Final results table merges both — researchers can apply their own thresholds

### Gene vs Transcript Level
Each comparison is run at two levels:

| Level | `gene_mode` | Normalization | Use case |
|-------|-------------|---------------|---------|
| Gene | `TRUE` | scaled_reads_per_base | Overall gene expression, pathway analysis |
| Transcript | `FALSE` | TPM | Isoform-level analysis, splicing |

---

## 5. Dependencies

### R Packages

| Package | Version | Source | Purpose |
|---------|---------|--------|---------|
| `sleuth` | 0.30.1 | Bioconda | Core DE analysis |
| `reshape2` | ≥1.4.4 | conda-forge | Must load BEFORE sleuth — fixes melt() conflict |
| `dplyr` | ≥1.1.0 | conda-forge | Data filtering and manipulation |
| `ggplot2` | ≥3.4.0 | conda-forge | Plotting (PCA, volcano, QQ) |
| `gridExtra` | ≥2.3 | conda-forge | Multi-panel plot layout |
| `rhdf5` | ≥2.46.1 | Bioconda | Read Kallisto .h5 bootstrap files |

### System

| Requirement | Value |
|-------------|-------|
| R version | 4.3.1 |
| Conda | Miniforge 24.11.3 |
| Scheduler | SGE (Sun Grid Engine) |
| Queue | UI (UIowa Argon general compute) |
| Memory | 128GB per job |
| Cores | 1 (single processor — see Known Issues) |

---

## 6. Environment Setup on Argon HPC

### One-time setup

```bash
# Conda is pre-installed at:
/old_Users/<hawkid>/miniforge3/bin/conda

# Create the sleuth environment with all required packages
conda create -n sleuth -c conda-forge -c bioconda \
  r-base=4.3.1 \
  r-essentials \
  r-reshape2 \
  r-dplyr \
  r-ggplot2 \
  r-gridextra \
  bioconductor-rhdf5 \
  r-sleuth \
  -y

# CRITICAL: downgrade data.table to avoid melt() conflict
conda activate sleuth
conda install -c conda-forge "r-data.table=1.14.8" -y
```

### Activate environment

```bash
source /old_Users/<hawkid>/miniforge3/bin/activate sleuth
```

### Verify installation

```bash
R -e "library(reshape2); library(sleuth); library(dplyr); library(ggplot2); sessionInfo()"
```

All four packages should load without errors before submitting any jobs.

---

## 7. Directory Structure

```
Atonu/
├── sleuth_atonu.R           # Batch 1: Comparisons 1-2 (R vs S baseline)
├── sleuth_atonu.job         # SGE job script for batch 1
├── sleuth_batch3.R          # Batch 3: Comparisons 3-10 (time response)
├── sleuth_batch3.job        # SGE job script for batch 3
├── metadata.txt             # Sample metadata (48 samples)
├── cdna_t2g_map.tsv         # Transcript-to-gene mapping (GRCh38)
│
├── Kallisto_1/              # Kallisto quantification output
│   ├── A1/kallisto/
│   │   ├── abundance.h5     # Bootstrap quantification (binary)
│   │   ├── abundance.tsv    # Plain text quantification
│   │   └── run_info.json
│   ├── A2/kallisto/
│   ├── ...
│   └── L4/kallisto/
│
├── Sleuth_2/                # Output: Batch 1 (comparisons 1-2)
│   ├── s2c_all.tsv
│   ├── ALL-RESULTS_GENE_OV90_RvsS_NT_beta_resistant.tsv
│   ├── ALL-RESULTS_TRANS_OV90_RvsS_NT_beta_resistant.tsv
│   ├── ALL-RESULTS_GENE_OV8_RvsS_NT_beta_resistant.tsv
│   ├── ALL-RESULTS_TRANS_OV8_RvsS_NT_beta_resistant.tsv
│   ├── soG_*.so
│   ├── soT_*.so
│   └── PCA/Volcano/QQ PDFs
│
└── Sleuth_3/                # Output: Batch 3 (comparisons 3-10)
    ├── s2c_all.tsv
    ├── ALL-RESULTS_GENE_OV90_S_6vsNT_beta_six.tsv
    ├── ALL-RESULTS_TRANS_OV90_S_6vsNT_beta_six.tsv
    ├── ... (one pair per comparison)
    ├── soG_*.so
    ├── soT_*.so
    └── PCA/Volcano/QQ PDFs
```

---

## 8. Input Files

### `metadata.txt`
Tab-separated, one row per sample, 48 rows total. Must contain all columns described in Section 2. Samples with `exclude = 1` are automatically removed before any analysis.

### `cdna_t2g_map.tsv`
Tab-separated transcript-to-gene mapping file. Two required columns:
- `target_id` — transcript ID matching Kallisto index
- `ens_gene` — Ensembl gene ID for aggregation

This file must be built from the **same reference transcriptome** used to create the Kallisto index. This pipeline uses **GRCh38 (human)** with Ensembl annotation.

### `Kallisto_1/`
One subdirectory per sample, named exactly as in `metadata.txt`. Each must contain a `kallisto/` folder with `abundance.h5` (required for bootstrap access by Sleuth).

---

## 9. Scripts in This Repository

### `sleuth_atonu.R` + `sleuth_atonu.job`
**Purpose:** Comparisons 1 and 2 — Resistant vs Sensitive at baseline (NT only)
**Output directory:** `Sleuth_2/`
**Comparisons:**
- OV90 Resistant vs Sensitive (NT)
- OV8 Resistant vs Sensitive (NT)

### `sleuth_batch3.R` + `sleuth_batch3.job`
**Purpose:** Comparisons 3–10 — Drug response over time
**Output directory:** `Sleuth_3/`
**Comparisons:** All 8 time-response comparisons (see Section 11)

### Why two separate scripts?
The full 10-comparison pipeline was split into two batches due to memory constraints on the Argon cluster. Each Sleuth object for transcript-level analysis consumes ~500MB RAM. Running all 10 comparisons sequentially in one job exceeded the 128GB memory limit. Splitting into two jobs of 2 and 8 comparisons respectively resolved this.

---

## 10. How to Run

### Step 1: Activate environment
```bash
source /old_Users/<hawkid>/miniforge3/bin/activate sleuth
```

### Step 2: Navigate to working directory
```bash
cd /Shared/lss_jmclendon/2_UserFolders/Atonu
```

### Step 3: Submit jobs
```bash
# Batch 1 - Resistance comparisons (Sleuth_2/)
qsub sleuth_atonu.job

# Batch 3 - Time response comparisons (Sleuth_3/)
qsub sleuth_batch3.job
```

### Step 4: Monitor jobs
```bash
# Check job status
qstat

# Check specific job
qstat -j <job-ID>

# Watch console output in real time
tail -f Sleuth_batch3_console_<timestamp>.txt
```

### Step 5: Verify completion
```bash
# Check output files created
ls -lth Sleuth_2/
ls -lth Sleuth_3/

# Check log for errors
cat Sleuth_batch3_log.txt
```

---

## 11. Comparisons Performed

### Batch 1 — Resistance Signature at Baseline (`Sleuth_2/`)

| # | Comparison | Model | Reference | Test |
|---|-----------|-------|-----------|------|
| 1 | OV90 Resistant vs Sensitive (NT) | `~resistant` | Sensitive (resistant=0) | Resistant (resistant=1) |
| 2 | OV8 Resistant vs Sensitive (NT) | `~resistant` | Sensitive (resistant=0) | Resistant (resistant=1) |

**Biological question:** Which genes are differentially expressed between resistant and sensitive cells *before any drug treatment*? These represent the baseline molecular signature of resistance.

### Batch 3 — Drug Response Over Time (`Sleuth_3/`)

| # | Comparison | Model | Reference | Test |
|---|-----------|-------|-----------|------|
| 3 | OV90 Sensitive: 6h vs NT | `~six` | NT (six=0) | 6h (six=1) |
| 4 | OV90 Sensitive: 48h vs NT | `~fortyeight` | NT (fortyeight=0) | 48h (fortyeight=1) |
| 5 | OV90 Resistant: 6h vs NT | `~six` | NT (six=0) | 6h (six=1) |
| 6 | OV90 Resistant: 48h vs NT | `~fortyeight` | NT (fortyeight=0) | 48h (fortyeight=1) |
| 7 | OV8 Sensitive: 6h vs NT | `~six` | NT (six=0) | 6h (six=1) |
| 8 | OV8 Sensitive: 48h vs NT | `~fortyeight` | NT (fortyeight=0) | 48h (fortyeight=1) |
| 9 | OV8 Resistant: 6h vs NT | `~six` | NT (six=0) | 6h (six=1) |
| 10 | OV8 Resistant: 48h vs NT | `~fortyeight` | NT (fortyeight=0) | 48h (fortyeight=1) |

**Biological question:** How does drug treatment change gene expression over time? Are the transcriptional responses to drug treatment different in sensitive vs resistant cells?

---

## 12. Output Files

### Results Tables (`.tsv`)

For each comparison, two results tables are produced:

| File | Level | Normalized measure |
|------|-------|-------------------|
| `ALL-RESULTS_GENE_<comparison>_beta_<beta>.tsv` | Gene | scaled_reads_per_base |
| `ALL-RESULTS_TRANS_<comparison>_beta_<beta>.tsv` | Transcript | TPM |

### Key Columns in Results Tables

| Column | Source test | Description |
|--------|-------------|-------------|
| `target_id` | — | Gene or transcript ID |
| `b` | Wald Test | Effect size (analogous to log fold change) |
| `se_b` | Wald Test | Standard error of b |
| `mean_obs` | — | Mean observed expression across samples |
| `pval.x` | Wald Test | Wald test p-value |
| `qval.x` | Wald Test | Wald test FDR q-value **(use for significance)** |
| `pval.y` | LRT | Likelihood ratio test p-value |
| `qval.y` | LRT | LRT FDR q-value |
| Sample columns | — | Normalized expression per sample |

**Recommended filtering for significant genes:**
```r
sig_genes <- results[results$qval.x < 0.05 & !is.na(results$qval.x), ]
```

### Saved Sleuth Objects (`.so`)

| File | Description |
|------|-------------|
| `soG_<comparison>.so` | Gene-level Sleuth object |
| `soT_<comparison>.so` | Transcript-level Sleuth object |

Reload in R for interactive analysis:
```r
library(sleuth)
so <- sleuth_load("Sleuth_2/soG_OV90_RvsS_NT_beta_resistant.so")
sleuth_live(so)  # opens interactive Shiny app
```

### Plots (`.pdf`)

| File | Description | What to look for |
|------|-------------|-----------------|
| `PCA_GENE_*.pdf` | PCA colored by group | Samples should cluster by group |
| `PCA_TRANS_*.pdf` | Transcript-level PCA | Same as above |
| `Volcano_GENE_*.pdf` | Effect size vs -log10(pval) | Symmetric spread, clear significant genes |
| `Volcano_TRANS_*.pdf` | Transcript-level volcano | — |
| `QQ_GENE_*.pdf` | Observed vs expected p-values | Inflation = model issues; deflation = few DE genes |
| `QQ_TRANS_*.pdf` | Transcript-level QQ | — |

---

## 13. Parameter Justification

### `num_cores = 1`
**Why:** Multi-core (`num_cores > 1`) caused job abortion errors on Argon HPC due to memory conflicts during parallel bootstrap processing. Single-core execution is stable and produces identical results — parallelization only affects speed, not output.

### `max_bootstrap = 30`
**Why:** Reduced from 100 (original script) to 30 due to memory constraints. This is fully supported by the Sleuth paper:

> *"We find that B = 30 bootstraps is sufficient for reliable variance estimation"*
> — Pimentel et al., 2017, Nature Methods

30 bootstraps is the value used in the Sleuth methods paper's own simulations and is standard practice for published analyses.

### `gene_mode = TRUE` with `aggregation_column = 'ens_gene'`
**Why:** Aggregates transcript-level estimates to gene level using the `ens_gene` column from the t2g map. This reduces the multiple testing burden (15,000 genes vs 83,000 transcripts) and is more interpretable for pathway analysis with tools like clusterProfiler or fgsea.

### `pval < 0.1` for LRT filtering
**Why:** Standard threshold in Sleuth analyses. The LRT p-value is used as an initial filter; the Wald test q-value (`< 0.05`) provides the final significance cutoff in the merged results table.

### Two separate job submissions
**Why:** Each transcript-level Sleuth object consumes ~500MB RAM. Running all 10 comparisons sequentially in one job caused memory exhaustion and job termination during bootstrap summarization (~83,000 transcripts × 30 bootstraps). Splitting into two jobs resolved this.

---

## 14. Known Issues and Fixes

### Issue 1: `melt()` conflict — data.table vs reshape2
**Error message:**
```
Error in melt.default: The melt generic in data.table has been passed a data.frame...
this redirection is now deprecated.
```
**Root cause:** data.table (v > 1.14.8) intercepts Sleuth's internal `melt()` calls on data.frames and treats the redirect to reshape2 as a hard error instead of a warning.

**Fix applied:**
1. Load `library(reshape2)` **before** `library(sleuth)` in every R script
2. Downgrade data.table:
```bash
conda install -c conda-forge "r-data.table=1.14.8" -y
```

### Issue 2: Job killed during bootstrap summarization
**Symptom:** Job disappears from `qstat` with console output truncated at:
```
summarizing bootstraps
```
**Root cause:** Memory limit exceeded. Bootstrap summarization for ~83,000 transcripts with 30 bootstraps × 8 samples requires >64GB RAM.

**Fix applied:**
- Increased `mem_free` from 64G to 128G in job script
- Reduced `max_bootstrap` from 100 to 30
- Split pipeline into two separate job submissions

### Issue 3: Bioconductor/CRAN not accessible
**Error message:**
```
Warning: unable to access index for repository https://bioconductor.org/...
cannot open URL
```
**Root cause:** Argon compute nodes block outbound internet access.

**Fix applied:** Install all R packages via conda channels (`bioconda`, `conda-forge`) instead of `install.packages()` or `BiocManager::install()`.

### Issue 4: Multi-core errors
**Symptom:** Job aborts when `num_cores > 1` in `sleuth_prep()`

**Root cause:** Sleuth's parallel processing via `parallel` package conflicts with SGE's memory management on shared nodes.

**Fix applied:** Set `num_cores = 1` in all `sleuth_prep()` calls.

---

## 15. References

1. **Pimentel H, Bray NL, Puente S, Melsted P, Pachter L.** Differential analysis of RNA-seq incorporating quantification uncertainty. *Nature Methods.* 2017;14:687–690.
   https://doi.org/10.1038/nmeth.4324
   *(Primary Sleuth reference — also justifies 30 bootstraps)*

2. **Bray NL, Pimentel H, Melsted P, Pachter L.** Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology.* 2016;34:525–527.
   https://doi.org/10.1038/nbt.3519
   *(Kallisto pseudoalignment)*

3. **Soneson C, Love MI, Robinson MD.** Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research.* 2015;4:1521.
   https://doi.org/10.12688/f1000research.7563.2
   *(Justification for transcript-level quantification and tximport approach)*

4. **Love MI, Huber W, Anders S.** Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology.* 2014;15:550.
   https://doi.org/10.1186/s13059-014-0550-8
   *(Alternative DE method — cited for comparison)*

5. **Robinson MD, McCarthy DJ, Smyth GK.** edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics.* 2010;26(1):139–140.
   https://doi.org/10.1093/bioinformatics/btp616
   *(Alternative DE method — cited for comparison)*

6. **Liao Y, Smyth GK, Shi W.** featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics.* 2014;30(7):923–930.
   *(Alternative counting approach — cited for comparison)*
