# Sleuth Differential Expression Analysis
### Ovarian Cancer Cell Lines — OV8 & OV90 (Sensitive vs Resistant)

**Original Script Author:** Jared McLendon, McLendon Lab, University of Iowa  
**Modified:** Atonu Chakrabortty, McLendon Lab, University of Iowa  
**Date:** April 05, 2026  
**Platform:** UIowa Argon HPC Cluster

---

## Table of Contents
1. [Scientific Background](#1-scientific-background)
2. [Experimental Design](#2-experimental-design)
3. [Analysis Overview](#3-analysis-overview)
4. [Dependencies](#4-dependencies)
5. [Environment Setup](#5-environment-setup)
6. [Directory Structure](#6-directory-structure)
7. [Input Files](#7-input-files)
8. [How to Run](#8-how-to-run)
9. [Comparisons Performed](#9-comparisons-performed)
10. [Output Files](#10-output-files)
11. [Parameter Justification](#11-parameter-justification)
12. [Known Issues and Fixes](#12-known-issues-and-fixes)
13. [Development Journey](#13-development-journey)
14. [References](#14-references)

---

## 1. Scientific Background

Ovarian cancer is one of the most lethal gynecological malignancies, largely due to the development of resistance to platinum-based chemotherapy. Understanding the transcriptomic differences between drug-sensitive and drug-resistant ovarian cancer cells is critical for identifying molecular mechanisms of resistance and potential therapeutic targets.

This pipeline performs RNA-seq differential expression analysis on two ovarian cancer cell lines:

- **OV8** — ovarian cancer cell line
- **OV90** — ovarian cancer cell line

For each cell line, sensitive and resistant variants were profiled at multiple treatment time points to capture both **baseline resistance signatures** and **dynamic transcriptional responses to drug treatment**.

RNA-seq reads were previously quantified using **Kallisto** (Bray et al., 2016), which uses pseudoalignment to rapidly estimate transcript-level abundances with bootstrap resampling for uncertainty quantification. This pipeline uses **Sleuth** (Pimentel et al., 2017) for downstream differential expression testing, which is specifically designed to work with Kallisto output and properly propagates quantification uncertainty into the statistical model.

---

## 2. Experimental Design

### Cell Lines and Conditions

| Group Label | Cell Line | Condition | Treatment | Samples |
|-------------|-----------|-----------|-----------|---------|
| OV8_S_NT | OV8 | Sensitive | No Treatment (0h) | A1–A4 |
| OV8_S_6 | OV8 | Sensitive | 6 hours | B1–B4 |
| OV8_S_48 | OV8 | Sensitive | 48 hours | E1–E4 |
| OV8_R_NT | OV8 | Resistant | No Treatment (0h) | C1–C4 |
| OV8_R_6 | OV8 | Resistant | 6 hours | D1–D4 |
| OV8_R_48 | OV8 | Resistant | 48 hours | F1–F4 |
| OV90_S_NT | OV90 | Sensitive | No Treatment (0h) | G1–G4 |
| OV90_S_6 | OV90 | Sensitive | 6 hours | H1–H4 |
| OV90_S_48 | OV90 | Sensitive | 48 hours | I1–I4 |
| OV90_R_NT | OV90 | Resistant | No Treatment (0h) | J1–J4 |
| OV90_R_6 | OV90 | Resistant | 6 hours | K1–K4 |
| OV90_R_48 | OV90 | Resistant | 48 hours | L1–L4 |

- **n = 4 biological replicates per group**
- **Total samples: 48**
- All samples with `exclude = 0` in metadata are included in analysis

### Metadata Columns

| Column | Description |
|--------|-------------|
| `sample` | Sample ID (A1–L4) |
| `group` | Group label |
| `ov8` | 1 = OV8 cell line |
| `ov90` | 1 = OV90 cell line |
| `resistant` | 1 = resistant condition |
| `sensitive` | 1 = sensitive condition |
| `zero` | 1 = no treatment (NT) |
| `six` | 1 = 6h treatment |
| `fortyeight` | 1 = 48h treatment |
| `exclude` | 1 = exclude sample from analysis |

---

## 3. Analysis Overview

### Why Sleuth?

When you run Kallisto, it doesn't just give you one estimate of transcript abundance — it runs bootstrap resampling internally to capture how uncertain each estimate is. That uncertainty is real and matters: two transcripts can have the same mean expression level but very different confidence around that mean depending on how well the reads aligned.

Standard DE tools like DESeq2 and edgeR only see a single count matrix. They have no idea about that uncertainty — it gets lost. Sleuth is designed specifically to take Kallisto's bootstrap output and use it. It models variance as two separate components:

```
Total variance = Biological variance + Technical variance (from bootstraps)
```

This means if a gene looks differentially expressed but the signal is mostly coming from technical quantification noise rather than genuine biology, Sleuth can flag that. It's a more honest model for pseudoalignment-based data.

### Statistical Framework

For each comparison, Sleuth fits two linear models per gene/transcript:

| Model | Formula | Meaning |
|-------|---------|---------|
| Full | `~ condition` | Expression depends on condition |
| Reduced | `~ 1` | Expression does not depend on condition (intercept only) |

Two complementary tests are run:

- **Likelihood Ratio Test (LRT):** Compares the full model to the reduced model. Asks — *does knowing the condition actually help explain the expression pattern?* Returns a p-value but no direction or magnitude.
- **Wald Test (WT):** Estimates the **b coefficient** for the condition variable — this is effectively the effect size (analogous to log fold change). Gives direction and magnitude of change, plus its own p-value and q-value.

Running both is intentional. The LRT is better for identifying genes where the model fit improves (broad signal), while the Wald test tells you how big the effect is and in which direction. Results from both tests are merged into a single output table. Genes significant in both are the most reliable hits.

### Analysis Levels

Each comparison is run at **two levels**:

1. **Gene level** (`gene_mode = TRUE`) — transcripts aggregated to gene level using `ens_gene` column from the t2g map. Uses `scaled_reads_per_base` as the normalized expression measure. Easier to interpret and more statistically powerful due to reduced multiple testing burden.
2. **Transcript level** (`gene_mode = FALSE`) — individual transcript results. Uses `TPM` as the normalized expression measure. Useful for isoform-level analysis.

### All-Samples Setup Block

Before any differential expression testing, the script runs a global setup on all 48 samples together (no filtering). This produces:
- Normalized expression matrices for all samples (gene and transcript level)
- Raw count matrices
- PCA plot of all samples colored by group
- Sample-level heatmap with clustering

This is useful for quality control — you can check whether samples cluster as expected by group, spot any outliers, and get an overview of the whole dataset before looking at individual comparisons.

---

## 4. Dependencies

### R Packages

| Package | Version | Source | Purpose |
|---------|---------|--------|---------|
| `sleuth` | 0.30.1 | Bioconda | Core DE analysis |
| `reshape2` | ≥1.4.4 | CRAN/conda-forge | Required to fix melt() conflict with data.table |
| `dplyr` | ≥1.0.0 | conda-forge | Data manipulation |
| `ggplot2` | ≥3.4.0 | conda-forge | Plotting |
| `gridExtra` | ≥2.3 | conda-forge | Multi-panel plot layout |

### System Requirements

- R version 4.3.1
- Miniforge/Conda
- SGE job scheduler (UIowa Argon HPC)
- **256GB RAM** (high-memory node, UI-HM queue)
- 4 CPU cores

---

## 5. Environment Setup

### One-time setup on Argon HPC

```bash
# Conda is available at:
/old_Users/<hawkid>/miniforge3/bin/conda

# Create the sleuth environment
conda create -n sleuth -c conda-forge -c bioconda \
  r-base=4.3.1 \
  r-essentials \
  r-devtools \
  r-ggplot2 \
  r-dplyr \
  r-gridextra \
  r-reshape2 \
  bioconductor-rhdf5 \
  r-sleuth \
  -y

# IMPORTANT: downgrade data.table to avoid melt() conflict with reshape2
conda activate sleuth
conda install -c conda-forge "r-data.table=1.14.8" -y
```

> **Note on biomaRt:** The original professor's script imported `library(biomaRt)` for gene annotation lookups. This was removed because: (1) biomaRt requires outbound internet access to query Ensembl, which Argon compute nodes block; (2) the t2g mapping file (`cdna_t2g_map.tsv`) already contains all the gene annotation needed — biomaRt was redundant. Do not attempt to install or use biomaRt on Argon compute nodes.

### Verify installation

```bash
source /old_Users/<hawkid>/miniforge3/bin/activate sleuth
R -e "library(sleuth); library(reshape2); library(dplyr); library(ggplot2); sessionInfo()"
```

---

## 6. Directory Structure

```
Atonu/
├── sleuth_run.R            # Parameterized R script — takes comp_index 1-10 as argument
├── sleuth_master.job       # SGE job script — loops through all 10 comparisons
├── metadata.txt            # Sample metadata table
├── cdna_t2g_map.tsv        # Transcript-to-gene mapping file (GRCh38)
├── Kallisto_1/             # Kallisto quantification output
│   ├── A1/
│   │   └── kallisto/
│   │       ├── abundance.h5      # bootstrap quantification (required for Sleuth)
│   │       ├── abundance.tsv
│   │       └── run_info.json
│   ├── A2/
│   │   └── kallisto/
│   ├── ...
│   └── L4/
│       └── kallisto/
├── Sleuth_4/               # Current output directory (created by script)
│   ├── s2c_all.tsv                          # Final sample table used
│   ├── s2mG_all_samples.tsv                 # All-samples gene expression (scaled_reads_per_base)
│   ├── s2mG_count_all_samples.tsv           # All-samples gene raw counts
│   ├── s2mT_TPM_all_samples.tsv             # All-samples transcript TPM
│   ├── s2mT_count_all_samples.tsv           # All-samples transcript raw counts
│   ├── PCA_soG_all_samples_group.pdf        # PCA of all samples (gene level)
│   ├── heatmap_soG_all_samples.pdf          # Sample heatmap
│   ├── PCA_soT_all_samples_group.pdf        # PCA of all samples (transcript level)
│   ├── ALL-RESULTS_GENE_samp_*_beta_*.tsv   # Gene-level DE results per comparison
│   ├── ALL-RESULTS_TRANS_samp_*_beta_*.tsv  # Transcript-level DE results per comparison
│   ├── soG_samp_*_beta_*.so                 # Saved gene-level Sleuth objects
│   ├── soT_samp_*_beta_*.so                 # Saved transcript-level Sleuth objects
│   ├── PCAplot_Group_*.pdf
│   ├── Volcano_GENE_*.pdf
│   ├── Volcano_Transcript_*.pdf
│   ├── QQ_GENE_*.pdf
│   └── QQ_Transcript_*.pdf
├── sleuth_comp1_console.txt    # Per-comparison R console output (for debugging)
├── sleuth_comp2_console.txt
├── ...
└── sleuth_master_log.txt       # SGE job log
```

Previous output directories (for reference):
- `Sleuth_2/` — Comparisons 1 & 2 completed (OV90_RvsS_NT, OV8_RvsS_NT)
- `Sleuth_3/` — Time response comparisons (failed due to memory — see Development Journey)

---

## 7. Input Files

### `metadata.txt`
Tab-separated file with one row per sample. Must contain columns:
`sample`, `group`, `ov8`, `ov90`, `resistant`, `sensitive`, `zero`, `six`, `fortyeight`, `exclude`

### `cdna_t2g_map.tsv`
Tab-separated transcript-to-gene mapping file with columns:
- `target_id` — Kallisto transcript ID
- `ens_gene` — Ensembl gene ID

This file must match the reference transcriptome used to build the Kallisto index. In this pipeline, the human reference genome **GRCh38** with Ensembl annotation was used.

### `Kallisto_1/`
Each sample subdirectory must contain a `kallisto/` folder with:
- `abundance.h5` — bootstrap quantification (binary, required for Sleuth)
- `abundance.tsv` — plain text quantification
- `run_info.json` — Kallisto run metadata

---

## 8. How to Run

### Submit the job

```bash
cd /Shared/lss_jmclendon/2_UserFolders/Atonu

# Copy updated scripts if coming from local machine
scp sleuth_run.R sleuth_master.job \
    HawkID@argon.hpc.uiowa.edu:/Shared/lss_jmclendon/2_UserFolders/Atonu/

# Submit
qsub sleuth_master.job
```

### What happens when you submit

1. SGE allocates a node on the **UI-HM** (high-memory) queue with 4 cores and 256GB RAM
2. The job loops `i` from 1 to 10
3. For each `i`, it calls:
   ```bash
   Rscript --vanilla sleuth_run.R $i Sleuth_4
   ```
4. Each comparison runs as a **completely fresh R process** — this is intentional. R does not release memory well after large objects, so running each comparison in its own process prevents memory accumulation across 10 comparisons
5. Console output for each comparison is captured to `sleuth_comp${i}_console.txt` for debugging
6. Comparison 1 additionally runs the all-samples setup block (PCA, heatmap, count matrices) before its DE analysis

### Monitor the job

```bash
# Check if running
qstat -u atonu-chakrabortty

# Check details including which compute node was assigned
qstat -j <job-id>

# Watch the master log in real time
tail -f /Shared/lss_jmclendon/2_UserFolders/Atonu/sleuth_master_log.txt

# Check progress of a specific comparison
cat sleuth_comp3_console.txt
```

### Job script (`sleuth_master.job`) — key settings

```bash
#!/bin/bash
#$ -q UI-HM          # High-memory queue (nodes with 1-2TB RAM)
#$ -pe smp 4         # 4 cores on one node
#$ -l mem_free=256G  # 256GB RAM allocation
#$ -l h_rt=200:00:00 # 200 hour wall time limit
```

> **Important:** Use the **full path** to Rscript. Do not rely on `conda activate` to set the PATH in SGE jobs — it does not work reliably in non-interactive shells. The full path is:
> `/old_Users/chakrabortty/miniforge3/envs/sleuth/bin/Rscript`

---

## 9. Comparisons Performed

### Category 1: Resistance Signature (Baseline)
Untreated samples only. Tests which genes are differentially expressed between resistant and sensitive cells **before any drug treatment**. This identifies the pre-existing transcriptomic signature of resistance.

| # | Comparison | Model | Samples Included |
|---|-----------|-------|---------|
| 1 | OV90 Resistant vs Sensitive (NT) | `~resistant` | OV90 zero==1 (J1–J4 vs G1–G4) |
| 2 | OV8 Resistant vs Sensitive (NT) | `~resistant` | OV8 zero==1 (C1–C4 vs A1–A4) |

### Category 2: Drug Response Over Time
Tests how gene expression changes in response to drug treatment **within** each group. Reference is always untreated (NT, 0h). This identifies dynamic transcriptional responses to treatment and whether sensitive vs resistant cells respond differently over time.

| # | Comparison | Model | Samples Included |
|---|-----------|-------|---------|
| 3 | OV90 Sensitive: 6h vs NT | `~six` | OV90 sensitive, zero==1 or six==1 |
| 4 | OV90 Sensitive: 48h vs NT | `~fortyeight` | OV90 sensitive, zero==1 or fortyeight==1 |
| 5 | OV90 Resistant: 6h vs NT | `~six` | OV90 resistant, zero==1 or six==1 |
| 6 | OV90 Resistant: 48h vs NT | `~fortyeight` | OV90 resistant, zero==1 or fortyeight==1 |
| 7 | OV8 Sensitive: 6h vs NT | `~six` | OV8 sensitive, zero==1 or six==1 |
| 8 | OV8 Sensitive: 48h vs NT | `~fortyeight` | OV8 sensitive, zero==1 or fortyeight==1 |
| 9 | OV8 Resistant: 6h vs NT | `~six` | OV8 resistant, zero==1 or six==1 |
| 10 | OV8 Resistant: 48h vs NT | `~fortyeight` | OV8 resistant, zero==1 or fortyeight==1 |

---

## 10. Output Files

All final output files are in `Sleuth_4/`.

### All-Samples QC Files (generated once, comparison 1)

| File | Description |
|------|-------------|
| `s2c_all.tsv` | The exact sample table passed to Sleuth — useful for verifying which samples were included |
| `s2mG_all_samples.tsv` | Gene-level normalized expression (scaled_reads_per_base) across all 48 samples |
| `s2mG_count_all_samples.tsv` | Gene-level raw estimated counts across all samples |
| `s2mT_TPM_all_samples.tsv` | Transcript-level TPM across all samples |
| `s2mT_count_all_samples.tsv` | Transcript-level raw estimated counts |
| `PCA_soG_all_samples_group.pdf` | PCA plot — all samples, gene level, colored by group |
| `heatmap_soG_all_samples.pdf` | Hierarchical clustering heatmap of all samples |
| `PCA_soT_all_samples_group.pdf` | PCA plot — all samples, transcript level |

### Differential Expression Results Tables

| File | Level | Contents |
|------|-------|---------|
| `ALL-RESULTS_GENE_samp_<comparison>_beta_<beta>.tsv` | Gene | Combined Wald + LRT results + per-sample normalized counts |
| `ALL-RESULTS_TRANS_samp_<comparison>_beta_<beta>.tsv` | Transcript | Combined Wald + LRT results + per-sample TPM |

**Key columns in results tables:**

| Column | Source | Description |
|--------|--------|-------------|
| `target_id` | Both | Gene or transcript ID |
| `b` | Wald Test | Effect size (log fold change equivalent) |
| `se_b` | Wald Test | Standard error of b |
| `pval.x` | Wald Test | Wald test p-value |
| `qval.x` | Wald Test | Wald test FDR-adjusted q-value |
| `pval.y` | LRT | Likelihood ratio test p-value |
| `qval.y` | LRT | LRT FDR-adjusted q-value |
| `mean_obs` | Both | Mean observed expression |
| Sample columns | Both | Per-sample normalized expression |

### Saved Sleuth Objects

| File | Description |
|------|-------------|
| `soG_samp_<comparison>_beta_<beta>.so` | Gene-level Sleuth object |
| `soT_samp_<comparison>_beta_<beta>.so` | Transcript-level Sleuth object |

These can be reloaded in R for interactive exploration:
```r
library(sleuth)
so <- sleuth_load("Sleuth_4/soG_samp_OV90_RvsS_NT_beta_resistant_.so")
sleuth_live(so)  # opens interactive Shiny browser
```

### Plots (per comparison)

| File | Description |
|------|-------------|
| `PCAplot_Group_<comparison>_beta_<beta>.pdf` | PCA of comparison-specific samples |
| `Volcano_GENE_<comparison>_beta_<beta>.pdf` | Volcano plot — gene level |
| `Volcano_Transcript_<comparison>_beta_<beta>_transcript.pdf` | Volcano plot — transcript level |
| `QQ_GENE_<comparison>_beta_<beta>.pdf` | QQ plot — gene level (check model fit) |
| `QQ_Transcript_<comparison>_beta_<beta>_transcript.pdf` | QQ plot — transcript level |

---

## 11. Parameter Justification

### `num_cores = 4`
Sleuth's `num_cores` parameter parallelizes the bootstrap resampling calculations. The job requests `smp 4` from SGE, allocating 4 cores on a single node. Setting `num_cores = 4` in R matches this allocation so all 4 cores are used for bootstrap summarization — the most computationally intensive step.

> **Earlier attempts used `num_cores = 1`:** During early runs on the standard `UI` queue, `num_cores > 1` caused job aborts. This was specific to the queue/node configuration at the time. On `UI-HM` nodes with `smp 4`, `num_cores = 4` works correctly and speeds up the bootstrap step.

### `max_bootstrap = 100`
The original professor's script used 100 bootstraps. This is the recommended value in the Sleuth paper and provides reliable variance estimation. It was temporarily reduced to 30 during early debugging to reduce memory load, but with the move to the UI-HM queue (256GB) it is restored to 100.

### Queue: `UI-HM`
The standard `UI` queue nodes have limited memory relative to the needs of this analysis. Sleuth holds bootstrap data for ~83,000 transcripts across 8–12 samples in memory simultaneously. At 100 bootstraps, this can exceed 128GB for the larger comparisons. The `UI-HM` queue nodes provide 1–2TB of RAM, making 256GB a comfortable allocation with headroom.

### `gene_mode = TRUE` with `aggregation_column = 'ens_gene'`
Aggregates transcript-level estimates to gene level using the t2g mapping. Reduces the multiple testing burden from ~83,000 transcripts to ~20,000 genes, and is more interpretable for pathway analysis downstream. Both gene and transcript level are run for completeness.

### `pval < 0.1` threshold for LRT
Standard threshold in Sleuth analyses. The LRT p-value filter is intentionally more lenient than 0.05 because the Wald test q-value (FDR-corrected) provides the final significance filter. Looking at `shared_results` (genes significant in both LRT and WT) is the most conservative approach.

### Running comparisons sequentially (not parallel)
It would be possible to submit 10 separate jobs in parallel. Instead, all 10 run sequentially within one job. The reason: each comparison uses up to 256GB. Running them in parallel would require 10× the memory (over 2TB), which could exhaust node resources. Sequential runs are safer and simpler to monitor. Runtime is long but the 200-hour wall time limit accommodates this.

---

## 12. Known Issues and Fixes

### Issue 1: `melt()` conflict between `data.table` and `reshape2`
**Error:**
```
Error in melt.default: The melt generic in data.table has been passed a data.frame
and will attempt to redirect to the relevant reshape2 method, however this
redirect is now deprecated and will be removed in a future version...
```
**Cause:** Sleuth uses `melt()` internally. In `data.table` versions > 1.14.8, the redirect from `data.table::melt()` to `reshape2::melt()` was made a hard error instead of a warning.

**Fix — two things required, both are necessary:**
1. Load `library(reshape2)` **before** `library(sleuth)` in the R script. This ensures `reshape2::melt` is registered in the namespace first.
2. Downgrade `data.table` to 1.14.8:
```bash
conda install -c conda-forge "r-data.table=1.14.8" -y
```
Doing only one of these is not sufficient. Both are needed.

---

### Issue 2: Job silently killed during bootstrap summarization
**Symptom:** Job disappears from `qstat` with no error message. Console output is truncated mid-analysis, usually at the line:
```
summarizing bootstraps
```
**Cause:** Memory limit exceeded. The node's OOM (out of memory) killer terminates the process without warning when it exceeds `mem_free`. At 100 bootstraps across ~83,000 transcripts, Sleuth holds very large matrices in memory. On the standard `UI` queue with 128GB allocation, this was consistently hit.

**Fix:** Move to the `UI-HM` queue with `mem_free=256G`. The high-memory nodes (1–2TB available) have headroom for the full 100 bootstraps.

> **Interim workaround used:** `max_bootstrap = 30` was used temporarily to complete early comparisons while debugging. This is statistically valid (Pimentel et al. 2017 note 30 is sufficient), but 100 is preferred and is now restored.

---

### Issue 3: Bioconductor packages unavailable during install
**Error:**
```
Warning: unable to access index for repository https://bioconductor.org/packages/...
```
**Cause:** Argon compute nodes do not have outbound internet access to CRAN or Bioconductor repositories.

**Fix:** Install all R packages via conda from the `bioconda` and `conda-forge` channels before submitting jobs. Never use `install.packages()` or `BiocManager::install()` from within an Argon job — it will fail silently or hang.

---

### Issue 4: biomaRt installation hanging indefinitely
**Symptom:** `BiocManager::install("biomaRt")` starts but never completes on the cluster.

**Cause:** Two problems: (1) internet blocked on compute nodes, and (2) biomaRt is not needed — the `cdna_t2g_map.tsv` file already contains all the gene annotation required. biomaRt was in the professor's original script for a different project context.

**Fix:** Removed `library(biomaRt)` entirely. Not needed and not installable on Argon.

---

### Issue 5: Jobs completing in under 1 second with no output
**Symptom:** Job runs but exits almost immediately. No output files are created. No error.

**Cause:** `conda activate sleuth` inside the SGE job script does not reliably update `$PATH` in non-interactive shell sessions. The `Rscript` command resolves to the wrong (system) R installation, which doesn't have sleuth installed and exits silently.

**Fix:** Use the **full absolute path** to Rscript from the conda environment:
```bash
/old_Users/chakrabortty/miniforge3/envs/sleuth/bin/Rscript --vanilla sleuth_run.R $i Sleuth_4
```
Do not rely on `which Rscript` or the PATH being set correctly in SGE jobs.

---

### Issue 6: `Error: object 'sleuth_run.R' not found`
**Symptom:** R exits immediately with:
```
Error: object 'sleuth_run.R' not found
```
**Cause:** The text `sleuth_run.R` was accidentally pasted as the first line of the file itself when creating it, so R tried to evaluate the filename as a variable expression.

**Diagnosis:** Check the raw first bytes of the file:
```bash
od -c sleuth_run.R | head -3
```
If the output starts with `s   l   e   u   t   h   _   r   u   n   .   R`, the filename is embedded as line 1.

**Fix:** Open the file in nano (`nano sleuth_run.R`), use Ctrl+Home to go to the very beginning, and delete the erroneous first line. Save with Ctrl+O, exit with Ctrl+X. Verify the file starts with `args <- commandArgs(...)`.

---

### Issue 7: Memory accumulation across multiple comparisons
**Symptom:** Early comparisons complete fine but later ones fail with memory errors, even with generous `mem_free`.

**Cause:** R does not aggressively return memory to the OS after `rm()`. When multiple comparisons run in the same R session (e.g., looped inside one R script), each comparison's objects accumulate in memory even after being removed.

**Fix:** Run each comparison as a **separate Rscript process** via the bash loop in `sleuth_master.job`. Each process starts with a clean memory slate. The `rm(); gc()` calls within each comparison provide additional cleanup within a single run.

---

### Issue 8: Cannot allocate vector / out of memory on login node
**Symptom:** Running the script interactively on the login node fails with:
```
Error: cannot allocate vector of size 125.5 MB
```
**Cause:** Login nodes on Argon have severely limited memory (a few GB at most). They are not meant for computation — only for file management, script editing, and job submission.

**Fix:** Always run computational work via `qsub`. Never run Sleuth on the login node. This error is expected and does not indicate a problem with the script itself.

---

## 13. Development Journey

This section documents how the pipeline evolved — what was tried, what broke, and why decisions were made the way they were. It's written partly as a record and partly as a guide for anyone who inherits this project and runs into the same problems.

---

### Starting point: the professor's script

The starting point was Jared McLendon's original Sleuth script, written for a different project. It was designed to run two binary comparisons (resistant vs sensitive, for OV90 and OV8) and was structured as a single, self-contained R script that you'd run once.

The structure was:

1. **All-samples block** — prep all 48 samples into one Sleuth object, generate count matrices, PCA, heatmap
2. **Comparison 1** — OV90 resistant vs sensitive (untreated)
3. **Comparison 2** — OV8 resistant vs sensitive (untreated)

This worked for the original two comparisons and results were saved in `Sleuth_2/`. The professor's analysis logic — specifically the combination of LRT + Wald test, the `sleuth_to_matrix` expression export, the `left_join` to merge results with expression data, and the column cleanup with `select(-starts_with(...))` — was kept intact throughout all modifications. That part was not touched.

---

### The new comparisons: 8 more to add

The original two comparisons covered the baseline resistance signature (untreated resistant vs sensitive). The next scientific question was the **time response**: how does gene expression change at 6h and 48h after treatment, separately for sensitive and resistant cells, and separately for OV8 and OV90? That's 2 cell lines × 2 conditions × 2 time points = 8 more comparisons.

The first instinct was to copy the script and manually change the filter, beta, and formula variables for each comparison. That would have produced 8 separate R scripts, each nearly identical. Instead, the comparisons were encoded as a list inside one parameterized script, and a single index argument (`comp_index`) was passed from the command line to select which one to run:

```r
args       <- commandArgs(trailingOnly = TRUE)
comp_index <- as.integer(args[1])
```

This meant one script handles all 10 comparisons. The bash job loops through `1 2 3 4 5 6 7 8 9 10` and calls the same script with a different index each time.

---

### First crash: the melt() error

The very first run of the extended script failed with a cryptic error about `melt()`. This took some time to track down. Sleuth uses `melt()` internally, and a newer version of `data.table` had changed how it handles `melt()` calls — what used to be a warning became a hard error. The fix required two things: loading `reshape2` before `sleuth` to register `reshape2::melt` in the namespace first, and downgrading `data.table` to version 1.14.8 via conda. Neither fix alone was sufficient.

---

### Second crash: jobs killed silently

After the melt fix, jobs started running. The early comparisons (1 and 2) completed successfully and produced output in `Sleuth_2/`. But when the extended script tried to run the time-response comparisons (comparisons 3–10), the job would disappear from `qstat` without any error message. The console output file was truncated right at the line `summarizing bootstraps`.

That line is the most memory-intensive point in Sleuth's workflow — it's where it processes all bootstrap resamples across all transcripts simultaneously. The node ran out of memory and the OS killed the process. No error, no warning — just gone.

The immediate workaround was to reduce `max_bootstrap` from 100 to 30. This is statistically defensible (the Sleuth methods paper uses 30 as the standard), and it got the jobs to complete. But 30 is not as good as 100, so the plan was always to find more memory and restore it.

---

### Discovering the high-memory queue

The Argon cluster has a standard `UI` queue and a `UI-HM` (high-memory) queue. The HM nodes have 1–2TB of RAM. To check if access was available:

```bash
qconf -su atonu-chakrabortty   # check account details
qstat -g c                     # check queue availability
qhost -q                       # check which nodes have which queues
```

The output confirmed access to `UI-HM`. Switching to `UI-HM` with `mem_free=256G` and `smp 4` (4 cores) gave the job enough room to run at full 100 bootstraps without being killed.

---

### The filename-in-the-file bug

One run failed with a particularly confusing error:
```
Error: object 'sleuth_run.R' not found
```

This looked like a variable name error, but the script didn't have any variable called `sleuth_run.R`. The actual cause was that when the file was created, the filename itself (`sleuth_run.R`) was accidentally pasted as the literal first line of the file. R read it, tried to evaluate it as an expression, couldn't find an object with that name, and crashed.

The way to diagnose this is:
```bash
od -c sleuth_run.R | head -3
```
If the first characters are `s l e u t h _ r u n . R`, the filename is embedded. Fixed by opening in nano and deleting the first line.

---

### The PATH problem in SGE jobs

Early in the project, some jobs would appear to start and then complete almost instantly — in under a second — with no output files created. No errors were logged.

The cause was that `conda activate sleuth` inside the bash job script doesn't reliably update `$PATH` in non-interactive SGE shells. The shell found `Rscript` from the system R installation, which didn't have sleuth installed. R loaded, found no packages, and exited cleanly.

The fix was to stop relying on `conda activate` for the PATH and instead use the full absolute path to Rscript:
```bash
/old_Users/chakrabortty/miniforge3/envs/sleuth/bin/Rscript
```
This is the recommended pattern for all conda-based tools in SGE job scripts on Argon.

---

### Memory accumulation: why each comparison is a fresh process

An early design had all 10 comparisons running inside a single R session, in a loop. This seemed efficient. The problem: after each comparison, even with `rm(soG); gc()`, R holds onto a lot of memory that it doesn't return to the OS. By comparison 5 or 6, the accumulated unreleased memory pushed the job over its limit.

The solution was to make the bash loop in `sleuth_master.job` call `Rscript` once per comparison, not loop inside R. Each call is an independent process with a clean memory slate:

```bash
for i in 1 2 3 4 5 6 7 8 9 10; do
  Rscript --vanilla sleuth_run.R $i Sleuth_4 > sleuth_comp${i}_console.txt 2>&1
done
```

Each comparison's output goes to its own console file, which makes debugging easier — if comparison 7 fails, you check `sleuth_comp7_console.txt` without having to dig through one giant log.

---

### Final configuration

After all iterations, the final working configuration is:

| Parameter | Value | Why |
|-----------|-------|-----|
| Queue | `UI-HM` | Needs high-memory node to avoid OOM kills |
| Cores | `smp 4` | Enables parallel bootstrap processing in R |
| Memory | `256G` | Sufficient for 100 bootstraps across ~83K transcripts |
| `num_cores` | `4` | Matches smp 4; parallelizes bootstrap summarization |
| `max_bootstrap` | `100` | Restored to original; reliable variance estimation |
| Per-comparison processes | Yes | Prevents memory accumulation |
| Full Rscript path | Yes | Avoids PATH resolution issues in SGE |

---

## 14. References

1. **Pimentel H, Bray NL, Puente S, Melsted P, Pachter L.** Differential analysis of RNA-seq incorporating quantification uncertainty. *Nature Methods.* 2017;14:687–690. https://doi.org/10.1038/nmeth.4324

2. **Bray NL, Pimentel H, Melsted P, Pachter L.** Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology.* 2016;34:525–527. https://doi.org/10.1038/nbt.3519

3. **Love MI, Huber W, Anders S.** Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology.* 2014;15:550. https://doi.org/10.1186/s13059-014-0550-8 *(cited for comparison of DE methods)*

4. **Robinson MD, McCarthy DJ, Smyth GK.** edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics.* 2010;26(1):139–140. *(cited for comparison of DE methods)*

5. **Soneson C, Love MI, Robinson MD.** Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research.* 2015;4:1521. https://doi.org/10.12688/f1000research.7563.2 *(justification for transcript-level quantification)*
