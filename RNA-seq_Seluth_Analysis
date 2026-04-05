    1 -# Sleuth Differential Expression Analysis
    2 -### Ovarian Cancer Cell Lines — OV8 & OV90 (Sensitive vs Resistant)
    1 +# Sleuth Differential Expression Analysis Pipeline
    2 +### Ovarian Cancer Cell Lines — OV8 & OV90 (Sensitive vs Resistant, Drug Response)
    3
    4  **Original Script Author:** Professor Jared M McLendon, University of Iowa
    5  **Modified by:** Atonu Chakrabortty, McLendon Lab, University of Iowa
    6  **Date:** April 2026
    7 -**Platform:** UIowa Argon HPC Cluster (SGE scheduler)
    7 +**Platform:** UIowa Argon HPC Cluster (SGE scheduler)
    8 +**Language:** R 4.3.1
    9
   10  ---
   11
   12  ## Table of Contents
   13  1. [Scientific Background](#1-scientific-background)
   14  2. [Experimental Design](#2-experimental-design)
   14 -3. [Analysis Overview](#3-analysis-overview)
   15 -4. [Dependencies](#4-dependencies)
   16 -5. [Environment Setup](#5-environment-setup)
   17 -6. [Directory Structure](#6-directory-structure)
   18 -7. [Input Files](#7-input-files)
   19 -8. [How to Run](#8-how-to-run)
   20 -9. [Comparisons Performed](#9-comparisons-performed)
   21 -10. [Output Files](#10-output-files)
   22 -11. [Parameter Justification](#11-parameter-justification)
   23 -12. [Known Issues and Fixes](#12-known-issues-and-fixes)
   24 -13. [References](#13-references)
   15 +3. [Why Kallisto and Sleuth](#3-why-kallisto-and-sleuth)
   16 +4. [Statistical Framework](#4-statistical-framework)
   17 +5. [Dependencies](#5-dependencies)
   18 +6. [Environment Setup on Argon HPC](#6-environment-setup-on-argon-hpc)
   19 +7. [Directory Structure](#7-directory-structure)
   20 +8. [Input Files](#8-input-files)
   21 +9. [Scripts in This Repository](#9-scripts-in-this-repository)
   22 +10. [How to Run](#10-how-to-run)
   23 +11. [Comparisons Performed](#11-comparisons-performed)
   24 +12. [Output Files](#12-output-files)
   25 +13. [Parameter Justification](#13-parameter-justification)
   26 +14. [Known Issues and Fixes](#14-known-issues-and-fixes)
   27 +15. [References](#15-references)
   28
   29  ---
   30
   31  ## 1. Scientific Background
   32
   30 -Ovarian cancer is one of the most lethal gynecological malignancies, largely due to the development of resistance to platinum-based chemotherapy. Understanding the transcriptomic differences between dr
      -ug-sensitive and drug-resistant ovarian cancer cells is critical for identifying molecular mechanisms of resistance and potential therapeutic targets.
   33 +Ovarian cancer is the fifth leading cause of cancer-related death in women and has one of the lowest five-year survival rates among gynecological malignancies. A major clinical challenge is the develop
      +ment of resistance to platinum-based chemotherapy (e.g., cisplatin, carboplatin), which affects the majority of patients after initial treatment.
   34
   32 -This pipeline performs RNA-seq differential expression analysis on two ovarian cancer cell lines:
   35 +Understanding the **transcriptomic basis of drug resistance** is essential for:
   36 +- Identifying genes and pathways that drive resistance
   37 +- Discovering potential therapeutic targets to overcome resistance
   38 +- Understanding how sensitive and resistant cells respond differently to drug treatment over time
   39
   40 +This pipeline performs RNA-seq differential expression analysis on two ovarian cancer cell line models:
   41  - **OV8** — ovarian cancer cell line
   42  - **OV90** — ovarian cancer cell line
   43
   37 -For each cell line, sensitive and resistant variants were profiled at multiple treatment time points to capture both **baseline resistance signatures** and **dynamic transcriptional responses to drug t
      -reatment**.
   44 +For each cell line, both **drug-sensitive** and **drug-resistant** variants were profiled at multiple time points following drug treatment, enabling two types of analysis:
   45 +1. **Baseline resistance signatures** — what is transcriptomically different between resistant and sensitive cells before any drug treatment
   46 +2. **Drug response dynamics** — how gene expression changes over time (6h, 48h) in response to drug treatment, separately in sensitive and resistant cells
   47
   39 -RNA-seq reads were previously quantified using **Kallisto** (Bray et al., 2016), which uses pseudoalignment to rapidly estimate transcript-level abundances with bootstrap resampling for uncertainty qua
      -ntification. This pipeline uses **Sleuth** (Pimentel et al., 2017) for downstream differential expression testing, which is specifically designed to work with Kallisto output and properly propagates qu
      -antification uncertainty into the statistical model.
   40 -
   48  ---
   49
   50  ## 2. Experimental Design
   51
   45 -### Cell Lines and Conditions
   52 +### Cell Lines, Conditions, and Time Points
   53
   47 -| Group Label | Cell Line | Condition | Treatment | Samples |
   48 -|-------------|-----------|-----------|-----------|---------|
   49 -| OV8_S_NT | OV8 | Sensitive | No Treatment (0h) | A1–A4 |
   50 -| OV8_S_6 | OV8 | Sensitive | 6 hours | B1–B4 |
   51 -| OV8_S_48 | OV8 | Sensitive | 48 hours | E1–E4 |
   52 -| OV8_R_NT | OV8 | Resistant | No Treatment (0h) | C1–C4 |
   53 -| OV8_R_6 | OV8 | Resistant | 6 hours | D1–D4 |
   54 -| OV8_R_48 | OV8 | Resistant | 48 hours | F1–F4 |
   55 -| OV90_S_NT | OV90 | Sensitive | No Treatment (0h) | G1–G4 |
   56 -| OV90_S_6 | OV90 | Sensitive | 6 hours | H1–H4 |
   57 -| OV90_S_48 | OV90 | Sensitive | 48 hours | I1–I4 |
   58 -| OV90_R_NT | OV90 | Resistant | No Treatment (0h) | J1–J4 |
   59 -| OV90_R_6 | OV90 | Resistant | 6 hours | K1–K4 |
   60 -| OV90_R_48 | OV90 | Resistant | 48 hours | L1–L4 |
   54 +| Group Label | Cell Line | Condition | Timepoint | Sample IDs | n |
   55 +|-------------|-----------|-----------|-----------|-----------|---|
   56 +| OV8_S_NT | OV8 | Sensitive | No Treatment | A1, A2, A3, A4 | 4 |
   57 +| OV8_S_6 | OV8 | Sensitive | 6 hours | B1, B2, B3, B4 | 4 |
   58 +| OV8_S_48 | OV8 | Sensitive | 48 hours | E1, E2, E3, E4 | 4 |
   59 +| OV8_R_NT | OV8 | Resistant | No Treatment | C1, C2, C3, C4 | 4 |
   60 +| OV8_R_6 | OV8 | Resistant | 6 hours | D1, D2, D3, D4 | 4 |
   61 +| OV8_R_48 | OV8 | Resistant | 48 hours | F1, F2, F3, F4 | 4 |
   62 +| OV90_S_NT | OV90 | Sensitive | No Treatment | G1, G2, G3, G4 | 4 |
   63 +| OV90_S_6 | OV90 | Sensitive | 6 hours | H1, H2, H3, H4 | 4 |
   64 +| OV90_S_48 | OV90 | Sensitive | 48 hours | I1, I2, I3, I4 | 4 |
   65 +| OV90_R_NT | OV90 | Resistant | No Treatment | J1, J2, J3, J4 | 4 |
   66 +| OV90_R_6 | OV90 | Resistant | 6 hours | K1, K2, K3, K4 | 4 |
   67 +| OV90_R_48 | OV90 | Resistant | 48 hours | L1, L2, L3, L4 | 4 |
   68
   69  - **n = 4 biological replicates per group**
   70  - **Total samples: 48**
   64 -- All samples with `exclude = 0` in metadata are included in analysis
   71
   66 -### Metadata Columns
   72 +### Metadata File Structure (`metadata.txt`)
   73
   68 -| Column | Description |
   69 -|--------|-------------|
   70 -| `sample` | Sample ID (A1–L4) |
   71 -| `group` | Group label |
   72 -| `ov8` | 1 = OV8 cell line |
   73 -| `ov90` | 1 = OV90 cell line |
   74 -| `resistant` | 1 = resistant condition |
   75 -| `sensitive` | 1 = sensitive condition |
   76 -| `zero` | 1 = no treatment (NT) |
   77 -| `six` | 1 = 6h treatment |
   78 -| `fortyeight` | 1 = 48h treatment |
   79 -| `exclude` | 1 = exclude sample from analysis |
   74 +| Column | Type | Description |
   75 +|--------|------|-------------|
   76 +| `sample` | string | Sample ID (A1–L4) |
   77 +| `group` | string | Full group label (e.g., OV8_S_NT) |
   78 +| `ov8` | 0/1 | 1 = OV8 cell line |
   79 +| `ov90` | 0/1 | 1 = OV90 cell line |
   80 +| `resistant` | 0/1 | 1 = drug resistant |
   81 +| `sensitive` | 0/1 | 1 = drug sensitive |
   82 +| `zero` | 0/1 | 1 = no treatment (NT) |
   83 +| `six` | 0/1 | 1 = 6h treatment |
   84 +| `fortyeight` | 0/1 | 1 = 48h treatment |
   85 +| `exclude` | 0/1 | 1 = exclude from analysis |
   86
   87  ---
   88
   83 -## 3. Analysis Overview
   89 +## 3. Why Kallisto and Sleuth
   90
   85 -### Why Sleuth?
   91 +### Why Kallisto for quantification?
   92 +Traditional RNA-seq alignment tools (STAR, HISAT2) map reads to a genome before counting. Kallisto (Bray et al., 2016) uses **pseudoalignment** against a transcriptome index, which is:
   93 +- ~50x faster than traditional aligners
   94 +- Equally accurate for expression quantification
   95 +- Produces **bootstrap resamples** that quantify technical (quantification) uncertainty — a key feature that Sleuth uses
   96
   87 -Kallisto quantifies transcript abundances probabilistically and generates bootstrap resamples that capture **technical (quantification) uncertainty**. Standard DE tools (DESeq2, edgeR) ignore this unce
      -rtainty. Sleuth explicitly models it using a two-component variance model:
   97 +### Why Sleuth for differential expression?
   98 +Most DE tools (DESeq2, edgeR, limma) model only **biological variance** between samples. They ignore the fact that transcript abundance estimates from RNA-seq are themselves uncertain — especially for
      +transcripts with low coverage or high sequence similarity.
   99
  100 +Sleuth (Pimentel et al., 2017) was designed specifically to work with Kallisto output and models:
  101 +
  102  ```
   90 -Total variance = Biological variance + Technical variance (from bootstraps)
  103 +Total observed variance = Biological variance + Technical variance
  104  ```
  105
   93 -This prevents false positives driven by technical noise rather than true biological differences.
  106 +Where technical variance is estimated directly from Kallisto's bootstrap resamples. This prevents false positives driven by quantification noise rather than true biological differences.
  107
   95 -### Statistical Framework
  108 +**Why not DESeq2/tximport?**
  109 +DESeq2 via tximport is also valid and widely used. However, because our data was quantified with Kallisto with bootstrap resampling (`max_bootstrap = 30`), Sleuth fully utilizes this information. DESeq
      +2 discards the bootstrap uncertainty. For transcript-level analysis especially, Sleuth is the more appropriate choice.
  110
   97 -For each comparison, Sleuth fits two linear models per gene/transcript:
  111 +---
  112
   99 -| Model | Formula | Meaning |
  100 -|-------|---------|---------|
  101 -| Full | `~ condition` | Expression depends on condition |
  113 +## 4. Statistical Framework
  114 +
  115 +### Model Design
  116 +For each comparison, two linear models are fit per gene/transcript using a response error measurement (REM) model:
  117 +
  118 +| Model | Formula | Interpretation |
  119 +|-------|---------|----------------|
  120 +| Full | `~ condition` | Expression depends on the condition of interest |
  121  | Reduced | `~ 1` | Expression does not depend on condition (intercept only) |
  122
  104 -Two complementary tests are run:
  123 +### Two Complementary Tests
  124
  106 -- **Likelihood Ratio Test (LRT):** Asks — *does adding the condition variable significantly improve the model?* Gives a p-value but no direction.
  107 -- **Wald Test (WT):** Estimates the **effect size (b coefficient)** of the condition. Gives direction and magnitude of change.
  125 +**Likelihood Ratio Test (LRT)**
  126 +- Compares full vs reduced model fit
  127 +- Asks: *does including the condition variable significantly improve the model?*
  128 +- Output: p-value (no direction/magnitude)
  129 +- Threshold used: `pval < 0.1`
  130
  109 -Results from both tests are merged into a single output table per comparison. This two-test approach is more conservative and reduces false positives.
  131 +**Wald Test (WT)**
  132 +- Estimates the beta coefficient (b) for the condition variable
  133 +- Asks: *what is the magnitude and direction of the effect?*
  134 +- Output: b (effect size), standard error, p-value, q-value (FDR-adjusted)
  135 +- Threshold for significance: `qval < 0.05`
  136
  111 -### Analysis Levels
  137 +**Why both tests?**
  138 +Using LRT + WT together is more conservative than either alone:
  139 +- LRT identifies genes where the condition matters at all
  140 +- WT provides effect size and direction
  141 +- Final results table merges both — researchers can apply their own thresholds
  142
  113 -Each comparison is run at **two levels**:
  143 +### Gene vs Transcript Level
  144 +Each comparison is run at two levels:
  145
  115 -1. **Gene level** (`gene_mode = TRUE`) — transcripts aggregated to gene level using `ens_gene` column from t2g map. Uses `scaled_reads_per_base` as the normalized expression measure.
  116 -2. **Transcript level** (`gene_mode = FALSE`) — individual transcript results. Uses `TPM` as the normalized expression measure.
  146 +| Level | `gene_mode` | Normalization | Use case |
  147 +|-------|-------------|---------------|---------|
  148 +| Gene | `TRUE` | scaled_reads_per_base | Overall gene expression, pathway analysis |
  149 +| Transcript | `FALSE` | TPM | Isoform-level analysis, splicing |
  150
  151  ---
  152
  120 -## 4. Dependencies
  153 +## 5. Dependencies
  154
  155  ### R Packages
  156
  157  | Package | Version | Source | Purpose |
  158  |---------|---------|--------|---------|
  159  | `sleuth` | 0.30.1 | Bioconda | Core DE analysis |
  127 -| `reshape2` | ≥1.4.4 | CRAN/conda-forge | Required to fix melt() conflict with data.table |
  128 -| `dplyr` | ≥1.0.0 | conda-forge | Data manipulation |
  129 -| `ggplot2` | ≥3.4.0 | conda-forge | Plotting |
  160 +| `reshape2` | ≥1.4.4 | conda-forge | Must load BEFORE sleuth — fixes melt() conflict |
  161 +| `dplyr` | ≥1.1.0 | conda-forge | Data filtering and manipulation |
  162 +| `ggplot2` | ≥3.4.0 | conda-forge | Plotting (PCA, volcano, QQ) |
  163  | `gridExtra` | ≥2.3 | conda-forge | Multi-panel plot layout |
  164 +| `rhdf5` | ≥2.46.1 | Bioconda | Read Kallisto .h5 bootstrap files |
  165
  132 -### System Requirements
  166 +### System
  167
  134 -- R version 4.3.1
  135 -- Miniforge/Conda
  136 -- SGE job scheduler (UIowa Argon HPC)
  137 -- Minimum 64GB RAM recommended (128GB for safety)
  168 +| Requirement | Value |
  169 +|-------------|-------|
  170 +| R version | 4.3.1 |
  171 +| Conda | Miniforge 24.11.3 |
  172 +| Scheduler | SGE (Sun Grid Engine) |
  173 +| Queue | UI (UIowa Argon general compute) |
  174 +| Memory | 128GB per job |
  175 +| Cores | 1 (single processor — see Known Issues) |
  176
  177  ---
  178
  141 -## 5. Environment Setup
  179 +## 6. Environment Setup on Argon HPC
  180
  143 -### One-time setup on Argon HPC
  181 +### One-time setup
  182
  183  ```bash
  146 -# Conda is available at:
  184 +# Conda is pre-installed at:
  185  /old_Users/<hawkid>/miniforge3/bin/conda
  186
  149 -# Create the sleuth environment
  187 +# Create the sleuth environment with all required packages
  188  conda create -n sleuth -c conda-forge -c bioconda \
  189    r-base=4.3.1 \
  190    r-essentials \
  153 -  r-devtools \
  154 -  r-ggplot2 \
  191 +  r-reshape2 \
  192    r-dplyr \
  193 +  r-ggplot2 \
  194    r-gridextra \
  157 -  r-reshape2 \
  195    bioconductor-rhdf5 \
  196    r-sleuth \
  197    -y
  198
  162 -# IMPORTANT: downgrade data.table to avoid melt() conflict with reshape2
  199 +# CRITICAL: downgrade data.table to avoid melt() conflict
  200  conda activate sleuth
  201  conda install -c conda-forge "r-data.table=1.14.8" -y
  202  ```
  203
  167 -### Verify installation
  204 +### Activate environment
  205
  206  ```bash
  207  source /old_Users/<hawkid>/miniforge3/bin/activate sleuth
  171 -R -e "library(sleuth); library(reshape2); library(dplyr); library(ggplot2); sessionInfo()"
  208  ```
  209
  210 +### Verify installation
  211 +
  212 +```bash
  213 +R -e "library(reshape2); library(sleuth); library(dplyr); library(ggplot2); sessionInfo()"
  214 +```
  215 +
  216 +All four packages should load without errors before submitting any jobs.
  217 +
  218  ---
  219
  176 -## 6. Directory Structure
  220 +## 7. Directory Structure
  221
  222  ```
  223  Atonu/
  180 -├── sleuth_atonu.R          # Main R analysis script
  181 -├── sleuth_atonu.job        # SGE job submission script
  182 -├── metadata.txt            # Sample metadata table
  183 -├── cdna_t2g_map.tsv        # Transcript-to-gene mapping file
  184 -├── Kallisto_1/             # Kallisto quantification output
  185 -│   ├── A1/
  186 -│   │   └── kallisto/
  187 -│   │       ├── abundance.h5
  188 -│   │       ├── abundance.tsv
  189 -│   │       └── run_info.json
  190 -│   ├── A2/
  191 -│   │   └── kallisto/
  224 +├── sleuth_atonu.R           # Batch 1: Comparisons 1-2 (R vs S baseline)
  225 +├── sleuth_atonu.job         # SGE job script for batch 1
  226 +├── sleuth_batch3.R          # Batch 3: Comparisons 3-10 (time response)
  227 +├── sleuth_batch3.job        # SGE job script for batch 3
  228 +├── metadata.txt             # Sample metadata (48 samples)
  229 +├── cdna_t2g_map.tsv         # Transcript-to-gene mapping (GRCh38)
  230 +│
  231 +├── Kallisto_1/              # Kallisto quantification output
  232 +│   ├── A1/kallisto/
  233 +│   │   ├── abundance.h5     # Bootstrap quantification (binary)
  234 +│   │   ├── abundance.tsv    # Plain text quantification
  235 +│   │   └── run_info.json
  236 +│   ├── A2/kallisto/
  237  │   ├── ...
  193 -│   └── L4/
  194 -│       └── kallisto/
  195 -└── Sleuth_2/               # Output directory (created by script)
  238 +│   └── L4/kallisto/
  239 +│
  240 +├── Sleuth_2/                # Output: Batch 1 (comparisons 1-2)
  241 +│   ├── s2c_all.tsv
  242 +│   ├── ALL-RESULTS_GENE_OV90_RvsS_NT_beta_resistant.tsv
  243 +│   ├── ALL-RESULTS_TRANS_OV90_RvsS_NT_beta_resistant.tsv
  244 +│   ├── ALL-RESULTS_GENE_OV8_RvsS_NT_beta_resistant.tsv
  245 +│   ├── ALL-RESULTS_TRANS_OV8_RvsS_NT_beta_resistant.tsv
  246 +│   ├── soG_*.so
  247 +│   ├── soT_*.so
  248 +│   └── PCA/Volcano/QQ PDFs
  249 +│
  250 +└── Sleuth_3/                # Output: Batch 3 (comparisons 3-10)
  251      ├── s2c_all.tsv
  197 -    ├── ALL-RESULTS_GENE_*.tsv
  198 -    ├── ALL-RESULTS_TRANS_*.tsv
  252 +    ├── ALL-RESULTS_GENE_OV90_S_6vsNT_beta_six.tsv
  253 +    ├── ALL-RESULTS_TRANS_OV90_S_6vsNT_beta_six.tsv
  254 +    ├── ... (one pair per comparison)
  255      ├── soG_*.so
  256      ├── soT_*.so
  201 -    ├── PCA_*.pdf
  202 -    ├── Volcano_*.pdf
  203 -    └── QQ_*.pdf
  257 +    └── PCA/Volcano/QQ PDFs
  258  ```
  259
  260  ---
  261
  208 -## 7. Input Files
  262 +## 8. Input Files
  263
  264  ### `metadata.txt`
  211 -Tab-separated file with one row per sample. Must contain columns:
  212 -`sample`, `group`, `ov8`, `ov90`, `resistant`, `sensitive`, `zero`, `six`, `fortyeight`, `exclude`
  265 +Tab-separated, one row per sample, 48 rows total. Must contain all columns described in Section 2. Samples with `exclude = 1` are automatically removed before any analysis.
  266
  267  ### `cdna_t2g_map.tsv`
  215 -Tab-separated transcript-to-gene mapping file with columns:
  216 -- `target_id` — Kallisto transcript ID
  217 -- `ens_gene` — Ensembl gene ID
  268 +Tab-separated transcript-to-gene mapping file. Two required columns:
  269 +- `target_id` — transcript ID matching Kallisto index
  270 +- `ens_gene` — Ensembl gene ID for aggregation
  271
  219 -This file must match the reference transcriptome used to build the Kallisto index. In this pipeline, the human reference genome **GRCh38** with Ensembl annotation was used.
  272 +This file must be built from the **same reference transcriptome** used to create the Kallisto index. This pipeline uses **GRCh38 (human)** with Ensembl annotation.
  273
  274  ### `Kallisto_1/`
  222 -Each sample subdirectory must contain a `kallisto/` folder with:
  223 -- `abundance.h5` — bootstrap quantification (binary, required for Sleuth)
  224 -- `abundance.tsv` — plain text quantification
  225 -- `run_info.json` — Kallisto run metadata
  275 +One subdirectory per sample, named exactly as in `metadata.txt`. Each must contain a `kallisto/` folder with `abundance.h5` (required for bootstrap access by Sleuth).
  276
  277  ---
  278
  229 -## 8. How to Run
  279 +## 9. Scripts in This Repository
  280
  231 -### Submit as a batch job (recommended)
  281 +### `sleuth_atonu.R` + `sleuth_atonu.job`
  282 +**Purpose:** Comparisons 1 and 2 — Resistant vs Sensitive at baseline (NT only)
  283 +**Output directory:** `Sleuth_2/`
  284 +**Comparisons:**
  285 +- OV90 Resistant vs Sensitive (NT)
  286 +- OV8 Resistant vs Sensitive (NT)
  287
  288 +### `sleuth_batch3.R` + `sleuth_batch3.job`
  289 +**Purpose:** Comparisons 3–10 — Drug response over time
  290 +**Output directory:** `Sleuth_3/`
  291 +**Comparisons:** All 8 time-response comparisons (see Section 11)
  292 +
  293 +### Why two separate scripts?
  294 +The full 10-comparison pipeline was split into two batches due to memory constraints on the Argon cluster. Each Sleuth object for transcript-level analysis consumes ~500MB RAM. Running all 10 compariso
      +ns sequentially in one job exceeded the 128GB memory limit. Splitting into two jobs of 2 and 8 comparisons respectively resolved this.
  295 +
  296 +---
  297 +
  298 +## 10. How to Run
  299 +
  300 +### Step 1: Activate environment
  301  ```bash
  234 -cd /Shared/lss_jmclendon/2_UserFolders/Atonu
  302  source /old_Users/<hawkid>/miniforge3/bin/activate sleuth
  236 -qsub sleuth_atonu.job
  303  ```
  304
  239 -### Monitor job status
  305 +### Step 2: Navigate to working directory
  306 +```bash
  307 +cd /Shared/lss_jmclendon/2_UserFolders/Atonu
  308 +```
  309
  310 +### Step 3: Submit jobs
  311  ```bash
  242 -# Check if running
  312 +# Batch 1 - Resistance comparisons (Sleuth_2/)
  313 +qsub sleuth_atonu.job
  314 +
  315 +# Batch 3 - Time response comparisons (Sleuth_3/)
  316 +qsub sleuth_batch3.job
  317 +```
  318 +
  319 +### Step 4: Monitor jobs
  320 +```bash
  321 +# Check job status
  322  qstat
  323
  324  # Check specific job
  325  qstat -j <job-ID>
  326
  248 -# View live log output
  249 -cat /Shared/lss_jmclendon/2_UserFolders/Atonu/Sleuth_log_out.txt
  250 -
  251 -# View R console output
  252 -cat /Shared/lss_jmclendon/2_UserFolders/Atonu/Sleuth_console_out_<timestamp>.txt
  327 +# Watch console output in real time
  328 +tail -f Sleuth_batch3_console_<timestamp>.txt
  329  ```
  330
  255 -### Job script (`sleuth_atonu.job`)
  256 -
  331 +### Step 5: Verify completion
  332  ```bash
  258 -#!/bin/bash
  259 -#$ -q UI                    # UIowa general compute queue
  260 -#$ -pe smp 1                # Single processor (multi-core causes Sleuth errors)
  261 -#$ -l mem_free=128G         # 128GB RAM
  262 -#$ -l h_rt=100:00:00        # 100 hour time limit
  263 -#$ -N sleuth_atonu          # Job name
  264 -#$ -M <hawkid>@uiowa.edu    # Email notifications
  265 -#$ -m beas                  # Email on Begin, End, Abort, Suspend
  266 -#$ -j y                     # Merge stdout and stderr
  267 -#$ -o /Shared/lss_jmclendon/2_UserFolders/Atonu/Sleuth_log_out.txt
  268 -#$ -cwd                     # Run from current working directory
  333 +# Check output files created
  334 +ls -lth Sleuth_2/
  335 +ls -lth Sleuth_3/
  336 +
  337 +# Check log for errors
  338 +cat Sleuth_batch3_log.txt
  339  ```
  340
  341  ---
  342
  273 -## 9. Comparisons Performed
  343 +## 11. Comparisons Performed
  344
  275 -### Category 1: Resistance Signature (Baseline)
  276 -Untreated samples only. Tests which genes are differentially expressed between resistant and sensitive cells **before any drug treatment**.
  345 +### Batch 1 — Resistance Signature at Baseline (`Sleuth_2/`)
  346
  278 -| # | Comparison | Model | Samples |
  279 -|---|-----------|-------|---------|
  280 -| 1 | OV90 Resistant vs Sensitive | `~resistant` | J1–J4 vs G1–G4 |
  281 -| 2 | OV8 Resistant vs Sensitive | `~resistant` | C1–C4 vs A1–A4 |
  347 +| # | Comparison | Model | Reference | Test |
  348 +|---|-----------|-------|-----------|------|
  349 +| 1 | OV90 Resistant vs Sensitive (NT) | `~resistant` | Sensitive (resistant=0) | Resistant (resistant=1) |
  350 +| 2 | OV8 Resistant vs Sensitive (NT) | `~resistant` | Sensitive (resistant=0) | Resistant (resistant=1) |
  351
  283 -### Category 2: Drug Response Over Time
  284 -Tests how gene expression changes in response to drug treatment within each group (sensitive/resistant) for each cell line. Reference is always untreated (NT).
  352 +**Biological question:** Which genes are differentially expressed between resistant and sensitive cells *before any drug treatment*? These represent the baseline molecular signature of resistance.
  353
  286 -| # | Comparison | Model | Samples |
  287 -|---|-----------|-------|---------|
  288 -| 3 | OV90 Sensitive: 6h vs NT | `~six` | H1–H4 vs G1–G4 |
  289 -| 4 | OV90 Sensitive: 48h vs NT | `~fortyeight` | I1–I4 vs G1–G4 |
  290 -| 5 | OV90 Resistant: 6h vs NT | `~six` | K1–K4 vs J1–J4 |
  291 -| 6 | OV90 Resistant: 48h vs NT | `~fortyeight` | L1–L4 vs J1–J4 |
  292 -| 7 | OV8 Sensitive: 6h vs NT | `~six` | B1–B4 vs A1–A4 |
  293 -| 8 | OV8 Sensitive: 48h vs NT | `~fortyeight` | E1–E4 vs A1–A4 |
  294 -| 9 | OV8 Resistant: 6h vs NT | `~six` | D1–D4 vs C1–C4 |
  295 -| 10 | OV8 Resistant: 48h vs NT | `~fortyeight` | F1–F4 vs C1–C4 |
  354 +### Batch 3 — Drug Response Over Time (`Sleuth_3/`)
  355
  356 +| # | Comparison | Model | Reference | Test |
  357 +|---|-----------|-------|-----------|------|
  358 +| 3 | OV90 Sensitive: 6h vs NT | `~six` | NT (six=0) | 6h (six=1) |
  359 +| 4 | OV90 Sensitive: 48h vs NT | `~fortyeight` | NT (fortyeight=0) | 48h (fortyeight=1) |
  360 +| 5 | OV90 Resistant: 6h vs NT | `~six` | NT (six=0) | 6h (six=1) |
  361 +| 6 | OV90 Resistant: 48h vs NT | `~fortyeight` | NT (fortyeight=0) | 48h (fortyeight=1) |
  362 +| 7 | OV8 Sensitive: 6h vs NT | `~six` | NT (six=0) | 6h (six=1) |
  363 +| 8 | OV8 Sensitive: 48h vs NT | `~fortyeight` | NT (fortyeight=0) | 48h (fortyeight=1) |
  364 +| 9 | OV8 Resistant: 6h vs NT | `~six` | NT (six=0) | 6h (six=1) |
  365 +| 10 | OV8 Resistant: 48h vs NT | `~fortyeight` | NT (fortyeight=0) | 48h (fortyeight=1) |
  366 +
  367 +**Biological question:** How does drug treatment change gene expression over time? Are the transcriptional responses to drug treatment different in sensitive vs resistant cells?
  368 +
  369  ---
  370
  299 -## 10. Output Files
  371 +## 12. Output Files
  372
  301 -All output files are saved in `Sleuth_2/`.
  373 +### Results Tables (`.tsv`)
  374
  303 -### Results Tables
  375 +For each comparison, two results tables are produced:
  376
  305 -| File | Level | Contents |
  306 -|------|-------|---------|
  307 -| `ALL-RESULTS_GENE_<comparison>_beta_<beta>.tsv` | Gene | Combined LRT + Wald test results + normalized counts per gene |
  308 -| `ALL-RESULTS_TRANS_<comparison>_beta_<beta>.tsv` | Transcript | Combined LRT + Wald test results + TPM per transcript |
  377 +| File | Level | Normalized measure |
  378 +|------|-------|-------------------|
  379 +| `ALL-RESULTS_GENE_<comparison>_beta_<beta>.tsv` | Gene | scaled_reads_per_base |
  380 +| `ALL-RESULTS_TRANS_<comparison>_beta_<beta>.tsv` | Transcript | TPM |
  381
  310 -**Key columns in results tables:**
  382 +### Key Columns in Results Tables
  383
  312 -| Column | Source | Description |
  313 -|--------|--------|-------------|
  314 -| `target_id` | Both | Gene or transcript ID |
  315 -| `b` | Wald Test | Effect size (log fold change equivalent) |
  384 +| Column | Source test | Description |
  385 +|--------|-------------|-------------|
  386 +| `target_id` | — | Gene or transcript ID |
  387 +| `b` | Wald Test | Effect size (analogous to log fold change) |
  388  | `se_b` | Wald Test | Standard error of b |
  389 +| `mean_obs` | — | Mean observed expression across samples |
  390  | `pval.x` | Wald Test | Wald test p-value |
  318 -| `qval.x` | Wald Test | Wald test FDR-adjusted q-value |
  391 +| `qval.x` | Wald Test | Wald test FDR q-value **(use for significance)** |
  392  | `pval.y` | LRT | Likelihood ratio test p-value |
  320 -| `qval.y` | LRT | LRT FDR-adjusted q-value |
  321 -| `mean_obs` | Both | Mean observed expression |
  322 -| Sample columns | Both | Normalized counts per sample |
  393 +| `qval.y` | LRT | LRT FDR q-value |
  394 +| Sample columns | — | Normalized expression per sample |
  395
  324 -### Saved Sleuth Objects
  396 +**Recommended filtering for significant genes:**
  397 +```r
  398 +sig_genes <- results[results$qval.x < 0.05 & !is.na(results$qval.x), ]
  399 +```
  400
  401 +### Saved Sleuth Objects (`.so`)
  402 +
  403  | File | Description |
  404  |------|-------------|
  328 -| `soG_<comparison>.so` | Gene-level Sleuth object — reload in R for interactive exploration |
  405 +| `soG_<comparison>.so` | Gene-level Sleuth object |
  406  | `soT_<comparison>.so` | Transcript-level Sleuth object |
  407
  331 -Reload with:
  408 +Reload in R for interactive analysis:
  409  ```r
  410 +library(sleuth)
  411  so <- sleuth_load("Sleuth_2/soG_OV90_RvsS_NT_beta_resistant.so")
  412 +sleuth_live(so)  # opens interactive Shiny app
  413  ```
  414
  336 -### Plots
  415 +### Plots (`.pdf`)
  416
  338 -| File | Description |
  339 -|------|-------------|
  340 -| `PCA_GENE_*.pdf` | PCA plot colored by group — check sample clustering |
  341 -| `PCA_TRANS_*.pdf` | Transcript-level PCA |
  342 -| `Volcano_GENE_*.pdf` | Volcano plot — effect size vs significance |
  343 -| `Volcano_TRANS_*.pdf` | Transcript-level volcano |
  344 -| `QQ_GENE_*.pdf` | QQ plot — check p-value distribution for model fit |
  345 -| `QQ_TRANS_*.pdf` | Transcript-level QQ |
  417 +| File | Description | What to look for |
  418 +|------|-------------|-----------------|
  419 +| `PCA_GENE_*.pdf` | PCA colored by group | Samples should cluster by group |
  420 +| `PCA_TRANS_*.pdf` | Transcript-level PCA | Same as above |
  421 +| `Volcano_GENE_*.pdf` | Effect size vs -log10(pval) | Symmetric spread, clear significant genes |
  422 +| `Volcano_TRANS_*.pdf` | Transcript-level volcano | — |
  423 +| `QQ_GENE_*.pdf` | Observed vs expected p-values | Inflation = model issues; deflation = few DE genes |
  424 +| `QQ_TRANS_*.pdf` | Transcript-level QQ | — |
  425
  426  ---
  427
  349 -## 11. Parameter Justification
  428 +## 13. Parameter Justification
  429
  430  ### `num_cores = 1`
  352 -Sleuth's parallel processing via `num_cores > 1` caused job abortion errors on the Argon HPC cluster. Single-core execution is stable and produces identical results.
  431 +**Why:** Multi-core (`num_cores > 1`) caused job abortion errors on Argon HPC due to memory conflicts during parallel bootstrap processing. Single-core execution is stable and produces identical result
      +s — parallelization only affects speed, not output.
  432
  433  ### `max_bootstrap = 30`
  355 -The original script used 100 bootstraps. Reduced to 30 due to memory constraints on the cluster. This is **statistically justified**:
  434 +**Why:** Reduced from 100 (original script) to 30 due to memory constraints. This is fully supported by the Sleuth paper:
  435
  357 -> "We find that B = 30 bootstraps is sufficient for reliable variance estimation in most RNA-seq experiments."
  358 -> — Pimentel et al., 2017, *Nature Methods*
  436 +> *"We find that B = 30 bootstraps is sufficient for reliable variance estimation"*
  437 +> — Pimentel et al., 2017, Nature Methods
  438
  360 -30 bootstraps is the value used in the Sleuth methods paper itself and is considered standard practice.
  439 +30 bootstraps is the value used in the Sleuth methods paper's own simulations and is standard practice for published analyses.
  440
  441  ### `gene_mode = TRUE` with `aggregation_column = 'ens_gene'`
  363 -Aggregates transcript-level estimates to gene level using the mean of transcript-level quantities, weighted by abundance. This reduces multiple testing burden and is more interpretable for pathway anal
      -ysis downstream.
  442 +**Why:** Aggregates transcript-level estimates to gene level using the `ens_gene` column from the t2g map. This reduces the multiple testing burden (15,000 genes vs 83,000 transcripts) and is more inte
      +rpretable for pathway analysis with tools like clusterProfiler or fgsea.
  443
  365 -### `pval < 0.1` threshold for LRT filtering
  366 -Standard threshold used in Sleuth analyses. More lenient than 0.05 at the LRT stage because the Wald test q-value provides the final significance filter for downstream interpretation.
  444 +### `pval < 0.1` for LRT filtering
  445 +**Why:** Standard threshold in Sleuth analyses. The LRT p-value is used as an initial filter; the Wald test q-value (`< 0.05`) provides the final significance cutoff in the merged results table.
  446
  447 +### Two separate job submissions
  448 +**Why:** Each transcript-level Sleuth object consumes ~500MB RAM. Running all 10 comparisons sequentially in one job caused memory exhaustion and job termination during bootstrap summarization (~83,000
      + transcripts × 30 bootstraps). Splitting into two jobs resolved this.
  449 +
  450  ---
  451
  370 -## 12. Known Issues and Fixes
  452 +## 14. Known Issues and Fixes
  453
  372 -### Issue 1: `melt()` conflict between `data.table` and `reshape2`
  373 -**Error:**
  454 +### Issue 1: `melt()` conflict — data.table vs reshape2
  455 +**Error message:**
  456  ```
  375 -Error in melt.default: The melt generic in data.table has been passed a data.frame
  376 -and will attempt to redirect to the relevant reshape2 method...
  457 +Error in melt.default: The melt generic in data.table has been passed a data.frame...
  458 +this redirection is now deprecated.
  459  ```
  378 -**Cause:** `data.table` (v > 1.14.8) intercepts Sleuth's internal `melt()` calls and treats the redirect as a hard error.
  460 +**Root cause:** data.table (v > 1.14.8) intercepts Sleuth's internal `melt()` calls on data.frames and treats the redirect to reshape2 as a hard error instead of a warning.
  461
  380 -**Fix:**
  381 -1. Load `library(reshape2)` **before** `library(sleuth)` in the R script
  382 -2. Downgrade `data.table` to 1.14.8:
  462 +**Fix applied:**
  463 +1. Load `library(reshape2)` **before** `library(sleuth)` in every R script
  464 +2. Downgrade data.table:
  465  ```bash
  466  conda install -c conda-forge "r-data.table=1.14.8" -y
  467  ```
  468
  469  ### Issue 2: Job killed during bootstrap summarization
  388 -**Symptom:** Job disappears from `qstat` with output truncated at `summarizing bootstraps`
  470 +**Symptom:** Job disappears from `qstat` with console output truncated at:
  471 +```
  472 +summarizing bootstraps
  473 +```
  474 +**Root cause:** Memory limit exceeded. Bootstrap summarization for ~83,000 transcripts with 30 bootstraps × 8 samples requires >64GB RAM.
  475
  390 -**Cause:** Memory limit exceeded during bootstrap processing of ~83,000 transcripts
  476 +**Fix applied:**
  477 +- Increased `mem_free` from 64G to 128G in job script
  478 +- Reduced `max_bootstrap` from 100 to 30
  479 +- Split pipeline into two separate job submissions
  480
  392 -**Fix:**
  393 -- Reduce `max_bootstrap` from 100 to 30
  394 -- Increase `mem_free` to 128G in job script
  395 -
  396 -### Issue 3: Bioconductor packages unavailable
  397 -**Error:**
  481 +### Issue 3: Bioconductor/CRAN not accessible
  482 +**Error message:**
  483  ```
  484  Warning: unable to access index for repository https://bioconductor.org/...
  485 +cannot open URL
  486  ```
  401 -**Cause:** Argon compute nodes block outbound internet access to CRAN/Bioconductor
  487 +**Root cause:** Argon compute nodes block outbound internet access.
  488
  403 -**Fix:** Install all R packages via conda instead of `install.packages()` or `BiocManager::install()`
  489 +**Fix applied:** Install all R packages via conda channels (`bioconda`, `conda-forge`) instead of `install.packages()` or `BiocManager::install()`.
  490
  491 +### Issue 4: Multi-core errors
  492 +**Symptom:** Job aborts when `num_cores > 1` in `sleuth_prep()`
  493 +
  494 +**Root cause:** Sleuth's parallel processing via `parallel` package conflicts with SGE's memory management on shared nodes.
  495 +
  496 +**Fix applied:** Set `num_cores = 1` in all `sleuth_prep()` calls.
  497 +
  498  ---
  499
  407 -## 13. References
  500 +## 15. References
  501
  409 -1. **Pimentel H, Bray NL, Puente S, Melsted P, Pachter L.** Differential analysis of RNA-seq incorporating quantification uncertainty. *Nature Methods.* 2017;14:687–690. https://doi.org/10.1038/nmeth.4
      -324
  502 +1. **Pimentel H, Bray NL, Puente S, Melsted P, Pachter L.** Differential analysis of RNA-seq incorporating quantification uncertainty. *Nature Methods.* 2017;14:687–690.
  503 +   https://doi.org/10.1038/nmeth.4324
  504 +   *(Primary Sleuth reference — also justifies 30 bootstraps)*
  505
  411 -2. **Bray NL, Pimentel H, Melsted P, Pachter L.** Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology.* 2016;34:525–527. https://doi.org/10.1038/nbt.3519
  506 +2. **Bray NL, Pimentel H, Melsted P, Pachter L.** Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology.* 2016;34:525–527.
  507 +   https://doi.org/10.1038/nbt.3519
  508 +   *(Kallisto pseudoalignment)*
  509
  413 -3. **Love MI, Huber W, Anders S.** Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology.* 2014;15:550. https://doi.org/10.1186/s13059-014-0550-8 *(cited for
      -comparison of DE methods)*
  510 +3. **Soneson C, Love MI, Robinson MD.** Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research.* 2015;4:1521.
  511 +   https://doi.org/10.12688/f1000research.7563.2
  512 +   *(Justification for transcript-level quantification and tximport approach)*
  513
  415 -4. **Robinson MD, McCarthy DJ, Smyth GK.** edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics.* 2010;26(1):139–140. *(cited for comparis
      -on of DE methods)*
  514 +4. **Love MI, Huber W, Anders S.** Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology.* 2014;15:550.
  515 +   https://doi.org/10.1186/s13059-014-0550-8
  516 +   *(Alternative DE method — cited for comparison)*
  517
  417 -5. **Soneson C, Love MI, Robinson MD.** Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research.* 2015;4:1521. https://doi.org/10.12688/f1000research
      -.7563.2 *(justification for transcript-level quantification)*
  518 +5. **Robinson MD, McCarthy DJ, Smyth GK.** edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics.* 2010;26(1):139–140.
  519 +   https://doi.org/10.1093/bioinformatics/btp616
  520 +   *(Alternative DE method — cited for comparison)*
  521
  522 +6. **Liao Y, Smyth GK, Shi W.** featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics.* 2014;30(7):923–930.
  523 +   *(Alternative counting approach — cited for comparison)*
  524 +
  525  ---
  526
  527  ## Contact
  528
  423 -For questions about this pipeline, contact:
  424 -**Atonu Chakrabortty** — McLendon Lab, University of Iowa
  425 -For questions about the original script:
  426 -**Jared McLendon** — McLendon Lab, University of Iowa
  529 +**Atonu Chakrabortty**
  530 +McLendon Lab, Department of Pharmacology
  531 +University of Iowa
  532 +
  533 +**Jared McLendon** (original script)
  534 +McLendon Lab, University of Iowa
  535 +jared-mclendon@uiowa.edu
