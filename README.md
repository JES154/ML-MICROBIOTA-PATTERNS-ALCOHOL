Gut Microbiome Analysis Pipeline
QIIME2 Processing · PICRUSt2 Functional Profiling · Diversity Analysis · Machine Learning Classification
---

📋 Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Stage 1 — Data Acquisition](#stage-1--data-acquisition)
- [Stage 2 — QIIME2 Processing](#stage-2--qiime2-processing)
- [Stage 3 — R Statistical Analysis](#stage-3--r-statistical-analysis)
- [Stage 4 — Functional Profiling & Machine Learning](#stage-4--functional-profiling--machine-learning)
- [Outputs Summary](#outputs-summary)
- [Citation](#citation)

---

## Quick Start

If you only want to reproduce the statistical analyses and figures, 
download the processed input files included in this repository, place their respective path and run:

MICROBIAL_ANALYSIS_USER.R
FUNCTIONAL_ANALYSIS_OPT.R


## Overview

This repository contains the complete bioinformatics pipeline used to characteriZe gut microbiome composition and predict host phenotype from 16S rRNA amplicon sequencing data. Raw FASTQ files are sourced from the **NCBI BioProject** database and processed through four sequential stages:

```
NCBI BioProject
      │
      ▼
 [SRA Toolkit]  ──────────────────────────────  Raw FASTQ files
      │
      ▼
 qiime_metadata.sh  ──────────────────────────  metadata.tsv
      │
      ▼
 FULL_PROCESS_QIIME_RAW.sh  ──────────────────  feature-table.tsv
                                                 taxonomy.tsv
      │
      ▼
 PICRUST2.sh  ────────────────────────────────  picrust2_out_pipeline/
      │
      ▼
 R: MICROBIAL_ANALYSIS_USER.R  ──────────────────────  Ridge plots, PCoA, Alpha diversity
      │
      ▼
 R: FUNCTIONAL_ANALYSIS_OPT.R  ───────────────────────  Top 5 enzymes & pathways
                                                 SVM / Random Forest / XGBoost
                                                 PCoA (enzymes and pathways)
```


## Repository Structure

```
.
├── scripts/
│   ├── qiime_metadata.sh               # Stage 1: metadata formatting
│   ├── FULL_PROCESS_QIIME_RAW.sh       # Stage 2: full QIIME2 pipeline
│   └── PICRUST2.sh                     # Stage 2: PICRUSt2 inference
├── R/
│   ├── MICROBIAL_ANALYSIS_USER.R       # Stage 3: composition & diversity
│   └── FUNCTIONAL_ANALYSIS_OPT.R       # Stage 4: functional ML profiling
├── README.md
└── LICENSE
```

---

## Requirements

### Bioinformatics (Linux / HPC)

| Tool | Version | Installation |
|------|---------|-------------|
| QIIME2 | 2024.2 
| PICRUSt2 | 2.5.2 
| SRA Toolkit | 3.0 
| silva-138-99-nb-classifier 

### R (4.2+)

All packages are auto-installed at the top of each R script
---
Note: The input files required for the R analyses have already been generated and uploaded to this repository. Therefore, Stages 1 and 2 can be omitted. They are included only to illustrate the complete preprocessing workflow. Users may proceed directly by downloading the provided files and running the analyses in RStudio.

## Stage 1 — Data Acquisition

Download raw FASTQ files from NCBI using the SRA Toolkit:

```bash
# Download
prefetch <SRR_accession>

# Convert to FASTQ
fasterq-dump <SRR_accession> --split-files -O ./raw_reads/
```

> Raw sequencing datasets were obtained from the NCBI BioProject database:

- PRJEB46824
- PRJEB16755
- PRJNA867698
- PRJNA517050

Processed tables used for statistical analyses are included in this repository.

---

## Stage 2 — QIIME2 Processing

qiime_metadata.sh and FULL_PROCESS_QIIME_RAW.sh shell scripts must be executed **in order** with the QIIME2 conda environment active:

```bash
conda activate qiime2-amplicon-2024.2
```

### 1. `qiime_metadata.sh`
Validates and formats sample metadata into a QIIME2-compatible mapping file.

```bash
bash scripts/qiime_metadata.sh
```

**Output:** `metadata.tsv`

---

### 2. `FULL_PROCESS_QIIME_RAW.sh`
Runs the complete QIIME2 amplicon pipeline:

- FASTQ import via manifest file
- DADA2 denoising and ASV inference
- Taxonomic classification (SILVA 138 99% OTUs)
- Export of feature table and taxonomy

```bash
bash scripts/FULL_PROCESS_QIIME_RAW.sh
```

> ⚠️ Update the `CLASSIFIER_PATH` variable inside the script to point to your local SILVA classifier file.

**relevant outputs:**
```
qiime2_output/exported/
├── feature-table.tsv    # ASV count matrix (samples × ASVs)
└── taxonomy.tsv         # Taxonomic lineage per ASV
```

---

### 3. `PICRUST2.sh`
Runs the full PICRUSt2 prediction pipeline on DADA2 representative sequences.

```bash
conda activate PICRUST2
```
```bash
bash scripts/PICRUST2.sh
```

**Output directory:** `picrust2_out_pipeline/`
```
picrust2_out_pipeline/
├── EC_metagenome_out/
│   └── pred_metagenome_unstrat.tsv    # Enzyme (EC number) abundances
└── pathways_out/
    └── path_abun_unstrat.tsv          # MetaCyc pathway abundances
```

---

## Stage 3 — R Statistical Analysis

Open `MICROBIAL_ANALYSIS_USER.R` and update the **configuration block** at the top:

```r
# ---- Experiment paths ----
exp1 <- list(otu="path/HC/feature-table_controls.tsv",
             tax="path/HC/taxonomy_controls.tsv",
             meta="path/HC/META_CONTROLS.tsv", name="HC")

exp2 <- list(otu="path/HD/feature-table_HD.tsv",
             tax="path/HD/taxonomy_HD.tsv",
             meta="path/HD/MD_60ALCOHOL.tsv", name="HD")

exp3 <- list(otu="path/VHD/feature-table_VHD.tsv",
             tax="path/VHD/taxonomy_VHD.tsv",
             meta="path/VHD/MD_118ALCOHOL.tsv", name="VHD")

)
```

Then run the script **line by line** or source it entirely. Outputs are saved automatically to `results/`.

### What this script produces

#### 🔵 Phylum Composition
- **Ridge plots** — kernel density distributions of relative phylum abundance (%) per group

#### 🟢 Alpha Diversity
Indices calculated: **Observed Species**, **Chao1**, **Shannon**, **Simpson**

- Boxplots with jittered points per group
- Wilcoxon rank-sum test results (two groups)

#### 🟠 Beta Diversity
- Bray-Curtis dissimilarity matrix
- PCoA ordination (axes 1–2 and 1–3) with 95% confidence ellipses
- **PERMANOVA** (group centroid separation, 999 permutations)
- **PERMDISP** (within-group dispersion homogeneity, 999 permutations)

---

## Stage 4 — Functional Profiling & Machine Learning

Open `FUNCTIONAL_ANALYSIS_OPT.R` and update the `picrust=` field inside **each** of the three `exp` list objects to point to your PICRUSt2 output directories:

```r
exp1 <- list(
  otu     = "path/HC/feature-table_controls.tsv",
  tax     = "path/HC/taxonomy_controls.tsv",
  meta    = "path/HC/META_CONTROLS.tsv",
  name    = "HC",
  picrust = "path/HC/picrust2_out_HC"       # <-- update this
)
exp2 <- list(
  otu     = "path/HD/feature-table_HD.tsv",
  tax     = "path/HD/taxonomy_HD.tsv",
  meta    = "path/HD/MD_60ALCOHOL.tsv",
  name    = "HD",
  picrust = "path/HD/picrust2_out_HD"       # <-- update this
)
exp3 <- list(
  otu     = "path/VHD/feature-table_VHD.tsv",
  tax     = "path/VHD/taxonomy_VHD.tsv",
  meta    = "path/VHD/MD_118ALCOHOL.tsv",
  name    = "VHD",
  picrust = "path/VHD/picrust2_out_VHD"    # <-- update this
)
```

The `name` field in each list defines the group label used in all figures and statistical outputs. Update it to match your experimental groups if needed.

Once paths are configured, run the full script. It executes the following stages in order:

```
1. Load PICRUSt2 Data
        │
        ▼
2. Combine Experiments (Pathways & Enzymes separately)
        │
        ▼
3. Machine Learning (RF / SVM / XGBoost)
   ├── 5-fold cross-validation, 3 repeats
   ├── Model comparison (Global Score)
   └── Feature importance extraction
        │
        ▼
4. Beta Diversity (PCoA + PERMANOVA)
        │
        ▼
5. Violin Plots (Top features, importance bar, Dunn stats)
        │
        ▼
6. Statistical Tests (Kruskal-Wallis + Dunn + η²)
        │
        ▼
7. Export Results (PNG figures + Excel tables)
```

---

### What this script produces

#### Feature Matrix Preparation

PICRUSt2 enzyme and pathway tables are transposed into samples × features matrices, merged with group labels, and normalized to relative abundances per sample. Pathways and enzymes are processed independently throughout the pipeline.

#### ML Classification

Three classifiers are trained and compared using `caret` with repeated k-fold cross-validation:

| Algorithm | Key settings |
|---|---|
| **Random Forest** | `method = "rf"`, feature importance enabled |
| **SVM (Radial)** | `method = "svmRadial"`, center + scale preprocessing, `tuneLength = 5` |
| **XGBoost** | `method = "xgbTree"`, `tuneLength = 3` |

```
Method:  Repeated k-fold cross-validation
Folds:   5
Repeats: 3
Metrics: Accuracy, Kappa, F1 Score, Sensitivity, Specificity
```

The best-performing algorithm is selected by **Global Score** (mean of all five metrics). Feature importances from the winning model are used in all downstream plots and tables.

#### Top 5 Features Ranked by Importance

The top 5 **MetaCyc pathways** and top 5 **Enzyme Commission numbers** are extracted from the best model's importance scores and carried forward into violin plots and statistical tests.

#### Figures produced

- `model_performance.png` — side-by-side bar plots comparing Accuracy, Kappa, F1, Sensitivity, Specificity, and Global Score for all three algorithms, faceted by Pathways vs. Enzymes
- `p_path.png` — PCoA (Bray-Curtis) ordination for pathway abundances, annotated with PERMANOVA R² and p-value
- `p_enz.png` — PCoA (Bray-Curtis) ordination for enzyme abundances, annotated with PERMANOVA R² and p-value
- `Figure_Combined_Pathways_Enzymes.png` — two-panel figure (A: pathways, B: enzymes) showing violin plots of the top 5 features per level, with ML importance bars and Dunn post-hoc significance brackets

#### Statistical tables produced

- `Summary_Pathways_Statistics.xlsx` — Kruskal-Wallis p-values and η² effect sizes for top pathway features
- `Summary_Enzymes_Statistics.xlsx` — Kruskal-Wallis p-values and η² effect sizes for top enzyme features
- `Posthoc_Dunn_Pathways.xlsx` — pairwise Dunn test results (Bonferroni-adjusted) for pathways
- `Posthoc_Dunn_Enzymes.xlsx` — pairwise Dunn test results (Bonferroni-adjusted) for enzymes

---

> 📄 Associated manuscript submitted to *Pattern Analysis and Applications* (Springer Nature).
---

<p align="center">
  <sub>Built with QIIME2 · PICRUSt2 · R · vegan · phyloseq</sub>
</p>




