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

The three shell scripts must be executed **in order** with the QIIME2 conda environment active:

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
exp1 <- list(
  otu  = "<path>/qiime2_output/exported/feature-table.tsv",
  tax  = "<path>/qiime2_output/exported/taxonomy.tsv",
  meta = "<path>/metadata.tsv"
)

GROUP_COL    <- "GROUP"                          # exact column name in metadata
GROUP_LEVELS <- c("HC", "HD", "VHD")            # control group first
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

Open `FUNCTIONAL_ANALYSIS_OPT.R` and set the PICRUSt2 output path:

```r
picrust2_path <- "<path>/picrust2_out_pipeline"
```

This script uses the same metadata grouping as Stage 3 (`GROUP_COL`, `GROUP_LEVELS`).

### What this script produces

#### Feature Matrix Preparation
PICRUSt2 enzyme and pathway tables are transposed to samples × features matrices, merged with group labels, filtered for zero-variance features, and normalised to relative abundances per sample.

#### ML Classification
Three classifiers are trained to discriminate `Healthy control` vs `60_AUD`:

| Classifier | Details |
|-----------|---------|
| **SVM** | Radial basis function kernel, tuned via cross-validation |
| **Random Forest** | 500 trees, out-of-bag error estimation |
| **XGBoost** | Gradient boosted trees with early stopping |

Training strategy: stratified 5-fold cross-validation, 10 repetitions.  
Performance metrics: Accuracy, AUC-ROC, Sensitivity, Specificity.

#### 🏆 Top 5 Features Ranked by Importance
Feature importance scores are extracted from the best-performing algorithm (by AUC-ROC). The **top 5 enzymes (EC numbers)** and **top 5 MetaCyc pathways** are reported, ranked by mean importance across CV folds.

#### Figures produced
- PCoA plots of enzyme and pathway abundance matrices
- Violin plots of top 5 discriminative enzymes and pathways by group
- Algorithm performance comparison (SVM vs RF vs XGBoost)

---

> Associated manuscript SUBMITTED to *Pattern Analysis and Applications* (Springer Nature).

---

<p align="center">
  <sub>Built with QIIME2 · PICRUSt2 · R · vegan · phyloseq</sub>
</p>
