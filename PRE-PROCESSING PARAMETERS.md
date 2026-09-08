## BioProject-Specific Quality Filtering

Because sequencing quality profiles differed across studies, 
trimming and truncation parameters were selected separately for each BioProject. 
Parameters were determined by inspecting the corresponding QIIME 2 `demux.qzv` quality visualizations.

Primer sequences were identified from the original publication associated with each BioProject and removed before denoising if still present. 
The table below summarizes the BioProject-specific preprocessing parameters used to remove adapters, primers, low-quality read regions, 
and reads that did not meet the selected quality criteria before ASV inference.

All quality-control decisions were based on sequencing-quality diagnostics and were made independently of alcohol-intake labels and 
machine-learning model performance. Samples with poor read quality or substantial contamination were excluded according to predefined 
technical quality-control criteria.

> **Note:** Quality filtering reduces technical noise and sequencing artifacts but does not eliminate zero values from the
resulting feature table. Zero counts are expected in 16S rRNA amplicon datasets because microbial feature tables are naturally sparse.
Subsequent filtering removed singleton ASVs and features with fewer than 10 total counts across the dataset.

## BioProject-Specific Quality Filtering

| BioProject | Group | Sequencing layout | Forward primer | Reverse primer | Trim-left F | Trunc-len F | Trim-left R | Trunc-len R | Max EE F | Max EE R | Samples removed | Exclusion reason |
|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| PRJNA517050 | VHD | Paired-end | `ACTCCTACGGGAGGCAGCAG` | `GGACTACHVGGGTWTCTAAT` | 0 | 251 | 0 | 240 | 2 | 2 | 0 | — |
| PRJEB46824 | HC | Paired-end | `CCTAYGGGRBGCASCAG` | `GGACTACNNGGGTATCTAAT` | 0 | 0 | 0 | 0 | 2 | 2 | 0 | — |
| PRJNA867698 | HD | Single-end | `CCTGTTCGATACCCGCACTTTCGAGCTTCAG` | NA | 0 | NA | NA | NA | NA| NA | 2 | Low read quality and substantial contamination |
| PRJEB16755 | HC | Paired-end | `CCTAYGGGRBGCASCAG` | `GGACTACNNGGGTATCTAAT` | 0 | 0 | 0 | 0 | 2 | 2 | 0 | — |
