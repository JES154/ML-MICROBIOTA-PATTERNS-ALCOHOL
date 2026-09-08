#!/bin/bash
set -euo pipefail

# =========================================================
# QIIME2 Pipeline - Auto-detect SE vs PE with manifest
# Import > Cutadapt > DADA2 > Taxonomy > Phylogeny > Diversity
# Author: Juan Esparza
# =========================================================

# # -----------------------------
# # CREATE OUTPUTS DESTINATION FOLDER
# # -----------------------------
CURRENT_DIR="$PWD"
OUT="$CURRENT_DIR/qiime2_output"
# mkdir -p "$OUT"
# echo "[INFO] Output directory: $OUT"

# # -----------------------------
# # 1. GENERATE QIIME2 MANIFEST
# # -----------------------------
if ls *_2.fastq 1>/dev/null 2>&1; then
    SEQ_TYPE="paired"
else
    SEQ_TYPE="single"
fi

echo "[INFO] Detected sequence type: $SEQ_TYPE"

if [[ "$SEQ_TYPE" == "paired" ]]; then

    MANIFEST="$OUT/manifest_paired.tsv"
    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST"

    for R1 in *_1.fastq; do
        SAMPLE=$(basename "$R1" _1.fastq)
        R2="${R1/_1.fastq/_2.fastq}"
        if [[ -f "$R2" ]]; then
            echo -e "$SAMPLE\t$CURRENT_DIR/$R1\t$CURRENT_DIR/$R2" >> "$MANIFEST"
        else
            echo "[WARN] R2 not found for $SAMPLE"
        fi
    done

else

    MANIFEST="$OUT/manifest_single.tsv"
    echo -e "sample-id\tabsolute-filepath" > "$MANIFEST"

    for R1 in *.fastq; do
        SAMPLE=$(basename "$R1" .fastq)
        echo -e "$SAMPLE\t$CURRENT_DIR/$R1" >> "$MANIFEST"
    done

fi

echo "[INFO] Manifest generated: $MANIFEST"
head "$MANIFEST"

# -----------------------------
# 2.  QIIME2 MANIFEST IMPORTATION
# -----------------------------
if [[ "$SEQ_TYPE" == "paired" ]]; then
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path "$MANIFEST" \
      --output-path "$OUT/paired-end-demux.qza" \
      --input-format PairedEndFastqManifestPhred33V2

    DEMUX="$OUT/paired-end-demux.qza"
else
    qiime tools import \
      --type 'SampleData[SequencesWithQuality]' \
      --input-path "$MANIFEST" \
      --output-path "$OUT/single-end-demux.qza" \
      --input-format SingleEndFastqManifestPhred33V2

    DEMUX="$OUT/single-end-demux.qza"
fi

# -----------------------------
# 3. SAMPLES QUALITY VISUALIZATION
# -----------------------------
qiime demux summarize \
  --i-data "$DEMUX" \
  --o-visualization "$OUT/demux.qzv"

# IDENTIFIQUEMOS SI HAY PRIMERS
# awk 'NR%4==2' SRR20993010.fastq | head -n 10  #SELECCIONAMOS LA MUESTRA QUE QUEREMOS OBSERVAR E IDENTIFICAMOS PATRONES AL INICIO Y FINAL
# awk 'NR%4==2' SRR8488424_2.fastq | head -n 10  #SELECCIONAMOS LA MUESTRA QUE QUEREMOS OBSERVAR E IDENTIFICAMOS PATRONES AL INICIO Y FINAL

# UNA VEZ OBTENIDOS (PUEDE SER CON IA, AGREGARLOS A CUTADAPT)

# -----------------------------
# 4. CUTADAPT
# -----------------------------
if [[ "$SEQ_TYPE" == "paired" ]]; then
    qiime cutadapt trim-paired \
      --i-demultiplexed-sequences "$DEMUX" \
      --p-front-f CCTAYGGGRBGCASCAG \
      --p-front-r GGACTACNNGGGTATCTAAT \
      --o-trimmed-sequences "$OUT/trimmed-seqs.qza" \
      --verbose
else
    qiime cutadapt trim-single \
      --i-demultiplexed-sequences "$DEMUX" \
      --p-front CCTGTTCGATACCCGCACTTTCGAGCTTCAG \
      --o-trimmed-sequences "$OUT/trimmed-seqs.qza" \
      --verbose
fi

qiime demux summarize \
  --i-data "$OUT/trimmed-seqs.qza" \
  --o-visualization "$OUT/trimmed-seqs.qzv"

TRIMMED="$OUT/trimmed-seqs.qza"

# -----------------------------
# 5. DADA2
# -----------------------------
if [[ "$SEQ_TYPE" == "paired" ]]; then
    qiime dada2 denoise-paired \
      --i-demultiplexed-seqs "$DEMUX" \
      --p-trunc-len-f 0 \
      --p-trunc-len-r 0 \
      --p-max-ee-f 2 \
      --p-max-ee-r 2 \
      --o-table "$OUT/table.qza" \
      --o-representative-sequences "$OUT/rep-seqs.qza" \
      --o-denoising-stats "$OUT/stats.qza" \
      --p-n-threads 8
else
    qiime dada2 denoise-single \
      --i-demultiplexed-seqs "$DEMUX" \
      --p-trunc-len 0 \
      --p-max-ee 2 \
      --o-table "$OUT/table.qza" \
      --o-representative-sequences "$OUT/rep-seqs.qza" \
      --o-denoising-stats "$OUT/stats.qza" \
      --p-n-threads 8
fi

# -----------------------------
# 6.  SINGLETONS FILTERING
# -----------------------------
qiime feature-table filter-features \
  --i-table "$OUT/table.qza" \
  --p-min-frequency 2 \
  --o-filtered-table "$OUT/table-no-singletons.qza"

qiime feature-table filter-seqs \
  --i-data "$OUT/rep-seqs.qza" \
  --i-table "$OUT/table-no-singletons.qza" \
  --o-filtered-data "$OUT/rep-seqs-no-singletons.qza"

TABLE="$OUT/table-no-singletons.qza"
REPSEQS="$OUT/rep-seqs-no-singletons.qza"

# -----------------------------
# 7. RESULTS VISUALIZATION
# -----------------------------
qiime metadata tabulate \
  --m-input-file "$OUT/stats.qza" \
  --o-visualization "$OUT/stats.qzv"

qiime feature-table summarize \
  --i-table "$TABLE" \
  --o-visualization "$OUT/table.qzv"

qiime feature-table tabulate-seqs \
  --i-data "$REPSEQS" \
  --o-visualization "$OUT/rep-seqs.qzv"

# -----------------------------
# 8. FEATURE TABLE EXPORT FOR R
# -----------------------------
qiime tools export \
  --input-path "$TABLE" \
  --output-path "$OUT/exported-feature-table"

R ANALYSIS
biom convert \
  -i "$OUT/exported-feature-table/feature-table.biom" \
  -o "$OUT/feature-table.tsv" \
  --to-tsv

METADATA="/mnt/c/BIOINFORMATICA/PRJNA867698/HD_SAMPLES/MD_HD.tsv"
# -----------------------------
# 10. TAXONOMIC CLASSIFICATION
# -----------------------------
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/c/BIOINFORMATICA/CODIGOS/silva-138-99-nb-classifier.qza \
  --i-reads "$REPSEQS" \
  --o-classification "$OUT/taxonomy.qza" \
  --p-n-jobs 1 \
  --p-reads-per-batch 1000

EXPORT TAXONOMY FOR R
qiime tools export \
  --input-path "$OUT/taxonomy.qza" \
  --output-path "$OUT/exported-taxonomy"

BAR PLOTS QIIME2
qiime taxa barplot \
  --i-table "$TABLE" \
  --i-taxonomy "$OUT/taxonomy.qza" \
  --m-metadata-file "$METADATA" \
  --o-visualization "$OUT/taxa-bar-plots.qzv"

# -----------------------------
# 12. DIVERSITY METRICS
# -----------------------------
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$OUT/rooted-tree.qza" \
  --i-table "$TABLE" \
  --p-sampling-depth 2500 \
  --m-metadata-file "$METADATA" \
  --output-dir "$OUT/core-metrics-results"

echo "[INFO] FINISHED PRE-PROCESSING. ALL OUTPUTS WERE EXPORTED TO: $OUT"

# =========================================================
# QIIME2 -> PICRUSt2 Pipeline + KO export
# =========================================================


PICRUST2_OUT="$OUT/picrust2_out"
PICRUST2_PIPELINE="$PICRUST2_OUT/picrust2_out_pipeline"
# mkdir -p "$PICRUST2_OUT"

# # -----------------------------
# # 12. Export ASV sequences and feature table for PICRUSt2
# # -----------------------------
# qiime tools export \
#   --input-path "$OUT/rep-seqs.qza" \
#   --output-path "$PICRUST2_OUT/rep_seqs_fasta"

# qiime tools export \
#   --input-path "$OUT/table.qza" \
#   --output-path "$PICRUST2_OUT/table_biom"

# biom convert \
#   -i "$PICRUST2_OUT/table_biom/feature-table.biom" \
#   -o "$PICRUST2_OUT/table_raw.tsv" \
#   --to-tsv

# -----------------------------
# 13. Run PICRUSt2
# -----------------------------
THREADS=4

picrust2_pipeline.py \
  -s "$PICRUST2_OUT/rep_seqs_fasta/dna-sequences.fasta" \
  -i "$OUT/exported-feature-table/feature-table.biom" \
  -o "$PICRUST2_PIPELINE" \
  -p "$THREADS" \
  --per_sequence_contrib \
  --stratified \
  --min_align 0.1 \
  --verbose

# # -----------------------------
# # 14. Export PICRUSt2 results for R
# # -----------------------------
# RUN INDIVIDUALLY
mkdir -p $PICRUST2_OUT/R_STUDIO
# EC tables
cp "$PICRUST2_PIPELINE/EC_metagenome_out/pred_metagenome_unstrat.tsv" "$PICRUST2_OUT/R_STUDIO/EC_unstrat.tsv"
cp "$PICRUST2_PIPELINE/EC_metagenome_out/pred_metagenome_strat.tsv"   "$PICRUST2_OUT/R_STUDIO/EC_strat.tsv"

# Pathway tables
cp "$PICRUST2_PIPELINE/pathways_out/path_abun_unstrat.tsv" "$PICRUST2_OUT/R_STUDIO/pathways_unstrat.tsv"
cp "$PICRUST2_PIPELINE/pathways_out/path_abun_strat.tsv"   "$PICRUST2_OUT/R_STUDIO/pathways_strat.tsv"

# KO tables
cp "$PICRUST2_PIPELINE/KO_metagenome_out/pred_metagenome_unstrat.tsv" "$PICRUST2_OUT/R_STUDIO/KO_unstrat.tsv"
cp "$PICRUST2_PIPELINE/KO_metagenome_out/pred_metagenome_strat.tsv"   "$PICRUST2_OUT/R_STUDIO/KO_strat.tsv"

echo "[INFO] PICRUSt2 complete. Files exported to: $PICRUST2_OUT