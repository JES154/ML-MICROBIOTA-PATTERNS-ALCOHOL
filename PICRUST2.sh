#!/bin/bash
set -e

# =========================================================
# 🧬 QIIME2 → PICRUSt2 Pipeline + KO export
# =========================================================

OUT_DIR="/mnt/c/BIOINFORMATICA_TESIS_RESPALDO/PRUEBA_QIIME2/qiime2_output"
PICRUST2_OUT="$OUT_DIR/picrust2_out"
PICRUST2_PIPELINE="$PICRUST2_OUT/picrust2_out_pipeline"

# -----------------------------
# 13. Run PICRUSt2
# -----------------------------
# --min_align 0.1: recommended for short reads (V3-V4)
picrust2_pipeline.py \
  -s "$PICRUST2_OUT/rep_seqs_fasta/dna-sequences.fasta" \
  -i "$OUT_DIR/exported/feature-table.biom" \
  -o "PICRUST2_PIPELINE" \
  -p 8 \
  --per_sequence_contrib \
  --stratified \
  --min_align 0.1 \
  --verbose
# ---------------------------------------------------------
# 3️⃣ Exportar resultados para R
# ---------------------------------------------------------
# Predicciones desglosadas y agregadas
cp "$PICRUST2_OUT/picrust2_out_pipeline_v4/EC_metagenome_out/pred_metagenome_unstrat.tsv" "$PICRUST2_OUT/EC_unstrat.tsv"
cp "$PICRUST2_OUT/picrust2_out_pipeline_v4/EC_metagenome_out/pred_metagenome_strat.tsv" "$PICRUST2_OUT/EC_strat.tsv"

cp "$PICRUST2_OUT/picrust2_out_pipeline_v4/pathways_out/path_abun_unstrat.tsv" "$PICRUST2_OUT/pathways_unstrat.tsv"
cp "$PICRUST2_OUT/picrust2_out_pipeline_v4/pathways_out/path_abun_strat.tsv" "$PICRUST2_OUT/pathways_strat.tsv"

# KO tables
cp "$PICRUST2_OUT/picrust2_out_pipeline_v4/KO_metagenome_out/pred_metagenome_unstrat.tsv" "$PICRUST2_OUT/KO_unstrat.tsv"
cp "$PICRUST2_OUT/picrust2_out_pipeline_v4/KO_metagenome_out/pred_metagenome_strat.tsv" "$PICRUST2_OUT/KO_strat.tsv"

echo "[✅ INFO] PICRUSt2 finalizado. Archivos exportados en $PICRUST2_OUT"
