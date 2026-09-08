#!/bin/bash
# =============================================================
# Convert Excel metadata (.xlsx) to QIIME 2–ready TSV (validated)
# Author: Juan Esparza (optimized for microbiome QIIME2 workflows)
# Requires: csvkit (in2csv), awk, tr, grep, qiime
# =============================================================

# ----------------------------
# Parameters
# ----------------------------
INPUT_XLSX="$1"           # Input Excel file
OUTPUT_TSV="$2"           # Output QIIME 2 metadata TSV
REMOVE_UNDETERMINED=true  # Set to true to remove 'Undetermined' samples

# ----------------------------
# Check input
# ----------------------------
if [[ -z "$INPUT_XLSX" || -z "$OUTPUT_TSV" ]]; then
    echo "Usage: $0 input.xlsx output.tsv"
    exit 1
fi

if ! command -v in2csv &> /dev/null; then
    echo "❌ Error: in2csv not found. Install csvkit first (pip install csvkit)."
    exit 1
fi

# ----------------------------
# Convert Excel to temporary CSV
# ----------------------------
TMP_CSV=$(mktemp)
echo "[INFO] Converting Excel to CSV..."
in2csv "$INPUT_XLSX" > "$TMP_CSV"

# ----------------------------
# Extract header and replace first column name
# ----------------------------
HEADER=$(head -n 1 "$TMP_CSV" | tr -d '\r')
HEADER=$(echo "$HEADER" | awk -F, 'BEGIN{OFS="\t"} {$1="#SampleID"; print}')

# ----------------------------
# Process body
# ----------------------------
echo "[INFO] Cleaning metadata and converting commas to tabs..."
if [ "$REMOVE_UNDETERMINED" = true ]; then
    tail -n +2 "$TMP_CSV" | grep -v -i "Undetermined" | awk -F, 'BEGIN{OFS="\t"} {print}' > "${OUTPUT_TSV}.body"
else
    tail -n +2 "$TMP_CSV" | awk -F, 'BEGIN{OFS="\t"} {print}' > "${OUTPUT_TSV}.body"
fi

# Remove extra quotes, spaces, and empty lines
sed -i 's/"//g; s/ *$//; /^\s*$/d' "${OUTPUT_TSV}.body"

# ----------------------------
# Combine header and body
# ----------------------------
echo -e "$HEADER" > "$OUTPUT_TSV"
cat "${OUTPUT_TSV}.body" >> "$OUTPUT_TSV"

# ----------------------------
# Validate file structure
# ----------------------------
echo "[INFO] Checking for tabs and proper header..."
if grep -q "," "$OUTPUT_TSV"; then
    echo "⚠️ Warning: Commas detected — converting to tabs."
    tr ',' '\t' < "$OUTPUT_TSV" > "${OUTPUT_TSV}.tmp" && mv "${OUTPUT_TSV}.tmp" "$OUTPUT_TSV"
fi

# ----------------------------
# QIIME 2 metadata validation
# ----------------------------
if command -v qiime &> /dev/null; then
    echo "[INFO] Validating with QIIME 2..."
    if qiime metadata tabulate --m-input-file "$OUTPUT_TSV" --o-visualization metadata_validation.qzv 2>/dev/null; then
        echo "✅ Metadata validated successfully with QIIME 2."
        echo "   Visualization: metadata_validation.qzv"
    else
        echo "⚠️ QIIME2 validation failed — check formatting manually."
    fi
fi

# ----------------------------
# Cleanup
# ----------------------------
rm -f "$TMP_CSV" "${OUTPUT_TSV}.body"

echo "✅ QIIME 2 metadata ready: $OUTPUT_TSV"



# COMO USAR: 
# ./qiime_metadata.sh METADATA_60_118.xlsx METADATA_60_118.tsv


#./qiime_metadata.sh /mnt/c/BIOINFORMATICA/PRJNA517050/VHD_SAMPLES/MD_vHD.xlsx /mnt/c/BIOINFORMATICA/PRJNA517050/VHD_SAMPLES/MD_VHD.tsv