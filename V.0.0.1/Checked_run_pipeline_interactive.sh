#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# MASTER PIPELINE
# 1. Config: Interactive setup (Reads, Map, DB, Organism).
# 2. Step 1: Generate initial sample.tsv (make_samples_tsv.sh).
# 3. Step 2: Run Kraken2 Analysis.
# 4. Step 3: Summarize Kraken2 Results.
# 5. Step 4: Extract Target Organism Reads.
# 6. Step 5: Generate HTML Charts.
# 7. Step 6: Stage Extracted Reads to 'input/' folder.
# 8. Step 7: Summarize Counts in 'input/'.
# 9. Step 8: Generate new sample.tsv for 'input/'.
# ==============================================================================

# --- Defaults ---
DEFAULT_DB="/data/db/kraken2/standard"
DEFAULT_ORG="Pseudomonas aeruginosa"
DEFAULT_MODE="exact"
SAMPLE_TSV_NAME="sample.tsv"
THREADS=42
PARALLEL=1
STAGING_DIR="input"  # Directory where extracted files will be moved

# --- SECTION 1: USER INPUTS ---
echo "========================================================"
echo "?? CONFIGURATION (Press Enter to accept defaults)"
echo "========================================================"

# 1. Ask for Input Directory (Required)
while [[ -z "${READS_DIR:-}" ]]; do
    read -rp "1. Enter the directory containing RAW FASTQ files: " READS_DIR
    if [[ ! -d "$READS_DIR" ]]; then
        echo "   ? Directory not found. Please try again."
        READS_DIR=""
    fi
done

# 2. Ask for CSV Mapping File (Required)
while [[ -z "${MAP_FILE:-}" ]]; do
    read -rp "2. Enter the path to your CSV mapping file: " MAP_FILE
    if [[ ! -f "$MAP_FILE" ]]; then
        echo "   ? File not found. Please try again."
        MAP_FILE=""
    fi
done

# 3. Ask for Kraken DB (Default provided)
read -rp "3. Enter Kraken2 DB path [default: $DEFAULT_DB]: " USER_DB
KRAKEN_DB="${USER_DB:-$DEFAULT_DB}"

# 4. Ask for Target Organism (Default provided)
read -rp "4. Enter Target Organism to extract [default: $DEFAULT_ORG]: " USER_ORG
TARGET_ORGANISM="${USER_ORG:-$DEFAULT_ORG}"

echo "========================================================"
echo "?? STARTING PIPELINE"
echo "   • Raw Reads: $READS_DIR"
echo "   • Mapping:   $MAP_FILE"
echo "   • DB:        $KRAKEN_DB"
echo "   • Target:    $TARGET_ORGANISM"
echo "   • Staging:   $STAGING_DIR/ (Extracted files will go here)"
echo "========================================================"

# --- ENSURE SCRIPTS ARE EXECUTABLE ---
chmod +x make_samples_tsv.sh \
         run_kraken2_identify_summary.sh \
         kraken2_summary.sh \
         extract_kraken2_by_name_auto_fastq.sh \
         stage_extracted_fastq.sh \
         summarize_input_counts.sh \
         make_input_tsv.sh

# ------------------------------------------------------------------------------
# STEP 1: Generate sample.tsv (Initial)
# ------------------------------------------------------------------------------
echo ">> [Step 1/8] Generating initial $SAMPLE_TSV_NAME..."

# Inputs piped: Dir -> Filename -> Default Group -> Use Map (y) -> Map File
printf "%s\n%s\nunknown\ny\n%s\n" "$READS_DIR" "$SAMPLE_TSV_NAME" "$MAP_FILE" | ./make_samples_tsv.sh

echo "? Step 1 Complete."
echo "--------------------------------------------------------"

# ------------------------------------------------------------------------------
# STEP 2: Run Kraken2 Classification
# ------------------------------------------------------------------------------
echo ">> [Step 2/8] Running Kraken2 identification..."

./run_kraken2_identify_summary.sh \
    --kraken_db "$KRAKEN_DB" \
    --mode reads \
    --reads_dir "$READS_DIR" \
    --threads "$THREADS" \
    --parallel "$PARALLEL" \
    --min_percent 0.1 \
    --percent_of total \
    --sample_tsv "$SAMPLE_TSV_NAME"

echo "? Step 2 Complete."
echo "--------------------------------------------------------"

# ------------------------------------------------------------------------------
# STEP 3: Summarize Kraken2 Results
# ------------------------------------------------------------------------------
echo ">> [Step 3/8] Generating Kraken2 Summary Report..."

./kraken2_summary.sh

echo "? Step 3 Complete."
echo "--------------------------------------------------------"

# ------------------------------------------------------------------------------
# STEP 4: Extract Specific Sequences
# ------------------------------------------------------------------------------
echo ">> [Step 4/8] Extracting sequences for: $TARGET_ORGANISM"

./extract_kraken2_by_name_auto_fastq.sh \
    "$TARGET_ORGANISM" \
    "$DEFAULT_MODE" \
    --kraken_dir results/tools/0_kraken2

echo "? Step 4 Complete."
echo "--------------------------------------------------------"

# ------------------------------------------------------------------------------
# STEP 5: Generate HTML Charts
# ------------------------------------------------------------------------------
echo ">> [Step 5/8] Generating Interactive HTML Charts..."

# Check if python is available
if command -v python3 &>/dev/null; then
    PYTHON_CMD="python3"
else
    PYTHON_CMD="python"
fi

$PYTHON_CMD make_genus_chart_html.py
$PYTHON_CMD make_organism_chart_html.py

echo "? Step 5 Complete: Charts created in results/Charts/"
echo "--------------------------------------------------------"

# ------------------------------------------------------------------------------
# STEP 6: Stage Extracted Files
# ------------------------------------------------------------------------------
echo ">> [Step 6/8] Staging extracted FASTQ files to '$STAGING_DIR/'..."

# Ensure input directory exists
mkdir -p "$STAGING_DIR"

./stage_extracted_fastq.sh \
    "$TARGET_ORGANISM" \
    "$STAGING_DIR/" \
    --kraken_dir results/tools/0_kraken2

echo "? Step 6 Complete."
echo "--------------------------------------------------------"

# ------------------------------------------------------------------------------
# STEP 7: Summarize Counts in Staging Directory
# ------------------------------------------------------------------------------
echo ">> [Step 7/8] Summarizing counts for files in '$STAGING_DIR/'..."

./summarize_input_counts.sh "$STAGING_DIR/"

echo "? Step 7 Complete: Check results/summary_input_count.xls"
echo "--------------------------------------------------------"

# ------------------------------------------------------------------------------
# STEP 8: Regenerate sample.tsv for Staging Directory
# ------------------------------------------------------------------------------
echo ">> [Step 8/8] Generating new sample.tsv for staged files..."

# Usage: ./make_input_tsv.sh [FASTQ_DIR] [DEFAULT_GROUP] [MAPPING_CSV]
# We use the Staging Directory as the new FASTQ source
./make_input_tsv.sh "$STAGING_DIR/" "unknown" "$MAP_FILE"

echo "? Step 8 Complete: New sample.tsv generated."
echo "========================================================"
echo "?? PIPELINE FINISHED SUCCESSFULLY"
echo "========================================================"