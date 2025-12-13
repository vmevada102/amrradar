#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
#  ONE HEALTH AMR MASTER PIPELINE (End-to-End)
#  Phase 1: Raw Data -> Kraken2 ID -> Extraction -> Staging
#  Phase 2: Assembly (AMR Radar) -> AMRFinder -> Statistics
#  Phase 3: Comparative Genomics -> Pangenome -> Mapping -> MLST
# ==============================================================================

# --- Auto-Detect System Resources ---
# Detect total CPUs and subtract 2 for system stability
CPU_TOTAL=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
AUTO_THREADS=$((CPU_TOTAL - 2))
[[ "$AUTO_THREADS" -lt 1 ]] && AUTO_THREADS=1

# --- Defaults ---
DEFAULT_DB="/data/db/kraken2/standard"
DEFAULT_ORG="Pseudomonas aeruginosa"
DEFAULT_MODE="exact"
SAMPLE_TSV_NAME="sample.tsv"
STAGING_DIR="input"
RESULTS_DIR="results"

# --- SECTION 1: CONFIGURATION ---
echo "========================================================"
echo "?? PIPELINE CONFIGURATION"
echo "   Detected System CPUs: $CPU_TOTAL"
echo "   Auto-set Threads:     $AUTO_THREADS"
echo "========================================================"

# 1. Input Directory
while [[ -z "${READS_DIR:-}" ]]; do
    read -rp "1. Enter directory containing RAW FASTQ files: " READS_DIR
    if [[ ! -d "$READS_DIR" ]]; then echo "   ? Directory not found."; READS_DIR=""; fi
done

# 2. CSV Mapping File
while [[ -z "${MAP_FILE:-}" ]]; do
    read -rp "2. Enter path to CSV mapping file: " MAP_FILE
    if [[ ! -f "$MAP_FILE" ]]; then echo "   ? File not found."; MAP_FILE=""; fi
done

# 3. Kraken DB
read -rp "3. Enter Kraken2 DB path [default: $DEFAULT_DB]: " USER_DB
KRAKEN_DB="${USER_DB:-$DEFAULT_DB}"

# 4. Target Organism
read -rp "4. Enter Target Organism to extract [default: $DEFAULT_ORG]: " USER_ORG
TARGET_ORGANISM="${USER_ORG:-$DEFAULT_ORG}"

# 5. Parallel Samples
read -rp "5. Enter number of parallel samples [default: 4]: " USER_PARALLEL
PARALLEL_SAMPLES="${USER_PARALLEL:-4}"

# 6. Reference Accession (New)
read -rp "6. Enter NCBI Accession Number for Reference (e.g. NC_002516.2): " ACCESSION_NUM
while [[ -z "$ACCESSION_NUM" ]]; do
    echo "   ? Accession number is required."
    read -rp "6. Enter NCBI Accession Number for Reference: " ACCESSION_NUM
done

echo "========================================================"
echo "?? STARTING PIPELINE"
echo "   • Input Reads:  $READS_DIR"
echo "   • Mapping CSV:  $MAP_FILE"
echo "   • Target Org:   $TARGET_ORGANISM"
echo "   • Accession:    $ACCESSION_NUM"
echo "   • Threads:      $AUTO_THREADS"
echo "   • Parallelism:  $PARALLEL_SAMPLES"
echo "========================================================"

# --- ENSURE EXECUTABLES ---
chmod +x make_samples_tsv.sh \
         run_kraken2_identify_summary.sh \
         kraken2_summary.sh \
         extract_kraken2_by_name_auto_fastq.sh \
         stage_extracted_fastq.sh \
         summarize_input_counts.sh \
         make_input_tsv.sh \
         One_Health-AMR-Radar.sh \
         One_Health-AMRfinder-module.sh \
         amr_statistics.sh \
         download_fasta.sh \
         run_step1_3_compare.sh \
         One_Health-AMR-Radar-prokka.sh \
         One_Health-AMR-Radar-pangenome.sh \
         run_bacterial_mapping.sh \
         run_mlst.sh

# ==============================================================================
# PHASE 1: IDENTIFICATION & EXTRACTION
# ==============================================================================

echo ">> [Phase 1] Generating initial $SAMPLE_TSV_NAME..."
# Pipe inputs: Dir -> Filename -> Default Group -> Use Map (y) -> Map File
printf "%s\n%s\nunknown\ny\n%s\n" "$READS_DIR" "$SAMPLE_TSV_NAME" "$MAP_FILE" | ./make_samples_tsv.sh

echo ">> [Phase 1] Running Kraken2 identification..."
./run_kraken2_identify_summary.sh \
    --kraken_db "$KRAKEN_DB" \
    --mode reads \
    --reads_dir "$READS_DIR" \
    --threads "$AUTO_THREADS" \
    --parallel 1 \
    --min_percent 0.1 \
    --percent_of total \
    --sample_tsv "$SAMPLE_TSV_NAME"

echo ">> [Phase 1] Summarizing Kraken2 Results..."
./kraken2_summary.sh

echo ">> [Phase 1] Extracting sequences for: $TARGET_ORGANISM"
./extract_kraken2_by_name_auto_fastq.sh \
    "$TARGET_ORGANISM" \
    "$DEFAULT_MODE" \
    --kraken_dir results/tools/0_kraken2

echo ">> [Phase 1] Generating Interactive HTML Charts..."
python3 make_genus_chart_html.py || echo "Warning: make_genus_chart_html.py failed"
python3 make_organism_chart_html.py || echo "Warning: make_organism_chart_html.py failed"

echo ">> [Phase 1] Staging extracted FASTQ files to '$STAGING_DIR/'..."
mkdir -p "$STAGING_DIR"
./stage_extracted_fastq.sh "$TARGET_ORGANISM" "$STAGING_DIR/" --kraken_dir results/tools/0_kraken2

echo ">> [Phase 1] Summarizing counts in '$STAGING_DIR/'..."
./summarize_input_counts.sh "$STAGING_DIR/"

echo ">> [Phase 1] Generating new $SAMPLE_TSV_NAME for staged files..."
# This sample.tsv now points to the EXTRACTED reads in input/
./make_input_tsv.sh "$STAGING_DIR/" "unknown" "$MAP_FILE"

# ==============================================================================
# PHASE 2: ASSEMBLY & AMR ANALYSIS
# ==============================================================================

echo "--------------------------------------------------------"
echo ">> [Phase 2] Starting AMR Analysis on Extracted Reads"
echo "--------------------------------------------------------"

# Step 7: Run AMR Radar (Assembly & QC)
echo ">> [Step 7] Running AMR Radar (Assembly & QC)..."
./One_Health-AMR-Radar.sh \
    -i "$SAMPLE_TSV_NAME" \
    -o "$RESULTS_DIR" \
    -t "$AUTO_THREADS" \
    -p "$PARALLEL_SAMPLES" \
    --resume

# Step 8: Run AMRfinder Module
echo ">> [Step 8] Running AMRfinder Module..."
./One_Health-AMRfinder-module.sh \
    --results "$RESULTS_DIR" \
    --threads "$AUTO_THREADS" \
    --parallel "$PARALLEL_SAMPLES" \
    --resume

# Step 9-11: Statistical Analysis & Network
echo ">> [Step 9] Running Statistical Analysis & Network Enrichment..."
# 9a. General Statistics
./amr_statistics.sh --run-python One_Health_AMR_analysis.py || echo "Warning: amr_statistics.sh failed"

# 9b. Network MGE Enrichment
python3 amr_network_mge_enrichment.py "$RESULTS_DIR/OneHealthAMR_AMRFinder_summary.xlsx"

# 9c. Extract Specific Gene Sequences
echo ">> [Step 9c] Extracting Specific Gene Sequences..."
python3 extract_gene_and_class_amr.py --dedup --mapping --flank 0 --translate

# Step 10: Comprehensive Analysis Pipeline
echo ">> [Step 10] Running Comprehensive Analysis Pipeline..."
python3 amr_analysis_pipeline.py

# Step 11: Statistical Charts (T-test / Chi2)
echo ">> [Step 11] Running Statistics Pipeline (Charts/Tests)..."
python3 amr_statistics_pipeline.py

# Step 12: Generate AMRFinder Reports
echo ">> [Step 12] Generating Final AMRFinder Charts..."
python3 generate_amr_reports_autoenv.py

# ==============================================================================
# PHASE 3: COMPARATIVE GENOMICS & PANGENOME (New Steps)
# ==============================================================================

echo "--------------------------------------------------------"
echo ">> [Phase 3] Starting Comparative Genomics & Pangenome"
echo "--------------------------------------------------------"

# Step 15: Download Reference
echo ">> [Step 15] Downloading Reference FASTA for $ACCESSION_NUM..."
./download_fasta.sh "$ACCESSION_NUM"

# Step 12 (Ref from prompt): Comparative genomics analysis
echo ">> [Step 12] Running Comparative Genomics (QUAST/FastANI/Mashtree)..."
./run_step1_3_compare.sh --threads "$AUTO_THREADS"

# Step 13 (Ref from prompt): Resistome comparison charts
echo ">> [Step 13] Generating FastANI Resistome Charts..."
python3 make_fastani_resistome_charts.py || echo "Warning: make_fastani_resistome_charts.py failed (check if resistome_extended exists)"

# Step 14 (Ref from prompt): Annotation and Pangenome
echo ">> [Step 14] Running Prokka Annotation..."
# Note: Using AUTO_THREADS for -c (threads per job) and PARALLEL_SAMPLES for -t (parallel jobs)
./One_Health-AMR-Radar-prokka.sh \
    --sample-file "$SAMPLE_TSV_NAME" \
    -r "$RESULTS_DIR" \
    -t "$PARALLEL_SAMPLES" \
    -c "$AUTO_THREADS"

echo ">> [Step 14] Running Pangenome Analysis..."
./One_Health-AMR-Radar-pangenome.sh \
    -r "$RESULTS_DIR" \
    -t "$AUTO_THREADS"

# Step 16: Mapping and SNP Detection
echo ">> [Step 16] Running Bacterial Mapping & SNP Detection..."
# Uses results/sequence.fasta (downloaded in Step 15) by default
./run_bacterial_mapping.sh

# Step 17: MLST
echo ">> [Step 17] Running MLST Typing..."
./run_mlst.sh

echo "========================================================"
echo "?? COMPLETE PIPELINE FINISHED SUCCESSFULLY"
echo "   Results stored in: $RESULTS_DIR"
echo "   Charts stored in:  $RESULTS_DIR/Charts"
echo "========================================================"
