#!/usr/bin/env bash
set -u  # Removed -e to handle specific errors manually without exiting immediately
set -o pipefail

# ==============================================================================
#  ONE HEALTH AMR MASTER PIPELINE (End-to-End)
#  Phase 1: Raw Data -> Kraken2 ID -> Extraction -> Staging
#  Phase 2: Assembly (AMR Radar) -> AMRFinder -> Statistics
#  Phase 3: Comparative Genomics -> Pangenome -> Mapping -> MLST
#
#  Usage: ./run_complete_pipeline.sh [-resume]
# ==============================================================================

# --- Resume Logic ---
RESUME_FLAG=false
if [[ "${1:-}" == "-resume" ]]; then
    RESUME_FLAG=true
    echo ">> [INFO] Resume mode ENABLED. Skipping completed steps."
fi

# --- Auto-Detect System Resources ---
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
LOG_DIR="logs"
mkdir -p "$LOG_DIR"

# --- Helper Function for Resume ---
check_done() {
    local step_marker="$LOG_DIR/$1.done"
    if [[ "$RESUME_FLAG" == "true" && -f "$step_marker" ]]; then
        echo ">> [SKIP] Step '$1' already completed."
        return 0 # True, step is done
    fi
    return 1 # False, step needs running
}

mark_done() {
    touch "$LOG_DIR/$1.done"
}

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

# 6. Reference Accession
read -rp "6. Enter NCBI Accession Number for Reference (e.g. NC_002516.2): " ACCESSION_NUM
while [[ -z "$ACCESSION_NUM" ]]; do
    echo "   ? Accession number is required."
    read -rp "6. Enter NCBI Accession Number for Reference: " ACCESSION_NUM
done

echo "========================================================"
echo "?? STARTING PIPELINE"
echo "========================================================"

# --- ENSURE EXECUTABLES ---
chmod +x *.sh 2>/dev/null || true

# ==============================================================================
# PHASE 1: IDENTIFICATION & EXTRACTION
# ==============================================================================

if ! check_done "phase1_kraken"; then
    echo ">> [Phase 1] Generating initial $SAMPLE_TSV_NAME..."
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

    ./kraken2_summary.sh
    mark_done "phase1_kraken"
fi

if ! check_done "phase1_extract"; then
    echo ">> [Phase 1] Extracting sequences for: $TARGET_ORGANISM"
    ./extract_kraken2_by_name_auto_fastq.sh \
        "$TARGET_ORGANISM" \
        "$DEFAULT_MODE" \
        --kraken_dir results/tools/0_kraken2
    
    # Generate charts
    python3 make_genus_chart_html.py || echo "Warning: make_genus_chart_html.py failed"
    python3 make_organism_chart_html.py || echo "Warning: make_organism_chart_html.py failed"
    mark_done "phase1_extract"
fi

if ! check_done "phase1_stage"; then
    echo ">> [Phase 1] Staging extracted FASTQ files to '$STAGING_DIR/'..."
    mkdir -p "$STAGING_DIR"
    ./stage_extracted_fastq.sh "$TARGET_ORGANISM" "$STAGING_DIR/" --kraken_dir results/tools/0_kraken2
    ./summarize_input_counts.sh "$STAGING_DIR/"
    
    # Generate NEW sample.tsv for staged files
    ./make_input_tsv.sh "$STAGING_DIR/" "unknown" "$MAP_FILE"
    mark_done "phase1_stage"
fi

# ==============================================================================
# PHASE 2: ASSEMBLY & AMR ANALYSIS
# ==============================================================================

if ! check_done "phase2_assembly"; then
    echo ">> [Step 7] Running AMR Radar (Assembly & QC)..."
    ./One_Health-AMR-Radar.sh \
        -i "$SAMPLE_TSV_NAME" \
        -o "$RESULTS_DIR" \
        -t "$AUTO_THREADS" \
        -p "$PARALLEL_SAMPLES" \
        --resume
    mark_done "phase2_assembly"
fi

if ! check_done "phase2_amrfinder"; then
    echo ">> [Step 8] Running AMRfinder Module..."
    ./One_Health-AMRfinder-module.sh \
        --results "$RESULTS_DIR" \
        --threads "$AUTO_THREADS" \
        --parallel "$PARALLEL_SAMPLES" \
        --resume
    mark_done "phase2_amrfinder"
fi

if ! check_done "phase2_stats"; then
    echo ">> [Step 9-11] Running Statistical Analysis..."
    ./amr_statistics.sh --run-python One_Health_AMR_analysis.py || echo "Warning: amr_statistics.sh failed"
    python3 amr_network_mge_enrichment.py "$RESULTS_DIR/OneHealthAMR_AMRFinder_summary.xlsx"
    python3 extract_gene_and_class_amr.py --dedup --mapping --flank 0 --translate
    python3 amr_analysis_pipeline.py
    python3 amr_statistics_pipeline.py
    python3 generate_amr_reports_autoenv.py
    mark_done "phase2_stats"
fi

# ==============================================================================
# PHASE 3: COMPARATIVE GENOMICS & PANGENOME
# ==============================================================================

if ! check_done "phase3_reference"; then
    echo ">> [Step 15] Downloading Reference FASTA..."
    ./download_fasta.sh "$ACCESSION_NUM"
    mark_done "phase3_reference"
fi

if ! check_done "phase3_compgen"; then
    echo ">> [Step 12] Running Comparative Genomics..."
    ./run_step1_3_compare.sh --threads "$AUTO_THREADS"
    python3 make_fastani_resistome_charts.py || echo "Warning: resistome charts failed"
    mark_done "phase3_compgen"
fi

if ! check_done "phase3_prokka"; then
    echo ">> [Step 14] Running Prokka & Pangenome..."
    ./One_Health-AMR-Radar-prokka.sh \
        --sample-file "$SAMPLE_TSV_NAME" \
        -r "$RESULTS_DIR" \
        -t "$PARALLEL_SAMPLES" \
        -c "$AUTO_THREADS"
    
    ./One_Health-AMR-Radar-pangenome.sh \
        -r "$RESULTS_DIR" \
        -t "$AUTO_THREADS"
    mark_done "phase3_prokka"
fi

# --- Step 16: Mapping (With Error Handling) ---
if ! check_done "phase3_mapping"; then
    echo ">> [Step 16] Running Bacterial Mapping & SNP Detection..."
    echo ">> Note: Errors in this step (e.g. length mismatch) will be logged but won't stop the pipeline."
    
    # We trap the exit code to ensure the script continues even if mapping fails
    set +e # Disable exit-on-error specifically for this block
    
    ./run_bacterial_mapping.sh -resume 2>&1 | tee "$LOG_DIR/mapping_step_error.log"
    MAPPING_EXIT_CODE=${PIPESTATUS[0]}
    
    set -e # Re-enable exit-on-error
    
    if [ $MAPPING_EXIT_CODE -ne 0 ]; then
        echo "??  [WARNING] Step 16 (Mapping) encountered errors. Check logs at: $LOG_DIR/mapping_step_error.log"
        echo "??  Proceeding to Step 17 (MLST)..."
    else
        echo "? Step 16 completed successfully."
        mark_done "phase3_mapping"
    fi
fi

# --- Step 17: MLST ---
if ! check_done "phase3_mlst"; then
    echo ">> [Step 17] Running MLST Typing..."
    ./run_mlst.sh
    mark_done "phase3_mlst"
fi

echo "========================================================"
echo "?? COMPLETE PIPELINE FINISHED SUCCESSFULLY"
echo "   Results stored in: $RESULTS_DIR"
echo "   Charts stored in:  $RESULTS_DIR/Charts"
echo "   Logs stored in:    $LOG_DIR"
echo "========================================================"
