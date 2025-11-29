#!/usr/bin/env bash

################################################################################
# AMR Radar Master Pipeline
# Automates all steps from raw data ? AMR profiling ? comparative genomics ? 
# biomarker discovery.
# Status will be printed step-by-step.
# This version asks for REQUIRED user inputs before starting.
################################################################################

set -euo pipefail
IFS=$'\n\t'

# --------------------- Ask for User Input ---------------------
echo "===== AMR-Radar Pipeline Setup ====="

read -rp "Enter Kraken2 DB Path (default: /data/db/kraken2/standard): " KRAKEN_DB
KRAKEN_DB="${KRAKEN_DB:-/data/db/kraken2/standard}"

read -rp "Enter THREADS (e.g., 40): " THREADS
read -rp "Enter PARALLEL samples (e.g., 10): " PARALLEL
read -rp "Enter ORGANISM to extract (e.g., Pseudomonas aeruginosa): " ORGANISM
read -rp "Enter PERCENT_OF (total OR classified): " PERCENT_OF
read -rp "Enter READS_DIR (default: raw): " READS_DIR
READS_DIR="${READS_DIR:-raw}"
read -rp "Enter INPUT_DIR (default: input): " INPUT_DIR
INPUT_DIR="${INPUT_DIR:-input}"

RESULTS_DIR="results"
KRAKEN_DIR="${RESULTS_DIR}/tools/0_kraken2"
ASSEMBLY_DIR="${RESULTS_DIR}/tools/3_assembly"
SAMPLE_TSV="sample.tsv"

# --------------------- Helpers ---------------------
section() { echo -e "\n============================"; echo "$1"; echo "============================"; }
log()     { echo "[INFO] $1"; }
error()   { echo "[ERROR] $1"; exit 1; }


# --------------------- Start -----------------------
section "Checking required files and environment"

REQUIRED_FILES=(
  "run_kraken2_identify_summary.sh"
  "kraken2_summary.sh"
  "extract_kraken2_by_name_auto_fastq.sh"
  "stage_extracted_fastq.sh"
  "summarize_input_counts.sh"
  "make_input_tsv.sh"
  "One_Health-AMR-Radar.sh"
  "One_Health-AMRfinder-module.sh"
  "extract_gene_and_class_amr.py"
  "amr_statistics.sh"
  "One_Health_AMR_analysis.py"
  "amr_network_mge_enrichment.py"
  "amr_analysis_pipeline.py"
  "run_resistome_extended.sh"
  "run_step1_3_compare.sh"
  "One_Health-AMR-Radar-prokka.sh"
  "One_Health-AMR-Radar-pangenome.sh"
  "run_biomarker_pipeline.sh"
  "make_genus_chart_html.py"
  "make_organism_chart_html.py"
)

for f in "${REQUIRED_FILES[@]}"; do
  [[ -e "$f" ]] || error "Missing required file: $f"
done
log "All required scripts found."

command -v python3 &>/dev/null || error "Python3 is required but not found."

[[ -f $SAMPLE_TSV ]] || error "sample.tsv missing! Create it with ./make_samples_tsv.sh"

mkdir -p "$RESULTS_DIR" "$INPUT_DIR"

# -------------------- FULL WORKFLOW --------------------------
section "Running Kraken2 (reads level)"
./run_kraken2_identify_summary.sh --kraken_db "$KRAKEN_DB" --mode reads --reads_dir "$READS_DIR" \
        --threads "$THREADS" --parallel "$PARALLEL" --min_percent 0.1 --percent_of "$PERCENT_OF"

section "Running Kraken2 (assembly level)"
./run_kraken2_identify_summary.sh --kraken_db "$KRAKEN_DB" --mode assembly --assembly_dir "$ASSEMBLY_DIR" \
        --threads "$THREADS" --parallel 1 --min_percent 0.1 --percent_of "$PERCENT_OF"

./kraken2_summary.sh

section "Extracting Organism-Specific Reads"
./extract_kraken2_by_name_auto_fastq.sh "$ORGANISM" exact --kraken_dir "$KRAKEN_DIR"
python3 make_genus_chart_html.py
python3 make_organism_chart_html.py

section "Staging FASTQs to input/"
./stage_extracted_fastq.sh "$ORGANISM" "$INPUT_DIR/" --kraken_dir "$KRAKEN_DIR"
./summarize_input_counts.sh "$INPUT_DIR/"
read -rp "Check results/summary_input_count.xls manually. Press ENTER to continue..."

section "Regenerate sample.tsv"
./make_input_tsv.sh "$INPUT_DIR"

section "Run AMR Radar main module"
./One_Health-AMR-Radar.sh -i "$SAMPLE_TSV" -o "$RESULTS_DIR" -t "$THREADS" -p "$PARALLEL"

section "Extract AMR Genes"
./One_Health-AMRfinder-module.sh
python3 extract_gene_and_class_amr.py --dedup --mapping --flank 0 --translate

section "Statistics"
./amr_statistics.sh --run-python One_Health_AMR_analysis.py
python3 amr_network_mge_enrichment.py results/OneHealthAMR_AMRFinder_summary.xlsx
python3 amr_analysis_pipeline.py

section "Resistome Extended"
./run_resistome_extended.sh -s "$SAMPLE_TSV" -o "$RESULTS_DIR" -t "$THREADS"

section "Comparative Genomics"
./run_step1_3_compare.sh

section "Resistome Comparison Charts"
python3 make_fastani_resistome_charts.py

section "Prokka Annotation + Pangenome"
./One_Health-AMR-Radar-prokka.sh --sample-file "$SAMPLE_TSV" -r "$RESULTS_DIR" -t "$THREADS" -c "$THREADS"
./One_Health-AMR-Radar-pangenome.sh -r "$RESULTS_DIR" -t "$THREADS"

section "Biomarker Discovery"
./run_biomarker_pipeline.sh

section "COMPLETED!"
log "Pipeline finished. Results in: $RESULTS_DIR/"
