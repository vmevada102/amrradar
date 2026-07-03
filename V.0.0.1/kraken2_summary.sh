#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
#  kraken2_summary.sh
# ------------------------------------------------------------------------------
#  Summarise all Kraken2 log files:
#    results/tools/0_kraken2/<sample>.kraken2.log
#
#  Output:
#    results/summary_kraken2.xls  (tab-separated, Excel-compatible)
#
#  Columns:
#    Sr No
#    Sample Name
#    No of sequences Processed
#    No of sequence classified
#    No of sequence unclassified
#    % classified Sequences
#    % unclassified Sequences
#
#  Usage:
#    ./kraken2_summary.sh
# ==============================================================================

TOOLS_DIR="results/tools/0_kraken2"
RESULTS_DIR="results"
OUTFILE="${RESULTS_DIR}/summary_kraken2.xls"

# --- checks ---
if [[ ! -d "$TOOLS_DIR" ]]; then
  echo "ERROR: directory '$TOOLS_DIR' not found." >&2
  exit 1
fi

mkdir -p "$RESULTS_DIR"

shopt -s nullglob
logs=( "${TOOLS_DIR}"/*.kraken2.log )
shopt -u nullglob

if (( ${#logs[@]} == 0 )); then
  echo "ERROR: No *.kraken2.log files found in '$TOOLS_DIR'" >&2
  exit 1
fi

# --- write header ---
{
  echo -e "Sr No\tSample Name\tNo of sequences Processed\tNo of sequence classified\tNo of sequence unclassified\t% classified Sequences\t% unclassified Sequences"

  sr=1
  for log in "${logs[@]}"; do
    base="$(basename "$log")"
    sample="${base%.kraken2.log}"

    # total sequences processed
    total_seq=$(awk '/sequences .*processed/ {print $1; exit}' "$log")

    # classified
    classified_seq=$(awk '/sequences classified/ {print $1; exit}' "$log")
    classified_pct=$(grep -m1 'sequences classified' "$log" \
                     | sed -n 's/.*(\([0-9.]\+\)%.*/\1/p')

    # unclassified
    unclassified_seq=$(awk '/sequences unclassified/ {print $1; exit}' "$log")
    unclassified_pct=$(grep -m1 'sequences unclassified' "$log" \
                       | sed -n 's/.*(\([0-9.]\+\)%.*/\1/p')

    # fallback to 0 if something is missing
    total_seq=${total_seq:-0}
    classified_seq=${classified_seq:-0}
    unclassified_seq=${unclassified_seq:-0}
    classified_pct=${classified_pct:-0}
    unclassified_pct=${unclassified_pct:-0}

    echo -e "${sr}\t${sample}\t${total_seq}\t${classified_seq}\t${unclassified_seq}\t${classified_pct}\t${unclassified_pct}"

    sr=$((sr + 1))
  done
} > "$OUTFILE"

echo "Kraken2 summary written to: $OUTFILE"
