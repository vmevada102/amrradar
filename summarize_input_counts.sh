#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# summarize_input_counts.sh
# Count reads in R1 and R2 FASTQ files and write a TSV (.xls) summary.
#
# Usage:
#   ./summarize_input_counts.sh <INPUT_DIR> [OUTPUT_XLS]
#
# Examples:
#   ./summarize_input_counts.sh input/
#   ./summarize_input_counts.sh input/ results/summary_input_count.xls
#
# Notes:
#   - Detects pairs: <sample>_R1.(fastq|fq)[.gz] and <sample>_R2.*
#   - Handles .gz or plain FASTQ files
#   - Output is a tab-separated file readable by Excel
#
# Options:
#   -h | --help   Show this help message and exit
# ------------------------------------------------------------

print_help() {
  grep '^# ' "$0" | sed 's/^# //'
}

if [[ $# -eq 0 ]]; then
  echo "Error: No input directory provided. Use -h for help."
  exit 1
fi

case "${1:-}" in
  -h|--help) print_help; exit 0 ;;
esac

INPUT_DIR="$1"
OUTFILE="${2:-results/summary_input_count.xls}"
mkdir -p "$(dirname "$OUTFILE")"

# Count reads = lines/4 (supports gz)
count_reads() {
  local f="$1"
  if [[ ! -s "$f" ]]; then
    echo 0
    return
  fi
  if [[ "$f" == *.gz ]]; then
    gzip -cd -- "$f" | awk 'END{print (NR>0?NR/4:0)}'
  else
    awk 'END{print (NR>0?NR/4:0)}' "$f"
  fi
}

# Header
echo -e "samplename\tr1\tr2" > "$OUTFILE"

# Find all R1 files safely
shopt -s nullglob
declare -a r1_files=()
for ext in fastq fq fastq.gz fq.gz; do
  for f in "$INPUT_DIR"/*_R1.$ext; do
    [[ -e "$f" ]] && r1_files+=("$f")
  done
done

if (( ${#r1_files[@]} == 0 )); then
  echo "No R1 FASTQ files found in: $INPUT_DIR"
  echo "? Summary written to: $OUTFILE"
  exit 0
fi

declare -A seen
for R1 in "${r1_files[@]}"; do
  base="$(basename "$R1")"
  sample="${base%_R1*}"     # everything before first _R1

  # Skip duplicates if multiple extensions of same sample exist
  [[ -n "${seen[$sample]+x}" ]] && continue
  seen[$sample]=1

  # Find matching R2 (any known extension)
  R2=""
  for ext in fastq.gz fq.gz fastq fq; do
    cand="$INPUT_DIR/${sample}_R2.$ext"
    if [[ -e "$cand" ]]; then
      R2="$cand"
      break
    fi
  done

  r1c="$(count_reads "$R1")"
  if [[ -z "$R2" ]]; then
    echo "WARNING: Missing R2 for sample '$sample'" >&2
    r2c=0
  else
    r2c="$(count_reads "$R2")"
  fi

  echo -e "${sample}\t${r1c}\t${r2c}" >> "$OUTFILE"
done

echo "? Summary written to: $OUTFILE"
