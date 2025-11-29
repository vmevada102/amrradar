#!/usr/bin/env bash
set -euo pipefail

# ===============================================================
#  make_input_tsv.sh   (FINAL VERSION AS PER REQUIREMENT)
# ---------------------------------------------------------------
#  Reads existing sample.tsv, backs it up ? sample0.tsv,
#  regenerates sample.tsv using FASTQ files in input/ directory,
#  copies group column from sample0.tsv (based on sample name).
#
#  Usage:
#    ./make_input_tsv.sh [FASTQ_DIR] [DEFAULT_GROUP]
#
#  Defaults:
#    FASTQ_DIR     = input
#    DEFAULT_GROUP = unknown
# ===============================================================

# --- Args & defaults ---
FASTQ_DIR="${1:-input}"
DEFAULT_GROUP="${2:-unknown}"

# Resolve paths
FASTQ_DIR="${FASTQ_DIR/#\~/$HOME}"
FASTQ_DIR="$(realpath "$FASTQ_DIR" 2>/dev/null || echo "$FASTQ_DIR")"

OUT_TSV="sample.tsv"
BAK_TSV="sample0.tsv"

echo "== make_input_tsv.sh =="
echo "FASTQ_DIR     : $FASTQ_DIR"
echo "INPUT file    : sample.tsv"
echo "BACKUP file   : sample0.tsv"
echo "OUTPUT file   : sample.tsv"
echo "DEFAULT_GROUP : $DEFAULT_GROUP"
echo "---------------------------------------"

# --- Check input directory ---
if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "ERROR: FASTQ directory '$FASTQ_DIR' does not exist!" >&2
  exit 1
fi

# --- Backup existing sample.tsv to sample0.tsv ---
if [[ -f "$OUT_TSV" ]]; then
  echo "Backing up existing sample.tsv -> sample0.tsv"
  cp -f "$OUT_TSV" "$BAK_TSV"
else
  echo "No existing sample.tsv found — creating new one."
fi

# --- Load group info from backup (if exists) ---
declare -A GROUP_MAP

if [[ -f "$BAK_TSV" ]]; then
  echo "Loading group mapping from sample0.tsv..."
  while IFS=$'\t' read -r s r1 r2 g; do
    [[ -z "$s" ]] && continue
    [[ "$s" == "sample" ]] && continue
    GROUP_MAP["$s"]="$g"
  done < "$BAK_TSV"
  echo "Loaded ${#GROUP_MAP[@]} group entries."
else
  echo "WARNING: No sample0.tsv found. All samples will get group: '$DEFAULT_GROUP'"
fi

# --- Create new sample.tsv ---
echo -e "sample\tr1\tr2\tgroup" > "$OUT_TSV"

tmpfile=$(mktemp)
trap 'rm -f "$tmpfile"' EXIT

find "$FASTQ_DIR" -type f \( -iname "*_R1*.fastq*" -o -iname "*_1*.fastq*" \) -print0 > "$tmpfile"

if [[ ! -s "$tmpfile" ]]; then
  echo "ERROR: No R1 FASTQ files found in $FASTQ_DIR" >&2
  exit 2
fi

declare -A ADDED

while IFS= read -r -d '' R1; do
  fname="$(basename "$R1")"

  sample="${fname%%_R1*}"
  sample="${sample%%_1*}"
  sample="${sample%%.fastq*}"
  sample="${sample%%.fq*}"

  R2=""
  for p in "_R2" "_2" "-2" ".2"; do
    cand="${R1/_R1/$p}"
    cand="${cand/_1/$p}"
    if [[ -f "$cand" ]]; then
      R2="$cand"
      break
    fi
  done

  if [[ -z "$R2" ]]; then
    echo "Warning: No R2 found for sample '$sample' (R1: $R1)" >&2
    continue
  fi

  group="$DEFAULT_GROUP"
  if [[ -n "${GROUP_MAP[$sample]+set}" ]]; then
    group="${GROUP_MAP[$sample]}"
  fi

  if [[ -z "${ADDED[$sample]:-}" ]]; then
    printf "%s\t%s\t%s\t%s\n" "$sample" "$R1" "$R2" "$group" >> "$OUT_TSV"
    ADDED["$sample"]=1
  fi

done < "$tmpfile"

echo "---------------------------------------"
echo "NEW sample.tsv created with $(($(wc -l < "$OUT_TSV") - 1)) samples"
echo "Backup file:  sample0.tsv"
echo "Default group used for missing entries: '$DEFAULT_GROUP'"
echo "DONE ?"
