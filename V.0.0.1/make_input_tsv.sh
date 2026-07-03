#!/usr/bin/env bash
set -euo pipefail

# ===============================================================
#  make_input_tsv.sh   (UPDATED FOR CSV GROUP MAPPING)
# ---------------------------------------------------------------
#  Reads FASTQ files in a directory, generates sample.tsv with:
#      sample   r1   r2   group
#
#  Group column is taken from a CSV mapping file:
#      sample_group_mapping.csv (columns: sample,group,...)
#
#  Usage:
#    ./make_input_tsv.sh [FASTQ_DIR] [DEFAULT_GROUP] [MAPPING_CSV]
#
#  Defaults:
#    FASTQ_DIR     = input
#    DEFAULT_GROUP = unknown
#    MAPPING_CSV   = sample_group_mapping.csv
# ===============================================================

# --- Args & defaults ---
FASTQ_DIR="${1:-input}"
DEFAULT_GROUP="${2:-unknown}"
MAPPING_CSV="${3:-sample_group_mapping.csv}"

# Resolve paths
FASTQ_DIR="${FASTQ_DIR/#\~/$HOME}"
FASTQ_DIR="$(realpath "$FASTQ_DIR" 2>/dev/null || echo "$FASTQ_DIR")"

OUT_TSV="sample.tsv"
BAK_TSV="sample0.tsv"

echo "== make_input_tsv.sh =="
echo "FASTQ_DIR      : $FASTQ_DIR"
echo "MAPPING_CSV    : $MAPPING_CSV"
echo "INPUT file     : sample.tsv (if exists, will be backed up)"
echo "BACKUP file    : sample0.tsv"
echo "OUTPUT file    : sample.tsv"
echo "DEFAULT_GROUP  : $DEFAULT_GROUP"
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

# --- Load group info from CSV mapping file ---
declare -A GROUP_MAP

if [[ -f "$MAPPING_CSV" ]]; then
  echo "Loading group mapping from CSV: $MAPPING_CSV"
  first=1
  while IFS=, read -r s g rest; do
    # Skip empty lines
    [[ -z "$s" ]] && continue

    # Skip header if present (e.g. 'sample,group')
    if (( first )); then
      first=0
      shdr="${s,,}"
      if [[ "$shdr" == "sample" ]]; then
        continue
      fi
    fi

    # Trim possible spaces
    s="${s%%[[:space:]]*}"
    g="${g%%[[:space:]]*}"

    [[ -z "$s" ]] && continue
    GROUP_MAP["$s"]="$g"
  done < "$MAPPING_CSV"
  echo "Loaded ${#GROUP_MAP[@]} group entries from CSV."
else
  echo "WARNING: Mapping CSV '$MAPPING_CSV' not found."
  echo "         All samples will use default group: '$DEFAULT_GROUP'"
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

  # Derive sample name
  sample="${fname%%_R1*}"
  sample="${sample%%_1*}"
  sample="${sample%%.fastq*}"
  sample="${sample%%.fq*}"

  # Find R2
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

  # Get group from map or use default
  group="$DEFAULT_GROUP"
  if [[ -n "${GROUP_MAP[$sample]+set}" ]]; then
    group="${GROUP_MAP[$sample]}"
  fi

  # Avoid duplicates
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
