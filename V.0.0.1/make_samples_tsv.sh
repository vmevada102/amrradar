#!/usr/bin/env bash
set -euo pipefail

# ===============================================================
#  make_samples_tsv.sh  (UPDATED FOR CSV GROUP MAPPING)
# ---------------------------------------------------------------
#  Interactive script to generate samples.tsv for paired-end FASTQ
#  Columns: sample    r1    r2    group
#
#  Group info can be provided via a CSV mapping file:
#      sample,group,...
# ===============================================================

echo "👋 Welcome to the samples.tsv generator!"
echo "---------------------------------------"

# --- Ask for FASTQ folder ---
read -rp "📁 Enter the folder path where your FASTQ files are located: " FASTQ_DIR
FASTQ_DIR="${FASTQ_DIR/#\~/$HOME}"   # expand ~
FASTQ_DIR="$(realpath "$FASTQ_DIR" 2>/dev/null || echo "$FASTQ_DIR")"

if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "❌ ERROR: Folder '$FASTQ_DIR' does not exist!"
  exit 1
fi

# --- Ask for output TSV file ---
read -rp "📝 Enter output TSV file name [default: samples.tsv]: " OUT_TSV
OUT_TSV="${OUT_TSV:-samples.tsv}"
OUT_TSV="$(realpath "$OUT_TSV")"

# --- Ask for default group ---
read -rp "👥 Enter default group name (e.g., disease / healthy / environment): " DEFAULT_GROUP
DEFAULT_GROUP="${DEFAULT_GROUP:-unknown}"

# --- Optionally use mapping file ---
read -rp "📄 Do you have a CSV mapping file with sample,group info? (y/n): " use_map
MAPPING_FILE=""
declare -A GROUP_MAP
if [[ "$use_map" =~ ^[Yy]$ ]]; then
  read -rp "➡️  Enter mapping CSV file path [e.g., sample_group_mapping.csv]: " MAPPING_FILE
  MAPPING_FILE="${MAPPING_FILE/#\~/$HOME}"
  if [[ ! -f "$MAPPING_FILE" ]]; then
    echo "❌ Mapping file '$MAPPING_FILE' not found!"
    exit 2
  fi
  echo "✅ Using mapping CSV file: $MAPPING_FILE"

  first=1
  while IFS=, read -r s g rest; do
    # Skip empty lines
    [[ -z "$s" ]] && continue

    # Skip header if present
    if (( first )); then
      first=0
      shdr="${s,,}"
      if [[ "$shdr" == "sample" ]]; then
        continue
      fi
    fi

    # Trim trailing spaces
    s="${s%%[[:space:]]*}"
    g="${g%%[[:space:]]*}"

    [[ -z "$s" ]] && continue
    GROUP_MAP["$s"]="$g"
  done < "$MAPPING_FILE"

  echo "ℹ️  Loaded ${#GROUP_MAP[@]} group entries from CSV."
fi

# --- Prepare output ---
mkdir -p "$(dirname "$OUT_TSV")"
echo -e "sample\tr1\tr2\tgroup" > "$OUT_TSV"

# --- Detect R1 FASTQs ---
tmpfile=$(mktemp)
trap 'rm -f "$tmpfile"' EXIT

find "$FASTQ_DIR" -type f \( -iname "*_R1*.fastq*" -o -iname "*_1*.fastq*" \) -print0 > "$tmpfile"

if [[ ! -s "$tmpfile" ]]; then
  echo "❌ No FASTQ files with R1 pattern found in $FASTQ_DIR"
  exit 3
fi

# --- Pair and write samples ---
declare -A ADDED
while IFS= read -r -d '' R1; do
  fname="$(basename "$R1")"

  sample="${fname%%_R1*}"
  sample="${sample%%_1*}"
  sample="${sample%%.fastq*}"
  sample="${sample%%.fq*}"

  R2=""
  for pattern in "_R2" "_2" "-2" ".2"; do
    cand="${R1/_R1/$pattern}"
    cand="${cand/_1/$pattern}"
    if [[ -f "$cand" ]]; then
      R2="$cand"
      break
    fi
  done

  if [[ -z "$R2" ]]; then
    echo "⚠️  Warning: No R2 found for sample $sample (R1: $R1)" >&2
    continue
  fi

  # Get group from CSV map or default
  group="${GROUP_MAP[$sample]:-$DEFAULT_GROUP}"

  if [[ -z "${ADDED[$sample]:-}" ]]; then
    echo -e "${sample}\t${R1}\t${R2}\t${group}" >> "$OUT_TSV"
    ADDED["$sample"]=1
  fi
done < "$tmpfile"

# --- Show result ---
echo ""
echo "✅ Created: $OUT_TSV"
echo "----------------------------------------------------"
column -t -s $'\t' "$OUT_TSV" | head -n 10 || true
echo "----------------------------------------------------"
echo "📊 Total samples detected: $(($(wc -l < "$OUT_TSV") - 1))"
echo "👥 Default group used when mapping missing: '$DEFAULT_GROUP'"
echo "Done! 🎉"
