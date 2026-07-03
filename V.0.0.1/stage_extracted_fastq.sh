#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# stage_extracted_fastq.sh
# ------------------------------------------------------------

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <organism_or_folder_label> <dest_input_dir> [--kraken_dir DIR] [--overwrite] [--dry-run]"
  exit 1
fi

LABEL_RAW="$1"; shift
DEST_DIR="$1"; shift || true

KRAKEN_DIR=""
OVERWRITE=0
DRYRUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --kraken_dir) KRAKEN_DIR="$2"; shift 2 ;;
    --overwrite)  OVERWRITE=1; shift ;;
    --dry-run)    DRYRUN=1; shift ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

# Determine Kraken folder prefix
if [[ -z "$KRAKEN_DIR" ]]; then
  KRAKEN_DIR=$(ls -d results/tools/*_kraken2 2>/dev/null | sort -V | tail -n1 || true)
fi
if [[ -z "$KRAKEN_DIR" || ! -d "$KRAKEN_DIR" ]]; then
  echo "ERROR: could not resolve kraken_dir"
  exit 1
fi

RUN_BASENAME=$(basename "$KRAKEN_DIR")
NUM_PREFIX="${RUN_BASENAME%%_kraken2}"
[[ "$NUM_PREFIX" == "$RUN_BASENAME" ]] && NUM_PREFIX="$RUN_BASENAME"

# Sanitization matching extractor behavior (SINGLE underscore)
sanitize_label() {
  local s="$1"
  s="${s//[ \/]/_}"      # replace space and slash with _
  echo "$s" | sed 's/[^A-Za-z0-9_.-]//g'
}

if [[ "$LABEL_RAW" =~ ^taxid_[0-9]+$ ]]; then
  ORGSAFE="$LABEL_RAW"
else
  ORGSAFE=$(sanitize_label "$LABEL_RAW")
fi

SRC_DIR="results/tools/${NUM_PREFIX}_kraken2_extracted/${ORGSAFE}"

if [[ ! -d "$SRC_DIR" ]]; then
  echo "ERROR: source folder not found: $SRC_DIR"
  echo "Available extracted folders:"
  find "results/tools" -maxdepth 2 -type d -name "*_kraken2_extracted" \
    | sed 's/^/  /'
  exit 1
fi

mkdir -p "$DEST_DIR"

copy_one() {
  local src="$1"
  local base=$(basename "$src")
  local out

  if [[ "$base" == *.fq.gz ]]; then
    out="${base%.fq.gz}.fastq.gz"
    cp "$src" "$DEST_DIR/$out"
  elif [[ "$base" == *.fq ]]; then
    out="${base%.fq}.fastq.gz"
    gzip -c "$src" > "$DEST_DIR/$out"
  elif [[ "$base" == *.fastq.gz ]]; then
    cp "$src" "$DEST_DIR/$base"
    return
  else
    return
  fi

  echo "  ? $DEST_DIR/$out"
}

echo "Staging files from: $SRC_DIR ? $DEST_DIR"
shopt -s nullglob
for f in "$SRC_DIR"/*; do
  copy_one "$f"
done

echo "? Done"
