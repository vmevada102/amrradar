#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------------------
# extract_kraken2_by_name_auto_fastq.sh  (FASTQ, paired; taxid-or-name aware)
#
# Extract paired FASTQ reads (with qualities) assigned by Kraken2 to a
# target organism, using per-sample .report (for taxids) and .kraken
# (for read->taxon). Works whether .kraken col3 is TAXID or NAME (--use-names).
#
# Features:
#   • Exact (single node) or clade (subtree) from the .report
#   • Match by normalized scientific name OR by --taxid <ID>
#   • Auto-detect .kraken mode (taxid vs name)
#   • AUTO-finds mates: <sample>.classified__1.fq/.classified__2.fq (and variants)
#   • Keeps qualities and pairing; supports .gz
#
# USAGE (examples)
#   Name exact:
#     ./extract_kraken2_by_name_auto_fastq.sh "Pseudomonas aeruginosa" exact \
#       --kraken_dir results/tools/0_kraken2
#
#   Taxid (recommended; genus=286, species=287):
#     ./extract_kraken2_by_name_auto_fastq.sh _ clade --taxid 286 \
#       --kraken_dir results/tools/0_kraken2
#
# OUTPUT
#   results/tools/<N>_kraken2_extracted/<OrganismSafe>/<sample>_R1.fq(.gz)
#   results/tools/<N>_kraken2_extracted/<OrganismSafe>/<sample>_R2.fq(.gz)
#   results/summary_extracted-<OrganismSafe>.tsv
# --------------------------------------------------------------------

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <organism_name_or_->_unused_with_taxid> <exact|clade> [--kraken_dir DIR] [--taxid ID]"
  exit 1
fi

TARGET_NAME="$1"      # If --taxid used and this is "_", folder becomes taxid_<ID>
MATCH_MODE="$2"       # exact|clade
shift 2

KRAKEN_DIR=""
TARGET_TAXID=""
SEARCH_DEPTH=3   # recursive search depth for classified mates

while [[ $# -gt 0 ]]; do
  case "$1" in
    --kraken_dir)    KRAKEN_DIR="$2"; shift 2 ;;
    --taxid)         TARGET_TAXID="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

# Guess kraken results dir if not provided
if [[ -z "$KRAKEN_DIR" ]]; then
  KRAKEN_DIR=$(ls -d results/tools/*_kraken2 2>/dev/null | sort -V | tail -n1 || true)
fi
if [[ -z "$KRAKEN_DIR" || ! -d "$KRAKEN_DIR" ]]; then
  echo "ERROR: Kraken2 results directory not found. Provide with --kraken_dir."
  exit 1
fi

# Output dirs / labels
ORG_LABEL="$TARGET_NAME"
if [[ -n "$TARGET_TAXID" && "$TARGET_NAME" == "_" ]]; then
  ORG_LABEL="taxid_${TARGET_TAXID}"
fi
ORGSAFE=$(echo "$ORG_LABEL" | tr ' /' '__' | tr -cd '[:alnum:]_.-')

RUN_BASENAME=$(basename "$KRAKEN_DIR")
NUM_PREFIX="${RUN_BASENAME%%_kraken2}"
[[ "$NUM_PREFIX" == "$RUN_BASENAME" ]] && NUM_PREFIX="${RUN_BASENAME}"
OUTDIR="results/tools/${NUM_PREFIX}_kraken2_extracted/${ORGSAFE}"
mkdir -p "$OUTDIR"

SUMMARY="results/summary_extracted-${ORGSAFE}.tsv"
if [[ ! -f "$SUMMARY" ]]; then
  echo -e "group\tsamplename\torganismname\tsamplefile\ttotalnoofreads" > "$SUMMARY"
fi

echo ">> Using kraken dir: $KRAKEN_DIR"
if [[ -n "$TARGET_TAXID" ]]; then
  echo ">> Target taxid: $TARGET_TAXID ($MATCH_MODE)"
else
  echo ">> Target: $TARGET_NAME ($MATCH_MODE)"
fi
echo ">> Output: $OUTDIR"

shopt -s nullglob globstar

# ----------------- Helpers -----------------

reader() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then gzip -cd -- "$f"; else cat -- "$f"; fi
}

# Normalize scientific names consistently (lowercase, strip (), collapse WS)
norm_name_awk='
  function norm(s,   t){
    t=tolower(s)
    gsub(/\([^)]*\)/,"",t)
    gsub(/[[:space:]]+/," ",t)
    sub(/^ /,"",t); sub(/ $/,"",t)
    return t
  }
'

# Collect target taxids from report by name or by explicit taxid; for clade mode,
# include descendants based on indentation. (Trim leading WS before dropping 5 cols.)
collect_taxids_from_report() {
  local report="$1" target_name="$2" mode="$3" target_taxid="$4"
  if [[ -n "$target_taxid" ]]; then
    awk -v tgt="$target_taxid" -v mode="$mode" '
      BEGIN{found=0; target_indent=-1}
      {
        taxid=$5
        if(!found){
          if(taxid==tgt){
            name=$0
            gsub(/\r/,"",name); sub(/^[ \t]+/,"",name)
            for(i=1;i<=5;i++){ sub(/^[^ \t]+[ \t]+/,"",name) }
            indent=0; while (substr(name, indent+1, 1)==" ") indent++
            target_indent=indent; found=1; print taxid
            if(mode=="exact") exit
          }
        } else {
          name=$0
          gsub(/\r/,"",name); sub(/^[ \t]+/,"",name)
          for(i=1;i<=5;i++){ sub(/^[^ \t]+[ \t]+/,"",name) }
          indent=0; while (substr(name, indent+1, 1)==" ") indent++
          if(indent<=target_indent) exit
          print taxid
        }
      }' "$report"
    return
  fi

  awk -v target="$target_name" -v mode="$mode" "$norm_name_awk"'  # inject norm()
    BEGIN{ t=norm(target); found=0; target_indent=-1 }
    {
      taxid=$5
      name=$0
      gsub(/\r/,"",name); sub(/^[ \t]+/,"",name)
      for(i=1;i<=5;i++){ sub(/^[^ \t]+[ \t]+/,"",name) }
      indent=0; while (substr(name, indent+1, 1)==" ") indent++
      n=norm(name)

      if(!found){
        ok = (mode=="exact" ? (n==t) : (n==t || index(n, t" ")==1))
        if(ok){
          found=1; target_indent=indent; print taxid
          if(mode=="exact") exit
        }
      } else {
        if(indent<=target_indent) exit
        print taxid
      }
    }' "$report"
}

# Build set of normalized names for taxids present in TAXSET (for name-mode kraken).
build_names_from_taxset() {
  local report="$1" taxset_file="$2" names_out="$3"
  awk "$norm_name_awk"' NR==FNR{ keep[$1]=1; next }
    {
      taxid=$5
      if(taxid in keep){
        name=$0
        gsub(/\r/,"",name); sub(/^[ \t]+/,"",name)
        for(i=1;i<=5;i++){ sub(/^[^ \t]+[ \t]+/,"",name) }
        print norm(name)
      }
    }' "$taxset_file" "$report" | sort -u > "$names_out"
}

# Detect whether .kraken 3rd field is numeric taxid or a scientific name.
detect_kraken_mode() {
  local kraken="$1"
  awk '
    BEGIN{num=0; str=0}
    NR<=200 && $1=="C" {
      if($3 ~ /^[0-9]+$/) num++; else str++;
    }
    END{
      if(num>str) print "taxid"; else print "name";
    }' "$kraken"
}

# Extract read IDs from .kraken for any of the target TAXIDs (taxid-mode).
ids_from_kraken_by_taxset_taxidmode() {
  local kraken="$1" taxset_file="$2" ids_out="$3"
  awk 'NR==FNR{tax[$1]=1; next}
       $1=="C" && $3 ~ /^[0-9]+$/ {
         rid=$2; taxid=$3;
         if (tax[taxid]) { gsub(/\/[12]$/,"",rid); print rid }
       }' "$taxset_file" "$kraken" | sort -u > "$ids_out"
}

# Extract read IDs from .kraken for any of the target NAMES (name-mode).
# Rebuild scientific name up to ")"/or before first "123|123"; normalize.
ids_from_kraken_by_names_namemode() {
  local kraken="$1" names_file="$2" ids_out="$3"
  awk '
    function norm(s,   t){
      t=tolower(s)
      gsub(/\([^)]*\)/,"",t)
      gsub(/[[:space:]]+/," ",t)
      sub(/^ /,"",t); sub(/ $/,"",t)
      return t
    }
    NR==FNR { ok[norm($0)]=1; next }

    $1=="C" {
      rid=$2
      name=$3
      if (NF>=4) {
        for (i=4; i<=NF; i++) {
          if ($i ~ /\)$/) { name = name " " $i; break }
          if ($i ~ /^[0-9]+\|[0-9]+$/) { break }
          name = name " " $i
        }
      }
      n = norm(name)
      if (n in ok) {
        gsub(/\/[12]$/,"",rid)
        print rid
      }
    }
  ' "$names_file" "$kraken" | sort -u > "$ids_out"
}

# Stream-extract FASTQ blocks for matching IDs
extract_fastq_by_ids() {
  local ids="$1" fq_in="$2" fq_out="$3"
  reader "$fq_in" | \
  awk 'BEGIN{ while((getline line < "'$ids'")>0){keep[line]=1} }
       {
         if (NR%4==1){
           id=$0; sub(/^@/,"",id); sub(/[ \t].*/,"",id); gsub(/\/[12]$/,"",id)
           keeprec = (id in keep)
         }
         if (keeprec) print
       }' > "$fq_out"
}

# Find best matching classified mate file automatically.
find_best_match() {
  local base_dir="$1" sample="$2" mate="$3"
  local -a found=()
  local pats=(
    "$base_dir/**/$sample.classified_${mate}.fq"
    "$base_dir/**/$sample.classified_${mate}.fastq"
    "$base_dir/**/$sample.classified_${mate}.fq.gz"
    "$base_dir/**/$sample.classified_${mate}.fastq.gz"

    "$base_dir/**/$sample.classified__${mate}.fq"
    "$base_dir/**/$sample.classified__${mate}.fastq"
    "$base_dir/**/$sample.classified__${mate}.fq.gz"
    "$base_dir/**/$sample.classified__${mate}.fastq.gz"

    "$base_dir/**/$sample"*".*classified*${mate}*.fq"
    "$base_dir/**/$sample"*".*classified*${mate}*.fastq"
    "$base_dir/**/$sample"*".*classified*${mate}*.fq.gz"
    "$base_dir/**/$sample"*".*classified*${mate}*.fastq.gz"
  )
  while IFS= read -r -d '' x; do
    local rel="${x#${base_dir}/}"; [[ "$x" != "$rel" ]] || continue
    local depth=$(( $(grep -o "/" <<<"$rel" | wc -l) + 1 ))
    (( depth <= SEARCH_DEPTH )) && found+=("$x")
  done < <(
    for pat in "${pats[@]}"; do
      for x in $pat; do [[ -e "$x" ]] && printf '%s\0' "$x"; done
    done
  )
  if (( ${#found[@]} == 0 )); then echo ""; return 0; fi
  local best="" bestsz=-1
  for f in "${found[@]}"; do
    local sz=""
    if [[ "$f" == *.gz ]]; then
      sz=$(gzip -l "$f" 2>/dev/null | awk 'NR==2{print $2}')
      [[ -z "$sz" ]] && sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
    else
      sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
    fi
    [[ -z "$sz" ]] && sz=0
    if (( sz > bestsz )); then best="$f"; bestsz=$sz; fi
  done
  echo "$best"
}

# ----------------- Main -----------------

for KRAKEN in "$KRAKEN_DIR"/*.kraken; do
  [[ -e "$KRAKEN" ]] || continue
  SAMPLE=$(basename "$KRAKEN" .kraken)
  REPORT="$KRAKEN_DIR/$SAMPLE.report"
  if [[ ! -f "$REPORT" ]]; then
    echo "WARN: report missing for $SAMPLE, skipping."
    continue
  fi

  echo ">> Sample: $SAMPLE"

  # Auto-find classified mates
  R1=$(find_best_match "$KRAKEN_DIR" "$SAMPLE" 1)
  R2=$(find_best_match "$KRAKEN_DIR" "$SAMPLE" 2)
  if [[ -z "$R1" || -z "$R2" ]]; then
    echo "   WARN: classified FASTQ not found for $SAMPLE, skipping."
    continue
  fi

  # Taxid set (by explicit taxid or by name)
  TAXTMP="$(mktemp)"
  collect_taxids_from_report "$REPORT" "$TARGET_NAME" "$MATCH_MODE" "$TARGET_TAXID" > "$TAXTMP"
  if [[ ! -s "$TAXTMP" ]]; then
    echo "   No taxids found for target in report; skipping $SAMPLE."
    rm -f "$TAXTMP"
    continue
  fi

  # Detect kraken mode
  MODE=$(detect_kraken_mode "$KRAKEN")

  IDSTMP="$(mktemp)"
  if [[ "$MODE" == "taxid" ]]; then
    ids_from_kraken_by_taxset_taxidmode "$KRAKEN" "$TAXTMP" "$IDSTMP"
  else
    NAMETMP="$(mktemp)"
    build_names_from_taxset "$REPORT" "$TAXTMP" "$NAMETMP"
    ids_from_kraken_by_names_namemode "$KRAKEN" "$NAMETMP" "$IDSTMP"
    rm -f "$NAMETMP"
  fi

  if [[ ! -s "$IDSTMP" ]]; then
    echo "   No reads for target in .kraken; skipping $SAMPLE."
    rm -f "$TAXTMP" "$IDSTMP"
    continue
  fi

  OUT_R1="$OUTDIR/${SAMPLE}_R1.fq"
  OUT_R2="$OUTDIR/${SAMPLE}_R2.fq"
  [[ "$R1" == *.gz ]] && OUT_R1="${OUT_R1}.gz"
  [[ "$R2" == *.gz ]] && OUT_R2="${OUT_R2}.gz"

  # Extract R1
  if [[ "$OUT_R1" == *.gz ]]; then
    TMP_R1="$(mktemp)"; extract_fastq_by_ids "$IDSTMP" "$R1" "$TMP_R1"; gzip -c "$TMP_R1" > "$OUT_R1"; rm -f "$TMP_R1"
  else
    extract_fastq_by_ids "$IDSTMP" "$R1" "$OUT_R1"
  fi
  # Extract R2
  if [[ "$OUT_R2" == *.gz ]]; then
    TMP_R2="$(mktemp)"; extract_fastq_by_ids "$IDSTMP" "$R2" "$TMP_R2"; gzip -c "$TMP_R2" > "$OUT_R2"; rm -f "$TMP_R2"
  else
    extract_fastq_by_ids "$IDSTMP" "$R2" "$OUT_R2"
  fi

  # Count reads (R1)
  if [[ "$OUT_R1" == *.gz ]]; then
    COUNT=$(gzip -cd "$OUT_R1" | awk 'END{print NR/4}')
  else
    COUNT=$(awk 'END{print NR/4}' "$OUT_R1")
  fi

  echo -e "${NUM_PREFIX}\t${SAMPLE}\t${ORG_LABEL}\t${SAMPLE}\t${COUNT}" >> "$SUMMARY"
  echo "   Extracted ${COUNT} read(s) for R1 (R2 similar)."
  echo "     R1: $OUT_R1"
  echo "     R2: $OUT_R2"

  rm -f "$TAXTMP" "$IDSTMP"
done

echo "Done. Summary: $SUMMARY"
