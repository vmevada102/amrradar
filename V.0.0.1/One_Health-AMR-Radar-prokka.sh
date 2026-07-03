#!/usr/bin/env bash
#
# One_Health-AMR-Radar-prokka.sh
#
# Combined Prokka runner + simple summary generator.
# - Runs Prokka on assemblies under results/tools/3_assembly/<sample>/
# - Saves outputs to results/tools/6_prokka/<sample>/
# - Generates summary results/annotation_prokka.xls (tab-separated)
#
set -euo pipefail
IFS=$'\n\t'

# ---------------- helpers ----------------
log(){ echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }
status(){ echo "[$(date +'%Y-%m-%d %H:%M:%S')] STATUS: $*"; }

usage(){
  cat <<EOF
Usage: One_Health-AMR-Radar-prokka.sh [options]

Options:
  -r DIR            results root (default: results)
  -t INT            number of parallel Prokka jobs (default: 2)
  -c INT            cpus per Prokka job (default: 4)
  --sample-file F   sample.tsv with columns: sample<TAB>r1<TAB>r2<TAB>group (optional)
  --clean-tmp       remove /tmp/prokka_* leftover dirs before running (use with care)
  --no-resume       do NOT resume; re-run all samples (default: resume enabled)
  -h                show help
EOF
}

# ---------------- defaults ----------------
RESULTS_ROOT="results"
PARALLEL_JOBS=10  # Changed from 2 to 10
PROKKA_CPUS=4     # Will be dynamically set later if not passed via -c
SAMPLE_FILE=""
CLEAN_TMP=false
RESUME=true

# ---------------- parse args ----------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r) RESULTS_ROOT="$2"; shift 2;;
    -t) PARALLEL_JOBS="$2"; shift 2;;
    -c) PROKKA_CPUS="$2"; shift 2;;
    --sample-file) SAMPLE_FILE="$2"; shift 2;;
    --clean-tmp) CLEAN_TMP=true; shift 1;;
    --no-resume) RESUME=false; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

# ---------------- validate environment ----------------
if [[ ! -d "$RESULTS_ROOT" ]]; then
  echo "ERROR: results root not found: $RESULTS_ROOT" >&2
  exit 1
fi

# --- NEW LOGIC: Set PROKKA_CPUS dynamically ---
# Only calculate if the user didn't override the default (4) via the -c flag
if [[ "$PROKKA_CPUS" -eq 4 ]]; then
  # Try to detect available CPUs
  if command -v nproc >/dev/null 2>&1; then
    TOTAL_CPUS=$(nproc)
  elif [[ -f /proc/cpuinfo ]]; then
    TOTAL_CPUS=$(grep -c ^processor /proc/cpuinfo)
  else
    TOTAL_CPUS=0
  fi
  
  # Calculate PROKKA_CPUS: Total available - 4
  if (( TOTAL_CPUS > 4 )); then
    PROKKA_CPUS=$((TOTAL_CPUS - 4))
    log "Dynamically set PROKKA_CPUS to $PROKKA_CPUS (Total CPUs: $TOTAL_CPUS, reserving 4)."
  elif (( TOTAL_CPUS > 0 )); then
    # If less than 5 CPUs, use 1 CPU
    PROKKA_CPUS=1
    log "Dynamically set PROKKA_CPUS to 1 (Total CPUs: $TOTAL_CPUS, insufficient to reserve 4)."
  else
    # Fallback to the default 4 if detection failed
    PROKKA_CPUS=4 
    log "Warning: Could not reliably detect Total CPUs. PROKKA_CPUS remains $PROKKA_CPUS."
  fi
fi
# --- END NEW LOGIC ---






TOOLS_ROOT="$RESULTS_ROOT/tools"

# --- MODIFIED: Dynamic directory finding ---

# 1. Find the Assembly Directory (ends with *assembly)
ASSEMBLY_DIR=$(find "$TOOLS_ROOT" -maxdepth 1 -type d -name "*_assembly" -print -quit || true)

if [[ -z "$ASSEMBLY_DIR" ]]; then
  echo "ERROR: Could not find any directory matching '*_assembly' under $TOOLS_ROOT" >&2
  exit 1
fi

# 2. Determine the next consecutive number (N) for the output directory
# Find the highest number prefix (e.g., '3' from '3_assembly') and increment it.
# Default to 1 if no numbered steps are found.
LAST_STEP_NUM=0
mapfile -t STEP_DIRS < <(find "$TOOLS_ROOT" -mindepth 1 -maxdepth 1 -type d -name "[0-9]*_*" -print 2>/dev/null | sort)

for dir in "${STEP_DIRS[@]}"; do
  # Extract the number from the directory name (e.g., '3' from '3_assembly')
  num_prefix=$(basename "$dir" | cut -d'_' -f1)
  if [[ "$num_prefix" =~ ^[0-9]+$ ]]; then
    if (( num_prefix > LAST_STEP_NUM )); then
      LAST_STEP_NUM=$num_prefix
    fi
  fi
done

NEXT_STEP_NUM=$((LAST_STEP_NUM + 1))

# 3. Define the new Output Directory
OUT_STEP_DIR="$TOOLS_ROOT/${NEXT_STEP_NUM}_prokka"
SUMMARY_OUT="$RESULTS_ROOT/annotation_prokka.xls"

log "Using Assembly Directory: $(basename "$ASSEMBLY_DIR")"
log "Creating Prokka Output Directory: $(basename "$OUT_STEP_DIR")"
# --- END MODIFIED ---

mkdir -p "$OUT_STEP_DIR"


# check prokka
if ! command -v prokka >/dev/null 2>&1; then
  echo "ERROR: prokka not found in PATH. Install via conda/mamba: mamba install -y -c conda-forge -c bioconda prokka" >&2
  exit 1
fi

# detect GNU parallel
USE_PARALLEL=false
if command -v parallel >/dev/null 2>&1; then
  USE_PARALLEL=true
  log "GNU parallel found -> will use it for $PARALLEL_JOBS parallel jobs."
else
  log "GNU parallel not found -> falling back to bash job semaphore with $PARALLEL_JOBS jobs."
fi

# optional clean tmp
if [[ "$CLEAN_TMP" == true ]]; then
  log "Cleaning /tmp/prokka_* directories (pattern) - make sure you are ok with this."
  ls -ld /tmp/prokka_* 2>/dev/null || true
  rm -rf /tmp/prokka_* 2>/dev/null || true
fi

# ---------------- collect assembly jobs ----------------
mapfile -t SAMPLE_DIRS < <(find "$ASSEMBLY_DIR" -mindepth 1 -maxdepth 1 -type d -print 2>/dev/null | sort)
if [[ ${#SAMPLE_DIRS[@]} -eq 0 ]]; then
  log "No sample directories under $ASSEMBLY_DIR. Nothing to do."
  exit 0
fi

declare -a JOBS
for sdir in "${SAMPLE_DIRS[@]}"; do
  sample=$(basename "$sdir")
  asm=""
  asm=$(find "$sdir" -maxdepth 1 -type f -name "*.assembly.fasta" -print -quit || true)
  if [[ -z "$asm" ]]; then
    asm=$(find "$sdir" -maxdepth 1 -type f -name "*.contigs.fasta" -print -quit || true)
  fi
  if [[ -z "$asm" ]]; then
    asm=$(find "$sdir" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" \) -print -quit || true)
  fi
  if [[ -n "$asm" ]]; then
    JOBS+=("$sample"$'\t'"$asm")
  else
    log "No assembly found in $sdir -> skipping"
  fi
done

TOTAL=${#JOBS[@]}
if (( TOTAL == 0 )); then
  log "No assemblies located. Exiting."
  exit 0
fi
log "Found $TOTAL assemblies to annotate."

# ---------------- sample file -> group mapping ----------------
declare -A SAMPLE_TO_GROUP
if [[ -n "$SAMPLE_FILE" ]]; then
  if [[ -f "$SAMPLE_FILE" ]]; then
    while IFS=$'\t' read -r s r1 r2 g || [[ -n "$s" ]]; do
      [[ "${s,,}" == "sample" ]] && continue
      [[ -z "$s" ]] && continue
      SAMPLE_TO_GROUP["$s"]="$g"
    done < "$SAMPLE_FILE"
    log "Loaded ${#SAMPLE_TO_GROUP[@]} group mappings from $SAMPLE_FILE"
  else
    log "Warning: sample file not found: $SAMPLE_FILE (Group column will be empty)"
  fi
else
  log "No sample file provided; Group column will be empty."
fi

# ---------------- prokka worker ----------------
run_prokka_job(){
  IFS=$'\t' read -r SAMPLE ASM <<< "$1"
  OUT_SAMPLE_DIR="$OUT_STEP_DIR/$SAMPLE"
  GFF_FILE="$OUT_SAMPLE_DIR/${SAMPLE}.gff"

  # resume logic
  if [[ "$RESUME" == "true" && -f "$OUT_SAMPLE_DIR/.done" && -f "$GFF_FILE" ]]; then
    status "SKIP prokka $SAMPLE (already .done)"
    return 0
  fi

  status "RUNNING prokka $SAMPLE (asm: $(basename "$ASM"))"

  # create unique TMP path that does NOT exist so Prokka will create it
  TMP_PARENT="/tmp"
  TMP_UNIQUE="prokka_${SAMPLE}_$(date +%s)_$RANDOM"
  TMPDIR="$TMP_PARENT/$TMP_UNIQUE"

  TMP_STDOUT="$TMP_PARENT/${TMP_UNIQUE}.prokka.stdout"
  TMP_STDERR="$TMP_PARENT/${TMP_UNIQUE}.prokka.stderr"

  # run prokka; Prokka creates TMPDIR
  if prokka --outdir "$TMPDIR" --prefix "$SAMPLE" --cpus "$PROKKA_CPUS" --kingdom Bacteria --addgenes --compliant "$ASM" > "$TMP_STDOUT" 2> "$TMP_STDERR"; then
    # atomically move into final location (or copy)
    rm -rf "$OUT_SAMPLE_DIR" 2>/dev/null || true
    if mv "$TMPDIR" "$OUT_SAMPLE_DIR" 2>/dev/null; then
      :
    else
      mkdir -p "$OUT_SAMPLE_DIR"
      cp -a "$TMPDIR"/* "$OUT_SAMPLE_DIR"/ || true
      rm -rf "$TMPDIR"
    fi
    [[ -f "$TMP_STDOUT" ]] && mv "$TMP_STDOUT" "$OUT_SAMPLE_DIR/prokka.stdout" || true
    [[ -f "$TMP_STDERR" ]] && mv "$TMP_STDERR" "$OUT_SAMPLE_DIR/prokka.stderr" || true

    # check for GFF and mark done
    if [[ -f "$OUT_SAMPLE_DIR/${SAMPLE}.gff" ]]; then
      touch "$OUT_SAMPLE_DIR/.done"
      status "DONE prokka $SAMPLE"
      return 0
    else
      status "FAILED prokka $SAMPLE (GFF missing after prokka)"
      return 2
    fi
  else
    mkdir -p "$OUT_SAMPLE_DIR"
    [[ -f "$TMP_STDOUT" ]] && mv "$TMP_STDOUT" "$OUT_SAMPLE_DIR/prokka.stdout" || true
    [[ -f "$TMP_STDERR" ]] && mv "$TMP_STDERR" "$OUT_SAMPLE_DIR/prokka.stderr" || true
    if [[ -d "$TMPDIR" ]]; then
      cp -a "$TMPDIR"/* "$OUT_SAMPLE_DIR"/ 2>/dev/null || true
      rm -rf "$TMPDIR"
    fi
    status "FAILED prokka $SAMPLE (see $OUT_SAMPLE_DIR/prokka.stderr)"
    return 1
  fi
}

export -f run_prokka_job
export -f log
export -f status
export PROKKA_CPUS
export OUT_STEP_DIR
export RESUME

# ---------------- run prokka jobs in parallel ----------------
if [[ "$USE_PARALLEL" == true ]]; then
  TMP_JOBFILE=$(mktemp /tmp/onehealth_prokka_jobs_XXXX)
  printf "%s\n" "${JOBS[@]}" > "$TMP_JOBFILE"
  log "Launching Prokka via GNU parallel with $PARALLEL_JOBS concurrent jobs..."
  parallel --jobs "$PARALLEL_JOBS" --halt now,fail=1 --joblog "$OUT_STEP_DIR/parallel.joblog" run_prokka_job :::: "$TMP_JOBFILE"
  rc=$?
  rm -f "$TMP_JOBFILE"
  if [[ $rc -ne 0 ]]; then
    log "One or more Prokka jobs returned non-zero (parallel exit code $rc). Check per-sample logs under $OUT_STEP_DIR"
  fi
else
  log "Launching Prokka with bash job control ($PARALLEL_JOBS concurrent jobs)..."
  running=0
  idx=0
  declare -a pids
  for job in "${JOBS[@]}"; do
    ((idx++))
    run_prokka_job "$job" &
    pids[${idx}]=$!
    ((running++))
    if (( running >= PARALLEL_JOBS )); then
      wait -n || true
      # recompute running jobs count
      running=0
      for pid in "${pids[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
          ((running++))
        fi
      done
    fi
  done
  wait
fi

# ---------------- build simple summary from per-sample .txt ----------------
log "Building summary $SUMMARY_OUT from per-sample .txt files (first .txt in each sample folder)."

# header
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "Group" "samplename" "file_used" "organism" "contigs" "bases" "CDS" "gene" "rRNA" "repeat_region" "tRNA" > "$SUMMARY_OUT"

for sample_dir in "$OUT_STEP_DIR"/*; do
  [[ -d "$sample_dir" ]] || continue
  sample=$(basename "$sample_dir")
  asm_path="$(find "$ASSEMBLY_DIR/$sample" -maxdepth 1 -type f \( -name "*.assembly.fasta" -o -name "*.contigs.fasta" -o -name "*.fasta" -o -name "*.fa" \) -print -quit || true)"

  txtfile="$(find "$sample_dir" -maxdepth 1 -type f -name "*.txt" -print -quit || true)"

  organism=""; contigs=""; bases=""; cds=""; gene=""; rRNA=""; repeat_region=""; trna=""; tmrna=""
  if [[ -z "$txtfile" ]]; then
    # no txt found
    status="MISSING_TXT"
  else
    # parse key: value lines robustly
    while IFS= read -r line; do
      [[ -z "$line" ]] && continue
      line_trim="$(echo "$line" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
      if echo "$line_trim" | grep -q ":"; then
        key=$(echo "$line_trim" | awk -F: '{print $1}' | tr '[:upper:]' '[:lower:]' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        val=$(echo "$line_trim" | awk -F: '{ $1=""; sub(/^[:space:]*/, ""); print }' )
        val="$(echo "$val" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
        case "$key" in
          organism) organism="$val" ;;
          contigs) contigs="$val" ;;
          bases) bases="$val" ;;
          cds) cds="$val" ;;
          gene) gene="$val" ;;
          rrna|r_rna) rRNA="$val" ;;
          repeat_region|repeat) repeat_region="$val" ;;
          trna|tRNA) trna="$val" ;;
          tmrna|tmRNA) tmrna="$val" ;;
          *) ;; 
        esac
      fi
    done < "$txtfile"
    status="OK"
  fi

  group="${SAMPLE_TO_GROUP[$sample]:-}"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$group" "$sample" "$asm_path" "$organism" "$contigs" "$bases" "$cds" "$gene" "$rRNA" "$repeat_region" "$trna" >> "$SUMMARY_OUT"
done

log "Wrote summary: $SUMMARY_OUT"
log "Done."
