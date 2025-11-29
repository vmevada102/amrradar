#!/usr/bin/env bash
#
# One_Health-AMR-Radar-pangenome.sh
# Pangenome runner (Panaroo or Roary) with resume, core alignment, iqtree (classic),
# Rtab output, and optional parallel postprocessing.
#
set -euo pipefail
IFS=$'\n\t'

log(){ echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }
status(){ echo "[$(date +'%Y-%m-%d %H:%M:%S')] STATUS: $*"; }

usage(){
  cat <<EOF
Usage: $0 [options]
Options:
  -r DIR              results root (default: results)
  --pan-tool TOOL     panaroo (default) | roary
  -t INT              threads (default: 4, but dynamically set to available CPUs - 4)
  --no-resume         force re-run even if results exist
  --parallel-post     enable GNU parallel for postprocessing (auto if installed)
  --no-parallel-post  disable parallel postprocessing
  -h                  help
EOF
}

# ---------------- defaults ----------------
RESULTS_ROOT="results"
PAN_TOOL="panaroo"
THREADS=4 # Will be dynamically set later if not passed via -t
RESUME=true
PARALLEL_POST="auto"

# ---------------- parse args ----------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r) RESULTS_ROOT="$2"; shift 2;;
    --pan-tool) PAN_TOOL="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    --no-resume) RESUME=false; shift 1;;
    --parallel-post) PARALLEL_POST="yes"; shift 1;;
    --no-parallel-post) PARALLEL_POST="no"; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

# ---------------- directories and validation ----------------
TOOLS_ROOT="$RESULTS_ROOT/tools"

# --- Dynamic THREADS/CPUs Calculation (Total CPUs - 4) ---
USER_THREAD_OVERRIDE=false
if [[ "$THREADS" -ne 4 ]]; then # Check if the user overrode the default 4
    USER_THREAD_OVERRIDE=true
fi

if [[ "$USER_THREAD_OVERRIDE" == false ]]; then
  # Try to detect available CPUs
  if command -v nproc >/dev/null 2>&1; then
    TOTAL_CPUS=$(nproc)
  elif [[ -f /proc/cpuinfo ]]; then
    TOTAL_CPUS=$(grep -c ^processor /proc/cpuinfo)
  else
    TOTAL_CPUS=0
  fi
  
  # Calculate THREADS: Total available - 4
  if (( TOTAL_CPUS > 4 )); then
    THREADS=$((TOTAL_CPUS - 4))
    log "Dynamically set THREADS to $THREADS (Total CPUs: $TOTAL_CPUS, reserving 4)."
  elif (( TOTAL_CPUS > 0 )); then
    # If less than 5 CPUs, use 1 thread
    THREADS=1
    log "Dynamically set THREADS to 1 (Total CPUs: $TOTAL_CPUS, insufficient to reserve 4)."
  else
    # Fallback to the default 4 if detection failed
    THREADS=4 
    log "Warning: Could not reliably detect Total CPUs. THREADS remains $THREADS."
  fi
fi

# 2. Find Prokka Directory (ends with *prokka)
PROK_DIR=$(find "$TOOLS_ROOT" -maxdepth 1 -type d -name "*prokka" -print -quit || true)

if [[ -z "$PROK_DIR" ]]; then
  echo "ERROR: Could not find any directory matching '*prokka' under $TOOLS_ROOT" >&2
  exit 1
fi
log "Using Prokka Input Directory: $(basename "$PROK_DIR")"

# 3. Determine the next consecutive number (N) for the output directory
# Find the highest number prefix (e.g., '6' from '6_prokka') and increment it.
# Default to 1 if no numbered steps are found.
LAST_STEP_NUM=0
mapfile -t STEP_DIRS < <(find "$TOOLS_ROOT" -mindepth 1 -maxdepth 1 -type d -name "[0-9]*_*" -print 2>/dev/null | sort)

for dir in "${STEP_DIRS[@]}"; do
  # Extract the number from the directory name (e.g., '6' from '6_prokka')
  num_prefix=$(basename "$dir" | cut -d'_' -f1)
  if [[ "$num_prefix" =~ ^[0-9]+$ ]]; then
    if (( num_prefix > LAST_STEP_NUM )); then
      LAST_STEP_NUM=$num_prefix
    fi
  fi
done

NEXT_STEP_NUM=$((LAST_STEP_NUM + 1))

# 4. Define the new Pangenome Output Directory
PANG_DIR="$TOOLS_ROOT/${NEXT_STEP_NUM}_pangenome"
log "Creating Pangenome Output Directory: $(basename "$PANG_DIR")"
# --- END NEW LOGIC ---

mkdir -p "$PANG_DIR"

# helper to list missing bins
need_bin(){
  local miss=()
  for prog in "$@"; do
    if ! command -v "$prog" >/dev/null 2>&1; then
      miss+=("$prog")
    fi
  done
  printf "%s " "${miss[@]:-}"
}

# choose required programs and install hints
if [[ "$PAN_TOOL" == "panaroo" ]]; then
  REQUIRED=(panaroo mafft iqtree python3)
  INSTALL_CMD="mamba install -y -c bioconda panaroo mafft python && conda install -y -c bioconda/label/cf201901 iqtree"
elif [[ "$PAN_TOOL" == "roary" ]]; then
  REQUIRED=(roary cd-hit mafft iqtree python3)
  INSTALL_CMD="mamba install -y -c bioconda roary cd-hit mafft python && conda install -y -c bioconda/label/cf201901 iqtree"
else
  echo "Invalid --pan-tool option: $PAN_TOOL" >&2
  exit 1
fi

MISSING=$(need_bin "${REQUIRED[@]}")
if [[ -n "${MISSING// }" ]]; then
  echo "ERROR: Missing required tools: $MISSING" >&2
  echo "Install with:"
  echo "  $INSTALL_CMD"
  exit 1
fi

# optional tools
FASTTREE_AVAILABLE=true
if ! command -v fasttree >/dev/null 2>&1; then
  FASTTREE_AVAILABLE=false
fi

PARALLEL_AVAILABLE=false
if command -v parallel >/dev/null 2>&1; then
  PARALLEL_AVAILABLE=true
fi

# decide parallel-post
if [[ "$PARALLEL_POST" == "auto" ]]; then
  if [[ "$PARALLEL_AVAILABLE" == true ]]; then
    PARALLEL_POST="yes"
  else
    PARALLEL_POST="no"
  fi
fi

log "Pan tool: $PAN_TOOL, threads: $THREADS, resume: $RESUME, parallel-post: $PARALLEL_POST"

# collect GFFs into array
mapfile -t GFF_PATHS < <(find "$PROK_DIR" -mindepth 2 -maxdepth 2 -type f -name "*.gff" -print | sort)
NUM_SAMPLES=${#GFF_PATHS[@]}
if (( NUM_SAMPLES == 0 )); then
  log "No GFF files found in $PROK_DIR. Run Prokka first."
  exit 0
fi
log "Found $NUM_SAMPLES GFF files (example: ${GFF_PATHS[0]})"

# resume check & run pangenome tool
DONE="$PANG_DIR/.done"
GPA="$PANG_DIR/gene_presence_absence.csv"

if [[ "$RESUME" == true && -f "$DONE" && -f "$GPA" ]]; then
  log "Previous run detected (resume enabled). Skipping pangenome step."
else
  if [[ "$RESUME" == false ]]; then
    log "--no-resume -> removing previous outputs."
    rm -rf "$PANG_DIR"/*
    mkdir -p "$PANG_DIR"
  fi

  if [[ "$PAN_TOOL" == "panaroo" ]]; then
    status "Running Panaroo with $THREADS threads..."
    # pass array elements as separate args
    panaroo -i "${GFF_PATHS[@]}" -o "$PANG_DIR" -t "$THREADS" --clean-mode moderate --aligner mafft
  else
    status "Running Roary with $THREADS threads..."
    roary -p "$THREADS" -f "$PANG_DIR" -e -n "${GFF_PATHS[@]}"
  fi

  touch "$DONE"
  log "Pangenome analysis finished."
fi

# detect outputs
CORE_ALIGN="$PANG_DIR/core_gene_alignment.aln"
if [[ ! -f "$CORE_ALIGN" ]]; then
  if [[ -f "$PANG_DIR/core_gene_alignment.fa" ]]; then CORE_ALIGN="$PANG_DIR/core_gene_alignment.fa"; fi
  if [[ -f "$PANG_DIR/core_gene_alignment.filtered.aln" ]]; then CORE_ALIGN="$PANG_DIR/core_gene_alignment.filtered.aln"; fi
fi

GPA="$PANG_DIR/gene_presence_absence.csv"
RTAB="$PANG_DIR/gene_presence_absence.Rtab"
IQTREE_PREFIX="$PANG_DIR/core_gene_alignment"
TREEFILE="$PANG_DIR/core_gene_alignment.treefile"

# define postprocessing functions
run_iqtree(){
  if [[ -z "$CORE_ALIGN" || ! -f "$CORE_ALIGN" ]]; then
    log "IQ-TREE: core alignment not found; skipping."
    return 0
  fi
  log "IQ-TREE: running on $CORE_ALIGN with $THREADS threads..."
  iqtree -s "$CORE_ALIGN" -m GTR+G -nt "$THREADS" -bb 1000 -pre "$IQTREE_PREFIX" || log "IQ-TREE failed (check logs)."
}

run_fasttree(){
  if [[ "$FASTTREE_AVAILABLE" == false ]]; then
    log "FastTree not installed; skipping."
    return 0
  fi
  if [[ -z "$CORE_ALIGN" || ! -f "$CORE_ALIGN" ]]; then
    log "FastTree: no core alignment; skipping."
    return 0
  fi
  log "FastTree: building quick tree..."
  fasttree -nt "$CORE_ALIGN" > "$PANG_DIR/core_gene_alignment.fasttree.nwk" 2> "$PANG_DIR/fasttree.log" || log "FastTree failed."
  log "FastTree: wrote core_gene_alignment.fasttree.nwk"
}

run_csv2rtab(){
  if [[ ! -f "$GPA" ]]; then
    log "CSV->Rtab: gene_presence_absence.csv missing; skipping."
    return 0
  fi
  log "CSV->Rtab: converting $GPA -> $RTAB"
  python3 - "$GPA" "$RTAB" "$NUM_SAMPLES" <<'PY'
import csv,sys
gpa=sys.argv[1]
rtab=sys.argv[2]
N=int(sys.argv[3])
with open(gpa,newline='',encoding='utf-8') as fh:
    reader=csv.reader(fh)
    header=next(reader)
    sample_cols=header[-N:]
    sample_start=len(header)-N
    with open(rtab,'w',encoding='utf-8') as out:
        out.write('gene\t' + '\t'.join(sample_cols) + '\n')
        for row in reader:
            gene=row[0]
            pres = ['1' if v.strip() else '0' for v in row[sample_start:sample_start+N]]
            out.write(gene + '\t' + '\t'.join(pres) + '\n')
print("Wrote", rtab)
PY
}

run_summary(){
  if [[ ! -f "$GPA" ]]; then
    log "Summary: gene_presence_absence.csv missing; skipping."
    return 0
  fi
  SUMMARY="$PANG_DIR/pangenome_summary.tsv"
  log "Summary: writing $SUMMARY"
  python3 - "$GPA" "$NUM_SAMPLES" "$SUMMARY" <<'PY'
import csv,sys,os
gpa=sys.argv[1]
N=int(sys.argv[2])
outp=sys.argv[3]
total=0
core=0
with open(gpa,newline='',encoding='utf-8') as fh:
    r=csv.reader(fh)
    header=next(r)
    sample_start=len(header)-N
    for row in r:
        total+=1
        pres=sum(1 for v in row[sample_start:sample_start+N] if v.strip())
        if pres==N:
            core+=1
with open(outp,'w',encoding='utf-8') as o:
    o.write(f"gene_presence_absence\t{gpa}\n")
    o.write(f"total_gene_clusters\t{total}\n")
    o.write(f"number_of_samples\t{N}\n")
    o.write(f"core_gene_count\t{core}\n")
    o.write(f"accessory_gene_count\t{total-core}\n")
    core_align = os.path.abspath(gpa).replace('.csv','_alignment?') if os.path.exists(gpa) else 'MISSING'
    o.write(f"core_alignment\t{core_align}\n")
    o.write(f"iqtree_tree\t{os.path.join(os.path.abspath(os.path.dirname(outp)),'core_gene_alignment.treefile')}\n")
    o.write(f"Rtab\t{os.path.join(os.path.abspath(os.path.dirname(outp)),'gene_presence_absence.Rtab')}\n")
print("Wrote summary to", outp)
PY
}

# run postprocessing in parallel or serial
log "Postprocessing tasks: CSV->RTAB, summary, FastTree (quick), IQ-TREE (full)."

if [[ "$PARALLEL_POST" == "yes" && "$PARALLEL_AVAILABLE" == true ]]; then
  log "GNU parallel detected: running postprocessing tasks in parallel."
  export -f run_iqtree run_fasttree run_csv2rtab run_summary log status
  # run tasks in parallel
  parallel ::: "run_csv2rtab" "run_summary" "run_fasttree" "run_iqtree"
else
  log "Running postprocessing serially."
  run_csv2rtab
  run_summary
  run_fasttree
  run_iqtree
fi

log "Pangenome step finished. Outputs in $PANG_DIR"
exit 0