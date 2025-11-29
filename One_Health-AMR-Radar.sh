#!/usr/bin/env bash
#
# One_Health-AMR-Radar.sh
# Parallel-capable stepwise pipeline WITHOUT Kraken2.
#
set -euo pipefail
IFS=$'\n\t'

# ----------------------
# Logging helper (defined early)
# ----------------------
log(){ echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }

# Defaults
SAMPLE_TSV="sample.tsv"
OUTDIR="results"
THREADS=8
PARALLEL=1
ASSEMBLER="spades"
TRIMMER="fastp"
RESUME="false"
MEM_OVERRIDE=""

print_help(){
  cat <<'HELP'
One Health-AMR-Radar - pipeline (parallel-capable) WITHOUT Kraken2

Usage:
  One_Health-AMR-Radar.sh -i sample.tsv -o results -t 8 -p 4

Options:
  -i FILE         sample TSV (tab-delimited) header: sample r1 r2 group
  -o DIR          output directory (default: results)
  -t INT          threads per tool (default: 8)
  -p INT          number of samples to process in parallel (default: 1)
  -a NAME         assembler: spades or shovill (default: spades)
  --trimmer NAME  trimmer: fastp or trimmomatic (default: fastp)
  --mem MB        override memory in MB (for assemblies)
  --resume        resume from previous run (skip steps with .done markers)
  -h              show help
Note: Kraken2 support has been removed from this script.
HELP
}

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) SAMPLE_TSV="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    -p|--parallel) PARALLEL="$2"; shift 2;;
    -a) ASSEMBLER="$2"; shift 2;;
    --trimmer) TRIMMER="$2"; shift 2;;
    --mem) MEM_OVERRIDE="$2"; shift 2;;
    --resume) RESUME="true"; shift 1;;
    -h) print_help; exit 0;;
    *) echo "Unknown option: $1"; print_help; exit 1;;
  esac
done

# sanity checks
re_is_int='^[0-9]+$'
if ! [[ "$PARALLEL" =~ $re_is_int && "$PARALLEL" -ge 1 ]]; then
  echo "ERROR: -p/--parallel must be integer >=1" >&2; exit 1
fi
if ! [[ "$THREADS" =~ $re_is_int && "$THREADS" -ge 1 ]]; then
  echo "ERROR: -t/--threads must be integer >=1" >&2; exit 1
fi
if [[ ! -f "$SAMPLE_TSV" ]]; then
  echo "ERROR: sample TSV not found: $SAMPLE_TSV" >&2; exit 1
fi

mkdir -p "$OUTDIR"

# Copy input sample file into results folder for traceability (timestamped)
if [[ -f "$SAMPLE_TSV" ]]; then
  cp -f "$SAMPLE_TSV" "$OUTDIR/sample_used_$(date +%Y%m%d_%H%M%S).tsv"
  log "Copied sample file to results folder for traceability"
fi

get_mem_mb(){
  if [[ -r /proc/meminfo ]]; then
    mem_total_kb=$(grep -i '^MemTotal:' /proc/meminfo | awk '{print $2}')
    echo $((mem_total_kb / 1024))
  else
    echo "4096"
  fi
}

if [[ -z "$MEM_OVERRIDE" ]]; then
  TOTAL_MEM_MB=$(get_mem_mb)
  USE_MEM_MB=$(( TOTAL_MEM_MB * 90 / 100 ))
else
  USE_MEM_MB=$MEM_OVERRIDE
fi

STEP_ROOT="$OUTDIR/tools"
STEP1="$STEP_ROOT/1_fastqc"
STEP2="$STEP_ROOT/2_trim"
STEP3="$STEP_ROOT/3_assembly"
STEP4="$STEP_ROOT/4_stats"
STEP5="$STEP_ROOT/5_checkm"

mkdir -p "$STEP1" "$STEP2" "$STEP3" "$STEP4" "$STEP5"

export THREADS PARALLEL USE_MEM_MB ASSEMBLER TRIMMER OUTDIR SAMPLE_TSV RESUME \
       STEP1 STEP2 STEP3 STEP4 STEP5

log "Starting One Health-AMR-Radar (Kraken2 removed)"
log "Samples file: $SAMPLE_TSV"
log "Outdir: $OUTDIR"
log "Threads per tool: $THREADS"
log "Parallel samples: $PARALLEL"
log "Memory (MB for assemblies): $USE_MEM_MB"
log "Assembler: $ASSEMBLER"
log "Trimmer: $TRIMMER"
log "Resume: $RESUME"
echo "---------------------------------------------------"

# Create helper parser script (parse_checkm_dict.py) under tools if not present
PARSER="$STEP_ROOT/parse_checkm_dict.py"
if [[ ! -f "$PARSER" ]]; then
  log "Creating helper parser: $PARSER"
  mkdir -p "$(dirname "$PARSER")"
  cat > "$PARSER" <<'PY'
#!/usr/bin/env python3
"""
parse_checkm_dict.py

Parses CheckM's bin_stats_ext.tsv when the second column is a Python-style dict.
Writes a TSV: bin\tcompleteness\tcontamination (with header).
Usage: parse_checkm_dict.py path/to/bin_stats_ext.tsv > parsed_checkm.tsv
"""
import sys, ast, io, re

if len(sys.argv) != 2:
    print("Usage: parse_checkm_dict.py bin_stats_ext.tsv", file=sys.stderr)
    sys.exit(2)

infile = sys.argv[1]
out = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
print("bin\tcompleteness\tcontamination", file=out)

with open(infile, 'r', encoding='utf-8') as fh:
    for raw in fh:
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t", 1)
        if len(parts) == 1:
            continue
        bin_name, dict_str = parts[0].strip(), parts[1].strip()
        comp = ""
        cont = ""
        try:
            d = ast.literal_eval(dict_str)
            if isinstance(d, dict):
                comp = d.get('Completeness') if 'Completeness' in d else d.get('completeness')
                cont = d.get('Contamination') if 'Contamination' in d else d.get('contamination')
        except Exception:
            # parsing failed - fall back to regex
            pass
        # regex fallback if values still missing
        if comp is None or comp == "":
            m = re.search(r"['\"]?Completeness['\"]?\s*:\s*([0-9]+(?:\.[0-9]+)?)", dict_str, re.IGNORECASE)
            comp = float(m.group(1)) if m else ""
        if cont is None or cont == "":
            m = re.search(r"['\"]?Contamination['\"]?\s*:\s*([0-9]+(?:\.[0-9]+)?)", dict_str, re.IGNORECASE)
            cont = float(m.group(1)) if m else ""
        # ensure numeric types print plainly
        print(f"{bin_name}\t{comp}\t{cont}", file=out)
PY
  chmod +x "$PARSER"
fi

# map header columns
header=$(head -n1 "$SAMPLE_TSV" | tr -d '\r')
IFS=$'\t' read -r -a cols <<< "$header"
declare -A colidx
for i in "${!cols[@]}"; do
  name="${cols[$i]}"
  colidx["$name"]=$((i+1))
done
required=(sample r1 r2 group)
for r in "${required[@]}"; do
  if [[ -z "${colidx[$r]:-}" ]]; then
    echo "ERROR: Required column '$r' missing from header" >&2; exit 1
  fi
done

mapfile -t SAMPLE_LINES < <(tail -n +2 "$SAMPLE_TSV")

wait_for_slot(){
  while true; do
    running=$(jobs -rp | wc -l)
    if [[ "$running" -lt "$PARALLEL" ]]; then break; fi
    sleep 1
  done
}

# compute & write read counts (prefers seqkit)
compute_and_write_counts(){
  local sample="$1"; local r1="$2"; local r2="$3"; local trim_dir="$4"; local trim_r1="$5"; local trim_r2="$6"
  local out_file="$trim_dir/read_counts.tsv"
  printf "raw_R1\traw_R2\ttrimmed_R1\ttrimmed_R2\n" > "$out_file"
  if command -v seqkit >/dev/null 2>&1; then
    raw1=0; raw2=0; trim1=0; trim2=0
    if [[ -f "$r1" ]]; then raw1=$(seqkit stats -a -T -j "$THREADS" "$r1" 2>/dev/null | awk -F'\t' '{print $4+0}'); fi
    if [[ -n "$r2" && -f "$r2" ]]; then raw2=$(seqkit stats -a -T -j "$THREADS" "$r2" 2>/dev/null | awk -F'\t' '{print $4+0}'); fi
    if [[ -f "$trim_r1" ]]; then trim1=$(seqkit stats -a -T -j "$THREADS" "$trim_r1" 2>/dev/null | awk -F'\t' '{print $4+0}'); fi
    if [[ -n "$r2" && -f "$trim_r2" ]]; then trim2=$(seqkit stats -a -T -j "$THREADS" "$trim_r2" 2>/dev/null | awk -F'\t' '{print $4+0}'); fi
    printf "%s\t%s\t%s\t%s\n" "${raw1}" "${raw2}" "${trim1}" "${trim2}" >> "$out_file"
    return 0
  fi
  if [[ -f "$trim_dir/fastp_report.json" ]]; then
    python3 - <<PY > "$out_file.tmp" 2>/dev/null
import json,os
p = os.path.join("$trim_dir","fastp_report.json")
try:
    j = json.load(open(p))
except Exception:
    print("0\t0\t0\t0"); raise SystemExit(0)
b = j.get("summary",{}).get("before_filtering",{})
a = j.get("summary",{}).get("after_filtering",{})
raw = int(b.get("total_reads",0))
aft = int(a.get("total_reads",0))
if "$r2":
    print(f"{raw}\t{raw}\t{aft}\t{aft}")
else:
    print(f"{raw}\t0\t{aft}\t0")
PY
    if [[ -f "$out_file.tmp" ]]; then mv "$out_file.tmp" "$out_file"; return 0; fi
  fi
  # fallback wc -l via gzip
  raw1=0; raw2=0; trim1=0; trim2=0
  if [[ -f "$r1" ]]; then raw1=$(( $(gzip -cd "$r1" 2>/dev/null | wc -l) / 4 )); fi
  if [[ -n "$r2" && -f "$r2" ]]; then raw2=$(( $(gzip -cd "$r2" 2>/dev/null | wc -l) / 4 )); fi
  if [[ -f "$trim_r1" ]]; then trim1=$(( $(gzip -cd "$trim_r1" 2>/dev/null | wc -l) / 4 )); fi
  if [[ -n "$r2" && -f "$trim_r2" ]]; then trim2=$(( $(gzip -cd "$trim_r2" 2>/dev/null | wc -l) / 4 )); fi
  printf "%s\t%s\t%s\t%s\n" "${raw1}" "${raw2}" "${trim1}" "${trim2}" >> "$out_file"
  return 0
}

# main per-sample processing
run_sample(){
  local line="$1"
  IFS=$'\t' read -r -a fields <<< "$line"
  local SAMPLEID="${fields[${colidx[sample]}-1]}"
  local R1="${fields[${colidx[r1]}-1]}"
  local R2="${fields[${colidx[r2]}-1]}"
  if [[ -z "$R2" || "$R2" == "NA" || "$R2" == "na" ]]; then R2=""; fi

  echo "[`date +'%Y-%m-%d %H:%M:%S'`] Tool: FastQC (sample: $SAMPLEID)"
  local S1_DIR="$STEP1/$SAMPLEID"; mkdir -p "$S1_DIR"
  if [[ "$RESUME" == "true" && -f "$S1_DIR/.done" ]]; then :; else
    if [[ -n "$R2" ]]; then
      fastqc -t "$THREADS" -o "$S1_DIR" "$R1" "$R2" > "$S1_DIR/fastqc.stdout" 2> "$S1_DIR/fastqc.stderr" || true
    else
      fastqc -t "$THREADS" -o "$S1_DIR" "$R1" > "$S1_DIR/fastqc.stdout" 2> "$S1_DIR/fastqc.stderr" || true
    fi
    touch "$S1_DIR/.done"
  fi

  echo "[`date +'%Y-%m-%d %H:%M:%S'`] Tool: Trimming ($TRIMMER) (sample: $SAMPLEID)"
  local S2_DIR="$STEP2/$SAMPLEID"; mkdir -p "$S2_DIR"
  local TRIM_R1="$S2_DIR/${SAMPLEID}_R1.trim.fastq.gz"
  local TRIM_R2="$S2_DIR/${SAMPLEID}_R2.trim.fastq.gz"
  if [[ "$RESUME" == "true" && -f "$S2_DIR/.done" ]]; then :; else
    if [[ "$TRIMMER" == "fastp" ]]; then
      if [[ -n "$R2" ]]; then
        fastp -i "$R1" -I "$R2" -o "$TRIM_R1" -O "$TRIM_R2" -w "$THREADS" --detect_adapter_for_pe \
          --html "$S2_DIR/fastp_report.html" --json "$S2_DIR/fastp_report.json" > "$S2_DIR/fastp.stdout" 2> "$S2_DIR/fastp.stderr" || true
      else
        fastp -i "$R1" -o "$TRIM_R1" -w "$THREADS" --html "$S2_DIR/fastp_report.html" --json "$S2_DIR/fastp_report.json" > "$S2_DIR/fastp.stdout" 2> "$S2_DIR/fastp.stderr" || true
      fi
    else
      if [[ -z "${TRIMMOMATIC_ADAPTERS:-}" ]]; then
        echo "ERROR: TRIMMOMATIC_ADAPTERS not set (sample: $SAMPLEID)" >&2; return 1
      fi
      if [[ -n "$R2" ]]; then
        trimmomatic PE -threads "$THREADS" "$R1" "$R2" \
          "$TRIM_R1" "$S2_DIR/${SAMPLEID}_R1.unpaired.fastq.gz" \
          "$TRIM_R2" "$S2_DIR/${SAMPLEID}_R2.unpaired.fastq.gz" \
          ILLUMINACLIP:${TRIMMOMATIC_ADAPTERS}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 > "$S2_DIR/trimmomatic.stdout" 2> "$S2_DIR/trimmomatic.stderr" || true
      else
        trimmomatic SE -threads "$THREADS" "$R1" "$TRIM_R1" SLIDINGWINDOW:4:20 MINLEN:50 > "$S2_DIR/trimmomatic.stdout" 2> "$S2_DIR/trimmomatic.stderr" || true
      fi
    fi

    compute_and_write_counts "$SAMPLEID" "$R1" "$R2" "$S2_DIR" "$TRIM_R1" "$TRIM_R2"
    touch "$S2_DIR/.done"
  fi

  echo "[`date +'%Y-%m-%d %H:%M:%S'`] Tool: Assembly ($ASSEMBLER) (sample: $SAMPLEID)"
  local S3_DIR="$STEP3/$SAMPLEID"; mkdir -p "$S3_DIR"
  local ASSEMBLY_FASTA="$S3_DIR/${SAMPLEID}.assembly.fasta"
  if [[ "$RESUME" == "true" && -f "$S3_DIR/.done" ]]; then :; else
    if [[ "$ASSEMBLER" == "spades" ]]; then
      MEM_GB=$(( USE_MEM_MB / 1024 )); [[ $MEM_GB -lt 1 ]] && MEM_GB=1
      if [[ -n "$R2" ]]; then
        spades.py -1 "$TRIM_R1" -2 "$TRIM_R2" -o "$S3_DIR/spades_out" --threads "$THREADS" -m "$MEM_GB" > "$S3_DIR/spades.stdout" 2> "$S3_DIR/spades.stderr" || true
      else
        spades.py --s1 "$TRIM_R1" -o "$S3_DIR/spades_out" --threads "$THREADS" -m "$MEM_GB" > "$S3_DIR/spades.stdout" 2> "$S3_DIR/spades.stderr" || true
      fi
      [[ -f "$S3_DIR/spades_out/contigs.fasta" ]] && cp "$S3_DIR/spades_out/contigs.fasta" "$ASSEMBLY_FASTA"
    elif [[ "$ASSEMBLER" == "shovill" ]]; then
      if [[ -n "$R2" ]]; then
        shovill --R1 "$TRIM_R1" --R2 "$TRIM_R2" --outdir "$S3_DIR/shovill_out" --threads "$THREADS" --ram "$USE_MEM_MB" > "$S3_DIR/shovill.stdout" 2> "$S3_DIR/shovill.stderr" || true
      else
        shovill --R1 "$TRIM_R1" --outdir "$S3_DIR/shovill_out" --threads "$THREADS" --ram "$USE_MEM_MB" > "$S3_DIR/shovill.stdout" 2> "$S3_DIR/shovill.stderr" || true
      fi
      [[ -f "$S3_DIR/shovill_out/contigs.fa" ]] && cp "$S3_DIR/shovill_out/contigs.fa" "$ASSEMBLY_FASTA"
    fi
    touch "$S3_DIR/.done"
  fi

  echo "[`date +'%Y-%m-%d %H:%M:%S'`] Tool: Stats (sample: $SAMPLEID)"
  local S4_DIR="$STEP4/$SAMPLEID"; mkdir -p "$S4_DIR"
  if [[ "$RESUME" == "true" && -f "$S4_DIR/.done" ]]; then :; else
    export ASSEMBLY_FASTA="$ASSEMBLY_FASTA"
    python3 - <<'PY' > "$S4_DIR/assembly_stats.tsv"
import os
from pathlib import Path
f = Path(os.environ.get("ASSEMBLY_FASTA",""))
if not f.exists():
    print("TotalBP\tContigs\tN50\tL50\tGC\tLargest"); raise SystemExit(0)
seqs=[]
for line in f.read_text().splitlines():
    if line.startswith(">"): seqs.append("")
    else: seqs[-1]+=line.strip()
lens=[len(s) for s in seqs if s]
if not lens:
    print("TotalBP\tContigs\tN50\tL50\tGC\tLargest"); raise SystemExit(0)
lens.sort(reverse=True)
T=sum(lens); cum=0; n50=l50=0
for i,l in enumerate(lens,1):
    cum+=l
    if cum*2>=T and not n50:
        n50,l50=l,i
GC=sum(s.count('G')+s.count('C')+s.count('g')+s.count('c') for s in seqs)/T*100
print("TotalBP\tContigs\tN50\tL50\tGC\tLargest")
print(f"{T}\t{len(lens)}\t{n50}\t{l50}\t{GC:.3f}\t{lens[0]}")
PY
    touch "$S4_DIR/.done"
  fi

  echo "[`date +'%Y-%m-%d %H:%M:%S'`] Tool: CheckM (sample: $SAMPLEID)"
  local S5_DIR="$STEP5/$SAMPLEID"; mkdir -p "$S5_DIR"
  if command -v checkm >/dev/null 2>&1; then
    if [[ "$RESUME" == "true" && -f "$S5_DIR/.done" ]]; then :; else
      if [[ -f "$ASSEMBLY_FASTA" ]]; then
        mkdir -p "$S5_DIR/input" "$S5_DIR/out"
        cp "$ASSEMBLY_FASTA" "$S5_DIR/input/contigs.fa" 2>/dev/null || true
        checkm lineage_wf -x fa "$S5_DIR/input" "$S5_DIR/out" -t "$THREADS" --pplacer_threads "$THREADS" > "$S5_DIR/checkm.stdout" 2> "$S5_DIR/checkm.stderr" || true
        if [[ -f "$S5_DIR/out/storage/bin_stats_ext.tsv" ]]; then
          # robust parsing via Python helper
          "$PARSER" "$S5_DIR/out/storage/bin_stats_ext.tsv" > "$S5_DIR/parsed_checkm.tsv" || true
        fi
      fi
      touch "$S5_DIR/.done"
    fi
  else
    echo "[`date +'%Y-%m-%d %H:%M:%S'`] Tool: CheckM (not available, skipped for $SAMPLEID)"
  fi

  log "Completed sample: $SAMPLEID"
}

# Launch sample jobs with concurrency control
for line in "${SAMPLE_LINES[@]}"; do
  [[ -z "${line//[[:space:]]/}" ]] && continue
  wait_for_slot
  run_sample "$line" &
done

wait

# Combine results into Excel summary
log "Combining results into Excel summary"
if [[ -f "generate_report.py" ]]; then
  python3 generate_report.py -i "$SAMPLE_TSV" -o "$OUTDIR/OneHealthAMR_summary.xlsx" -r "$OUTDIR" || log "generate_report.py returned non-zero"
  log "Summary written to $OUTDIR/OneHealthAMR_summary.xlsx"
else
  log "generate_report.py not found - skipping Excel generation"
fi

log "Pipeline finished at: $(date)"
