#!/usr/bin/env bash
set -euo pipefail

# run_kma_all.sh
# Usage: run_kma_all.sh -s sample.tsv -p 8 -t 4 [--resume]
# sample.tsv format: sample<TAB>R1.fastq.gz<TAB>R2.fastq.gz
# or provide -i input_dir with fastq files named SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz

usage() {
  cat <<'USAGE'
Usage: run_kma_all.sh -s sample.tsv -p <parallel_jobs> -t <threads_per_job> [--resume]
Requires: conda env with kma (or kma installed in PATH)
Outputs to results/tools/N+1_resistome/
USAGE
  exit 1
}

SAMPLE_TSV=""
PAR_JOBS=4
THREADS=4
RESUME=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s) SAMPLE_TSV="$2"; shift 2;;
    -p) PAR_JOBS="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    --resume) RESUME=true; shift ;;
    -h|--help) usage ;;
    *) echo "Unknown arg $1"; usage ;;
  esac
done

if [[ -z "$SAMPLE_TSV" ]]; then echo "Missing -s sample.tsv"; usage; fi

# check tools
for tool in kma bwa samtools awk sort; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: required tool '$tool' not found in PATH. Install via conda (e.g. conda install -c bioconda kma bwa samtools)"; exit 2
  fi
done

# find next output dir: results/tools/N+1_resistome
mkdir -p results/tools
Nmax=0
for d in results/tools/*_resistome results/tools/*; do
  [[ -d "$d" ]] || continue
  # find any directory matching number_resistome; extract number if present at start
  base=$(basename "$d")
  if [[ "$base" =~ ^([0-9]+)_resistome$ ]]; then
    n=${BASH_REMATCH[1]}
    (( n > Nmax )) && Nmax=$n
  fi
done
OUTDIR="results/tools/$((Nmax+1))_resistome"
mkdir -p "$OUTDIR"
echo "Outputs to: $OUTDIR"

# require CARD fasta (nucleotide) present in conda abricate DB or user-specified env
CARD_FA="/data/home/amr/miniconda3/envs/resistome1/db/card/sequences"
if [[ ! -f "$CARD_FA" ]]; then
  echo "ERROR: CARD fasta not found at $CARD_FA. Provide CARD fasta or set CARD_FA in script."
  exit 3
fi

# build KMA index (only once)
KMA_DB="$OUTDIR/card_kma_db"
if [[ ! -f "${KMA_DB}.index" && ! -f "${KMA_DB}.seq" ]]; then
  echo "Building KMA DB in $KMA_DB (this may take a while)"
  kma_index -i "$CARD_FA" -o "$KMA_DB" || { echo "kma_index failed"; exit 4; }
fi

# iterate samples in parallel
cmds=()
while IFS=$'\t' read -r sample r1 r2; do
  [[ -z "$sample" || "$sample" =~ ^# ]] && continue
  s_out="$OUTDIR/$sample"
  mkdir -p "$s_out"
  # skip if resume and output exists
  if $RESUME && [[ -f "$s_out/${sample}_kma.res" ]]; then
    echo "Skipping $sample (exists)"; continue
  fi
  # generate command
  cmds+=( "kma -ipe $r1 $r2 -o $s_out/${sample}_kma -t_db $KMA_DB -t ${THREADS} &> $s_out/${sample}_kma.log || true; \
           awk 'NR>1{print}' $s_out/${sample}_kma.res > $s_out/${sample}_kma_hits.tsv || true" )
done < "$SAMPLE_TSV"

# run commands in parallel
echo "Running ${#cmds[@]} jobs with $PAR_JOBS parallel jobs"
printf "%s\n" "${cmds[@]}" | xargs -P "$PAR_JOBS" -I CMD bash -c CMD

# Consolidate hits -> project-level matrix (presence/reads)
echo -e "sample\tgene\tcoverage\tdepth" > "$OUTDIR/amr_kma_hits.tsv"
for d in "$OUTDIR"/*/; do
  sname=$(basename "$d")
  if [[ -f "$d/${sname}_kma_hits.tsv" ]]; then
    awk -v S="$sname" 'NR>1{print S"\t"$0}' "$d/${sname}_kma_hits.tsv" >> "$OUTDIR/amr_kma_hits.tsv"
  fi
done

echo "KMA mapping complete. Per-sample results in $OUTDIR/* ; combined: $OUTDIR/amr_kma_hits.tsv"
