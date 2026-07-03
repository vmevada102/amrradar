#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Kraken2 runner — READS or ASSEMBLY mode, auto-parallel across samples
# ==============================================================================

# -------------------- CPU & defaults --------------------
detect_cpus() {
  if command -v nproc >/dev/null 2>&1; then nproc
  elif command -v getconf >/dev/null 2>&1; then getconf _NPROCESSORS_ONLN || echo 1
  else echo 1; fi
}
CPU_TOTAL="$(detect_cpus)"; CPU_TOTAL=$(( CPU_TOTAL>0 ? CPU_TOTAL : 1 ))

# "auto" defaults ? spread CPUs fairly across parallel samples
THREADS="auto"            # per-sample threads
PARALLEL="auto"           # number of samples to run concurrently

KRAKEN_DB="${KRAKEN2_DB:-}"
MIN_PERCENT=0.1
PCT_OF="total"            # {total|classified}
RESUME=0
MODE="reads"              # {reads|assembly}
READS_DIR="raw"
ASSEMBLY_DIR=""
ASSEMBLY_GLOB='results/tools/*_assembly/*/*.assembly.fasta'

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SAMPLE_TSV="$SCRIPT_DIR/sample.tsv"

TOOLS_ROOT="results/tools"
OUTDIR="$TOOLS_ROOT/0_kraken2"

IDENTIFIED_SUM="results/Identified_organisms.tsv"
IDENTIFIED_LONG="results/Identified_organisms_long.tsv"
IDENTIFIED_GENUS="results/Identified_organisms_genus.tsv"
GROUP_SUM="results/Identified_organisms_group_summary.tsv"

usage() {
  cat <<'USAGE'
Usage:
  run_kraken2_identify_summary.sh --kraken_db <DB_DIR> [options]

Required:
  --kraken_db PATH             Kraken2 DB directory (or set KRAKEN2_DB)

Modes & inputs:
  --mode {reads|assembly}      Default: reads
    reads:    sample.tsv (sample  r1  r2  group) + raw/ FASTQs
    assembly: per-sample assembly FASTA (auto-discover); r1/r2 ignored

Parallelism (auto is recommended):
  --threads N|auto             Threads PER SAMPLE (default: auto)
  --parallel N|auto            Max samples in parallel (default: auto = min(#samples, CPUs))

Other options:
  --min_percent P              Taxa percent threshold (default: 0.1)
  --percent_of MODE            {total|classified} (default: total)
  --sample_tsv PATH            Path to sample.tsv (default: alongside script)
  --reads_dir DIR              Where to look for reads (default: raw)
  --assembly_dir DIR           Optional root to search assemblies
  --assembly_glob GLOB         Glob for assemblies (default: results/tools/*_assembly/*/*.assembly.fasta)
  -r, --resume                 Skip samples already completed
  -h, --help                   This help

Outputs:
  Per-sample in results/tools/0_kraken2/
  Summaries in results/: Identified_organisms.tsv / _long.tsv / _genus.tsv / _group_summary.tsv
USAGE
}

# -------------------- Parse args --------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --kraken_db) KRAKEN_DB="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --parallel) PARALLEL="$2"; shift 2;;
    --min_percent) MIN_PERCENT="$2"; shift 2;;
    --percent_of) PCT_OF="$2"; shift 2;;
    --mode) MODE="$2"; shift 2;;
    --sample_tsv) SAMPLE_TSV="$2"; shift 2;;
    --reads_dir) READS_DIR="$2"; shift 2;;
    --assembly_dir) ASSEMBLY_DIR="$2"; shift 2;;
    --assembly_glob) ASSEMBLY_GLOB="$2"; shift 2;;
    -r|--resume) RESUME=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

case "$PCT_OF" in total|classified) :;; *) echo "ERROR: --percent_of must be total|classified"; exit 1;; esac
case "$MODE" in reads|assembly) :;; *) echo "ERROR: --mode must be reads|assembly"; exit 1;; esac

# -------------------- Checks --------------------
MISSING=()
command -v kraken2 >/dev/null 2>&1 || MISSING+=("kraken2 not in PATH")
[[ -z "$KRAKEN_DB" || ! -d "$KRAKEN_DB" ]] && MISSING+=("kraken2 DB (--kraken_db) not found")
[[ ! -f "$SAMPLE_TSV" ]] && MISSING+=("sample.tsv not found at $SAMPLE_TSV")
if [[ "$MODE" == "reads" && ! -d "$READS_DIR" ]]; then MISSING+=("--reads_dir '$READS_DIR' does not exist"); fi
if [[ ${#MISSING[@]} -gt 0 ]]; then
  echo "Missing requirements:"; for i in "${MISSING[@]}"; do echo " - $i"; done; exit 2
fi

# -------------------- Load samples --------------------
declare -a samples
declare -A s_r1 s_r2 s_grp
while IFS=$'\t' read -r sample r1 r2 grp rest; do
  [[ -z "${sample:-}" ]] && continue
  [[ "$sample" =~ ^[Ss]ample$ ]] && continue
  samples+=("$sample")
  s_r1["$sample"]="${r1:-}"
  s_r2["$sample"]="${r2:-}"
  s_grp["$sample"]="${grp:-NA}"
done < "$SAMPLE_TSV"
NS=${#samples[@]}
(( NS == 0 )) && { echo "No samples in $SAMPLE_TSV"; exit 1; }

# -------------------- Decide parallel plan (auto) --------------------
is_int() { [[ "$1" =~ ^[0-9]+$ ]]; }

if [[ "$PARALLEL" == "auto" ]]; then
  # run all samples in parallel but never exceed CPU count
  if (( NS <= CPU_TOTAL )); then PARALLEL="$NS"; else PARALLEL="$CPU_TOTAL"; fi
elif ! is_int "$PARALLEL" || (( PARALLEL < 1 )); then
  echo "ERROR: --parallel must be an integer or 'auto'"; exit 1
fi

if [[ "$THREADS" == "auto" ]]; then
  # distribute cores fairly across parallel samples; keep at least 1/thread
  # keep one safety core if many samples: floor((CPU_TOTAL)/PARALLEL)
  calc=$(( CPU_TOTAL / PARALLEL ))
  (( calc < 1 )) && calc=1
  THREADS="$calc"
elif ! is_int "$THREADS" || (( THREADS < 1 )); then
  echo "ERROR: --threads must be an integer or 'auto'"; exit 1
fi

echo "CPUs detected: ${CPU_TOTAL}"
echo "Samples found : ${NS}"
echo "Mode: $MODE | Threads per sample: $THREADS | Parallel samples: $PARALLEL"
echo "Kraken2 DB: $KRAKEN_DB"
echo "sample.tsv : $SAMPLE_TSV"
[[ "$MODE" == "reads"    ]] && echo "Reads directory: $READS_DIR"
[[ "$MODE" == "assembly" ]] && { echo "Assembly dir: ${ASSEMBLY_DIR:-<auto via glob>}"; echo "Assembly glob: $ASSEMBLY_GLOB"; }
echo "Min %: $MIN_PERCENT | Percent-of: $PCT_OF | Resume: $RESUME"
echo

# -------------------- Prep --------------------
mkdir -p "$OUTDIR"
TMPDIR="$OUTDIR/tmp"
mkdir -p "$TMPDIR/sum" "$TMPDIR/long" "$TMPDIR/genus"
mkdir -p "$(dirname "$IDENTIFIED_SUM")"

# -------------------- Helpers --------------------
sanitize_or_empty() {
  local p="$1"
  [[ -z "$p" ]] && { echo ""; return 0; }
  [[ "$p" == *"..."* || "$p" == *"…"* ]] && { echo ""; return 0; }
  [[ -f "$p" ]] && { echo "$p"; return 0; }
  # if basename exists in READS_DIR, use it
  local b; b="$(basename "$p")"
  [[ -f "$READS_DIR/$b" ]] && { echo "$READS_DIR/$b"; return 0; }
  echo ""
}

resolve_read() {
  local p="$1"
  p="$(sanitize_or_empty "$p")"
  echo "$p"
}

fastq_metrics() {
  local fq1="$1" fq2="$2"; local total=0 bases=0
  _proc() { local f="$1"; { [[ "$f" =~ \.gz$ ]] && gzip -cd "$f" || cat "$f"; } \
    | awk 'NR%4==2{c++; n+=length($0)} END{printf "%d\t%d\n", (c+0),(n+0)}'; }
  [[ -f "$fq1" ]] && { read -r c1 n1 < <(_proc "$fq1"); total=$((total+c1)); bases=$((bases+n1)); }
  [[ -f "$fq2" ]] && { read -r c2 n2 < <(_proc "$fq2"); total=$((total+c2)); bases=$((bases+n2)); }
  echo -e "${total}\t${bases}"
}

find_assembly_for_sample() {
  local sample="$1"
  if [[ -n "$ASSEMBLY_DIR" && -d "$ASSEMBLY_DIR" ]]; then
    shopt -s nullglob
    local hits=( "$ASSEMBLY_DIR"/**/*"$sample"*".fasta" "$ASSEMBLY_DIR"/**/*"$sample"*".fa" )
    shopt -u nullglob
    [[ ${#hits[@]} -gt 0 ]] && { echo "${hits[0]}"; return 0; }
  fi
  shopt -s nullglob
  local ghits=( $ASSEMBLY_GLOB )
  shopt -u nullglob
  for f in "${ghits[@]}"; do [[ "$(basename "$f")" == *"$sample"* ]] && { echo "$f"; return 0; }; done
  echo ""
}

fasta_metrics() {
  local fa="$1"
  awk '
    BEGIN{FS=""; seqlen=0; total=0; gc=0; nseq=0}
    /^>/{
      if(seqlen>0){lens[++nseq]=seqlen; total+=seqlen}
      seqlen=0; next
    }
    { for(i=1;i<=NF;i++){ ch=toupper($i); if(ch!="\n" && ch!="\r"){ if(ch=="G"||ch=="C")gc++; if(ch!=">")seqlen++ } } }
    END{
      if(seqlen>0){lens[++nseq]=seqlen; total+=seqlen}
      for(i=1;i<=nseq;i++) for(j=i+1;j<=nseq;j++) if(lens[j]>lens[i]){t=lens[i];lens[i]=lens[j];lens[j]=t}
      target50=total*0.5; target90=total*0.9; run=0; n50=0;l50=0;n90=0;l90=0
      for(i=1;i<=nseq;i++){ run+=lens[i]; if(n50==0 && run>=target50){n50=lens[i]; l50=i} if(n90==0 && run>=target90){n90=lens[i]; l90=i} }
      gpct=(total==0?0:gc/total*100.0)
      printf "%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", nseq,total,(n50?n50:0),(l50?l50:0),(n90?n90:0),(l90?l90:0),gpct
    }' "$fa"
}

# -------------------- Worker --------------------
process_one_sample() {
  local sample="$1"
  echo "Processing sample: $sample"
  local grp="${s_grp[$sample]:-NA}"

  local out_report="$OUTDIR/${sample}.report"
  local out_kraken="$OUTDIR/${sample}.kraken"
  local done_marker="$OUTDIR/${sample}.done"
  local logf="$OUTDIR/${sample}.kraken2.log"

  local sum_out="$TMPDIR/sum/${sample}.tsv"
  local long_out="$TMPDIR/long/${sample}.tsv"
  local genus_out="$TMPDIR/genus/${sample}.tsv"

  local classified_count=0 unclassified_count=0 total_seq=0 total_bases=0
  local N50="NA" L50="NA" N90="NA" L90="NA" GCpc="NA"

  if (( RESUME == 1 )) && [[ -s "$out_report" && -s "$out_kraken" && -f "$done_marker" ]]; then
    echo "  (resume) already done."
  else
    if [[ "$MODE" == "reads" ]]; then
      local r1_in r2_in
      r1_in="$(resolve_read "${s_r1[$sample]}")"
      r2_in="$(resolve_read "${s_r2[$sample]}")"

      # If unresolved (empty/bad paths), auto-discover by sample id
      if [[ -z "$r1_in" || -z "$r2_in" ]]; then
        shopt -s nullglob
        local cand_r1=( "$READS_DIR"/*"${sample}"*_R1*.fastq.gz "$READS_DIR"/*"${sample}"*_R1*.fastq )
        local cand_r2=( "$READS_DIR"/*"${sample}"*_R2*.fastq.gz "$READS_DIR"/*"${sample}"*_R2*.fastq )
        shopt -u nullglob
        [[ -z "$r1_in" ]] && r1_in="${cand_r1[0]:-}"
        [[ -z "$r2_in" ]] && r2_in="${cand_r2[0]:-}"
      fi

      if [[ -z "$r1_in" || -z "$r2_in" || ! -f "$r1_in" || ! -f "$r2_in" ]]; then
        echo "  WARNING: reads not found for $sample (looked in $READS_DIR)."
        echo -e "${grp}\t${sample}\t0\t0\tNA\tNA\tNA\tNA\tNA\tunclassified\t0\t0\t0.00\t" > "$sum_out"
        : > "$long_out"; : > "$genus_out"
        return 0
      fi

      read -r total_reads total_bases < <(fastq_metrics "$r1_in" "$r2_in")
      local class_out="$OUTDIR/${sample}.classified_#.fq"
      local unclass_out="$OUTDIR/${sample}.unclassified_#.fq"

      echo "  Running kraken2 on reads..."
      if ! kraken2 --db "$KRAKEN_DB" \
          --threads "$THREADS" \
          --use-names \
          --paired \
          --report "$out_report" \
          --output "$out_kraken" \
          --classified-out "$class_out" \
          --unclassified-out "$unclass_out" \
          "$r1_in" "$r2_in" >"$logf" 2>&1; then
        echo "  ERROR: kraken2 failed for $sample (see $logf)."
        echo -e "${grp}\t${sample}\t${total_reads}\t${total_bases}\tNA\tNA\tNA\tNA\tNA\tunclassified\t0\t${total_reads}\t0.00\t" > "$sum_out"
        : > "$long_out"; : > "$genus_out"
        touch "$done_marker"
        return 0
      fi

      classified_count=$(grep -c '^C' "$out_kraken" || true)
      unclassified_count=$(grep -c '^U' "$out_kraken" || true)
      total_seq=$(( classified_count + unclassified_count ))
      (( total_seq == 0 )) && total_seq=$total_reads

    else
      local asm
      asm="$(find_assembly_for_sample "$sample")"
      if [[ -z "$asm" || ! -f "$asm" ]]; then
        echo "  WARNING: no assembly found for $sample."
        echo -e "${grp}\t${sample}\t0\t0\tNA\tNA\tNA\tNA\tNA\tunclassified\t0\t0\t0.00\t" > "$sum_out"
        : > "$long_out"; : > "$genus_out"
        return 0
      fi

      read -r contigs total_bases N50 L50 N90 L90 GCpc < <(fasta_metrics "$asm")
      total_seq="$contigs"

      local class_out="$OUTDIR/${sample}.classified.fa"
      local unclass_out="$OUTDIR/${sample}.unclassified.fa"

      echo "  Running kraken2 on assembly..."
      if ! kraken2 --db "$KRAKEN_DB" \
          --threads "$THREADS" \
          --use-names \
          --report "$out_report" \
          --output "$out_kraken" \
          --classified-out "$class_out" \
          --unclassified-out "$unclass_out" \
          "$asm" >"$logf" 2>&1; then
        echo "  ERROR: kraken2 failed for $sample (see $logf)."
      fi

      classified_count=$(grep -c '^C' "$out_kraken" 2>/dev/null || echo 0)
      unclassified_count=$(grep -c '^U' "$out_kraken" 2>/dev/null || echo 0)
    fi
  fi

  local denom_total=$(( classified_count + unclassified_count ))
  local denom_for_pct
  if [[ "$PCT_OF" == "total" ]]; then denom_for_pct="$denom_total"; else denom_for_pct="$classified_count"; fi
  (( denom_for_pct < 1 )) && denom_for_pct=1

  # Parse report ? species/genus tables
  local sp_tmp="$OUTDIR/${sample}.tax_sp.tsv"
  local ge_tmp="$OUTDIR/${sample}.tax_ge.tsv"
  : > "$sp_tmp"; : > "$ge_tmp"
  if [[ -s "$out_report" ]]; then
    awk -v minp="$MIN_PERCENT" -F'\t' '($4=="S" && $2>0 && $1+0>=minp){n=$6; gsub(/^[ \t]+|[ \t]+$/,"",n); printf "%s\t%s\t%s\tS\n",$1,$2,n}' "$out_report" > "$sp_tmp" || true
    awk -v minp="$MIN_PERCENT" -F'\t' '($4=="G" && $2>0 && $1+0>=minp){n=$6; gsub(/^[ \t]+|[ \t]+$/,"",n); printf "%s\t%s\t%s\tG\n",$1,$2,n}' "$out_report" > "$ge_tmp" || true
  fi

  declare -A org_clade org_rank
  local have_species=0
  if [[ -s "$sp_tmp" ]]; then
    have_species=1
    while IFS=$'\t' read -r pct clade name rank; do
      org_clade["$name"]=$(( ${org_clade["$name"]:-0} + ${clade%.*} )); org_rank["$name"]="S"
    done < "$sp_tmp"
  fi
  if (( have_species == 0 )) && [[ -s "$ge_tmp" ]]; then
    while IFS=$'\t' read -r pct clade name rank; do
      org_clade["$name"]=$(( ${org_clade["$name"]:-0} + ${clade%.*} )); org_rank["$name"]="G"
    done < "$ge_tmp"
  fi

  local detail_parts=()
  local identified="unclassified"
  if (( ${#org_clade[@]} > 0 )); then
    while IFS=$'\t' read -r k c; do
      local r="${org_rank[$k]:-NA}"
      local p2; p2=$(awk -v num="$c" -v den="$denom_for_pct" 'BEGIN{printf "%.2f", (den==0?0:num/den*100)}')
      echo -e "${grp}\t${sample}\t${r}\t${k}\t${c}\t${p2}\t${PCT_OF}" >> "$long_out"
      detail_parts+=( "${k} (${c}, ${p2}%)" )
    done < <(for k in "${!org_clade[@]}"; do echo -e "$k\t${org_clade[$k]}"; done | sort -k2,2nr)
    identified=$(for row in $(for k in "${!org_clade[@]}"; do echo -e "${org_clade[$k]}\t$k"; done | sort -k1,1nr | awk -F'\t' '{print $2}'); do printf "%s," "$row"; done | sed 's/,$//')
  fi

  # Genus table
  if [[ -s "$ge_tmp" ]]; then
    declare -A g_clade
    while IFS=$'\t' read -r pct clade name rank; do g_clade["$name"]=$(( ${g_clade["$name"]:-0} + ${clade%.*} )); done < "$ge_tmp"
    if (( ${#g_clade[@]} > 0 )); then
      while IFS=$'\t' read -r g c; do
        local gp; gp=$(awk -v num="$c" -v den="$denom_for_pct" 'BEGIN{printf "%.2f", (den==0?0:num/den*100)}')
        echo -e "${grp}\t${sample}\t${g}\t${c}\t${gp}\t${PCT_OF}" >> "$genus_out"
      done < <(for g in "${!g_clade[@]}"; do echo -e "$g\t${g_clade[$g]}"; done | sort -k2,2nr)
    fi
  fi

  # Percent classified + final summary row
  local pct_classified
  pct_classified=$(awk -v c="$classified_count" -v t="$denom_total" 'BEGIN{printf "%.2f", (t==0?0:c/t*100)}')

  if [[ "$MODE" == "reads" ]]; then
    echo -e "${grp}\t${sample}\t${total_seq}\t${total_bases}\tNA\tNA\tNA\tNA\tNA\t${identified}\t${classified_count}\t${unclassified_count}\t${pct_classified}\t${detail_parts[*]}" > "$sum_out"
  else
    echo -e "${grp}\t${sample}\t${total_seq}\t${total_bases}\t${N50}\t${L50}\t${N90}\t${L90}\t${GCpc}\t${identified}\t${classified_count}\t${unclassified_count}\t${pct_classified}\t${detail_parts[*]}" > "$sum_out"
  fi

  touch "$done_marker"
  echo "  => ${sample} done (C:${classified_count} U:${unclassified_count} %C:${pct_classified})"
}

# -------------------- Launch parallel jobs (robust PID array; no 'jobs', no 'wait -n') --------------------
pids=()
for sample in "${samples[@]}"; do
  ( process_one_sample "$sample" ) &           # start in subshell
  pids+=("$!")                                  # remember PID

  # If we reached the parallel limit, wait for the oldest PID to finish
  if (( ${#pids[@]} >= PARALLEL )); then
    pid="${pids[0]}"
    wait "$pid" || true                         # never crash pipeline if one sample fails
    # drop the finished pid from the head of the array
    pids=("${pids[@]:1}")
  fi
done

# Wait for any remaining background PIDs
for pid in "${pids[@]}"; do
  wait "$pid" || true
done


# -------------------- Merge per-sample outputs --------------------
echo -e "group\tsample\ttotal_sequences_in_assembly\ttotal_bases\tN50\tL50\tN90\tL90\tGC_percent\tidentified_organism\tclassified_count\tunclassified_count\tpercent_classified\tidentified_organism_detail" > "$IDENTIFIED_SUM"
cat "$TMPDIR/sum/"*.tsv 2>/dev/null >> "$IDENTIFIED_SUM" || true

echo -e "group\tsample\trank\torganism\tclade_reads\tpercent\tpercent_basis" > "$IDENTIFIED_LONG"
cat "$TMPDIR/long/"*.tsv 2>/dev/null >> "$IDENTIFIED_LONG" || true

echo -e "group\tsample\tgenus\tclade_reads\tpercent\tpercent_basis" > "$IDENTIFIED_GENUS"
cat "$TMPDIR/genus/"*.tsv 2>/dev/null >> "$IDENTIFIED_GENUS" || true

# -------------------- Group-level summary --------------------
tmp_species="$(mktemp)"
awk -F'\t' 'NR>1 && $10!~/^unclassified$/ { g=$1; split($10,a,","); for(i=1;i<=length(a);i++){ s=a[i]; gsub(/^[ \t]+|[ \t]+$/,"",s); if(s!="") print g"\t"s } }' "$IDENTIFIED_SUM" > "$tmp_species" || true
echo -e "group\ttotal_samples\tunclassified_samples\tpercent_unclassified\ttop_species" > "$GROUP_SUM"
if [[ -s "$IDENTIFIED_SUM" ]]; then
  mapfile -t groups < <(awk -F'\t' 'NR>1{print $1}' "$IDENTIFIED_SUM" | sort -u)
  for g in "${groups[@]}"; do
    total=$(awk -F'\t' -v G="$g" 'NR>1 && $1==G{c++}END{print c+0}' "$IDENTIFIED_SUM")
    unclassified=$(awk -F'\t' -v G="$g" 'NR>1 && $1==G && $10=="unclassified"{c++}END{print c+0}' "$IDENTIFIED_SUM")
    pct=$(awk -v t="${total:-0}" -v u="${unclassified:-0}" 'BEGIN{printf "%.2f", (t==0?0:u/t*100)}')
    top=$(awk -F'\t' -v G="$g" '$1==G{print $2}' "$tmp_species" | sort | uniq -c | sort -k1,1nr | head -n 5 | awk '{$1=$1; c=$1; $1=""; sub(/^ /,""); printf "%s:%d,", $0, c }' | sed 's/,$//')
    echo -e "${g}\t${total}\t${unclassified}\t${pct}\t${top}" >> "$GROUP_SUM"
  done
fi
rm -f "$tmp_species"

echo
echo "Summaries:"
echo " - $IDENTIFIED_SUM"
echo " - $IDENTIFIED_LONG"
echo " - $IDENTIFIED_GENUS"
echo " - $GROUP_SUM"
echo "Per-sample outputs: $OUTDIR"
