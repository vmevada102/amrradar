#!/usr/bin/env bash
set -euo pipefail

# run_resistome_extended.sh
# Integrated pipeline: ABRICATE + AMRFinder -> mob_recon -> extract AMR refs
# -> map reads (bwa) -> KMA -> (optional DIAMOND) -> TPM + presence/absence.

# Usage:
#   ./run_resistome_extended.sh -s sample.tsv [-o results] [-t threads] [--resume]
#                               [--card-db /path/to/card_nucl.fa]
#                               [--card-prot /path/to/card_proteins.fasta]

SAMPLE_TSV=""
OUT_ROOT="results"
THREADS=""
RESUME=false
CARD_DB="/data/home/amr/miniconda3/envs/resistome1/db/card/sequences"
CARD_PROT=""
WORKDIR_PREFIX="results/tools"
DMND=""

usage() {
  cat <<USAGE
Usage: $0 -s sample.tsv [-o results] [-t threads] [--resume]
          [--card-db PATH] [--card-prot PATH]

sample.tsv columns (minimum 3):
  sample<TAB>R1.fastq.gz<TAB>R2.fastq.gz[<TAB>other columns ignored]

Assembly is ALWAYS auto-detected as:
  results/tools/*_assembly/<sample>/<sample>.assembly.fasta

Outputs go to:
  results/tools/N+1_resistome_extended
where N is max numeric prefix among results/tools/N_*
USAGE
  exit 1
}

# -------------------- parse args --------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s) SAMPLE_TSV="$2"; shift 2 ;;
    -o) OUT_ROOT="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    --resume) RESUME=true; shift ;;
    --card-db) CARD_DB="$2"; shift 2 ;;
    --card-prot) CARD_PROT="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown arg: $1"; usage ;;
  esac
done

if [[ -z "$SAMPLE_TSV" ]]; then echo "ERROR: -s sample.tsv required"; usage; fi
if [[ ! -f "$SAMPLE_TSV" ]]; then echo "ERROR: sample.tsv not found: $SAMPLE_TSV"; exit 2; fi

# -------------------- CPU + threads --------------------
CPU_TOTAL=$(nproc 2>/dev/null || echo 4)
if [[ -z "$THREADS" ]]; then
  THREADS=$((CPU_TOTAL - 2))
  (( THREADS < 1 )) && THREADS=1
fi

echo "[INFO] Detected CPU cores : $CPU_TOTAL"
echo "[INFO] Threads per job    : $THREADS"

# -------------------- tool check --------------------
auto_install_tools() {
  local missing_tools=("$@")
  echo "The following required tools are missing:"
  printf '  - %s\n' "${missing_tools[@]}"

  if ! command -v conda >/dev/null 2>&1; then
    echo "conda not found; please install these tools manually and rerun."
    exit 3
  fi

  read -rp "Install them now with conda? [y/N] " ans
  if [[ ! "$ans" =~ ^[Yy]$ ]]; then
    echo "Aborting. Install the tools manually and rerun."
    exit 3
  fi

  conda install -y -c bioconda -c conda-forge "${missing_tools[@]}" || {
    echo "conda install failed; fix manually and rerun."; exit 3; }

  local still_missing=()
  for t in "${missing_tools[@]}"; do
    command -v "$t" >/dev/null 2>&1 || still_missing+=("$t")
  done
  if (( ${#still_missing[@]} )); then
    echo "Still missing after install:"
    printf '  - %s\n' "${still_missing[@]}"
    exit 3
  fi
}

REQUIRED_TOOLS=(abricate amrfinder mob_recon bwa samtools kma diamond seqkit python3 awk sort)
MISSING=()
for t in "${REQUIRED_TOOLS[@]}"; do
  command -v "$t" >/dev/null 2>&1 || MISSING+=("$t")
done
(( ${#MISSING[@]} )) && auto_install_tools "${MISSING[@]}"

# -------------------- CARD nucleotide DB --------------------
ensure_card_db() {
  if [[ -f "$CARD_DB" ]]; then
    echo "[INFO] Using CARD nucleotide DB at: $CARD_DB"
    return
  fi
  echo "[WARN] CARD nucleotide DB not found at $CARD_DB"

  if command -v abricate >/dev/null 2>&1; then
    if [[ -n "${CONDA_PREFIX:-}" && -f "$CONDA_PREFIX/db/abricate/card/sequences" ]]; then
      CARD_DB="$CONDA_PREFIX/db/abricate/card/sequences"
      echo "[INFO] Found CARD DB at: $CARD_DB"
      return
    fi
  fi

  echo "I can run 'abricate --setupdb' to download CARD and other DBs."
  read -rp "Run 'abricate --setupdb' now? [y/N] " ans
  if [[ ! "$ans" =~ ^[Yy]$ ]]; then
    echo "Download the CARD DB manually and rerun with --card-db"; exit 4;
  fi

  abricate --setupdb || { echo "abricate --setupdb failed"; exit 4; }

  if [[ -n "${CONDA_PREFIX:-}" && -f "$CONDA_PREFIX/db/abricate/card/sequences" ]]; then
    CARD_DB="$CONDA_PREFIX/db/abricate/card/sequences"
  fi

  if [[ ! -f "$CARD_DB" ]]; then
    echo "CARD DB still not found; locate card/sequences and use --card-db"; exit 4;
  fi
  echo "[INFO] Using CARD nucleotide DB at: $CARD_DB"
}
ensure_card_db

# -------------------- choose OUTDIR: N+1_resistome_extended --------------------
mkdir -p "$WORKDIR_PREFIX"
Nmax=0
for d in "$WORKDIR_PREFIX"/*; do
  [[ -d "$d" ]] || continue
  b=$(basename "$d")
  if [[ "$b" =~ ^([0-9]+)_ ]]; then
    n=${BASH_REMATCH[1]}
    (( n > Nmax )) && Nmax=$n
  fi
done
OUTDIR="$WORKDIR_PREFIX/$((Nmax+1))_resistome_extended"
mkdir -p "$OUTDIR"
echo "[INFO] Outputs to: $OUTDIR"

# -------------------- CARD protein / DIAMOND DB --------------------
ensure_card_prot() {
  if [[ -n "$CARD_PROT" ]]; then
    if [[ -f "$CARD_PROT" ]]; then
      echo "[INFO] Using CARD protein FASTA at: $CARD_PROT"
      return
    else
      echo "[WARN] CARD protein FASTA not found at $CARD_PROT"
    fi
  fi

  if [[ -z "$CARD_PROT" && -n "${CONDA_PREFIX:-}" ]]; then
    for cand in "$CONDA_PREFIX/db/card/protein_fasta" \
                "$CONDA_PREFIX/db/card/protein_fasta.fasta" \
                "$CONDA_PREFIX/db/card/protein_fasta.fa"; do
      if [[ -f "$cand" ]]; then
        CARD_PROT="$cand"
        echo "[INFO] Found CARD protein FASTA at: $CARD_PROT"
        return
      fi
    done
  fi

  echo "[WARN] No CARD protein FASTA (--card-prot) specified or found."
  echo "I can try to download it and build a DIAMOND DB."
  read -rp "Download CARD protein FASTA and build DIAMOND DB now? [y/N] " ans
  if [[ ! "$ans" =~ ^[Yy]$ ]]; then
    echo "[INFO] DIAMOND step will be skipped."
    CARD_PROT=""
    return
  fi

  local url="https://card.mcmaster.ca/latest/protein_fasta"
  CARD_PROT="$OUTDIR/card_proteins.fasta"
  echo "[INFO] Downloading CARD protein FASTA to $CARD_PROT ..."
  if command -v wget >/dev/null 2>&1; then
    wget -O "$CARD_PROT" "$url" > "$OUTDIR/card_prot_download.log" 2>&1 || {
      echo "[WARN] Download failed; skipping DIAMOND (see card_prot_download.log)"; CARD_PROT=""; return; }
  elif command -v curl >/dev/null 2>&1; then
    curl -L "$url" -o "$CARD_PROT" > "$OUTDIR/card_prot_download.log" 2>&1 || {
      echo "[WARN] Download failed; skipping DIAMOND (see card_prot_download.log)"; CARD_PROT=""; return; }
  else
    echo "[WARN] Neither wget nor curl available; skipping DIAMOND."
    CARD_PROT=""
    return
  fi

  if [[ ! -s "$CARD_PROT" ]]; then
    echo "[WARN] Downloaded CARD protein FASTA is empty; skipping DIAMOND."
    CARD_PROT=""
    return
  fi
  echo "[INFO] Downloaded CARD protein FASTA to: $CARD_PROT"
}

ensure_card_prot

# KMA DB for CARD
KMA_DB="$OUTDIR/card_kma_db"
if [[ ! -f "${KMA_DB}.index" && ! -f "${KMA_DB}.seq" ]]; then
  echo "[INFO] Building KMA DB (CARD index) ..."
  kma_index -i "$CARD_DB" -o "$KMA_DB" > "$OUTDIR/kma_index.log" 2>&1 || {
    echo "kma_index failed, see $OUTDIR/kma_index.log"; exit 5; }
  echo "[INFO] KMA DB built (log: $OUTDIR/kma_index.log)"
fi

# DIAMOND DB
if [[ -n "$CARD_PROT" ]]; then
  DMND="$OUTDIR/card_proteins.dmnd"
  if [[ ! -f "$DMND" ]]; then
    echo "[INFO] Building DIAMOND DB at $DMND ..."
    if diamond makedb --in "$CARD_PROT" -d "$DMND" > "$OUTDIR/diamond_makedb.log" 2>&1; then
      echo "[INFO] DIAMOND DB built (log: $OUTDIR/diamond_makedb.log)"
    else
      echo "[WARN] diamond makedb failed; DIAMOND step will be skipped (see diamond_makedb.log)"
      DMND=""
    fi
  fi
else
  DMND=""
fi

mkdir -p script

# -------------------- helper scripts (auto-created) --------------------
if [[ ! -f script/compute_amr_quant.py ]]; then
  cat > script/compute_amr_quant.py <<'PY'
#!/usr/bin/env python3
import sys, pandas as pd
if len(sys.argv) < 3:
    print("Usage: compute_amr_quant.py <combined_counts.tsv> <out_tpm.tsv>")
    sys.exit(1)
inf, outf = sys.argv[1], sys.argv[2]
df = pd.read_csv(inf, sep='\t')
df['length_kb'] = df['length'].replace(0,1)/1000.0
df['RPK'] = df['mapped_reads'] / df['length_kb']
rows = []
for s, g in df.groupby('sample'):
    g = g.copy()
    sum_rpk = g['RPK'].sum()
    if sum_rpk <= 0:
        g['TPM'] = 0.0
    else:
        g['TPM'] = (g['RPK'] / sum_rpk) * 1e6
    total_mapped = g['mapped_reads'].sum()
    if total_mapped <= 0:
        g['RPKM'] = 0.0
    else:
        g['RPKM'] = g['RPK'] / (total_mapped/1e6)
    rows.append(g)
res = pd.concat(rows, ignore_index=True)
res.to_csv(outf, sep='\t', index=False)
print("Wrote", outf)
PY
  chmod +x script/compute_amr_quant.py
fi

if [[ ! -f script/parse_mob_plasmids.py ]]; then
  cat > script/parse_mob_plasmids.py <<'PY'
#!/usr/bin/env python3
import sys, os
d = sys.argv[1] if len(sys.argv) > 1 else "."
out = set()
pfn = os.path.join(d, "plasmids.fna")
if os.path.isfile(pfn):
    with open(pfn) as fh:
        for line in fh:
            if line.startswith(">"):
                out.add(line[1:].split()[0])
rpt = os.path.join(d, "mob_recon_report.txt")
if os.path.isfile(rpt):
    with open(rpt) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            for p in parts:
                if p.startswith("NODE_") or p.startswith("contig_") or ":" in p:
                    out.add(p.strip())
for x in sorted(out):
    print(x)
PY
  chmod +x script/parse_mob_plasmids.py
fi

# -------------------- sample processing --------------------
process_sample() {
  local idx="$1" total="$2" sample="$3" r1="$4" r2="$5"
  local sample_out="$OUTDIR/$sample"
  mkdir -p "$sample_out"

  echo "=== START [$idx/$total] $sample ==="

  if $RESUME && [[ -f "$sample_out/${sample}.pipeline.OK" ]]; then
    echo "[$sample] Resume: pipeline.OK exists, skipping."
    return 0
  fi

  # Auto-detect assembly ONLY from results/tools/*_assembly/<sample>/<sample>.assembly.fasta
  local assembly=""
  local adir
  adir=$(ls -d results/tools/*_assembly/"$sample" 2>/dev/null | head -n1 || true)
  if [[ -n "$adir" && -f "$adir/$sample.assembly.fasta" ]]; then
    assembly="$adir/$sample.assembly.fasta"
  fi

  if [[ -n "$assembly" ]]; then
    echo "[$sample] Using assembly: $assembly"
  else
    echo "[$sample] No assembly found at results/tools/*_assembly/$sample/$sample.assembly.fasta"
  fi

  # 1) ABRICATE
  if [[ -n "$assembly" ]]; then
    echo "[$sample] ABRICATE..."
    abricate --db card --minid 80 --mincov 80 "$assembly" > "$sample_out/abricate_card.tab" 2> "$sample_out/abricate_card.log" || true
    abricate --db resfinder --minid 80 --mincov 80 "$assembly" > "$sample_out/abricate_resfinder.tab" 2> "$sample_out/abricate_resfinder.log" || true
    abricate --db plasmidfinder --minid 80 --mincov 80 "$assembly" > "$sample_out/abricate_plasmidfinder.tab" 2> "$sample_out/abricate_plasmidfinder.log" || true
  else
    echo "[$sample] No assembly; skipping ABRICATE and AMRFinder and mob_recon."
  fi

  # 2) AMRFinder
  if [[ -n "$assembly" ]]; then
    echo "[$sample] AMRFinder..."
    amrfinder -n "$assembly" -o "$sample_out/amrfinder.tsv" --plus --threads "$THREADS" 2> "$sample_out/amrfinder.log" || true
  fi

  # 3) mob_recon
  local mob_dir="$sample_out/mob_recon"
  if [[ -n "$assembly" ]]; then
    echo "[$sample] mob_recon..."
    mob_recon -i "$assembly" -o "$mob_dir" -t "$THREADS" > "$sample_out/mob_recon.stdout" 2> "$sample_out/mob_recon.stderr" || true
  fi

  # 4) plasmid contigs list
  if [[ -d "$mob_dir" ]]; then
    python3 script/parse_mob_plasmids.py "$mob_dir" > "$sample_out/${sample}.mob_plasmid_contigs.txt" || true
  fi

  # 5) Build per-sample AMR reference FASTA from CARD
  local amr_refs="$sample_out/${sample}_amr_refs.fasta"
  : > "$amr_refs"
  local idlist="$sample_out/${sample}.abricate_ids.txt"
  : > "$idlist"

  if [[ -f "$sample_out/abricate_card.tab" ]]; then
    tail -n +2 "$sample_out/abricate_card.tab" | awk 'NF{print $1}' >> "$idlist" || true
  fi
  if [[ -f "$sample_out/amrfinder.tsv" ]]; then
    tail -n +2 "$sample_out/amrfinder.tsv" | awk 'NF{print $6}' >> "$idlist" || true
  fi
  if [[ -s "$idlist" ]]; then
    sort -u "$idlist" -o "$idlist"
    echo "[$sample] Extracting AMR refs from CARD..."
    seqkit grep -f "$idlist" "$CARD_DB" -o "$amr_refs" 2> "$sample_out/seqkit.log" || true
  fi

  # 6) Map reads to AMR refs (bwa)
  local bam="$sample_out/${sample}.amr.bam"
  if [[ -s "$amr_refs" ]]; then
    echo "[$sample] bwa mem to AMR refs..."
    bwa index "$amr_refs" >/dev/null 2>&1 || true
    bwa mem -t "$THREADS" "$amr_refs" "$r1" "$r2" 2> "$sample_out/bwa_amr.log" | \
      samtools view -bS - | samtools sort -o "$bam" || true
    if [[ -f "$bam" ]]; then
      samtools index "$bam" || true
      samtools idxstats "$bam" > "$sample_out/${sample}.amr_idxstats.tsv" || true
    fi
  else
    echo "[$sample] No per-sample AMR refs produced."
  fi

  # 7) KMA mapping to CARD DB
  echo "[$sample] KMA mapping to CARD DB..."
  kma -ipe "$r1" "$r2" -o "$sample_out/${sample}_kma" -t_db "$KMA_DB" -t "$THREADS" &> "$sample_out/${sample}_kma.log" || true
  if [[ -f "$sample_out/${sample}_kma.res" ]]; then
    tail -n +2 "$sample_out/${sample}_kma.res" > "$sample_out/${sample}_kma_hits.tsv" || true
  fi

  # 8) DIAMOND blastx (optional)
  if [[ -n "$DMND" && -f "$DMND" ]]; then
    echo "[$sample] DIAMOND blastx (R1)..."
    diamond blastx -p "$THREADS" --query "$r1" --db "$DMND" -a "$sample_out/${sample}_R1" -sensitive &> "$sample_out/${sample}_diamond_R1.log" || true
    diamond view -a "$sample_out/${sample}_R1.daa" -o "$sample_out/${sample}_R1.m8" || true
    echo "[$sample] DIAMOND blastx (R2)..."
    diamond blastx -p "$THREADS" --query "$r2" --db "$DMND" -a "$sample_out/${sample}_R2" -sensitive &> "$sample_out/${sample}_diamond_R2.log" || true
    diamond view -a "$sample_out/${sample}_R2.daa" -o "$sample_out/${sample}_R2.m8" || true
  fi

  # 9) Counts summary
  local cntf="$sample_out/${sample}.counts.tsv"
  : > "$cntf"
  if [[ -f "$sample_out/${sample}.amr_idxstats.tsv" ]]; then
    awk -v S="$sample" 'BEGIN{OFS="\t"} {print S,$1,$2,$3}' "$sample_out/${sample}.amr_idxstats.tsv" >> "$cntf" || true
  fi
  if [[ -f "$sample_out/${sample}_kma_hits.tsv" ]]; then
    awk -v S="$sample" 'BEGIN{OFS="\t"} {g=$1; len=$2+0; cov=$5+0; mapped=int(cov*len/100); print S,g,len,mapped}' \
      "$sample_out/${sample}_kma_hits.tsv" >> "$cntf" || true
  fi

  # 10) Annotate ABRICATE hits as plasmid/chromosome
  if [[ -f "$sample_out/${sample}.mob_plasmid_contigs.txt" && -f "$sample_out/abricate_card.tab" ]]; then
    awk -F'\t' -v OFS='\t' \
        'FNR==NR{p[$1]=1; next} /^#/ {next} {loc = ($2 in p) ? "plasmid" : "chromosome"; print $0,loc}' \
        "$sample_out/${sample}.mob_plasmid_contigs.txt" "$sample_out/abricate_card.tab" \
        > "$sample_out/abricate_card.annotated.tab" || true
  fi

  touch "$sample_out/${sample}.pipeline.OK"
  echo "=== DONE  [$idx/$total] $sample ==="
}

# -------------------- main loop over samples --------------------
TOTAL_SAMPLES=$(grep -v '^#' "$SAMPLE_TSV" | awk 'NF>0' | wc -l)
echo "[INFO] Total samples to process: $TOTAL_SAMPLES"

idx=0
while IFS=$'\t' read -r sample r1 r2 rest; do
  # skip empty lines and comments
  [[ -z "${sample:-}" || "$sample" =~ ^# ]] && continue
  # skip header
  if [[ "$sample" == "sample" || "$sample" == "Sample" ]]; then
    continue
  fi
  idx=$((idx+1))
  process_sample "$idx" "$TOTAL_SAMPLES" "$sample" "$r1" "$r2"
done < "$SAMPLE_TSV"

# -------------------- combine results --------------------
COMBINED="$OUTDIR/amr_combined_counts.tsv"
echo -e "sample\tgene\tlength\tmapped_reads" > "$COMBINED"
for d in "$OUTDIR"/*/; do
  [[ -d "$d" ]] || continue
  s=$(basename "$d")
  f="$d/${s}.counts.tsv"
  [[ -f "$f" ]] && cat "$f" >> "$COMBINED"
done

if [[ -s "$COMBINED" ]]; then
  python3 script/compute_amr_quant.py "$COMBINED" "$OUTDIR/amr_quant_tpm.tsv"
else
  echo "[WARN] No combined counts; skipping TPM."
fi

python3 - <<PY
import pandas as pd, os
outdir="$OUTDIR"
combined=os.path.join(outdir,"amr_quant_tpm.tsv")
rows=[]
if os.path.exists(combined):
    df=pd.read_csv(combined,sep='\t')
    df['present'] = ((df['mapped_reads']>0) | (df['TPM']>0))
    for (s,g),gdf in df.groupby(['sample','gene']):
        rows.append((s,g,int(gdf['present'].iloc[0])))
if rows:
    pd.DataFrame(rows, columns=['sample','gene','present']).to_csv(
        os.path.join(outdir,"presence_absence.tsv"), sep='\t', index=False
    )
    print("Wrote presence_absence.tsv")
else:
    print("No rows for presence/absence.")
PY

echo "[INFO] Pipeline complete."
echo "  - per-sample dirs: $OUTDIR/<sample>/"
echo "  - combined counts: $OUTDIR/amr_combined_counts.tsv"
echo "  - quantified TPM:  $OUTDIR/amr_quant_tpm.tsv (if any)"
echo "  - presence/absence: $OUTDIR/presence_absence.tsv (if any)"
