#!/usr/bin/env bash
# run_step1_3_compare.sh (finalized robustness for empty-newick and QUAST input)
set -euo pipefail

# ------------------------------
# AUTO THREAD DETECTION
# ------------------------------
# Detect total CPU cores and use (N - 2), but at least 1
CPU_TOTAL=$(nproc 2>/dev/null || echo 4)
threads=$((CPU_TOTAL - 2))
if [[ $threads -lt 1 ]]; then
  threads=1
fi

force=0
strict_visuals=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads|-t) threads="$2"; shift 2;;
    --force) force=1; shift;;
    --strict-visuals) strict_visuals=1; shift;;
    -h|--help) sed -n '1,260p' "$0"; exit 0;;
    *) echo "[WARN] Unknown arg: $1" >&2; shift;;
  esac
done

ts() { date +"%Y-%m-%d %H:%M:%S"; }

mkdir -p results/tools
next_N() {
  shopt -s nullglob
  local maxN=0
  local d base num
  for d in results/tools/*/ ; do
    base=$(basename "$d")
    if [[ "$base" =~ ^([0-9]+)_ ]]; then
      num="${BASH_REMATCH[1]}"
      (( num > maxN )) && maxN=$num
    fi
  done
  echo $((maxN + 1))
}

sanitize_label() {
  local s="$1"
  s="${s// /_}"; s="${s//,/_}"; s="${s//(/_}"; s="${s//)/_}"; s="${s//:/_}"; s="${s//;/_}"
  echo "$s"
}

N=$(next_N)
OUTDIR="results/tools/${N}_comparativegenomics"
mkdir -p "$OUTDIR"
LOG="$OUTDIR/pipeline.log"
exec > >(tee -a "$LOG") 2>&1

echo "[$(ts)] Output directory: $OUTDIR"
echo "[$(ts)] CPU total: $CPU_TOTAL, using threads: $threads"

REQUIRED_TOOLS=(quast.py fastANI mashtree)
MISSING_REQ=()
for t in "${REQUIRED_TOOLS[@]}"; do
  if ! command -v "$t" >/dev/null 2>&1; then
    MISSING_REQ+=("$t")
  fi
done
if (( ${#MISSING_REQ[@]} > 0 )); then
  echo "[$(ts)] [ERROR] Missing REQUIRED tools:"
  printf '  - %s\n' "${MISSING_REQ[@]}"
  echo "  mamba install -c bioconda -c conda-forge ${MISSING_REQ[*]}"
  exit 2
fi

HAS_MASH=0; command -v mash >/dev/null 2>&1 && HAS_MASH=1 || true
HAS_RAPIDNJ=0; command -v rapidnj >/dev/null 2>&1 && HAS_RAPIDNJ=1 || true
HAS_FASTME=0; command -v fastme >/dev/null 2>&1 && HAS_FASTME=1 || true

echo "[$(ts)] Searching assemblies: results/tools/*_assembly/*/*.assembly.fasta"
shopt -s nullglob

ASSEMBLIES=( results/tools/*_assembly/*/*.assembly.fasta )
FOUND=${#ASSEMBLIES[@]}
echo "[$(ts)] Found ${FOUND} source assemblies"

# List of all candidate assemblies (absolute paths)
> "$OUTDIR/assembly_files.list"

# Validate directly on original paths
BADLIST="$OUTDIR/assemblies_excluded.txt"; > "$BADLIST"
KEPTLIST="$OUTDIR/assemblies_kept.txt"; > "$KEPTLIST"
VALID_FILES=()

for f in "${ASSEMBLIES[@]}"; do
  # get absolute path (in case script is run from other directories)
  abs=$(readlink -f "$f" 2>/dev/null || echo "$f")
  sample=$(basename "$(dirname "$abs")")
  safe=$(sanitize_label "$sample")

  echo "$abs" >> "$OUTDIR/assembly_files.list"

  # file must exist and be >0
  if [[ ! -s "$abs" ]]; then
    echo "${safe}.fasta [empty file]" >> "$BADLIST"; continue
  fi
  if ! grep -m1 -q '^>' "$abs"; then
    echo "${safe}.fasta [no FASTA header]" >> "$BADLIST"; continue
  fi
  if ! awk 'BEGIN{ok=0} /^>/ {next} {if ($0 ~ /[ACGTNacgtn]/) ok=1} END{exit ok?0:1}' "$abs"; then
    echo "${safe}.fasta [no sequence lines]" >> "$BADLIST"; continue
  fi

  VALID_FILES+=( "$abs" )
  echo "${safe}.fasta" >> "$KEPTLIST"
done

VALID_COUNT=${#VALID_FILES[@]}
if [[ -s "$BADLIST" ]]; then
  echo "[$(ts)] [WARN] Excluding invalid FASTA files:"
  sed 's/^/  - /' "$BADLIST"
fi

echo "[$(ts)] Using ${VALID_COUNT} valid assemblies"
if (( VALID_COUNT == 0 )); then
  echo "[$(ts)] [ERROR] No valid assemblies remain after validation."
  echo "Please inspect $OUTDIR/assemblies_excluded.txt"
  exit 3
fi

# Debug preview of first few valid files
echo "[$(ts)] Valid examples:"
for x in "${VALID_FILES[@]:0:5}"; do echo "  - $(basename "$x")"; done

# --- QUAST ---
mkdir -p "$OUTDIR/quast"
QUAST_DONE_FLAG="$OUTDIR/.quast.done"
if [[ $force -eq 1 || ! -f "$QUAST_DONE_FLAG" ]]; then
  echo "[$(ts)] [QUAST] Running QUAST on ${VALID_COUNT} assemblies..."
  set +e
  quast.py -o "$OUTDIR/quast" -t "$threads" "${VALID_FILES[@]}" 2>&1
  rc=$?
  set -e
  if [[ $rc -ne 0 ]]; then
    echo "[$(ts)] [ERROR] QUAST failed with exit code $rc"
    echo "Debug: first valid file content (head -n 3):"
    head -n 3 "${VALID_FILES[0]}"
    exit 4
  fi
  touch "$QUAST_DONE_FLAG"
else
  echo "[$(ts)] [QUAST] Skipping (already done). Use --force to rerun."
fi

# --- FastANI ---
mkdir -p "$OUTDIR/fastani"
FASTANI_OUT="$OUTDIR/fastani/fastani.tsv"
FASTANI_DONE_FLAG="$OUTDIR/.fastani.done"
if [[ $force -eq 1 || ! -f "$FASTANI_DONE_FLAG" ]]; then
  echo "[$(ts)] [FastANI] Computing all-vs-all ANI on ${VALID_COUNT} assemblies..."
  printf "%s\n" "${VALID_FILES[@]}" > "$OUTDIR/fastani/assemblies.lst"
  set +e
  fastANI --ql "$OUTDIR/fastani/assemblies.lst" \
          --rl "$OUTDIR/fastani/assemblies.lst" \
          -o "$FASTANI_OUT" -t "$threads" 2>&1
  rc=$?
  set -e
  if [[ $rc -ne 0 ]]; then
    echo "[$(ts)] [ERROR] FastANI failed with exit code $rc"
    exit 5
  fi
  touch "$FASTANI_DONE_FLAG"
else
  echo "[$(ts)] [FastANI] Skipping (already done). Use --force to rerun."
fi

# --- Tree building ---
mkdir -p "$OUTDIR/mashtree"
MASHTREE_NWK="$OUTDIR/mashtree/mashtree.nwk"
MASHTREE_DONE_FLAG="$OUTDIR/.mashtree.done"
MASHTREE_LOG="$OUTDIR/mashtree/mashtree.log"

run_mashtree() {
  echo "[$(ts)] [Mashtree] Building quick tree..."
  set +e
  # Use the correct mashtree option: --numcpus (not --num-cpus)
  mashtree --numcpus "$threads" --outtree "$MASHTREE_NWK" "${VALID_FILES[@]}" > "$MASHTREE_LOG" 2>&1
  rc=$?
  set -e
  return $rc
}

if [[ $force -eq 1 || ! -f "$MASHTREE_DONE_FLAG" ]]; then
  if (( VALID_COUNT == 1 )); then
    one=$(basename "${VALID_FILES[0]}"); one="${one%.fasta}"
    echo "${one};" > "$MASHTREE_NWK"
    echo "[$(ts)] [Mashtree] Only one assembly ? trivial tree."
    touch "$MASHTREE_DONE_FLAG"
  else
    if run_mashtree; then
      touch "$MASHTREE_DONE_FLAG"
    else
      rc=$?
      echo "[$(ts)] [ERROR] Mashtree failed with exit code $rc"
      echo "---------- mashtree.log (last 200 lines) ----------"
      tail -n 200 "$MASHTREE_LOG" || true
      echo "---------------------------------------------------"
      echo "[$(ts)] [Fallback] Building tree via Mash + rapidnj/fastme"
      if ! command -v mash >/dev/null 2>&1; then
        echo "[$(ts)] [ERROR] mash not available. Install: mamba install -c bioconda mash"
        exit 6
      fi
      mash sketch -o "$OUTDIR/mashtree/all" "${VALID_FILES[@]}" >/dev/null 2>&1
      mash dist "$OUTDIR/mashtree/all.msh" "$OUTDIR/mashtree/all.msh" > "$OUTDIR/mashtree/mash_pairs.tsv"
      # Convert Mash pairs to PHYLIP using awk (names are stems)
      awk '
      BEGIN{OFS="\t"}
      {a=$1; b=$2; d=$3;
       gsub(/^.*\//,"",a); gsub(/\.fasta$/,"",a);
       gsub(/^.*\//,"",b); gsub(/\.fasta$/,"",b);
       names[a]=1; names[b]=1; pair[a,b]=d; pair[b,a]=d
      }
      END{
        # create stable order from input taxa list (VALID_FILES)
      }' /dev/null > /dev/null

      # Build taxa order from VALID_FILES
      : > "$OUTDIR/mashtree/taxa.list"
      for f in "${VALID_FILES[@]}"; do
        bn=$(basename "$f"); echo "${bn%.fasta}" >> "$OUTDIR/mashtree/taxa.list"
      done

      # Use Python if available for deterministic PHYLIP creation
      if command -v python3 >/dev/null 2>&1; then
        python3 - "$OUTDIR/mashtree/mash_pairs.tsv" "$OUTDIR/mashtree/taxa.list" "$OUTDIR/mashtree/mash_dist.phylip" <<'PY'
import sys, math
from pathlib import Path

pairs = Path(sys.argv[1]).read_text().splitlines()
taxa  = Path(sys.argv[2]).read_text().splitlines()
idx = {t:i for i,t in enumerate(taxa)}
n = len(taxa)
M = [[0.0 if i==j else math.nan for j in range(n)] for i in range(n)]
for line in pairs:
    parts = line.strip().split()
    if len(parts) < 3: continue
    a = Path(parts[0]).stem
    b = Path(parts[1]).stem
    try:
        d = float(parts[2])
    except Exception:
        continue
    if a in idx and b in idx:
        i, j = idx[a], idx[b]
        M[i][j] = d; M[j][i] = d
# Fill missing with small distance
for i in range(n):
    for j in range(n):
        if M[i][j] != M[i][j]:  # NaN
            M[i][j] = 0.001 if i!=j else 0.0
with open(sys.argv[3], "w") as out:
    out.write(f"{n}\n")
    for i in range(n):
        row = [f"{M[i][k]:.6f}" for k in range(i+1)]
        out.write(f"{taxa[i]}\t" + "\t".join(row) + "\n")
PY
      else
        # Pure awk fallback with arbitrary order
        awk '
        BEGIN{OFS="\t"}
        FNR==NR{bn=$0; order[++i]=bn; next}
        {
          a=$1; b=$2; d=$3;
          gsub(/^.*\//,"",a); gsub(/\.fasta$/,"",a);
          gsub(/^.*\//,"",b); gsub(/\.fasta$/,"",b);
          pair[a,b]=d; pair[b,a]=d; names[a]=1; names[b]=1
        }
        END{
          n=i; print n;
          for (r=1; r<=n; r++){
            name=order[r]; line=name;
            for (c=1; c<=r; c++){
              other=order[c];
              if (name==other) val=0; else {
                key1=name SUBSEP other; key2=other SUBSEP name;
                if (key1 in pair) val=pair[key1]; else if (key2 in pair) val=pair[key2]; else val=0.001
              }
              line=line OFS val
            }
            print line
          }
        }' "$OUTDIR/mashtree/taxa.list" "$OUTDIR/mashtree/mash_pairs.tsv" > "$OUTDIR/mashtree/mash_dist.phylip"
      fi

      if command -v rapidnj >/dev/null 2>/dev/null; then
        rapidnj -i pd -o t "$OUTDIR/mashtree/mash_dist.phylip" > "$MASHTREE_NWK"
      elif command -v fastme >/dev/null 2>/dev/null; then
        fastme -i "$OUTDIR/mashtree/mash_dist.phylip" -O "$MASHTREE_NWK" >/dev/null 2>&1
      else
        echo "[$(ts)] [ERROR] Need rapidnj or fastme. Install one: mamba install -c bioconda rapidnj  (or) fastme"
        exit 6
      fi
      touch "$MASHTREE_DONE_FLAG"
    fi
  fi
else
  echo "[$(ts)] [Mashtree] Skipping (already done). Use --force to rerun."
fi

# Visualization
TREE_PNG="$OUTDIR/mashtree/mashtree.png"
VIS_DONE_FLAG="$OUTDIR/.visual.done"
if [[ $force -eq 1 || ! -f "$VIS_DONE_FLAG" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    set +e
    python3 - <<PY
import sys
from pathlib import Path
try:
    from Bio import Phylo
    import matplotlib.pyplot as plt
except Exception:
    sys.exit(0)
nwk = Path("$MASHTREE_NWK")
png = Path("$TREE_PNG")
if not nwk.exists() or nwk.stat().st_size < 3:
    sys.exit(0)
tree = Phylo.read(str(nwk), "newick")
fig = plt.figure(figsize=(10, 16))
ax = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, do_show=False, axes=ax)
fig.tight_layout()
fig.savefig(str(png), dpi=200)
print(f"[VIS] Wrote {png}")
PY
    rc=$?
    set -e
    if [[ $rc -ne 0 ]]; then
      echo "[$(ts)] [WARN] Visualization non-zero exit ($rc)."
      if [[ $strict_visuals -eq 1 ]]; then exit 7; fi
    else
      touch "$VIS_DONE_FLAG"
    fi
  else
    echo "[$(ts)] [INFO] python3 not available; skipping PNG."
  fi
else
  echo "[$(ts)] [VIS] Skipping (already done). Use --force to rerun."
fi

# Index
INDEX="$OUTDIR/index.html"
if [[ $force -eq 1 || ! -f "$INDEX" ]]; then
  {
    echo "<!doctype html><meta charset='utf-8'><title>Comparative Genomics (Steps 1–3)</title>"
    echo "<h1>Comparative Genomics (Steps 1–3)</h1>"
    echo "<ul>"
    echo "<li><a href='quast/report.html'>QUAST report</a></li>"
    if [[ -f "$FASTANI_OUT" ]]; then echo "<li><a href='fastani/fastani.tsv'>FastANI matrix</a></li>"; fi
    if [[ -f "$MASHTREE_NWK" ]]; then echo "<li><a href='mashtree/mashtree.nwk'>Tree (Newick)</a></li>"; fi
    if [[ -f "$TREE_PNG" ]]; then echo "<li><a href='mashtree/mashtree.png'>Tree (PNG)</a></li>"; fi
    echo "</ul>"
    echo "<p>Log: <a href='pipeline.log'>pipeline.log</a></p>"
    echo "<p>Mashtree log: <a href='mashtree/mashtree.log'>mashtree.log</a></p>"
  } > "$INDEX"
fi

echo "[$(ts)] Done. Outputs in: $OUTDIR"
exit 0
