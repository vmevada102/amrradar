#!/usr/bin/env bash
# =====================================================================
#  One_Health-AMRfinder-module.sh  (patched v2)
#  Unified AMR analysis pipeline:
#   1) Run AMRFinderPlus on all assemblies
#   2) Combine all .amrfinder.tsv files -> adds 'sample' column
#   3) Merge 'group' info from sample.tsv (adds 'group' before 'sample')
#   4) Generate Excel + heatmaps (sample & group); zero-hit samples included
# =====================================================================

set -euo pipefail

SCRIPT_DIR="$(dirname "$(realpath "$0")")"
RESULTS="results"
THREADS=40
PARALLEL=20
RESUME=false
PYTHON_BIN=$(command -v python3 || echo python3)
AMRFINDER_BIN=$(command -v amrfinder || true)

# ------------------ ARGUMENT PARSING ------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results)   RESULTS="$2"; shift 2;;
    --threads)   THREADS="$2"; shift 2;;
    --parallel)  PARALLEL="$2"; shift 2;;
    --amrfinder) AMRFINDER_BIN="$2"; shift 2;;
    --resume)    RESUME=true; shift;;
    *)           echo "[WARN] Ignoring unknown arg: $1"; shift;;
  esac
done

if [[ -z "$AMRFINDER_BIN" ]]; then
  echo "ERROR: amrfinder not found. Use --amrfinder /path/to/amrfinder"
  exit 1
fi

TOOLS_DIR="${RESULTS%/}/tools"
mkdir -p "$TOOLS_DIR"

# ------------------ STEP 1: Determine next step number ------------------
next_step_number() {
  local base="$1"; local max=0
  while IFS= read -r -d $'\0' d; do
    bn=$(basename "$d")
    if [[ "$bn" =~ ^([0-9]+)_ ]]; then
      n=${BASH_REMATCH[1]}; n=$((10#$n))
      (( n > max )) && max=$n
    fi
  done < <(find "$base" -maxdepth 1 -mindepth 1 -type d -print0 2>/dev/null)
  echo $((max + 1))
}

STEP_NUM=$(next_step_number "$TOOLS_DIR")
AMR_DIR="$TOOLS_DIR/${STEP_NUM}_AMRFinder"
mkdir -p "$AMR_DIR"

echo "?? Results directory: $AMR_DIR"
echo "?? Both sample- and group-level heatmaps will be generated."

# ------------------ STEP 2: Run AMRFinder ------------------
echo "?? Searching for assemblies..."
mapfile -t ASSEMBLIES < <(find "$TOOLS_DIR" -type f -name "*.assembly.fasta" | sort)
if [[ ${#ASSEMBLIES[@]} -eq 0 ]]; then
  echo "No assemblies found under $TOOLS_DIR/*_assembly"
  exit 2
fi

for asm in "${ASSEMBLIES[@]}"; do
  sample=$(basename "$asm" | cut -d. -f1)
  outdir="$AMR_DIR/$sample"
  mkdir -p "$outdir"
  outtsv="$outdir/${sample}.amrfinder.tsv"
  if [[ "$RESUME" == true && -s "$outtsv" ]]; then
    echo "[$sample] Skipping (resume mode)"
    continue
  fi
  echo "[$sample] Running AMRFinder..."
  "$AMRFINDER_BIN" -n "$asm" -o "$outtsv" --threads "$THREADS" || echo "[$sample] warning: amrfinder failed"
done
echo "? All AMRFinder jobs finished."

# ------------------ STEP 3: Combine TSV Files (add 'sample' column) ------------------
echo "?? Combining *.amrfinder.tsv files..."
COMBINED_TSV="$AMR_DIR/combined_amrfinder.tsv"

if [[ ! -s "$COMBINED_TSV" ]]; then
  : > "$COMBINED_TSV"
  header_written=false
  shopt -s nullglob
  for f in "$AMR_DIR"/*/*.amrfinder.tsv; do
    sample="$(basename "$f")"
    sample="${sample%.amrfinder.tsv}"
    if [[ "$header_written" = false ]]; then
      # write header once, prefixed with 'sample'
      { printf "sample\t"; head -n1 "$f"; } >> "$COMBINED_TSV"
      header_written=true
    fi
    # append rows with sample as first column; skip header lines
    awk -v s="$sample" 'BEGIN{OFS="\t"} NR>1 {print s, $0}' "$f" >> "$COMBINED_TSV"
  done
  shopt -u nullglob
  if [[ ! -s "$COMBINED_TSV" ]]; then
    echo "[ERROR] No AMRFinder TSV content found to combine."; exit 1
  fi
fi

# Also keep a copy at top-level results (useful for downstream tools)
GLOBAL_COMBINED="${RESULTS%/}/combined_amrfinder.tsv"
cp -f "$COMBINED_TSV" "$GLOBAL_COMBINED"

# ------------------ STEP 4: Merge 'group' from sample.tsv (if present) ------------------
LOCAL_SAMPLE_FILE="${SCRIPT_DIR}/sample.tsv"
if [[ -f "$LOCAL_SAMPLE_FILE" ]]; then
  echo "?? Merging group column from $LOCAL_SAMPLE_FILE into combined tables..."
  "$PYTHON_BIN" - "$COMBINED_TSV" "$GLOBAL_COMBINED" "$LOCAL_SAMPLE_FILE" <<'PY'
import sys, pandas as pd
tool_combined, global_combined, sample_path = sys.argv[1:]

def merge_with_group(infile, samplefile):
    comb = pd.read_csv(infile, sep="\t", dtype=str)
    meta = pd.read_csv(samplefile, sep="\t", dtype=str)
    meta.columns = [c.lower() for c in meta.columns]
    comb.columns = [c.lower() for c in comb.columns]
    if "sample" not in meta.columns or "group" not in meta.columns:
        sys.exit("ERROR: sample.tsv must have 'sample' and 'group' columns.")
    meta_sel = meta[["sample", "group"]].drop_duplicates()
    merged = comb.merge(meta_sel, on="sample", how="left")
    # put 'group' before 'sample'
    cols = merged.columns.tolist()
    if "group" in cols and "sample" in cols:
        cols = ["group","sample"] + [c for c in cols if c not in ("group","sample")]
    merged = merged[cols]
    merged.to_csv(infile, sep="\t", index=False)
    print(f"? Added group column to: {infile}")

merge_with_group(tool_combined, sample_path)
merge_with_group(global_combined, sample_path)
PY
else
  echo "?? sample.tsv not found in ${SCRIPT_DIR}. Skipping group merge."
fi

# ------------------ STEP 5: Generate Excel (+ zero-hit) & Heatmaps ------------------
SUMMARY_XLSX="${RESULTS%/}/OneHealthAMR_AMRFinder_summary.xlsx"

"$PYTHON_BIN" - "$GLOBAL_COMBINED" "$SUMMARY_XLSX" "$AMR_DIR" <<'PY'
import sys, os, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, csv, tempfile
from matplotlib import cm, colors

inpath, out_xlsx, outdir = sys.argv[1:]
os.makedirs(outdir, exist_ok=True)
print("[INFO] Generating blue-white heatmaps with dynamic text contrast...")

# ---- Helpers ----
def clean(infile):
    tmp=tempfile.mkstemp()[1]
    with open(infile,"r",errors="replace") as fin, open(tmp,"w",newline="") as fout:
        w=csv.writer(fout,delimiter="\t")
        keep_header=False
        for row in csv.reader(fin,delimiter="\t"):
            if not row: continue
            if not keep_header:
                w.writerow(row); keep_header=True; continue
            # Drop spurious lines from concatenation
            if row[0].strip().lower() in ["group","sample"]: continue
            if "Protein id" in row or "Contig id" in row: continue
            w.writerow(row)
    return tmp

def build(df, key_col, cat_col, binary=False):
    df[key_col]=df[key_col].astype(str).str.strip()
    df[cat_col]=df[cat_col].fillna("Unknown").astype(str).str.strip()
    mat=pd.crosstab(df[key_col], df[cat_col])
    if binary:
        mat=(mat>0).astype(int)
    return mat.sort_index(axis=0).sort_index(axis=1)

def plot(mat, fname, title, ylabel, xlabel, annotate_counts, mode="counts"):
    samples=list(mat.index)
    cats=list(mat.columns)
    data=mat.T.values.astype(float)
    width=max(6,0.3*len(samples))
    height=max(6,0.25*len(cats))
    fig,ax=plt.subplots(figsize=(width,height))
    cmap=cm.get_cmap("Blues")
    vmin=0.0
    vmax=float(np.nanmax(data)) if np.isfinite(np.nanmax(data)) else 1.0
    if vmax == 0.0: vmax = 1.0   # keep color scale stable for all-zero panels
    norm=colors.Normalize(vmin=vmin, vmax=vmax)
    im=ax.imshow(data,aspect="auto",interpolation="nearest",cmap=cmap,norm=norm)
    if annotate_counts:
        for i in range(len(cats)):
            for j in range(len(samples)):
                val=data[i,j]
                if np.isnan(val) or val==0: continue
                rgba=cmap(norm(val))
                brightness=(rgba[0]*0.299+rgba[1]*0.587+rgba[2]*0.114)
                text_color="black" if brightness>0.5 else "white"
                ax.text(j,i,int(val),ha="center",va="center",color=text_color,fontsize=7,fontweight="bold")
    ax.set_xticks(np.arange(len(samples)))
    ax.set_xticklabels(samples,rotation=90,fontsize=8)
    ax.set_yticks(np.arange(len(cats)))
    ax.set_yticklabels(cats,fontsize=8)
    ax.set_xlabel(xlabel,fontsize=12)
    ax.set_ylabel(ylabel,fontsize=12)
    ax.set_title(title,fontsize=13)
    label="Gene count per category" if mode=="counts" else "Presence (1=yes, 0=no)"
    cbar=plt.colorbar(im,ax=ax,fraction=0.046,pad=0.04)
    cbar.set_label(label)
    for spine in ax.spines.values(): spine.set_visible(False)
    ax.tick_params(length=0)
    plt.tight_layout()
    plt.savefig(fname,dpi=300,facecolor="white")
    plt.close(fig)
    print(f"? Heatmap written: {fname}")

# ---- Load Data ----
tmp=clean(inpath)
df=pd.read_csv(tmp,sep="\t",dtype=str,low_memory=False)

# Column mapping
lower_map={c.lower():c for c in df.columns}
group_col=lower_map.get("group")
sample_col=lower_map.get("sample")
if sample_col is None:
    raise SystemExit("ERROR: Combined table must contain a 'sample' column.")
class_col=lower_map.get("class") or lower_map.get("name") or df.columns[2]
subclass_col=lower_map.get("subclass")
type_col=None
for t in ["type","element_type","feature","product"]:
    if t in lower_map: type_col=lower_map[t]; break
categories={"class":class_col}
if subclass_col: categories["subclass"]=subclass_col
if type_col: categories["type"]=type_col

# ---- Determine ALL processed samples (including zero-hit)
all_samples = sorted([d for d in os.listdir(outdir) if os.path.isdir(os.path.join(outdir, d))])
samples_with_hits = sorted(pd.unique(df[sample_col])) if sample_col in df.columns else []
zero_hit_samples = sorted(set(all_samples) - set(samples_with_hits))
print(f"[INFO] Samples with hits: {len(samples_with_hits)}; zero-hit: {len(zero_hit_samples)}; total: {len(all_samples)}")

# ---- Excel with zero-hit awareness ----
with pd.ExcelWriter(out_xlsx,engine="openpyxl") as w:
    # Raw combined as-is
    df.to_excel(w,"combined_raw",index=False)
    # Explicit list of zero-hit samples
    zdf = pd.DataFrame({"sample": list(zero_hit_samples), "genes_detected": 0})
    zdf.to_excel(w,"zero_hit_samples",index=False)

    # Build all matrices and reindex rows to include every processed sample
    for name,col in categories.items():
        count_mat = build(df, sample_col, col, binary=False).reindex(all_samples, fill_value=0)
        binary_mat = build(df, sample_col, col, binary=True).reindex(all_samples, fill_value=0)
        count_mat.to_excel(w,f"counts_{name}_by_sample")
        binary_mat.to_excel(w,f"binary_{name}_by_sample")

        if group_col:
            count_mat_g = build(df, group_col, col, binary=False)
            binary_mat_g = build(df, group_col, col, binary=True)
            count_mat_g.to_excel(w,f"counts_{name}_by_group")
            binary_mat_g.to_excel(w,f"binary_{name}_by_group")
print(f"? Excel summary written (zeros included): {out_xlsx}")

# ---- Heatmaps (zeros included for sample-level) ----
for name,col in categories.items():
    count_mat = build(df, sample_col, col, binary=False).reindex(all_samples, fill_value=0)
    binary_mat = build(df, sample_col, col, binary=True).reindex(all_samples, fill_value=0)

    plot(count_mat, os.path.join(outdir, f"heatmap_{name}_counts_by_sample.png"),
         f"AMRFinderPlus: {name} counts (by sample)",
         f"Antibiotic {name.capitalize()}", "Samples", True, "counts")
    plot(binary_mat, os.path.join(outdir, f"heatmap_{name}_binary_by_sample.png"),
         f"AMRFinderPlus: {name} presence (by sample)",
         f"Antibiotic {name.capitalize()}", "Samples", False, "binary")

    if group_col:
        count_mat_g = build(df, group_col, col, binary=False)
        binary_mat_g = build(df, group_col, col, binary=True)
        plot(count_mat_g, os.path.join(outdir, f"heatmap_{name}_counts_by_group.png"),
             f"AMRFinderPlus: {name} counts (by group)",
             f"Antibiotic {name.capitalize()}", "Groups", True, "counts")
        plot(binary_mat_g, os.path.join(outdir, f"heatmap_{name}_binary_by_group.png"),
             f"AMRFinderPlus: {name} presence (by group)",
             f"Antibiotic {name.capitalize()}", "Groups", False, "binary")
PY

echo "? Pipeline completed successfully."
echo "?? Excel summary: $SUMMARY_XLSX"
echo "?? Combined TSVs and blue-white heatmaps saved in: $AMR_DIR"
