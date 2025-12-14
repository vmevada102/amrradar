#!/bin/bash
set -euo pipefail

# ==================================================
# MLST Automation with QC, Confidence & Reporting
# ==================================================

echo "[INFO] Starting MLST pipeline..."

# --------------------------------------------------
# Step 1: Conda environment check
# --------------------------------------------------
if [[ -z "${CONDA_DEFAULT_ENV:-}" ]]; then
    echo "[ERROR] No conda environment active."
    exit 1
fi

echo "[INFO] Conda environment: $CONDA_DEFAULT_ENV"
echo "[INFO] Python: $(which python3)"

# --------------------------------------------------
# Step 2: Required packages (skip if installed)
# --------------------------------------------------
REQUIRED_PACKAGES=(mlst csvkit pandas matplotlib seaborn networkx openpyxl reportlab)

for pkg in "${REQUIRED_PACKAGES[@]}"; do
    if conda list "$pkg" >/dev/null 2>&1; then
        echo "[INFO] Package '$pkg' already installed — skipping."
    else
        echo "[INFO] Installing missing package: $pkg"
        conda install -y -c bioconda -c conda-forge "$pkg"
    fi
done

# --------------------------------------------------
# Step 3: Output directory
# --------------------------------------------------
maxN=$(ls -d results/tools/* 2>/dev/null | grep -oE '/[0-9]+' | tr -d '/' | sort -n | tail -1)
nextN=$(( ${maxN:-0} + 1 ))
OUTPUT_DIR="results/tools/${nextN}_MLST"
mkdir -p "$OUTPUT_DIR"

CSV_FILE="$OUTPUT_DIR/mlst_results.csv"
TSV_FILE="$OUTPUT_DIR/mlst_results.tsv"
XLSX_FILE="$OUTPUT_DIR/mlst_results.xlsx"
CMD_FILE="$OUTPUT_DIR/command.txt"

# --------------------------------------------------
# Step 4: Provenance
# --------------------------------------------------
{
    echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "User: $(whoami)"
    echo "Host: $(hostname)"
    echo "Conda env: $CONDA_DEFAULT_ENV"
    echo "Command: $0 $*"
} > "$CMD_FILE"

# --------------------------------------------------
# Step 5: Run MLST
# --------------------------------------------------
echo "Sample,Scheme,ST,Alleles" > "$CSV_FILE"

shopt -s nullglob
ASSEMBLIES=(results/tools/*_assembly/*/*.assembly.fasta)

if [[ ${#ASSEMBLIES[@]} -eq 0 ]]; then
    echo "[ERROR] No assembly FASTA files found."
    exit 1
fi

for fasta in "${ASSEMBLIES[@]}"; do
    sample=$(basename "$fasta" .assembly.fasta)
    echo "[INFO] Processing $sample"
    mlst "$fasta" | awk -v s="$sample" '{print s","$0}' >> "$CSV_FILE"
done

csvformat -T "$CSV_FILE" > "$TSV_FILE"

# --------------------------------------------------
# Step 6: Python QC, confidence, dashboard & PDF
# --------------------------------------------------
python3 <<EOF
import os
import pandas as pd
from statistics import mean

OUT="${OUTPUT_DIR}"

df=pd.read_csv("${CSV_FILE}")
df.to_excel("${XLSX_FILE}",index=False)

amr=pd.read_excel("results/OneHealthAMR_AMRFinder_summary.xlsx",sheet_name="combined_raw")

def species(x):
    if isinstance(x,str):
        p=x.split()
        return " ".join(p[:2]) if len(p)>=2 else None
    return None

amr["organism"]=amr["closest reference name"].apply(species)
spdf=amr[["sample","organism"]].dropna().drop_duplicates().set_index("sample")

VALID={
    "paeruginosa":"Pseudomonas aeruginosa",
    "ecoli":"Escherichia coli",
    "abaumannii":"Acinetobacter baumannii",
    "kpneumoniae":"Klebsiella pneumoniae"
}
REV={v:k for k,v in VALID.items()}

records=[]
scores=[]
qc="PASS"

for i,r in df.iterrows():
    sample=r["Sample"]
    scheme=str(r["Scheme"]).lower()
    sps=set(spdf.loc[spdf.index==sample,"organism"]) if sample in spdf.index else set()
    exp=VALID.get(scheme)

    status="OK"
    action="None"
    score=100

    if not sps:
        score=50
    elif exp and exp in sps and len(sps)==1:
        score=100
    elif exp and exp in sps:
        score=80
    elif exp and exp not in sps:
        qc="WARN"
        corr=None
        for s in sps:
            if s in REV:
                corr=REV[s]; break
        if corr:
            df.at[i,"Scheme"]=corr
            status="AUTO-CORRECTED"
            action=f"Scheme corrected to {corr}"
            score=70
        else:
            status="MISMATCH"
            action="Scheme-species mismatch"
            score=40

    scores.append(score)

    records.append({
        "sample":sample,
        "original_scheme":scheme,
        "final_scheme":df.at[i,"Scheme"],
        "species_inferred":"; ".join(sorted(sps)),
        "status":status,
        "confidence_score":score,
        "action":action
    })

val=pd.DataFrame(records)
val.to_csv(os.path.join(OUT,"mlst_species_validation.csv"),index=False)

with open(os.path.join(OUT,"QC_STATUS.txt"),"w") as f:
    f.write(qc)

avg=mean(scores)
overall="HIGH" if avg>=85 else "MODERATE" if avg>=65 else "LOW"

# ---- Dashboard ----
def color(r):
    return ["background-color:#ffcccc"]*len(r) if r["status"]!="OK" else ["background-color:#ccffcc"]*len(r)

html=f"""
<h1>MLST Dashboard</h1>
<h2>QC Status: <span style="color:{'red' if qc=='WARN' else 'green'}">{qc}</span></h2>
<h2>Overall Confidence: {overall} (Mean score = {avg:.1f})</h2>

<h3>MLST Results</h3>
{df.to_html(index=False)}

<h3>Validation & Confidence</h3>
{val.style.apply(color,axis=1).to_html()}
"""

with open(os.path.join(OUT,"dashboard.html"),"w") as f:
    f.write(html)

# ---- PDF ----
try:
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Table
    from reportlab.lib.styles import getSampleStyleSheet
    pdf=os.path.join(OUT,"MLST_QC_Report.pdf")
    styles=getSampleStyleSheet()
    doc=SimpleDocTemplate(pdf)
    elems=[
        Paragraph("MLST QC & Confidence Report",styles["Title"]),
        Paragraph(f"QC Status: {qc}",styles["Normal"]),
        Paragraph(f"Overall Confidence: {overall} (Mean = {avg:.1f})",styles["Normal"]),
        Table([val.columns.tolist()]+val.values.tolist())
    ]
    doc.build(elems)
except Exception as e:
    print("[WARN] PDF generation skipped:",e)

# ---- AMR dashboard badge ----
amr_html="results/tools/AMRNetwork/dashboard.html"
if os.path.exists(amr_html):
    with open(amr_html) as f: h=f.read()
    badge=f"<div style='padding:10px;border:2px solid {'red' if qc=='WARN' else 'green'};'>MLST QC: {qc} | Confidence: {overall}</div>"
    if "MLST QC" not in h:
        with open(amr_html+".bak","w") as f: f.write(h)
        with open(amr_html,"w") as f: f.write(badge+h)

print("[INFO] MLST QC Status:",qc)
print("[INFO] Overall MLST Confidence:",overall)
EOF

# --------------------------------------------------
# Final messages
# --------------------------------------------------
QC_STATUS=$(cat "$OUTPUT_DIR/QC_STATUS.txt")
echo "[INFO] MLST pipeline completed."
echo "[INFO] QC Status: $QC_STATUS"
echo "[INFO] MLST results generated in directory: $OUTPUT_DIR"
echo "[INFO] Dashboard, validation CSV, and PDF report generated."
