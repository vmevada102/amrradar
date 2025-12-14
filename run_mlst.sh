#!/bin/bash
set -euo pipefail

# ================================
# MLST Automation & Analysis Script
# ================================

echo "[INFO] Starting MLST pipeline..."

# -------------------------------
# Step 1: Check Conda Environment
# -------------------------------
if [[ -z "${CONDA_DEFAULT_ENV:-}" ]]; then
    echo "[ERROR] No conda environment active. Please activate the AMR conda environment."
    exit 1
fi

echo "[INFO] Using Conda environment: $CONDA_DEFAULT_ENV"
echo "[INFO] Python: $(which python3)"

# -------------------------------
# Step 2: Check required conda packages (skip if present)
# -------------------------------
echo "[INFO] Checking required conda packages..."

REQUIRED_PACKAGES=(
    mlst
    csvkit
    pandas
    matplotlib
    seaborn
    networkx
    openpyxl
)

for pkg in "${REQUIRED_PACKAGES[@]}"; do
    if conda list "$pkg" >/dev/null 2>&1; then
        echo "[INFO] Package '$pkg' already installed — skipping."
    else
        echo "[INFO] Package '$pkg' not found — installing."
        conda install -y -c bioconda -c conda-forge "$pkg"
    fi
done


# -------------------------------
# Step 3: Determine Output Folder
# -------------------------------
echo "[INFO] Determining next MLST results directory..."

maxN=$(ls -d results/tools/* 2>/dev/null | grep -oE '/[0-9]+' | tr -d '/' | sort -n | tail -1)

if [[ -z "$maxN" ]]; then
    nextN=1
else
    nextN=$((maxN + 1))
fi

OUTPUT_DIR="results/tools/${nextN}_MLST"
mkdir -p "$OUTPUT_DIR"

CSV_FILE="$OUTPUT_DIR/mlst_results.csv"
TSV_FILE="$OUTPUT_DIR/mlst_results.tsv"
XLSX_FILE="$OUTPUT_DIR/mlst_results.xlsx"
CMD_FILE="$OUTPUT_DIR/command.txt"

echo "[INFO] Results directory: $OUTPUT_DIR"

# -------------------------------
# Step 4: Record Provenance
# -------------------------------
{
    echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "User: $(whoami)"
    echo "Host: $(hostname)"
    echo "Conda environment: $CONDA_DEFAULT_ENV"
    echo "Python: $(which python3)"
    echo "Command: $0 $*"
} > "$CMD_FILE"

# -------------------------------
# Step 5: Run MLST
# -------------------------------
echo "[INFO] Running MLST typing..."
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
    result=$(mlst "$fasta")
    echo "$sample,$result" >> "$CSV_FILE"
done

# -------------------------------
# Step 6: CSV → TSV
# -------------------------------
csvformat -T "$CSV_FILE" > "$TSV_FILE"

# -------------------------------
# Step 7: Python Analysis + Plots + Dashboard
# -------------------------------
python3 <<EOF
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx

OUTPUT_DIR = "${OUTPUT_DIR}"
CSV_FILE = "${CSV_FILE}"
XLSX_FILE = "${XLSX_FILE}"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ------------------------
# Load data
# ------------------------
df = pd.read_csv(CSV_FILE)
df.to_excel(XLSX_FILE, index=False)

# ------------------------
# Summary statistics
# ------------------------
summary_file = os.path.join(OUTPUT_DIR, "summary.txt")
with open(summary_file, "w") as f:
    f.write("MLST Summary Statistics\n")
    f.write("=======================\n")
    f.write(f"Total samples: {len(df)}\n")
    f.write(f"Unique schemes: {df['Scheme'].nunique()}\n")
    f.write(f"Unique STs: {df['ST'].nunique()}\n\n")
    f.write("ST distribution:\n")
    f.write(str(df['ST'].value_counts()))

# ------------------------
# ST Distribution Plot
# ------------------------
plt.figure(figsize=(10,6))
sns.countplot(x="ST", data=df, order=df['ST'].value_counts().index)
plt.xticks(rotation=90)
plt.title("Sequence Type Distribution")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "st_distribution.png"))
plt.close()

# ------------------------
# Organism summary
# ------------------------
organism_file = os.path.join(OUTPUT_DIR, "organism_summary.txt")
with open(organism_file, "w") as f:
    f.write("Organism Detection Summary\n")
    f.write("==========================\n\n")
    for _, row in df.iterrows():
        f.write(f"{row['Sample']} → {row['Scheme']} (ST{row['ST']})\n")

# ------------------------
# Minimum Spanning Tree
# ------------------------
allele_cols = [c for c in df.columns if c.startswith("allele")]
mst_png = os.path.join(OUTPUT_DIR, "mst.png")

if allele_cols:
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_node(row['Sample'], st=row['ST'])

    profiles = df[allele_cols].astype(str).values.tolist()
    samples = df['Sample'].tolist()

    for i in range(len(samples)):
        for j in range(i+1, len(samples)):
            dist = sum(a != b for a, b in zip(profiles[i], profiles[j]))
            G.add_edge(samples[i], samples[j], weight=dist)

    pos = nx.spring_layout(G, seed=42)
    plt.figure(figsize=(8,8))
    nx.draw(G, pos, with_labels=True, node_color="lightblue", font_size=8)
    plt.title("Minimum Spanning Tree (Allelic Profiles)")
    plt.savefig(mst_png)
    plt.close()

# ------------------------
# HTML Dashboard
# ------------------------
with open(summary_file) as f:
    summary_text = f.read()
with open(organism_file) as f:
    organism_text = f.read()

table_html = df.to_html(index=False)

dashboard_file = os.path.join(OUTPUT_DIR, "dashboard.html")

html = f"""
<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>MLST Dashboard</title>
<style>
body {{ font-family: Arial; margin: 20px; }}
h1 {{ color: #2c3e50; }}
table {{ border-collapse: collapse; width: 100%; }}
th, td {{ border: 1px solid #ccc; padding: 6px; }}
img {{ max-width: 900px; margin: 20px 0; }}
</style>
</head>
<body>
<h1>MLST Analysis Dashboard</h1>

<h2>Summary</h2>
<pre>{summary_text}</pre>

<h2>Organism Detection</h2>
<pre>{organism_text}</pre>

<h2>MLST Results</h2>
{table_html}

<h2>ST Distribution</h2>
<img src="st_distribution.png">

<h2>Minimum Spanning Tree</h2>
<img src="mst.png">

</body>
</html>
"""

with open(dashboard_file, "w") as f:
    f.write(html)

print("[INFO] Downstream analyses complete.")
print("[INFO] Dashboard generated:", dashboard_file)
EOF

# -------------------------------
# Step 8: Final Report
# -------------------------------
echo "[INFO] MLST analysis completed successfully."
echo "[INFO] Outputs:"
echo " - $CSV_FILE"
echo " - $TSV_FILE"
echo " - $XLSX_FILE"
echo " - $OUTPUT_DIR/summary.txt"
echo " - $OUTPUT_DIR/organism_summary.txt"
echo " - $OUTPUT_DIR/st_distribution.png"
echo " - $OUTPUT_DIR/mst.png (if available)"
echo " - $OUTPUT_DIR/dashboard.html"
echo " - $CMD_FILE"
