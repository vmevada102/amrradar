#!/bin/bash

# ================================
# MLST Automation & Analysis Script
# ================================
# Requirements:
#   - Conda environment active
#   - mlst (https://github.com/tseemann/mlst)
#   - csvkit (CSV/TSV handling)
#   - pandas, matplotlib, seaborn, networkx (Python analysis/plots)
# ================================

# Step 1: Check and install required tools
echo "[INFO] Checking required conda packages..."
REQUIRED_PACKAGES=("mlst" "csvkit" "pandas" "matplotlib" "seaborn" "networkx")

for pkg in "${REQUIRED_PACKAGES[@]}"; do
    if ! conda list | grep -q "^$pkg"; then
        echo "[INFO] Installing missing package: $pkg"
        conda install -y -c bioconda $pkg || conda install -y -c conda-forge $pkg
    else
        echo "[INFO] Package $pkg already installed."
    fi
done

# Step 2: Determine next output directory
echo "[INFO] Determining next MLST results directory..."
maxN=$(ls -d results/tools/*/ 2>/dev/null | grep -oE '[0-9]+_' | sed 's/_//' | sort -n | tail -1)

if [ -z "$maxN" ]; then
    nextN=1
else
    nextN=$((maxN + 1))
fi

OUTPUT_DIR="results/tools/${nextN}_mlst"
mkdir -p "$OUTPUT_DIR"

CSV_FILE="$OUTPUT_DIR/mlst_results.csv"
TSV_FILE="$OUTPUT_DIR/mlst_results.tsv"
XLSX_FILE="$OUTPUT_DIR/mlst_results.xlsx"
CMD_FILE="$OUTPUT_DIR/command.txt"

echo "[INFO] Results will be saved in: $OUTPUT_DIR"

# Step 3: Store the command used to run this script with timestamp and conda env
{
    echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Conda Environment: ${CONDA_DEFAULT_ENV:-unknown}"
    echo "Command: $0 $@"
    echo "Host: $(hostname)"
    echo "User: $(whoami)"
} > "$CMD_FILE"

# Step 4: Run MLST on all assemblies
echo "[INFO] Running MLST typing..."
echo "Sample,Scheme,ST,Alleles" > "$CSV_FILE"

for fasta in results/tools/*_assembly/*/*.assembly.fasta; do
    sample=$(basename "$fasta" .assembly.fasta)
    echo "[INFO] Processing $sample"
    result=$(mlst "$fasta")
    echo "$sample,$result" >> "$CSV_FILE"
done

# Step 5: Convert CSV ? TSV
csvformat -T "$CSV_FILE" > "$TSV_FILE"

# Step 6: Convert CSV ? Excel and run downstream analyses
python3 - <<EOF
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx

OUTPUT_DIR = "${OUTPUT_DIR}"
CSV_FILE = "${CSV_FILE}"
XLSX_FILE = "${XLSX_FILE}"

# Load MLST results
df = pd.read_csv(CSV_FILE)

# Save Excel
df.to_excel(XLSX_FILE, index=False)

# ================================
# Downstream Analyses
# ================================

# 1. Summary statistics
summary_file = os.path.join(OUTPUT_DIR, "summary.txt")
with open(summary_file, "w") as f:
    f.write("MLST Summary Statistics\\n")
    f.write("=======================\\n")
    f.write(f"Total samples: {len(df)}\\n")
    f.write(f"Unique schemes: {df['Scheme'].nunique()}\\n")
    f.write(f"Unique STs: {df['ST'].nunique()}\\n")
    f.write("\\nST distribution:\\n")
    f.write(str(df['ST'].value_counts()))

# 2. Diversity plot
plt.figure(figsize=(10,6))
sns.countplot(x="ST", data=df, order=df['ST'].value_counts().index)
plt.xticks(rotation=90)
plt.title("Sequence Type Distribution")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "st_distribution.png"))

# 3. Minimum Spanning Tree
allele_cols = [c for c in df.columns if c.startswith("allele")]
if allele_cols:
    G = nx.Graph()
    for i, row in df.iterrows():
        G.add_node(row['Sample'], st=row['ST'])
    samples = df['Sample'].tolist()
    profiles = df[allele_cols].astype(str).values.tolist()
    for i in range(len(samples)):
        for j in range(i+1, len(samples)):
            dist = sum(a!=b for a,b in zip(profiles[i], profiles[j]))
            G.add_edge(samples[i], samples[j], weight=dist)
    pos = nx.spring_layout(G)
    plt.figure(figsize=(8,8))
    nx.draw(G, pos, with_labels=True, node_color="lightblue", font_size=8)
    plt.title("Minimum Spanning Tree (Allelic Profiles)")
    plt.savefig(os.path.join(OUTPUT_DIR, "mst.png"))

# 4. Organism summary
organism_summary_file = os.path.join(OUTPUT_DIR, "organism_summary.txt")
with open(organism_summary_file, "w") as f:
    f.write("Organism Detection Summary\\n")
    f.write("==========================\\n\\n")
    for i, row in df.iterrows():
        organism = row['Scheme']
        st = row['ST']
        sample = row['Sample']
        f.write(f"{sample} ? Organism: {organism}, ST{st}\\n")

print("[INFO] Downstream analyses complete.")
EOF



# 5. Combined HTML dashboard
dashboard_file = os.path.join(OUTPUT_DIR, "dashboard.html")

html = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>MLST Dashboard</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
        table, th, td { border: 1px solid #ccc; }
        th, td { padding: 8px; text-align: left; }
        img { max-width: 800px; margin: 20px 0; }
        .section { margin-bottom: 40px; }
    </style>
</head>
<body>
    <h1>MLST Analysis Dashboard</h1>

    <div class="section">
        <h2>Summary Statistics</h2>
        <pre>{summary}</pre>
    </div>

    <div class="section">
        <h2>Organism Detection Summary</h2>
        <pre>{organism_summary}</pre>
    </div>

    <div class="section">
        <h2>MLST Results Table</h2>
        {table}
    </div>

    <div class="section">
        <h2>Sequence Type Distribution</h2>
        <img src="st_distribution.png" alt="ST Distribution">
    </div>

    <div class="section">
        <h2>Minimum Spanning Tree</h2>
        <img src="mst.png" alt="MST">
    </div>
</body>
</html>
"""

# Load summaries
with open(os.path.join(OUTPUT_DIR, "summary.txt")) as f:
    summary_text = f.read()
with open(os.path.join(OUTPUT_DIR, "organism_summary.txt")) as f:
    organism_text = f.read()

# Convert MLST results to HTML table
table_html = df.to_html(index=False)

# Write dashboard
with open(dashboard_file, "w") as f:
    f.write(html.format(summary=summary_text,
                        organism_summary=organism_text,
                        table=table_html))

print(f"[INFO] Dashboard generated: {dashboard_file}")


echo "[INFO] All results saved in $OUTPUT_DIR:"
echo " - $CSV_FILE"
echo " - $TSV_FILE"
echo " - $XLSX_FILE"
echo " - $CMD_FILE"
echo " - $OUTPUT_DIR/summary.txt"
echo " - $OUTPUT_DIR/st_distribution.png"
echo " - $OUTPUT_DIR/mst.png (if allelic profiles available)"
echo " - $OUTPUT_DIR/organism_summary.txt"
