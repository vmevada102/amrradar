#!/usr/bin/env python3
"""
make_fastani_resistome_charts.py

Combine FastANI and resistome TPM results to:
  1) Build an ANI-based clustermap with clade assignments
  2) Build per-sample resistome TPM heatmap
  3) Build clade-mean resistome TPM heatmap

Expected inputs:

  FastANI:
    results/tools/*_comparativegenomics/fastani/fastani.tsv

  Resistome TPM:
    results/tools/*_resistome_extended/amr_quant_tpm.tsv

Outputs (written to results/Charts/):

  - ANI_clustermap_with_clades.png
  - clade_assignments.tsv
  - Resistome_TPM_heatmap_samples.png
  - Resistome_TPM_heatmap_clades.png
"""

import os
import sys
import subprocess

# ---------------------------------------------------------------
# 1. Check required Python packages and install via conda if needed
# ---------------------------------------------------------------
REQUIRED_PKGS = ["pandas", "numpy", "scipy", "matplotlib", "seaborn"]


def check_and_install(packages):
    missing = []
    for pkg in packages:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)

    if not missing:
        print("[INFO] All required Python packages are available.")
        return

    print("\n[WARN] Missing Python packages detected:", ", ".join(missing))

    # Check if conda is available
    try:
        subprocess.run(["conda", "--version"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        print("[ERROR] 'conda' command not found in PATH.")
        print("        Please install the missing packages manually in your environment:")
        print("        ", " ".join(missing))
        sys.exit(1)

    ans = input("Install them now using conda? [y/N]: ").strip().lower()
    if ans != "y":
        print("[ERROR] Cannot continue without required packages. Exiting.")
        sys.exit(1)

    print("[INFO] Installing missing packages via conda:", ", ".join(missing))
    cmd = ["conda", "install", "-y", "-c", "conda-forge"] + missing
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        print("[ERROR] Conda installation failed:", e)
        sys.exit(1)

    print("[INFO] Installation complete. Restarting script with same arguments...")
    os.execv(sys.executable, [sys.executable] + sys.argv)


# Run package check before importing heavy libs
check_and_install(REQUIRED_PKGS)

# ---------------------------------------------------------------
# 2. Now safely import scientific stack
# ---------------------------------------------------------------
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy as sch
from scipy.spatial import distance as ssd

sns.set(context="notebook", style="white")

# ---------------------------------------------------------------------
# Helper: ensure directory for charts
# ---------------------------------------------------------------------
CHART_DIR = os.path.join("results", "Charts")
os.makedirs(CHART_DIR, exist_ok=True)
print("[INFO] Output charts will be written to:", CHART_DIR)

# ---------------------------------------------------------------------
# 3. Locate comparative genomics (FastANI) run
# ---------------------------------------------------------------------
print("[INFO] Searching for comparative genomics (FastANI) output...")
cg_dirs = sorted(glob.glob("results/tools/*_comparativegenomics"))

if not cg_dirs:
    raise SystemExit("[ERROR] No results/tools/*_comparativegenomics directory found")

cg_run = cg_dirs[-1]  # most recent
print(f"[INFO] Using comparative genomics run: {cg_run}")

FASTANI_FILE = os.path.join(cg_run, "fastani", "fastani.tsv")
if not os.path.isfile(FASTANI_FILE):
    raise SystemExit(f"[ERROR] FastANI file not found: {FASTANI_FILE}")

print(f"[INFO] Reading FastANI from: {FASTANI_FILE}")
with open(FASTANI_FILE) as fh:
    first_line = fh.readline().strip().split("\t")

# Typical fastANI columns: query, ref, ANI, matches, fragments
if first_line[0].lower() in ("query", "reference", "ref"):
    df_ani = pd.read_csv(FASTANI_FILE, sep="\t")
else:
    df_ani = pd.read_csv(
        FASTANI_FILE,
        sep="\t",
        header=None,
        names=["query", "ref", "ani", "matches", "fragments"],
    )


def sample_from_path(p: str) -> str:
    """Extract sample ID from file path by removing directory + fasta extensions."""
    base = os.path.basename(p)
    for ext in (".fasta", ".fa", ".fna", ".fas"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base


df_ani["query_sample"] = df_ani["query"].apply(sample_from_path)
df_ani["ref_sample"]   = df_ani["ref"].apply(sample_from_path)

samples = sorted(set(df_ani["query_sample"]).union(df_ani["ref_sample"]))
print(f"[INFO] Number of samples in FastANI: {len(samples)}")

# Build symmetric ANI matrix
ani_mat = pd.DataFrame(
    np.nan, index=samples, columns=samples, dtype=float
)

for _, row in df_ani.iterrows():
    qs = row["query_sample"]
    rs = row["ref_sample"]
    ani = row["ani"]
    if qs in ani_mat.index and rs in ani_mat.columns:
        ani_mat.loc[qs, rs] = ani
        ani_mat.loc[rs, qs] = ani

# Fill diagonal with 100 and any missing with min observed or 80
for s in samples:
    ani_mat.loc[s, s] = 100.0

if ani_mat.isna().values.any():
    observed = ani_mat.values[~np.isnan(ani_mat.values)]
    fill_value = float(observed.min()) if observed.size > 0 else 80.0
    ani_mat = ani_mat.fillna(fill_value)

# ---------------------------------------------------------------------
# 4. Cluster by ANI and define clades
# ---------------------------------------------------------------------
print("[INFO] Performing hierarchical clustering on ANI...")

# Distance = 100 - ANI
dist_array = ssd.squareform(100.0 - ani_mat.values)
Z = sch.linkage(dist_array, method="average")

# Define clades using ANI >= 95% threshold => distance <= 5
clade_labels = sch.fcluster(Z, t=5.0, criterion="distance")

clade_assignments = pd.DataFrame(
    {"sample": samples, "clade": clade_labels.astype(int)}
).sort_values(["clade", "sample"])

clade_tsv = os.path.join(CHART_DIR, "clade_assignments.tsv")
clade_assignments.to_csv(clade_tsv, sep="\t", index=False)
print(f"[INFO] Wrote clade assignments: {clade_tsv}")

# ---------------------------------------------------------------------
# 5. Plot ANI clustermap with clades
# ---------------------------------------------------------------------
print("[INFO] Building ANI clustermap...")

# Reorder matrix according to dendrogram leaves
dendro = sch.dendrogram(Z, no_plot=True, labels=samples)
ordered_samples = dendro["ivl"]
ani_ordered = ani_mat.loc[ordered_samples, ordered_samples]

g = sns.clustermap(
    ani_mat,
    row_linkage=Z,
    col_linkage=Z,
    cmap="viridis",
    linewidths=0.0,
    xticklabels=False,
    yticklabels=False,
)
g.fig.suptitle("ANI clustermap", y=1.02)

ani_plot = os.path.join(CHART_DIR, "ANI_clustermap_with_clades.png")
g.savefig(ani_plot, dpi=300, bbox_inches="tight")
plt.close("all")
print(f"[INFO] Saved ANI clustermap to: {ani_plot}")

# ---------------------------------------------------------------------
# 6. Locate resistome_extended results
# ---------------------------------------------------------------------
print("[INFO] Searching for resistome_extended output...")
resistome_dirs = sorted(glob.glob("results/tools/*_resistome_extended"))

if not resistome_dirs:
    raise SystemExit("[ERROR] No results/tools/*_resistome_extended directory found")

resistome_run = resistome_dirs[-1]  # most recent
print(f"[INFO] Using resistome run: {resistome_run}")

TPM_FILE = os.path.join(resistome_run, "amr_quant_tpm.tsv")
if not os.path.isfile(TPM_FILE):
    raise SystemExit(f"[ERROR] Cannot find TPM file: {TPM_FILE}")

print(f"[INFO] Reading TPM table from: {TPM_FILE}")
df_tpm = pd.read_csv(TPM_FILE, sep="\t")

required_cols = {"sample", "gene", "TPM"}
missing_cols = required_cols - set(df_tpm.columns)
if missing_cols:
    raise SystemExit(
        f"[ERROR] TPM file {TPM_FILE} is missing required columns: {missing_cols}"
    )

# ---------------------------------------------------------------------
# 7. Build per-sample TPM matrix (genes x samples)
# ---------------------------------------------------------------------
print("[INFO] Building TPM matrix (genes x samples)...")
tpm_matrix = df_tpm.pivot_table(
    index="gene", columns="sample", values="TPM", aggfunc="sum"
).fillna(0.0)

# Restrict to samples present in ANI (intersection)
common_samples = sorted(set(tpm_matrix.columns) & set(samples))
if not common_samples:
    raise SystemExit(
        "[ERROR] No overlapping samples between FastANI and TPM tables. "
        "Check sample naming consistency."
    )

# Reorder columns to follow ANI dendrogram order
ordered_common_samples = [s for s in ordered_samples if s in common_samples]
tpm_matrix = tpm_matrix[ordered_common_samples]

# Log-transform TPM for visualization
tpm_log = np.log10(tpm_matrix + 1.0)

# ---------------------------------------------------------------------
# 8. Per-sample resistome heatmap
# ---------------------------------------------------------------------
print("[INFO] Plotting per-sample resistome TPM heatmap...")
g2 = sns.clustermap(
    tpm_log,
    row_cluster=True,
    col_cluster=False,  # keep ANI order on columns
    cmap="mako",
    xticklabels=False,
    yticklabels=False,
)
g2.fig.suptitle("Resistome TPM (log10(TPM+1)) - samples", y=1.02)

sample_heat = os.path.join(CHART_DIR, "Resistome_TPM_heatmap_samples.png")
g2.savefig(sample_heat, dpi=300, bbox_inches="tight")
plt.close("all")
print(f"[INFO] Saved sample-level resistome heatmap to: {sample_heat}")

# ---------------------------------------------------------------------
# 9. Clade-level resistome heatmap (mean TPM per clade)
# ---------------------------------------------------------------------
print("[INFO] Aggregating TPM by clade...")
sample_to_clade = dict(
    zip(clade_assignments["sample"], clade_assignments["clade"])
)

# Keep only samples present in resistome TPM
sample_to_clade = {
    s: c for s, c in sample_to_clade.items() if s in tpm_matrix.columns
}
if not sample_to_clade:
    raise SystemExit(
        "[ERROR] No overlapping samples between clade assignments and TPM matrix."
    )

df_tpm_long = tpm_matrix.stack().reset_index()
df_tpm_long.columns = ["gene", "sample", "TPM"]
df_tpm_long["clade"] = df_tpm_long["sample"].map(sample_to_clade)
df_tpm_long = df_tpm_long.dropna(subset=["clade"])
df_tpm_long["clade"] = df_tpm_long["clade"].astype(int)

# Mean TPM per gene per clade
clade_matrix = df_tpm_long.pivot_table(
    index="gene", columns="clade", values="TPM", aggfunc="mean"
).fillna(0.0)

# Log-transform
clade_log = np.log10(clade_matrix + 1.0)

# Order clades numerically
clades_sorted = sorted(clade_log.columns)
clade_log = clade_log[clades_sorted]

print("[INFO] Plotting clade-level resistome TPM heatmap...")
g3 = sns.clustermap(
    clade_log,
    row_cluster=True,
    col_cluster=True,
    cmap="rocket_r",
    xticklabels=True,
    yticklabels=False,
)
g3.fig.suptitle("Resistome TPM (log10(TPM+1)) - clade means", y=1.02)

clade_heat = os.path.join(CHART_DIR, "Resistome_TPM_heatmap_clades.png")
g3.savefig(clade_heat, dpi=300, bbox_inches="tight")
plt.close("all")
print(f"[INFO] Saved clade-level resistome heatmap to: {clade_heat}")

print("[INFO] All charts written to:", CHART_DIR)
