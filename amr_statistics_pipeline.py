#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
AMRFinder statistics & report generator (publication-ready).

- Reads sample information from sample.tsv (sample, r1, r2, group).
- Finds all AMRFinder result files:
    results/tools/*_AMRFinder/<sample>/<sample>.amrfinder.tsv
- Merges them into one table with columns: group, sample, <AMRFinder columns>
- Computes summary tables and statistical tests (group-wise).
- Generates:
    * Gene presence heatmap (samples x genes)
    * Barplots (counts and percentages) for:
        type, subtype, class, subclass, method
    * Group-wise stacked bars (counts) and 100% stacked bars (percentages)
    * Category-level heatmaps (group x category, percentage)
    * Bubble charts (Krona-like overview, group x category)
    * PCA on gene presence/absence
    * Welch t-tests on total AMR hits between groups
- Writes everything into results/tools/N+1_AMR_Statistics
  (where N is the maximum existing N_* under results/tools).
- Creates:
    * AMR_statistics_report.html
    * interpretation.txt (text explanation of all analyses)

Requirements (install in your conda env):
  pandas, numpy, matplotlib, seaborn, scipy
"""

import argparse
import glob
import os
import re
from pathlib import Path
from typing import List, Optional, Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, fisher_exact, ttest_ind


# === CONFIG: adjust these if your AMRFinder columns are named differently ===
AMR_COLUMNS = {
    "gene_symbol": "gene_symbol",     # or "Gene symbol" / "Name"
    "type": "element_type",           # e.g. "Element type" / "type"
    "subtype": "element_subtype",     # e.g. "Element subtype" / "subtype"
    "class": "class",                 # e.g. "Class"
    "subclass": "subclass",           # e.g. "Subclass"
    "method": "method"                # e.g. "Method"
}
ID_COLUMN = "sequence_name"          # fallback ID for loci; change if needed

CATEGORY_LEVELS = ["type", "subtype", "class", "subclass", "method"]


def set_publication_style():
    """
    Global matplotlib style for publication-ready figures.
    Small, consistent fonts + vector-friendly settings.
    """
    mpl.rcParams.update({
        "figure.dpi": 300,
        "savefig.dpi": 600,
        "font.size": 8,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "legend.title_fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.autolayout": True,
    })
    sns.set_style("whitegrid")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine AMRFinder results and generate statistics and HTML report."
    )
    parser.add_argument(
        "--samples",
        default="sample.tsv",
        help="Sample sheet with columns: sample, r1, r2, group (default: sample.tsv)",
    )
    parser.add_argument(
        "--results-dir",
        default="results/tools",
        help="Base results tools directory (default: results/tools)",
    )
    parser.add_argument(
        "--amr-pattern",
        default="*_AMRFinder/*/*.amrfinder.tsv",
        help="Glob pattern under results-dir to find AMR TSVs "
             "(default: '*_AMRFinder/*/*.amrfinder.tsv')",
    )
    parser.add_argument(
        "--report-name",
        default="AMR_statistics_report.html",
        help="Name of HTML report file (inside new step dir)",
    )
    return parser.parse_args()


def find_next_step_dir(results_dir: str, suffix: str = "AMR_Statistics") -> Path:
    """
    Scan results_dir for N_* directories and create (N+1)_<suffix>.
    Also creates plots/ and tables/ subdirectories.
    """
    results_path = Path(results_dir)
    results_path.mkdir(parents=True, exist_ok=True)

    max_n = 0
    for d in results_path.iterdir():
        if d.is_dir() and re.match(r"^\d+_", d.name):
            try:
                n = int(d.name.split("_", 1)[0])
                if n > max_n:
                    max_n = n
            except ValueError:
                continue

    next_n = max_n + 1
    step_dir = results_path / f"{next_n}_{suffix}"
    (step_dir / "plots").mkdir(parents=True, exist_ok=True)
    (step_dir / "tables").mkdir(parents=True, exist_ok=True)

    return step_dir


def load_samples(sample_path: str) -> pd.DataFrame:
    """
    Load sample.tsv and check required columns.
    """
    df = pd.read_csv(sample_path, sep="\t", dtype=str)
    required = {"sample", "group"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError("sample.tsv must contain columns: {}. Missing: {}".format(required, missing))
    return df


def find_amr_files(results_dir: str, amr_pattern: str) -> List[Path]:
    pattern = str(Path(results_dir) / amr_pattern)
    files = [Path(p) for p in glob.glob(pattern)]
    if not files:
        raise FileNotFoundError("No AMRFinder TSV files found with pattern: {}".format(pattern))
    return files


def safe_get_col(df: pd.DataFrame, logical_name: str) -> Optional[str]:
    """
    Return actual column name for a logical AMR column (type, class, etc) if present.

    Matching strategy:
    1) Use AMR_COLUMNS mapping directly.
    2) Case-insensitive exact match.
    3) Case-insensitive using logical_name itself.
    4) Canonical form: remove spaces, underscores, hyphens and lowercase.
    """
    cols = list(df.columns)

    wanted = AMR_COLUMNS.get(logical_name)

    # direct mapping
    if wanted and wanted in cols:
        return wanted

    lower_map = {c.lower(): c for c in cols}

    if wanted and wanted.lower() in lower_map:
        return lower_map[wanted.lower()]

    if logical_name in cols:
        return logical_name
    if logical_name.lower() in lower_map:
        return lower_map[logical_name.lower()]

    def canon(s: str) -> str:
        return re.sub(r"[^a-z0-9]", "", s.lower())

    canon_map = {canon(c): c for c in cols}

    for cand in (wanted, logical_name):
        if not cand:
            continue
        key = canon(cand)
        if key in canon_map:
            return canon_map[key]

    print("[WARN] Could not find column for logical name '{}' in {}".format(logical_name, cols))
    return None


def merge_amr_files(amr_files: List[Path], sample_df: pd.DataFrame) -> pd.DataFrame:
    sample_to_group: Dict[str, str] = dict(zip(sample_df["sample"], sample_df["group"]))
    all_rows = []

    for f in amr_files:
        sample_name = f.stem  # "<sample>.amrfinder"
        if sample_name.endswith(".amrfinder"):
            sample_name = sample_name.replace(".amrfinder", "")

        group = sample_to_group.get(sample_name, "NA")

        df = pd.read_csv(f, sep="\t", dtype=str)
        df.insert(0, "sample", sample_name)
        df.insert(0, "group", group)
        all_rows.append(df)

    if not all_rows:
        raise ValueError("No AMR rows collected from AMRFinder TSV files.")

    merged = pd.concat(all_rows, ignore_index=True)
    merged = merged.fillna("")
    return merged


def write_table(df: pd.DataFrame, out_dir: Path, basename: str) -> Tuple[Path, Path]:
    tsv = out_dir / "tables" / f"{basename}.tsv"
    csv = out_dir / "tables" / f"{basename}.csv"
    df.to_csv(tsv, sep="\t", index=False)
    df.to_csv(csv, index=False)
    return tsv, csv


def plot_gene_presence_heatmap(merged: pd.DataFrame, out_dir: Path) -> Tuple[Path, Path]:
    gene_col = safe_get_col(merged, "gene_symbol") or ID_COLUMN
    if gene_col not in merged.columns:
        raise ValueError(
            "Neither gene_symbol nor {} found in merged table. Cannot make gene heatmap.".format(ID_COLUMN)
        )

    presence = (
        merged.groupby(["sample", gene_col])
        .size()
        .reset_index(name="count")
    )
    presence["present"] = 1

    mat = presence.pivot_table(
        index="sample",
        columns=gene_col,
        values="present",
        fill_value=0,
        aggfunc="max",
    )

    plt.figure(figsize=(max(6, mat.shape[1] * 0.3), max(4, mat.shape[0] * 0.3)))
    sns.heatmap(mat, cmap="viridis", cbar_kws={"label": "Presence (1/0)"})
    plt.title("AMR gene presence/absence across samples")
    plt.xlabel("Gene")
    plt.ylabel("Sample")
    png = out_dir / "plots" / "gene_presence_heatmap.png"
    pdf = out_dir / "plots" / "gene_presence_heatmap.pdf"
    plt.tight_layout()
    plt.savefig(png, dpi=600)
    plt.savefig(pdf)
    plt.close()

    return png, pdf


def plot_category_bars(
    merged: pd.DataFrame, out_dir: Path, level: str
) -> Tuple[Optional[Path], Optional[Path], Optional[Path], Optional[Path], Optional[Path]]:
    """
    For a given category level (type, subtype, class, subclass, method):
    - Plot overall barplot (counts).
    - Plot overall barplot (percentages).
    - Plot group-wise stacked bar (counts).
    - Plot group-wise 100% stacked bar (percentages).
    - Write count/percent table.

    Returns:
      (overall_counts_png,
       overall_percent_png,
       group_counts_png,
       group_percent_png,
       main_table_path)
    """
    col = safe_get_col(merged, level)
    if col is None:
        return None, None, None, None, None

    total_samples = merged["sample"].nunique()

    counts_all = (
        merged.groupby(col)["sample"]
        .nunique()
        .reset_index(name="n_samples_with_category")
    )
    counts_all["percent_samples"] = (
        counts_all["n_samples_with_category"] * 100.0 / float(total_samples)
    )
    counts_all = counts_all.sort_values("n_samples_with_category", ascending=False)

    counts_group = (
        merged.groupby(["group", col])["sample"]
        .nunique()
        .reset_index(name="n_samples_with_category")
    )
    group_sizes = (
        merged[["sample", "group"]]
        .drop_duplicates()
        .groupby("group")["sample"]
        .nunique()
        .rename("n_samples_in_group")
        .reset_index()
    )
    counts_group = counts_group.merge(group_sizes, on="group", how="left")
    counts_group["percent_samples"] = (
        counts_group["n_samples_with_category"] * 100.0
        / counts_group["n_samples_in_group"].astype(float)
    )

    # Overall barplot (counts)
    plt.figure(figsize=(10, 4 + len(counts_all) * 0.15))
    sns.barplot(data=counts_all, x="n_samples_with_category", y=col)
    plt.title(f"{level.capitalize()} - sample counts", fontsize=8)
    plt.xlabel("Count", fontsize=8)
    plt.ylabel(level.capitalize(), fontsize=8)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.tight_layout()
    png_all_counts = out_dir / "plots" / f"{level}_overall_bar_counts.png"
    pdf_all_counts = out_dir / "plots" / f"{level}_overall_bar_counts.pdf"
    plt.savefig(png_all_counts, dpi=600)
    plt.savefig(pdf_all_counts)
    plt.close()

    # Overall barplot (percentages)
    plt.figure(figsize=(10, 4 + len(counts_all) * 0.15))
    sns.barplot(data=counts_all, x="percent_samples", y=col)
    plt.title(f"{level.capitalize()} - sample percentages", fontsize=8)
    plt.xlabel("Percentage (%)", fontsize=8)
    plt.ylabel(level.capitalize(), fontsize=8)
    plt.xlim(0, 100)
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    plt.tight_layout()
    png_all_percent = out_dir / "plots" / f"{level}_overall_bar_percent.png"
    pdf_all_percent = out_dir / "plots" / f"{level}_overall_bar_percent.pdf"
    plt.savefig(png_all_percent, dpi=600)
    plt.savefig(pdf_all_percent)
    plt.close()

    # Group-wise stacked bar (counts)
    pivot_counts = counts_group.pivot_table(
        index="group",
        columns=col,
        values="n_samples_with_category",
        fill_value=0,
    )
    plt.figure(figsize=(max(10, pivot_counts.shape[1] * 0.3), 6))
    ax = pivot_counts.plot(kind="bar", stacked=True)
    plt.title(f"{level.capitalize()} - stacked counts", fontsize=8)
    plt.ylabel("Count", fontsize=8)
    plt.xlabel("Group", fontsize=8)
    plt.xticks(rotation=45, ha="right", fontsize=7)
    plt.yticks(fontsize=7)
    plt.legend(
        title=level.capitalize(),
        fontsize=6,
        title_fontsize=7,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )
    plt.tight_layout()
    png_group_counts = out_dir / "plots" / f"{level}_group_stacked_bar_counts.png"
    pdf_group_counts = out_dir / "plots" / f"{level}_group_stacked_bar_counts.pdf"
    plt.savefig(png_group_counts, dpi=600)
    plt.savefig(pdf_group_counts)
    plt.close()

    # Group-wise stacked bar (percentages; 100% stacked)
    pivot_percent = counts_group.pivot_table(
        index="group",
        columns=col,
        values="percent_samples",
        fill_value=0,
    )
    plt.figure(figsize=(max(10, pivot_percent.shape[1] * 0.3), 6))
    ax = pivot_percent.plot(kind="bar", stacked=True)
    plt.title(f"{level.capitalize()} - stacked percentages", fontsize=8)
    plt.ylabel("Percentage (%)", fontsize=8)
    plt.xlabel("Group", fontsize=8)
    plt.ylim(0, 100)
    plt.xticks(rotation=45, ha="right", fontsize=7)
    plt.yticks(fontsize=7)
    plt.legend(
        title=level.capitalize(),
        fontsize=6,
        title_fontsize=7,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )
    plt.tight_layout()
    png_group_percent = out_dir / "plots" / f"{level}_group_stacked_bar_percent.png"
    pdf_group_percent = out_dir / "plots" / f"{level}_group_stacked_bar_percent.pdf"
    plt.savefig(png_group_percent, dpi=600)
    plt.savefig(pdf_group_percent)
    plt.close()

    all_tab = out_dir / "tables" / f"{level}_overall_counts.tsv"
    grp_tab = out_dir / "tables" / f"{level}_group_counts.tsv"
    counts_all.to_csv(all_tab, sep="\t", index=False)
    counts_group.to_csv(grp_tab, sep="\t", index=False)

    return png_all_counts, png_all_percent, png_group_counts, png_group_percent, all_tab


def plot_category_heatmap(
    merged: pd.DataFrame, out_dir: Path, level: str
) -> Optional[Path]:
    """
    Heatmap of group x category showing percentage of samples in each group
    that have at least one hit in that category.
    """
    col = safe_get_col(merged, level)
    if col is None:
        return None

    counts_group = (
        merged.groupby(["group", col])["sample"]
        .nunique()
        .reset_index(name="n_samples_with_category")
    )
    group_sizes = (
        merged[["sample", "group"]]
        .drop_duplicates()
        .groupby("group")["sample"]
        .nunique()
        .rename("n_samples_in_group")
        .reset_index()
    )
    counts_group = counts_group.merge(group_sizes, on="group", how="left")
    counts_group["percent_samples"] = (
        counts_group["n_samples_with_category"] * 100.0
        / counts_group["n_samples_in_group"].astype(float)
    )

    mat = counts_group.pivot_table(
        index="group",
        columns=col,
        values="percent_samples",
        fill_value=0,
    )

    plt.figure(figsize=(max(8, mat.shape[1] * 0.3), 4 + mat.shape[0] * 0.5))
    sns.heatmap(mat, cmap="viridis", cbar_kws={"label": "Percentage of samples"})
    plt.title(f"{level.capitalize()} - group-wise percentage heatmap")
    plt.xlabel(level.capitalize())
    plt.ylabel("Group")
    png = out_dir / "plots" / f"{level}_group_heatmap.png"
    pdf = out_dir / "plots" / f"{level}_group_heatmap.pdf"
    plt.tight_layout()
    plt.savefig(png, dpi=600)
    plt.savefig(pdf)
    plt.close()

    return png


def plot_bubble_chart(
    merged: pd.DataFrame, out_dir: Path, level: str, top_n: int = 25
) -> Optional[Path]:
    """
    Bubble chart for group x category:
    - x-axis: group
    - y-axis: category (top_n by global prevalence)
    - bubble size: percentage of samples in that group with at least one hit
    """
    col = safe_get_col(merged, level)
    if col is None:
        return None

    total_samples = merged["sample"].nunique()

    global_counts = (
        merged.groupby(col)["sample"]
        .nunique()
        .reset_index(name="n_samples_with_category")
    )
    global_counts["percent_samples"] = (
        global_counts["n_samples_with_category"] * 100.0 / float(total_samples)
    )
    top_cats = (
        global_counts.sort_values("percent_samples", ascending=False)
        .head(top_n)[col]
        .tolist()
    )

    counts_group = (
        merged[merged[col].isin(top_cats)]
        .groupby(["group", col])["sample"]
        .nunique()
        .reset_index(name="n_samples_with_category")
    )
    group_sizes = (
        merged[["sample", "group"]]
        .drop_duplicates()
        .groupby("group")["sample"]
        .nunique()
        .rename("n_samples_in_group")
        .reset_index()
    )
    counts_group = counts_group.merge(group_sizes, on="group", how="left")
    counts_group["percent_samples"] = (
        counts_group["n_samples_with_category"] * 100.0
        / counts_group["n_samples_in_group"].astype(float)
    )

    groups = sorted(counts_group["group"].unique().tolist())
    cat_order = top_cats

    group_to_x = {g: i for i, g in enumerate(groups)}
    cat_to_y = {c: i for i, c in enumerate(cat_order)}

    counts_group["x"] = counts_group["group"].map(group_to_x)
    counts_group["y"] = counts_group[col].map(cat_to_y)

    plt.figure(figsize=(max(8, len(groups) * 1.5), max(6, len(cat_order) * 0.3)))
    sizes = counts_group["percent_samples"].values
    size_scale = 30.0
    plt.scatter(
        counts_group["x"],
        counts_group["y"],
        s=(sizes + 1.0) * size_scale,
        alpha=0.7,
        edgecolor="k",
    )
    plt.xticks(range(len(groups)), groups, rotation=45, ha="right")
    plt.yticks(range(len(cat_order)), cat_order)
    plt.xlabel("Group")
    plt.ylabel(level.capitalize())
    plt.title(f"{level.capitalize()} - group-wise bubble chart (top {top_n} categories)")
    png = out_dir / "plots" / f"{level}_group_bubble.png"
    pdf = out_dir / "plots" / f"{level}_group_bubble.pdf"
    plt.tight_layout()
    plt.savefig(png, dpi=600)
    plt.savefig(pdf)
    plt.close()

    return png


def run_association_tests(merged: pd.DataFrame, level: str) -> Optional[pd.DataFrame]:
    """
    For each category at a given level, test association between category presence
    and group using Fisher (2 groups) or Chi-squared (>2 groups) + BH FDR.
    """
    col = safe_get_col(merged, level)
    if col is None:
        return None

    pres = (
        merged.groupby(["sample", "group", col])
        .size()
        .reset_index(name="count")
    )
    pres["present"] = 1
    pres = pres[["sample", "group", col, "present"]].drop_duplicates()

    samples = merged[["sample", "group"]].drop_duplicates()
    categories = pres[col].drop_duplicates()
    samples["key"] = 1
    cats_df = categories.to_frame(name=col)
    cats_df["key"] = 1
    grid = samples.merge(cats_df, on="key").drop("key", axis=1)

    grid = grid.merge(pres, on=["sample", "group", col], how="left")
    grid["present"] = grid["present"].fillna(0).astype(int)

    results = []
    groups = grid["group"].unique()
    n_groups = len(groups)

    for cat in categories:
        sub = grid[grid[col] == cat]
        contingency = pd.crosstab(sub["group"], sub["present"])
        if contingency.shape[1] < 2:
            continue

        if n_groups == 2:
            try:
                oddsratio, pval = fisher_exact(contingency.values)
            except Exception:
                oddsratio, pval = np.nan, np.nan
            stat = oddsratio
            test_name = "Fisher_exact"
        else:
            chi2, pval, dof, expected = chi2_contingency(contingency.values)
            stat = chi2
            test_name = "Chi2"

        results.append(
            {
                level: cat,
                "test": test_name,
                "statistic": stat,
                "pvalue": pval,
            }
        )

    if not results:
        return None

    res_df = pd.DataFrame(results)
    res_df = res_df.sort_values("pvalue").reset_index(drop=True)
    m = len(res_df)
    res_df["rank"] = np.arange(1, m + 1)
    res_df["qvalue_BH"] = res_df["pvalue"] * m / res_df["rank"]
    res_df["qvalue_BH"] = res_df["qvalue_BH"].clip(upper=1.0)

    return res_df


def run_pca_on_genes(
    merged: pd.DataFrame, out_dir: Path, sample_df: pd.DataFrame
) -> Dict[str, object]:
    """
    Run PCA on gene presence/absence matrix (samples x genes).
    Returns dict with PC table, explained variance, and plot.
    """
    gene_col = safe_get_col(merged, "gene_symbol") or ID_COLUMN
    if gene_col not in merged.columns:
        return {}

    presence = (
        merged.groupby(["sample", gene_col])
        .size()
        .reset_index(name="count")
    )
    presence["present"] = 1

    mat = presence.pivot_table(
        index="sample",
        columns=gene_col,
        values="present",
        fill_value=0,
        aggfunc="max",
    )

    groups = (
        sample_df[["sample", "group"]]
        .drop_duplicates()
        .set_index("sample")
        .reindex(mat.index)
    )
    group_labels = groups["group"].fillna("NA").tolist()

    X = mat.values.astype(float)
    if X.shape[1] < 2 or X.shape[0] < 2:
        return {}

    X_mean = X.mean(axis=0, keepdims=True)
    X_centered = X - X_mean

    U, S, VT = np.linalg.svd(X_centered, full_matrices=False)
    PCs = U * S

    var_explained = (S ** 2) / (X.shape[0] - 1)
    var_explained_ratio = var_explained / var_explained.sum()

    pc_df = pd.DataFrame(
        {
            "sample": mat.index,
            "group": group_labels,
            "PC1": PCs[:, 0],
            "PC2": PCs[:, 1] if PCs.shape[1] > 1 else 0.0,
        }
    )
    pc_table_path = out_dir / "tables" / "pca_gene_presence.tsv"
    pc_df.to_csv(pc_table_path, sep="\t", index=False)

    plt.figure(figsize=(8, 6))
    unique_groups = sorted(set(group_labels))
    group_to_idx = {g: i for i, g in enumerate(unique_groups)}
    colors = [group_to_idx[g] for g in group_labels]
    plt.scatter(pc_df["PC1"], pc_df["PC2"], c=colors, alpha=0.8)
    handles = []
    for g, idx in group_to_idx.items():
        handles.append(plt.Line2D([], [], marker="o", linestyle="", label=g))
    plt.legend(handles=handles, title="Group", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.xlabel("PC1 ({:.1f}% variance)".format(var_explained_ratio[0] * 100.0))
    if len(var_explained_ratio) > 1:
        plt.ylabel("PC2 ({:.1f}% variance)".format(var_explained_ratio[1] * 100.0))
    else:
        plt.ylabel("PC2")
    plt.title("PCA of AMR gene presence/absence")
    pca_png = out_dir / "plots" / "pca_gene_presence.png"
    pca_pdf = out_dir / "plots" / "pca_gene_presence.pdf"
    plt.tight_layout()
    plt.savefig(pca_png, dpi=600)
    plt.savefig(pca_pdf)
    plt.close()

    return {
        "pc_table_path": pc_table_path,
        "explained_variance_ratio": var_explained_ratio.tolist(),
        "n_genes": mat.shape[1],
        "plot_path": pca_png,
    }


def run_overall_ttest_gene_counts(
    merged: pd.DataFrame, out_dir: Path, sample_df: pd.DataFrame
) -> Dict[str, object]:
    """
    Compute total AMR counts per sample and perform Welch t-tests between groups.
    """
    counts = (
        merged.groupby("sample")
        .size()
        .reset_index(name="total_amr_hits")
    )
    counts = counts.merge(
        sample_df[["sample", "group"]].drop_duplicates(),
        on="sample",
        how="left",
    )

    groups = counts["group"].dropna().unique()
    if len(groups) < 2:
        return {}

    results = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            g1, g2 = groups[i], groups[j]
            x = counts.loc[counts["group"] == g1, "total_amr_hits"].astype(float)
            y = counts.loc[counts["group"] == g2, "total_amr_hits"].astype(float)
            if len(x) < 2 or len(y) < 2:
                continue
            t_stat, p_val = ttest_ind(x, y, equal_var=False)
            results.append(
                {
                    "group1": g1,
                    "group2": g2,
                    "n_group1": len(x),
                    "n_group2": len(y),
                    "mean_hits_group1": x.mean(),
                    "mean_hits_group2": y.mean(),
                    "t_statistic": t_stat,
                    "pvalue": p_val,
                }
            )

    if not results:
        return {}

    res_df = pd.DataFrame(results)
    ttest_path = out_dir / "tables" / "ttest_total_amr_hits_between_groups.tsv"
    res_df.to_csv(ttest_path, sep="\t", index=False)

    return {
        "table_path": ttest_path,
        "n_comparisons": len(res_df),
    }


def write_html_report(
    step_dir: Path,
    merged: pd.DataFrame,
    heatmap_paths: Tuple[Path, Path],
    bar_info: Dict[str, Tuple[Optional[Path], Optional[Path], Optional[Path], Optional[Path], Optional[Path]]],
    cat_heatmaps: Dict[str, Optional[Path]],
    bubble_info: Dict[str, Optional[Path]],
    stat_tables: Dict[str, Tuple[pd.DataFrame, Path]],
    pca_info: Dict[str, object],
    ttest_info: Dict[str, object],
    sample_df: pd.DataFrame,
    report_name: str,
) -> Path:
    """
    Generate HTML report summarizing results, linking tables and plots.
    """
    report = step_dir / report_name

    def rel(p: Optional[Path]) -> str:
        if p is None:
            return ""
        return os.path.relpath(p, start=step_dir)

    n_samples = merged["sample"].nunique()
    n_groups = merged["group"].nunique()
    gene_col = safe_get_col(merged, "gene_symbol") or ID_COLUMN
    if gene_col in merged.columns:
        n_genes = merged[gene_col].nunique()
    else:
        n_genes = 0

    with open(report, "w") as fh:
        fh.write("<html><head><meta charset='utf-8'>\n")
        fh.write("<title>AMRFinder Statistics Report</title>\n")
        fh.write(
            "<style>"
            "body{font-family:Arial, sans-serif; margin:40px;}"
            "img{max-width:100%;}"
            "table{border-collapse:collapse;}"
            "th,td{border:1px solid #ccc; padding:4px 8px;}"
            "</style>\n"
        )
        fh.write("</head><body>\n")

        fh.write("<h1>AMRFinder Statistics Report</h1>\n")

        fh.write("<h2>Summary</h2>\n")
        fh.write("<ul>\n")
        fh.write("<li>Number of samples: <b>{}</b></li>\n".format(n_samples))
        fh.write("<li>Number of groups: <b>{}</b></li>\n".format(n_groups))
        fh.write("<li>Number of unique AMR genes ({}): <b>{}</b></li>\n".format(gene_col, n_genes))
        fh.write("</ul>\n")

        fh.write("<h3>Sample / Group overview</h3>\n")
        fh.write("<table><tr><th>Sample</th><th>Group</th></tr>\n")
        for _, row in sample_df[["sample", "group"]].drop_duplicates().iterrows():
            fh.write(
                "<tr><td>{}</td><td>{}</td></tr>\n".format(row["sample"], row["group"])
            )
        fh.write("</table>\n")

        fh.write("<h2>Merged AMRFinder results</h2>\n")
        fh.write("<p>Combined table with all samples and groups:</p>\n")
        fh.write("<ul>\n")
        fh.write("<li><a href='tables/amrfinder_merged.tsv'>amrfinder_merged.tsv (TSV)</a></li>\n")
        fh.write("<li><a href='tables/amrfinder_merged.csv'>amrfinder_merged.csv (CSV)</a></li>\n")
        fh.write("</ul>\n")

        fh.write("<h2>Gene presence/absence heatmap</h2>\n")
        if heatmap_paths:
            fh.write(
                "<p><img src='{}' alt='Gene presence heatmap'></p>\n".format(
                    rel(heatmap_paths[0])
                )
            )
        else:
            fh.write("<p><i>Gene-level heatmap could not be generated.</i></p>\n")

        fh.write("<h2>Category-wise summaries</h2>\n")
        for level, info in bar_info.items():
            fh.write("<h3>{}</h3>\n".format(level.capitalize()))
            if info is None:
                fh.write("<p><i>Column for {} not found in AMRFinder output.</i></p>\n".format(level))
                continue
            png_all_counts, png_all_percent, png_group_counts, png_group_percent, table_path = info
            fh.write("<p>Overall and group-wise distributions:</p>\n")
            if png_all_percent:
                fh.write(
                    "<p><b>Overall barplot (percentage of samples):</b><br><img src='{}'></p>\n".format(
                        rel(png_all_percent)
                    )
                )
            if png_all_counts:
                fh.write(
                    "<p><b>Overall barplot (number of samples):</b><br><img src='{}'></p>\n".format(
                        rel(png_all_counts)
                    )
                )
            if png_group_percent:
                fh.write(
                    "<p><b>Group-wise 100% stacked barplot (percentage):</b><br><img src='{}'></p>\n".format(
                        rel(png_group_percent)
                    )
                )
            if png_group_counts:
                fh.write(
                    "<p><b>Group-wise stacked barplot (counts):</b><br><img src='{}'></p>\n".format(
                        rel(png_group_counts)
                    )
                )
            if level in cat_heatmaps and cat_heatmaps[level] is not None:
                fh.write(
                    "<p><b>Group-wise percentage heatmap:</b><br><img src='{}'></p>\n".format(
                        rel(cat_heatmaps[level])
                    )
                )
            if level in bubble_info and bubble_info[level] is not None:
                fh.write(
                    "<p><b>Bubble chart (top categories):</b><br><img src='{}'></p>\n".format(
                        rel(bubble_info[level])
                    )
                )
            if table_path:
                fh.write(
                    "<p>Counts and percentages table: <a href='{}'>{}</a></p>\n".format(
                        rel(table_path), Path(table_path).name
                    )
                )

        fh.write("<h2>Group-wise association tests (Chi-square / Fisher)</h2>\n")
        if not stat_tables:
            fh.write("<p><i>No association tests could be computed.</i></p>\n")
        else:
            fh.write(
                "<p>For each category level we tested association with group using "
                "Fisher's exact test (2 groups) or Chi-squared test (> 2 groups). "
                "P-values were corrected using the Benjamini-Hochberg method.</p>\n"
            )
            fh.write("<ul>\n")
            for level, (_df, path) in stat_tables.items():
                fh.write(
                    "<li>{} statistics: <a href='{}'>{}</a></li>\n".format(
                        level.capitalize(), rel(path), Path(path).name
                    )
                )
            fh.write("</ul>\n")

        fh.write("<h2>PCA of AMR gene presence/absence</h2>\n")
        if pca_info:
            if pca_info.get("plot_path"):
                fh.write(
                    "<p><img src='{}' alt='PCA plot'></p>\n".format(
                        rel(pca_info["plot_path"])
                    )
                )
            evr = pca_info.get("explained_variance_ratio", [])
            if evr:
                fh.write(
                    "<p>PC1 explains {:.1f}% of variance; PC2 explains {:.1f}% of variance.</p>\n".format(
                        evr[0] * 100.0, (evr[1] * 100.0 if len(evr) > 1 else 0.0)
                    )
                )
            if pca_info.get("pc_table_path"):
                fh.write(
                    "<p>PCA coordinates table: <a href='{}'>{}</a></p>\n".format(
                        rel(pca_info["pc_table_path"]),
                        Path(pca_info["pc_table_path"]).name,
                    )
                )
        else:
            fh.write("<p><i>PCA could not be computed (not enough genes or samples).</i></p>\n")

        fh.write("<h2>T-tests on total AMR burden</h2>\n")
        if ttest_info and ttest_info.get("table_path"):
            fh.write(
                "<p>T-tests (Welch) comparing total AMR hits per sample between groups. "
                "Results table: <a href='{}'>{}</a></p>\n".format(
                    rel(ttest_info["table_path"]),
                    Path(ttest_info["table_path"]).name,
                )
            )
        else:
            fh.write("<p><i>T-tests could not be computed (not enough groups or samples).</i></p>\n")

        fh.write("<hr>\n")
        fh.write("<p><i>Generated by amr_statistics_pipeline.py</i></p>\n")
        fh.write("</body></html>\n")

    return report


def write_interpretation_txt(
    step_dir: Path,
    merged: pd.DataFrame,
    sample_df: pd.DataFrame,
    stat_tables: Dict[str, Tuple[pd.DataFrame, Path]],
    pca_info: Dict[str, object],
    ttest_info: Dict[str, object],
) -> Path:
    """
    Write a plain-text interpretation guide for the main plots and statistics.
    """
    interp_path = step_dir / "interpretation.txt"

    n_samples = merged["sample"].nunique()
    n_groups = merged["group"].nunique()
    gene_col = safe_get_col(merged, "gene_symbol") or ID_COLUMN
    n_genes = merged[gene_col].nunique() if gene_col in merged.columns else 0

    with open(interp_path, "w") as fh:
        fh.write("AMR STATISTICS INTERPRETATION GUIDE\n")
        fh.write("===================================\n\n")
        fh.write("Dataset overview\n")
        fh.write("----------------\n")
        fh.write("Number of samples: {}\n".format(n_samples))
        fh.write("Number of groups: {}\n".format(n_groups))
        fh.write("Number of unique AMR genes ({}): {}\n".format(gene_col, n_genes))
        fh.write("\n")

        fh.write("1. Gene presence/absence heatmap\n")
        fh.write("--------------------------------\n")
        fh.write(
            "The gene presence/absence heatmap shows samples on the y-axis and AMR genes on "
            "the x-axis. A filled cell indicates that at least one hit for that gene was "
            "detected in the corresponding sample.\n"
        )
        fh.write(
            "- Horizontal patterns (across many samples) indicate genes that are widely "
            "distributed in the dataset.\n"
        )
        fh.write(
            "- Vertical clusters of genes may indicate co-occurring resistance determinants.\n"
        )
        fh.write(
            "- Sample-level clusters (visually similar rows) can reflect similar resistome "
            "profiles, often corresponding to the same group or ecological niche.\n\n"
        )

        fh.write("2. Category-level barplots and stacked bars\n")
        fh.write("-------------------------------------------\n")
        fh.write(
            "For each category level (type, subtype, class, subclass, method), the script "
            "generates:\n"
        )
        fh.write(
            "  a) Overall barplots:\n"
            "     - Counts: number of samples with at least one hit per category.\n"
            "     - Percentages: percentage of samples with at least one hit per category.\n"
        )
        fh.write(
            "  b) Group-wise stacked barplots:\n"
            "     - Counts: absolute number of samples in each group with that category.\n"
            "     - 100% stacked: percentages where each group bar sums to 100%, so "
            "you can compare composition between groups even if group sizes differ.\n"
        )
        fh.write(
            "Categories with high percentages are common across the dataset; those enriched "
            "in a particular group often drive resistome differences.\n\n"
        )

        fh.write("3. Category-level heatmaps\n")
        fh.write("--------------------------\n")
        fh.write(
            "For each category level, the group-by-category heatmap shows the percentage of "
            "samples in each group in which that category is present.\n"
        )
        fh.write(
            "- Rows correspond to groups, columns to categories; warmer colors indicate a "
            "higher percentage of positive samples.\n"
        )
        fh.write(
            "- This is useful for visually identifying categories that are strongly associated "
            "with specific groups.\n\n"
        )

        fh.write("4. Bubble charts (Krona-like overview)\n")
        fh.write("--------------------------------------\n")
        fh.write(
            "The bubble charts summarize group-wise patterns for the top categories at each level. "
            "On these plots:\n"
        )
        fh.write(
            "- The x-axis corresponds to groups and the y-axis to categories (for example, AMR class).\n"
        )
        fh.write(
            "- The size of each bubble is proportional to the percentage of samples in that group "
            "where the category is present.\n"
        )
        fh.write(
            "- Large bubbles highlight dominant categories within groups; comparing bubble sizes "
            "across groups indicates which categories drive differences between groups.\n\n"
        )

        fh.write("5. Statistical association tests (Chi-square / Fisher)\n")
        fh.write("------------------------------------------------------\n")
        if stat_tables:
            fh.write(
                "For each category level, the script performs statistical tests for association "
                "between category presence (present vs absent) and group.\n"
            )
            fh.write(
                "- If there are exactly two groups, Fisher's exact test is used.\n"
            )
            fh.write(
                "- If there are more than two groups, a Chi-squared test is applied.\n"
            )
            fh.write(
                "- P-values are corrected using the Benjamini-Hochberg procedure to control the "
                "false discovery rate (q-values).\n"
            )
            fh.write(
                "- Categories with q-value < 0.05 can be interpreted as significantly associated "
                "with group identity, under the usual assumptions of the tests.\n\n"
            )
            fh.write(
                "Inspect the *_group_association_stats.tsv files in the tables/ directory. Focus "
                "on categories with the smallest q-values and cross-reference them with the "
                "barplots and heatmaps to understand which groups show enrichment.\n\n"
            )
        else:
            fh.write(
                "No valid association tables were generated (for example, due to insufficient "
                "variation or too few groups). In such cases, visual inspection of the barplots "
                "and heatmaps may still reveal qualitative patterns.\n\n"
            )

        fh.write("6. PCA of gene presence/absence\n")
        fh.write("--------------------------------\n")
        if pca_info:
            evr = pca_info.get("explained_variance_ratio", [])
            if evr:
                fh.write(
                    "Principal component analysis (PCA) was run on the binary matrix of gene "
                    "presence/absence across samples.\n"
                )
                fh.write(
                    "PC1 explains approximately {:.1f}% of the variance; PC2 explains approximately "
                    "{:.1f}%.\n".format(
                        evr[0] * 100.0, (evr[1] * 100.0 if len(evr) > 1 else 0.0)
                    )
                )
            else:
                fh.write(
                    "PCA was attempted, but explained variance could not be computed reliably.\n"
                )
            fh.write(
                "In the PCA scatterplot, each point is a sample, colored by group. Separation "
                "between groups in PC space suggests that overall resistome composition differs "
                "systematically between groups.\n"
            )
            fh.write(
                "- Tight clustering of samples from the same group indicates homogeneous resistomes.\n"
            )
            fh.write(
                "- Overlap between groups implies that resistome profiles are more similar or "
                "that grouping does not strongly follow AMR patterns.\n\n"
            )
        else:
            fh.write(
                "PCA could not be computed (for example, due to too few genes or samples). If "
                "needed, consider filtering for more variable genes or combining additional data.\n\n"
            )

        fh.write("7. T-tests on total AMR burden\n")
        fh.write("--------------------------------\n")
        if ttest_info and ttest_info.get("table_path"):
            fh.write(
                "The script computes the total number of AMR hits per sample and compares these "
                "distributions between groups using Welch's t-test.\n"
            )
            fh.write(
                "- A significant p-value (for example, < 0.05) indicates that the mean AMR burden differs "
                "between the two groups being compared.\n"
            )
            fh.write(
                "- Inspect the t-test results table (ttest_total_amr_hits_between_groups.tsv) to see "
                "which group pairs most strongly differ in total AMR counts.\n\n"
            )
        else:
            fh.write(
                "T-tests on total AMR burden were not performed or did not yield valid results "
                "(for example, due to too few samples per group).\n\n"
            )

        fh.write("8. Practical notes for publication\n")
        fh.write("----------------------------------\n")
        fh.write(
            "- The PNG files in the plots/ directory are suitable for quick review, while the PDF "
            "files are vector-based and suitable for journal submission.\n"
        )
        fh.write(
            "- When preparing final figures, you may crop or combine panels, but retain axis labels "
            "and legends so that plots remain interpretable.\n"
        )
        fh.write(
            "- Report which statistical tests were used (Chi-square, Fisher's exact, Welch's t-test) "
            "and how multiple testing was handled (Benjamini-Hochberg FDR).\n"
        )
        fh.write(
            "- When interpreting significance, always consider biological plausibility and sample size, "
            "not only q-values.\n"
        )

    return interp_path


def main():
    set_publication_style()

    args = parse_args()

    # 1) set up output dir
    step_dir = find_next_step_dir(args.results_dir, suffix="AMR_Statistics")
    print("[INFO] Writing outputs to: {}".format(step_dir))

    # 2) read sample sheet
    sample_df = load_samples(args.samples)

    # 3) find AMR files
    amr_files = find_amr_files(args.results_dir, args.amr_pattern)
    print("[INFO] Found {} AMRFinder TSV files".format(len(amr_files)))

    # 4) merge
    merged = merge_amr_files(amr_files, sample_df)
    merged_tsv, merged_csv = write_table(merged, step_dir, "amrfinder_merged")
    print("[INFO] Merged table written to: {}".format(merged_tsv))

    # 5) gene presence heatmap
    try:
        heatmap_paths = plot_gene_presence_heatmap(merged, step_dir)
    except Exception as e:
        print("[WARN] Could not generate gene presence heatmap: {}".format(e))
        heatmap_paths = ()

    # 6) category-level plots
    bar_info: Dict[str, Tuple[Optional[Path], Optional[Path], Optional[Path], Optional[Path], Optional[Path]]] = {}
    cat_heatmaps: Dict[str, Optional[Path]] = {}
    bubble_info: Dict[str, Optional[Path]] = {}

    for level in CATEGORY_LEVELS:
        print("[INFO] Generating barplots for {}".format(level))
        bar_info[level] = plot_category_bars(merged, step_dir, level)

        print("[INFO] Generating heatmap for {}".format(level))
        cat_heatmaps[level] = plot_category_heatmap(merged, step_dir, level)

        print("[INFO] Generating bubble chart for {}".format(level))
        bubble_info[level] = plot_bubble_chart(merged, step_dir, level)

    # 7) statistics: association tests
    stat_tables: Dict[str, Tuple[pd.DataFrame, Path]] = {}
    for level in CATEGORY_LEVELS:
        print("[INFO] Running association tests for {}".format(level))
        stats_df = run_association_tests(merged, level)
        if stats_df is not None:
            out_path = step_dir / "tables" / f"{level}_group_association_stats.tsv"
            stats_df.to_csv(out_path, sep="\t", index=False)
            stat_tables[level] = (stats_df, out_path)
            print("  -> {}".format(out_path))

    # 8) PCA
    print("[INFO] Running PCA on gene presence/absence")
    pca_info = run_pca_on_genes(merged, step_dir, sample_df)

    # 9) Overall t-tests on total AMR hits
    print("[INFO] Running t-tests on total AMR hits between groups")
    ttest_info = run_overall_ttest_gene_counts(merged, step_dir, sample_df)

    # 10) HTML report
    report = write_html_report(
        step_dir,
        merged,
        heatmap_paths,
        bar_info,
        cat_heatmaps,
        bubble_info,
        stat_tables,
        pca_info,
        ttest_info,
        sample_df,
        args.report_name,
    )
    print("[INFO] HTML report written to: {}".format(report))

    # 11) Interpretation text
    interp_path = write_interpretation_txt(
        step_dir, merged, sample_df, stat_tables, pca_info, ttest_info
    )
    print("[INFO] Interpretation text written to: {}".format(interp_path))


if __name__ == "__main__":
    main()
