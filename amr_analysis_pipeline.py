#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AMR analysis pipeline for results/OneHealthAMR_AMRFinder_summary.xlsx.

Features:
- Uses sample.tsv (if present) to set sample -> group mapping.
- Reads AMRFinder combined summary from results/OneHealthAMR_AMRFinder_summary.xlsx.
- Checks required Python packages; if missing, prints a note and attempts to install them with pip.
- Scans results/tools/N_* to find the maximum N and creates results/tools/(N+1)_amr_analysis_pipeline.
- Performs:
    Stage 1 - Core AMR Gene Summary
    Stage 2 - Group-Wise Comparison
    Stage 3 - Co-occurrence / Gene Interaction Networks
              + GraphML outputs for Cytoscape
              + Node attributes
              + Core-periphery (k-core) annotation
              + Per-group networks and comparison
    Stage 4 - Virulence and AMR Link (if VIRULENCE type exists)
    Stage 5 - Candidate Novel AMR Mutations (heuristic)
    Stage 6 - Machine Learning Without Phenotype (predict group from AMR genes)
- Writes:
    - CSV tables
    - PNG charts
    - GraphML files for Cytoscape
    - stage6_rf_interpretation.txt
    - figure_captions_and_methods_template.txt
    - all_stages_interpretation.txt
    - instructions_for_cytoscape.txt
    - run_log_and_file_structure.txt
"""

import sys
import subprocess
import importlib
from pathlib import Path
import itertools
import re

# ----------------------------------------------------------------------
# 0. Ensure required packages
# ----------------------------------------------------------------------
REQUIRED_PACKAGES = [
    "pandas",
    "numpy",
    "matplotlib",
    "seaborn",
    "networkx",
    "sklearn",
    "scipy",
]


def ensure_packages(packages):
    """
    Try importing each package. If ImportError, print a note,
    then attempt to install via pip, and import again.
    """
    for pkg in packages:
        try:
            importlib.import_module(pkg)
        except ImportError:
            pip_name = pkg if pkg != "sklearn" else "scikit-learn"
            print(
                "[INFO] Python package '{}' not found in this environment.\n"
                "       Attempting to install it automatically via pip ({})."
                .format(pkg, pip_name)
            )
            try:
                subprocess.check_call(
                    [sys.executable, "-m", "pip", "install", pip_name]
                )
                print("[INFO] Successfully installed '{}'.".format(pip_name))
            except Exception as e:
                print(
                    "[WARNING] Automatic installation of '{}' failed:\n    {}\n"
                    "          Please install it manually (for example: "
                    "'conda install {}' or 'pip install {}') and re-run the script."
                    .format(pip_name, e, pip_name, pip_name)
                )
            try:
                importlib.import_module(pkg)
            except ImportError:
                print(
                    "[ERROR] Still cannot import '{}' after installation attempt; "
                    "some functionality may fail.".format(pkg)
                )


ensure_packages(REQUIRED_PACKAGES)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn.metrics import classification_report, confusion_matrix
from scipy.stats import chi2_contingency

# ----------------------------------------------------------------------
# 1. Paths and dynamic results directory
# ----------------------------------------------------------------------
BASE_RESULTS = Path("results")
AMR_SUMMARY_XLSX = BASE_RESULTS / "OneHealthAMR_AMRFinder_summary.xlsx"
SAMPLE_TSV = Path("sample.tsv")

TOOLS_DIR = BASE_RESULTS / "tools"
TOOLS_DIR.mkdir(parents=True, exist_ok=True)


def get_next_tools_dir(base_tools_dir, suffix="amr_analysis_pipeline"):
    """
    Scan base_tools_dir for directories named like N_*,
    find max N, and create (N+1)_<suffix>.
    If none exist, N starts at 0.
    """
    max_n = -1
    for child in base_tools_dir.iterdir():
        if child.is_dir():
            name = child.name
            parts = name.split("_", 1)
            if len(parts) >= 2 and parts[0].isdigit():
                n = int(parts[0])
                if n > max_n:
                    max_n = n
    next_n = max_n + 1
    outdir = base_tools_dir / f"{next_n}_{suffix}"
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


OUTDIR = get_next_tools_dir(TOOLS_DIR, "amr_analysis_pipeline")
print("[INFO] Results will be written to:", OUTDIR.resolve())

sns.set(style="whitegrid")

# Analysis parameters
COMBINED_SHEET = "combined_raw"
TOP_N_GENES_FOR_PLOTS = 20
MIN_COOC_SAMPLE_COUNT = 2
TRAIN_TEST_SPLIT = 0.3

# Network visualization parameters
MAX_NETWORK_NODES = 30   # draw only top N most connected genes
DEGREE_LABEL_MIN = 2     # label nodes only if degree >= this


def sanitize_group_name(group):
    """
    Make a group name safe for filenames: keep alphanumeric and replace others with underscore.
    """
    if pd.isna(group):
        return "NA"
    g = str(group)
    g = re.sub(r"[^A-Za-z0-9]+", "_", g)
    g = g.strip("_")
    if not g:
        g = "group"
    return g

# ----------------------------------------------------------------------
# 2. Load data and integrate sample.tsv
# ----------------------------------------------------------------------
def load_sample_groups_from_tsv(sample_tsv_path):
    """
    If sample.tsv exists, read it and return mapping sample -> group (Series).
    Expects columns: sample, group.
    """
    if not sample_tsv_path.exists():
        print("[INFO] sample.tsv not found; will rely on 'group' column in AMR summary.")
        return None

    df_samp = pd.read_csv(sample_tsv_path, sep="\t")
    if "sample" not in df_samp.columns or "group" not in df_samp.columns:
        print("[WARNING] sample.tsv does not contain 'sample' and 'group' columns. Ignoring sample.tsv.")
        return None

    sg = df_samp[["sample", "group"]].drop_duplicates().set_index("sample")["group"]
    print("[INFO] Loaded {} sample group entries from sample.tsv.".format(len(sg)))
    return sg


def load_combined_amr(path_xlsx, sheet_name, sample_group_override):
    """
    Load combined AMRFinder output and, if sample_group_override is provided,
    override the group column using sample.tsv mapping.
    """
    if not path_xlsx.exists():
        raise FileNotFoundError("AMR summary file not found: {}".format(path_xlsx))

    df = pd.read_excel(path_xlsx, sheet_name=sheet_name)
    print(
        "[INFO] Loaded {} AMR hits from {} [sheet='{}']."
        .format(len(df), path_xlsx, sheet_name)
    )

    if "sample" not in df.columns:
        raise ValueError("Expected 'sample' column in combined AMR sheet, but it is missing.")

    if sample_group_override is not None:
        df["group_from_tsv"] = df["sample"].map(sample_group_override)
        if "group" in df.columns:
            df["group_original"] = df["group"]
            df["group"] = df["group_from_tsv"].combine_first(df["group_original"])
        else:
            df["group"] = df["group_from_tsv"]
        print("[INFO] Group column has been updated using sample.tsv (where available).")
    else:
        if "group" not in df.columns:
            raise ValueError("No 'group' information found (sample.tsv missing and 'group' column absent).")

    return df

# ----------------------------------------------------------------------
# 3. Gene presence matrix
# ----------------------------------------------------------------------
def make_gene_presence_matrix(df_combined):
    """
    Returns:
      gene_matrix: samples x genes (0/1 for presence of element symbol)
      sample_groups: group label per sample (Series)
    """
    sample_group = (
        df_combined[["sample", "group"]]
        .drop_duplicates()
        .set_index("sample")["group"]
    )

    df_presence = (
        df_combined
        .groupby(["sample", "element symbol"])
        .size()
        .reset_index(name="count")
    )
    df_presence["binary"] = 1
    gene_matrix = (
        df_presence
        .pivot_table(index="sample", columns="element symbol", values="binary", fill_value=0)
        .astype(int)
    )

    sample_groups = sample_group.reindex(gene_matrix.index)
    return gene_matrix, sample_groups

# ----------------------------------------------------------------------
# 4. Stage 1 - Core AMR Gene Summary
# ----------------------------------------------------------------------
def stage1_core_summary(df_combined, outdir):
    print("[Stage 1] Core AMR Gene Summary")

    gene_counts = df_combined["element symbol"].value_counts().rename("count")
    gene_counts.to_csv(outdir / "stage1_gene_counts.csv")

    class_counts = df_combined["class"].value_counts().rename("count")
    class_counts.to_csv(outdir / "stage1_class_counts.csv")

    top_genes = gene_counts.head(TOP_N_GENES_FOR_PLOTS)
    plt.figure(figsize=(10, 6))
    sns.barplot(x=top_genes.values, y=top_genes.index)
    plt.xlabel("Number of hits")
    plt.ylabel("AMR gene (element symbol)")
    plt.title("Top {} AMR genes by count".format(TOP_N_GENES_FOR_PLOTS))
    plt.tight_layout()
    plt.savefig(outdir / "stage1_top_genes_barplot.png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 5))
    sns.barplot(x=class_counts.index, y=class_counts.values)
    plt.xticks(rotation=45, ha="right")
    plt.xlabel("AMR class")
    plt.ylabel("Number of hits")
    plt.title("AMR class distribution across all isolates")
    plt.tight_layout()
    plt.savefig(outdir / "stage1_class_barplot.png", dpi=300)
    plt.close()

# ----------------------------------------------------------------------
# 5. Stage 2 - Group-Wise Comparison
# ----------------------------------------------------------------------
def _enrichment_test(table):
    """
    Perform simple chi-square test for each feature across groups.
    table: features x groups (counts of samples with the feature)
    Returns DataFrame with p-values.
    """
    results = []
    for feature, row in table.iterrows():
        counts = row.values
        if counts.sum() == 0:
            continue
        if np.count_nonzero(counts) < 2:
            results.append({"feature": feature, "p_value_ch2": np.nan})
            continue
        try:
            chi2, p, dof, _ = chi2_contingency([counts])
        except Exception:
            p = np.nan
        results.append({"feature": feature, "p_value_ch2": p})
    if not results:
        return pd.DataFrame(columns=["feature", "p_value_ch2"]).set_index("feature")
    return pd.DataFrame(results).set_index("feature")


def stage2_group_comparison(df_combined, outdir):
    print("[Stage 2] Group-Wise Comparison")

    group_sample_counts = (
        df_combined[["group", "sample"]].drop_duplicates()
        .groupby("group").size().rename("n_samples")
    )

    gene_group = (
        df_combined
        .groupby(["element symbol", "group"])["sample"]
        .nunique()
        .unstack(fill_value=0)
    )
    gene_group.to_csv(outdir / "stage2_gene_counts_by_group.csv")

    gene_group_prop = gene_group.divide(group_sample_counts, axis=1)
    gene_group_prop.to_csv(outdir / "stage2_gene_prop_by_group.csv")

    enrichment = _enrichment_test(gene_group)
    enrichment.to_csv(outdir / "stage2_gene_enrichment_pvalues.csv")

    if gene_group_prop.shape[0] > 0:
        gene_var = gene_group_prop.var(axis=1)
        top_var_genes = gene_var.sort_values(ascending=False).head(TOP_N_GENES_FOR_PLOTS).index

        plt.figure(figsize=(14, 10))
        sns.heatmap(
            gene_group_prop.loc[top_var_genes],
            annot=False,
            cmap="viridis",
            cbar_kws={"shrink": 0.7}
        )
        plt.title("Gene presence proportion by group (top variable genes)")
        plt.xlabel("Group")
        plt.ylabel("AMR gene")
        plt.xticks(rotation=45, ha="right", fontsize=9)
        plt.yticks(fontsize=8)
        plt.tight_layout()
        plt.savefig(outdir / "stage2_gene_group_heatmap.png", dpi=300, bbox_inches="tight")
        plt.close()

    class_group = (
        df_combined
        .groupby(["class", "group"])["sample"]
        .nunique()
        .unstack(fill_value=0)
    )
    class_group.to_csv(outdir / "stage2_class_counts_by_group.csv")
    class_group_prop = class_group.divide(group_sample_counts, axis=1)
    class_group_prop.to_csv(outdir / "stage2_class_prop_by_group.csv")

    n_classes, n_groups = class_group_prop.shape
    fig_w = max(8.0, 0.7 * n_groups + 4.0)
    fig_h = max(6.0, 0.5 * n_classes + 3.0)
    total_cells = n_classes * n_groups
    do_annot = total_cells <= 80

    plt.figure(figsize=(fig_w, fig_h))
    sns.heatmap(
        class_group_prop,
        annot=do_annot,
        fmt=".2f" if do_annot else "",
        cmap="magma",
        cbar_kws={"shrink": 0.8}
    )
    plt.title("AMR class presence proportion by group")
    plt.xlabel("Group")
    plt.ylabel("AMR class")
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=9)
    plt.tight_layout()
    plt.savefig(outdir / "stage2_class_group_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close()

# ----------------------------------------------------------------------
# 6. Stage 3 - Co-occurrence / Gene Interaction Networks
#     + GraphML, node attributes, core-periphery, per-group networks
# ----------------------------------------------------------------------
def compute_cooccurrence_edges(gene_matrix, min_cooc):
    """
    Compute co-occurrence edges (full matrix) for a given gene_matrix.
    Returns DataFrame with columns: gene1, gene2, cooccurrence_count.
    """
    genes = gene_matrix.columns
    coocc_edges = []

    for g1, g2 in itertools.combinations(genes, 2):
        v1 = gene_matrix[g1].values
        v2 = gene_matrix[g2].values
        co = np.sum((v1 == 1) & (v2 == 1))
        if co >= min_cooc:
            coocc_edges.append((g1, g2, co))

    edges_df = pd.DataFrame(coocc_edges, columns=["gene1", "gene2", "cooccurrence_count"])
    return edges_df


def stage3_cooccurrence_network(df_combined, gene_matrix, sample_groups, outdir):
    print("[Stage 3] Co-occurrence / Gene Interaction Networks")

    # 6.1 global co-occurrence
    edges_df = compute_cooccurrence_edges(gene_matrix, MIN_COOC_SAMPLE_COUNT)
    edges_df.to_csv(outdir / "stage3_cooccurrence_edges.csv", index=False)

    if edges_df.empty:
        print("[Stage 3] No gene pairs met the co-occurrence threshold; network not plotted.")
        return

    # 6.2 node attributes from AMR table
    # most frequent class per gene
    try:
        gene_class = (
            df_combined.groupby("element symbol")["class"]
            .agg(lambda x: x.value_counts().index[0])
        )
    except Exception:
        gene_class = pd.Series(dtype=str)

    # number of distinct isolates per gene
    gene_isolates = (
        df_combined.groupby("element symbol")["sample"]
        .nunique()
    )

    # 6.3 build full graph
    G = nx.Graph()
    for _, row in edges_df.iterrows():
        g1 = row["gene1"]
        g2 = row["gene2"]
        w = row["cooccurrence_count"]
        G.add_edge(g1, g2, weight=float(w))

    # 6.4 degree and weighted degree
    deg_unweighted = dict(G.degree())
    deg_weighted = dict(G.degree(weight="weight"))

    # 6.5 k-core (core-periphery)
    try:
        core_num = nx.core_number(G)
        max_core = max(core_num.values()) if core_num else 0
    except Exception:
        core_num = {n: 0 for n in G.nodes()}
        max_core = 0

    is_core = {n: int(core_num.get(n, 0) == max_core and max_core > 0) for n in G.nodes()}

    # 6.6 assign node attributes
    nx.set_node_attributes(G, {g: str(gene_class.get(g, "NA")) for g in G.nodes()}, "amr_class")
    nx.set_node_attributes(G, {g: int(gene_isolates.get(g, 0)) for g in G.nodes()}, "n_isolates")
    nx.set_node_attributes(G, {g: int(deg_unweighted.get(g, 0)) for g in G.nodes()}, "degree")
    nx.set_node_attributes(G, {g: float(deg_weighted.get(g, 0.0)) for g in G.nodes()}, "weighted_degree")
    nx.set_node_attributes(G, {g: int(core_num.get(g, 0)) for g in G.nodes()}, "core_number")
    nx.set_node_attributes(G, is_core, "is_core")

    # 6.7 export full graph for Cytoscape
    full_graphml = outdir / "stage3_cooccurrence_network_full.graphml"
    nx.write_graphml(G, full_graphml)
    print("[Stage 3] Wrote full co-occurrence network GraphML to:", full_graphml)

    # 6.8 subgraph of top nodes by weighted degree for plotting
    deg = deg_weighted
    sorted_nodes = sorted(deg, key=deg.get, reverse=True)
    top_nodes = sorted_nodes[:MAX_NETWORK_NODES]
    G_sub = G.subgraph(top_nodes).copy()

    if G_sub.number_of_nodes() < 2:
        print("[Stage 3] Too few nodes in subgraph for visualization.")
    else:
        # adjacency heatmap for subgraph
        nodes_order = sorted(G_sub.nodes())
        adj_matrix = nx.to_numpy_array(G_sub, nodelist=nodes_order, weight="weight")
        adj_df = pd.DataFrame(adj_matrix, index=nodes_order, columns=nodes_order)
        adj_df.to_csv(outdir / "stage3_cooccurrence_subgraph_matrix.csv")

        plt.figure(figsize=(12, 10))
        sns.heatmap(
            adj_df,
            cmap="Reds",
            square=True,
            cbar_kws={"shrink": 0.7}
        )
        plt.title("Co-occurrence counts for top {} AMR genes (subgraph)".format(len(nodes_order)))
        plt.xticks(rotation=90, fontsize=7)
        plt.yticks(rotation=0, fontsize=7)
        plt.tight_layout()
        plt.savefig(outdir / "stage3_cooccurrence_heatmap_subgraph.png", dpi=300)
        plt.close()

        # layout and network plot
        pos = nx.spring_layout(G_sub, k=0.5, seed=42)
        deg_sub = dict(G_sub.degree(weight="weight"))
        node_sizes = [200 + 80 * deg_sub[n] for n in G_sub.nodes()]
        edge_widths = [0.5 + 0.8 * d.get("weight", 1.0) for _, _, d in G_sub.edges(data=True)]

        plt.figure(figsize=(12, 12))
        nx.draw_networkx_edges(G_sub, pos, width=edge_widths, alpha=0.5)
        nx.draw_networkx_nodes(
            G_sub,
            pos,
            node_size=node_sizes,
            node_color="lightblue",
            edgecolors="black",
            linewidths=0.5,
            alpha=0.9,
        )
        labels = {n: n for n in G_sub.nodes() if deg_sub[n] >= DEGREE_LABEL_MIN}
        nx.draw_networkx_labels(G_sub, pos, labels=labels, font_size=8)

        plt.title(
            "AMR gene co-occurrence network (top {} genes by weighted degree, min co-occ >= {})"
            .format(len(G_sub.nodes()), MIN_COOC_SAMPLE_COUNT)
        )
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(outdir / "stage3_cooccurrence_network.png", dpi=300)
        plt.close()

        # write subgraph GraphML
        sub_graphml = outdir / "stage3_cooccurrence_network_subgraph.graphml"
        nx.write_graphml(G_sub, sub_graphml)
        print("[Stage 3] Wrote subgraph co-occurrence network GraphML to:", sub_graphml)

    # 6.9 per-group networks for comparison
    group_edge_rows = []
    group_values = sorted(sample_groups.dropna().unique())
    for grp in group_values:
        grp_samples = sample_groups[sample_groups == grp].index
        if len(grp_samples) < 2:
            continue
        gm_grp = gene_matrix.loc[grp_samples]
        edges_grp = compute_cooccurrence_edges(gm_grp, MIN_COOC_SAMPLE_COUNT)
        if edges_grp.empty:
            continue
        edges_grp["group"] = grp

        # save per-group edges
        safe_grp = sanitize_group_name(grp)
        edges_grp.to_csv(
            outdir / f"stage3_group_{safe_grp}_cooccurrence_edges.csv",
            index=False
        )

        # build per-group graph and export GraphML
        Gg = nx.Graph()
        for _, row in edges_grp.iterrows():
            g1 = row["gene1"]
            g2 = row["gene2"]
            w = row["cooccurrence_count"]
            Gg.add_edge(g1, g2, weight=float(w))
        # reuse node attributes from global G if present
        for n in Gg.nodes():
            if n in G.nodes:
                Gg.nodes[n].update(G.nodes[n])

        grp_graphml = outdir / f"stage3_group_{safe_grp}_network.graphml"
        nx.write_graphml(Gg, grp_graphml)

        group_edge_rows.append(edges_grp)

    if group_edge_rows:
        df_long = pd.concat(group_edge_rows, ignore_index=True)
        df_long.to_csv(outdir / "stage3_group_edges_long.csv", index=False)

        # wide format: rows gene1,gene2; columns groups; values cooccurrence_count
        df_wide = (
            df_long
            .pivot_table(
                index=["gene1", "gene2"],
                columns="group",
                values="cooccurrence_count",
                fill_value=0
            )
        )
        df_wide.to_csv(outdir / "stage3_group_edges_wide.csv")
        print("[Stage 3] Wrote per-group co-occurrence tables for comparison.")
    else:
        print("[Stage 3] No per-group co-occurrence edges above threshold were found.")

# ----------------------------------------------------------------------
# 7. Stage 4 - Virulence and AMR Link
# ----------------------------------------------------------------------
def stage4_virulence_amr_link(df_combined, gene_matrix, outdir):
    print("[Stage 4] Virulence and AMR Link")

    if "type" not in df_combined.columns:
        print("[Stage 4] 'type' column missing; skipping virulence analysis.")
        return

    type_counts = df_combined["type"].value_counts(dropna=True)
    has_vir = any(str(t).upper() == "VIRULENCE" for t in type_counts.index)
    if not has_vir:
        print("[Stage 4] No VIRULENCE entries in 'type'; skipping virulence analysis.")
        return

    df_tmp = df_combined.copy()
    df_tmp["type_norm"] = df_tmp["type"].astype(str).str.upper()
    vir_samples = (
        df_tmp[df_tmp["type_norm"] == "VIRULENCE"]["sample"]
        .drop_duplicates()
    )

    vir_flag = gene_matrix.index.isin(vir_samples)
    amr_gene_counts_per_sample = gene_matrix.sum(axis=1)
    vir_status = pd.Series(
        np.where(vir_flag, "Virulence+", "Virulence-"),
        index=gene_matrix.index,
        name="vir_status",
    )

    df_plot = pd.DataFrame({
        "amr_gene_count": amr_gene_counts_per_sample,
        "vir_status": vir_status,
    })
    df_plot.to_csv(outdir / "stage4_amr_burden_by_virulence.csv", index_label="sample")

    plt.figure(figsize=(6, 5))
    sns.boxplot(x="vir_status", y="amr_gene_count", data=df_plot)
    sns.stripplot(x="vir_status", y="amr_gene_count", data=df_plot, color="black", alpha=0.6)
    plt.ylabel("Number of AMR genes per isolate")
    plt.xlabel("")
    plt.title("AMR gene burden in virulence+ vs virulence- isolates")
    plt.tight_layout()
    plt.savefig(outdir / "stage4_amr_burden_boxplot.png", dpi=300)
    plt.close()

    df_plot.groupby("vir_status")["amr_gene_count"].describe().to_csv(
        outdir / "stage4_amr_burden_summary_stats.csv"
    )

# ----------------------------------------------------------------------
# 8. Stage 5 - Candidate Novel AMR Mutations (heuristic)
# ----------------------------------------------------------------------
def stage5_candidate_novel_mutations(df_combined, outdir):
    print("[Stage 5] Candidate Novel AMR Mutations (heuristic)")

    id_col = "% identity to reference"
    cov_col = "% coverage of reference"
    if id_col not in df_combined.columns or cov_col not in df_combined.columns:
        print("[Stage 5] Identity/coverage columns missing; skipping candidate variant detection.")
        return

    candidates = df_combined[
        (df_combined[id_col] < 98.0) | (df_combined[cov_col] < 95.0)
    ].copy()

    if candidates.empty:
        print("[Stage 5] No candidate novel hits by thresholds (<98% identity or <95% coverage).")
        return

    keep_cols = [
        "group", "sample", "element symbol", "element name",
        "class", "subclass", id_col, cov_col, "scope",
    ]
    keep_cols = [c for c in keep_cols if c in candidates.columns]
    candidates = candidates[keep_cols]
    candidates.to_csv(outdir / "stage5_candidate_novel_amr_hits.csv", index=False)

    cand_gene_counts = candidates["element symbol"].value_counts().rename("candidate_count")
    cand_gene_counts.to_csv(outdir / "stage5_candidate_counts_by_gene.csv")

    top_cand_genes = cand_gene_counts.head(TOP_N_GENES_FOR_PLOTS)
    plt.figure(figsize=(10, 6))
    sns.barplot(x=top_cand_genes.values, y=top_cand_genes.index)
    plt.xlabel("Number of candidate variant hits")
    plt.ylabel("AMR gene")
    plt.title("Top AMR genes with candidate novel variants")
    plt.tight_layout()
    plt.savefig(outdir / "stage5_candidate_genes_barplot.png", dpi=300)
    plt.close()

# ----------------------------------------------------------------------
# 9. Stage 6 - Machine Learning Without Phenotype (predict group)
# ----------------------------------------------------------------------
def stage6_write_interpretation(report_dict, top_feats, outdir):
    """
    Create a human-readable interpretation of the random forest results.
    """
    out_path = outdir / "stage6_rf_interpretation.txt"

    accuracy = report_dict.get("accuracy", None)
    macro = report_dict.get("macro avg", {})
    weighted = report_dict.get("weighted avg", {})

    lines = []
    lines.append("Random forest classification of groups from AMR gene presence")
    lines.append("================================================================")
    lines.append("")

    if accuracy is not None:
        lines.append("Overall accuracy on the held-out test set: {:.3f}".format(accuracy))
    if macro:
        lines.append(
            "Macro-averaged precision: {:.3f}, recall: {:.3f}, F1-score: {:.3f}".format(
                macro.get("precision", float("nan")),
                macro.get("recall", float("nan")),
                macro.get("f1-score", float("nan")),
            )
        )
    if weighted:
        lines.append(
            "Weighted-averaged precision: {:.3f}, recall: {:.3f}, F1-score: {:.3f}".format(
                weighted.get("precision", float("nan")),
                weighted.get("recall", float("nan")),
                weighted.get("f1-score", float("nan")),
            )
        )
    lines.append("")

    if accuracy is not None:
        if accuracy < 0.5:
            lines.append(
                "Interpretation: The classifier shows limited discriminative power; AMR gene profiles do not "
                "separate the groups strongly in this dataset, or the dataset is too small or unbalanced."
            )
        elif accuracy < 0.75:
            lines.append(
                "Interpretation: The classifier captures a moderate amount of signal in AMR gene profiles. "
                "Groups are partially separable, but there is substantial overlap in resistome composition."
            )
        else:
            lines.append(
                "Interpretation: The classifier achieves relatively high accuracy, suggesting that AMR gene "
                "profiles differ systematically between groups and carry useful information for predicting group."
            )
        lines.append("")

    lines.append("Class-wise performance (per group):")
    for label, metrics in report_dict.items():
        if label in ["accuracy", "macro avg", "weighted avg"]:
            continue
        if not isinstance(metrics, dict):
            continue
        lines.append(
            "  - Group '{}': precision {:.3f}, recall {:.3f}, F1-score {:.3f}, support {}".format(
                label,
                metrics.get("precision", float("nan")),
                metrics.get("recall", float("nan")),
                metrics.get("f1-score", float("nan")),
                metrics.get("support", 0),
            )
        )
    lines.append("")

    if top_feats is not None and not top_feats.empty:
        lines.append("Top AMR genes contributing to group discrimination (by feature importance):")
        for gene, imp in top_feats.items():
            lines.append("  - {} (importance {:.4f})".format(gene, imp))
        lines.append("")
        lines.append(
            "These genes appear to be the main determinants that differ between groups based on the current model "
            "and dataset. In a manuscript, they can be discussed as candidate group-specific AMR markers, with "
            "caution that feature importance reflects association and model-based importance, not necessarily causality."
        )
    else:
        lines.append("No feature importance information was available for interpretation.")

    with out_path.open("w") as f:
        for line in lines:
            f.write(line + "\n")

    print("[Stage 6] Wrote interpretation to:", out_path)


def stage6_ml_without_phenotype(gene_matrix, sample_groups, outdir):
    print("[Stage 6] Machine Learning Without Phenotype (predict group)")

    mask = sample_groups.notna()
    X = gene_matrix.loc[mask]
    y = sample_groups.loc[mask]

    vc = y.value_counts()
    print("[Stage 6] Group sizes before filtering:")
    print(vc.to_string())

    keep_classes = vc[vc >= 2].index
    if len(keep_classes) < 2:
        print(
            "[Stage 6] Not enough groups with >= 2 isolates "
            "(need at least 2 groups). Skipping machine-learning stage."
        )
        return

    mask2 = y.isin(keep_classes)
    X = X.loc[mask2]
    y = y.loc[mask2]

    vc2 = y.value_counts()
    print("[Stage 6] Group sizes after filtering:")
    print(vc2.to_string())

    if vc2.min() < 2:
        print(
            "[Stage 6] After filtering, at least one group still has < 2 samples; "
            "skipping machine-learning stage."
        )
        return

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=TRAIN_TEST_SPLIT,
        stratify=y,
        random_state=42,
    )

    clf = RandomForestClassifier(
        n_estimators=500,
        random_state=42,
        n_jobs=-1,
        class_weight="balanced",
    )
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)

    report = classification_report(y_test, y_pred, output_dict=True)
    pd.DataFrame(report).transpose().to_csv(outdir / "stage6_rf_classification_report.csv")

    cm = confusion_matrix(y_test, y_pred, labels=sorted(y.unique()))
    cm_df = pd.DataFrame(cm, index=sorted(y.unique()), columns=sorted(y.unique()))
    cm_df.to_csv(outdir / "stage6_confusion_matrix.csv")

    plt.figure(figsize=(9, 7))
    sns.heatmap(
        cm_df,
        annot=True,
        fmt="d",
        cmap="Blues",
        annot_kws={"size": 11}
    )
    plt.xlabel("Predicted group")
    plt.ylabel("True group")
    plt.title("Random forest - group prediction from AMR genes")
    plt.tight_layout()
    plt.savefig(outdir / "stage6_confusion_matrix_heatmap.png", dpi=300)
    plt.close()

    importances = pd.Series(clf.feature_importances_, index=X.columns)
    importances_sorted = importances.sort_values(ascending=False)

    importances_sorted.to_csv(
        outdir / "stage6_feature_importances_all_genes.csv",
        header=["importance"]
    )

    top_feats = importances_sorted.head(TOP_N_GENES_FOR_PLOTS)
    top_feats.to_csv(
        outdir / "stage6_top_feature_importances.csv",
        header=["importance"]
    )

    plt.figure(figsize=(9, 8))
    sns.barplot(x=top_feats.values, y=top_feats.index)
    plt.xlabel("Feature importance")
    plt.ylabel("AMR gene")
    plt.title("Top AMR genes contributing to group prediction")
    plt.tight_layout()
    plt.savefig(outdir / "stage6_feature_importances_barplot.png", dpi=300)
    plt.close()

    if X.shape[1] >= 2:
        pca = PCA(n_components=2, random_state=42)
        X_pca = pca.fit_transform(X)
        df_pca = pd.DataFrame(X_pca, columns=["PC1", "PC2"], index=X.index)
        df_pca["group"] = y
        df_pca.to_csv(outdir / "stage6_pca_coordinates.csv")

        plt.figure(figsize=(9, 7))
        sns.scatterplot(
            data=df_pca,
            x="PC1",
            y="PC2",
            hue="group",
            s=50,
            alpha=0.7,
            edgecolor="none"
        )
        plt.title("PCA of AMR gene presence profiles coloured by group")
        plt.tight_layout()
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
        plt.savefig(outdir / "stage6_pca_scatter.png", dpi=300, bbox_inches="tight")
        plt.close()

    stage6_write_interpretation(report, top_feats, outdir)

# ----------------------------------------------------------------------
# 10. Figure captions and methods template
# ----------------------------------------------------------------------
def write_figure_captions_and_methods(df_combined, outdir):
    groups = sorted(df_combined["group"].dropna().unique())
    classes = sorted(df_combined["class"].dropna().unique())
    n_samples = df_combined["sample"].nunique()

    text_path = outdir / "figure_captions_and_methods_template.txt"
    with text_path.open("w") as f:
        f.write("=== Suggested Figure Captions ===\n\n")

        f.write("Figure 1. Distribution of AMR classes across all isolates.\n")
        f.write(
            "Bar plot showing the total number of AMRFinder hits per antimicrobial class "
            "({}) among {} sequenced isolates.\n\n"
            .format(", ".join(classes), n_samples)
        )

        f.write("Figure 2. Heatmap of variable AMR genes by host or source group.\n")
        f.write(
            "Heatmap representing the proportion of isolates carrying each AMR determinant among groups "
            "({}). Rows correspond to genes and columns to groups.\n\n"
            .format(", ".join(groups))
        )

        f.write("Figure 3. AMR gene co-occurrence network and subgraph heatmap.\n")
        f.write(
            "Network representation of the most highly co-occurring AMR determinants across isolates. "
            "Nodes represent genes and edges are weighted by the number of isolates in which two genes "
            "co-occur. A complementary heatmap shows the co-occurrence matrix for the top genes by "
            "weighted degree. GraphML exports are provided for Cytoscape-based exploration.\n\n"
        )

        f.write("Figure 4. AMR gene burden in virulence-positive versus virulence-negative isolates.\n")
        f.write(
            "Boxplot comparing the number of distinct AMR genes per isolate between virulence-positive and "
            "virulence-negative genomes, as defined by AMRFinder 'type' annotation. Overlaid points represent "
            "individual isolates.\n\n"
        )

        f.write("Figure 5. Candidate novel AMR variants based on reduced identity or coverage.\n")
        f.write(
            "Bar plot summarizing AMR determinants with putative novel alleles, defined as hits exhibiting "
            "less than 98 percent identity or less than 95 percent coverage to their closest reference sequence.\n\n"
        )

        f.write("Figure 6. Random forest classification of host or source group using AMR gene presence.\n")
        f.write(
            "Performance of a random forest classifier trained to distinguish host or source groups using binary "
            "AMR gene presence as features, including confusion matrix, feature importance ranking, and a "
            "two-dimensional PCA projection of AMR profiles.\n\n"
        )

        f.write("=== Methods Section Template ===\n\n")

        f.write("Whole-genome sequencing and AMR gene detection.\n")
        f.write(
            "Short-read whole-genome sequencing data were processed using standard quality control and assembly "
            "procedures (details to be added here). Antimicrobial resistance determinants were identified using "
            "AMRFinderPlus on assembled contigs, and the resulting combined summary table was parsed in Python "
            "to extract gene-level information, AMR classes, and functional annotations.\n\n"
        )

        f.write("Grouping of isolates and integration of metadata.\n")
        f.write(
            "Isolates were assigned to host or source groups based on a sample manifest file (sample.tsv) and the "
            "group column in the AMRFinder summary. The final analysis included {} isolates distributed across "
            "the following groups: {}.\n\n"
            .format(n_samples, ", ".join(groups))
        )

        f.write("Core AMR gene summary and group-wise comparison.\n")
        f.write(
            "For each AMR determinant, the total number of hits and the number of distinct isolates carrying the "
            "gene were computed. AMR determinants were further classified into classes using AMRFinder annotations. "
            "Group-wise proportions of gene presence were derived by normalizing gene counts by the number of "
            "isolates per group. Statistical evidence for differential distribution across groups was explored "
            "using chi-square tests applied to contingency tables of gene presence by group.\n\n"
        )

        f.write("Co-occurrence analysis and network construction.\n")
        f.write(
            "Binary presence-absence matrices of AMR genes per isolate were used to quantify gene-gene co-occurrence. "
            "Co-occurrence counts corresponded to the number of isolates carrying both genes. Gene pairs exceeding "
            "a minimum co-occurrence threshold were represented as edges in an undirected network, with nodes "
            "representing AMR genes and edge weights reflecting co-occurrence counts. To aid interpretation, the "
            "network visualization focused on the most highly connected genes, and GraphML exports were generated "
            "for interactive exploration in Cytoscape.\n\n"
        )

        f.write("Detection of candidate novel AMR variants.\n")
        f.write(
            "Putative novel AMR variants were identified heuristically from AMRFinderPlus outputs by flagging hits "
            "with percent identity below 98 percent or query coverage below 95 percent relative to the closest "
            "reference sequence. These candidate variants were summarized at the gene level and stratified by "
            "host or source group to prioritize loci for further validation.\n\n"
        )

        f.write("Machine-learning analysis of AMR profiles.\n")
        f.write(
            "To assess whether AMR gene profiles carry a host or source-specific signal, a random forest classifier "
            "was trained to predict the group label from binary AMR gene presence. The dataset was split into training "
            "and test sets using stratified sampling, and classifier performance was evaluated using a held-out test "
            "set (classification report and confusion matrix). Feature importance values were used to rank AMR genes "
            "according to their contribution to group discrimination, and a full gene-level importance table was "
            "exported for downstream interpretation. A principal component analysis (PCA) of gene presence profiles "
            "was used to visualize clustering of isolates by group.\n\n"
        )

    print("[INFO] Wrote figure captions and methods template to:", text_path)

# ----------------------------------------------------------------------
# 11. All-stages interpretation text
# ----------------------------------------------------------------------
def write_all_stages_interpretation(df_combined, outdir):
    """
    Generate a single text file summarizing the interpretation for stages 1-6.
    """
    path_out = outdir / "all_stages_interpretation.txt"
    n_samples = df_combined["sample"].nunique()
    groups = sorted(df_combined["group"].dropna().unique())
    classes = sorted(df_combined["class"].dropna().unique())
    n_genes = df_combined["element symbol"].nunique()

    lines = []
    lines.append("Interpretation of AMR analysis pipeline (Stages 1-6)")
    lines.append("===================================================")
    lines.append("")
    lines.append("Dataset overview:")
    lines.append(
        "The analysis included {} isolates, grouped into {} categories ({}), "
        "with a total of {} distinct AMR determinants detected based on AMRFinder annotations."
        .format(n_samples, len(groups), ", ".join(groups) if groups else "no groups available", n_genes)
    )
    lines.append(
        "The main AMR classes captured in this dataset are: {}."
        .format(", ".join(classes) if classes else "not available")
    )
    lines.append("")

    # Stage 1
    lines.append("Stage 1 - Core AMR gene and class summary")
    try:
        gene_counts = pd.read_csv(outdir / "stage1_gene_counts.csv", index_col=0)["count"]
        class_counts = pd.read_csv(outdir / "stage1_class_counts.csv", index_col=0)["count"]
        top_genes = gene_counts.sort_values(ascending=False).head(5)
        top_classes = class_counts.sort_values(ascending=False).head(5)

        lines.append(
            "Stage 1 quantified how frequently each AMR gene and AMR class was detected across the collection. "
            "The barplot 'stage1_top_genes_barplot.png' highlights the most commonly detected AMR determinants, "
            "while 'stage1_class_barplot.png' summarizes the distribution of AMR classes."
        )
        lines.append(
            "The most frequent AMR genes in this dataset include: {}."
            .format(", ".join(["{} (n={})".format(g, int(c)) for g, c in top_genes.items()]))
        )
        lines.append(
            "At the class level, the dominant AMR categories are: {}."
            .format(", ".join(["{} (n={})".format(cn, int(c)) for cn, c in top_classes.items()]))
        )
    except Exception:
        lines.append(
            "Stage 1 results could not be summarized from the output files. "
            "Please check that 'stage1_gene_counts.csv' and 'stage1_class_counts.csv' were generated correctly."
        )
    lines.append("")

    # Stage 2
    lines.append("Stage 2 - Group-wise comparison of AMR profiles")
    try:
        class_group_prop = pd.read_csv(outdir / "stage2_class_prop_by_group.csv", index_col=0)
        lines.append(
            "Stage 2 compared AMR gene and class distributions across groups. The heatmap "
            "'stage2_gene_group_heatmap.png' focuses on the most variable genes across groups, and "
            "'stage2_class_group_heatmap.png' shows how frequently each AMR class is observed in each group."
        )
        class_var = class_group_prop.var(axis=1)
        top_var_classes = class_var.sort_values(ascending=False).head(5)
        lines.append(
            "Across groups, some AMR classes show clear differences in prevalence. The classes with the largest "
            "between-group variation include: {}."
            .format(", ".join(["{} (variance {:.3f})".format(cn, v) for cn, v in top_var_classes.items()]))
        )
    except Exception:
        lines.append(
            "Stage 2 heatmaps and proportions could not be summarized from the output files. "
            "Please check that 'stage2_gene_prop_by_group.csv' and 'stage2_class_prop_by_group.csv' were generated."
        )
    lines.append("")

    # Stage 3
    lines.append("Stage 3 - Co-occurrence and interaction between AMR genes")
    try:
        edges_df = pd.read_csv(outdir / "stage3_cooccurrence_edges.csv")
        n_pairs = edges_df.shape[0]
        genes_in_pairs = set(edges_df["gene1"]).union(set(edges_df["gene2"]))
        lines.append(
            "Stage 3 explored how AMR genes co-occur within the same isolates. Gene pairs that co-occur in at least "
            "{} isolates were retained, resulting in {} qualifying pairs involving {} distinct genes. The file "
            "'stage3_cooccurrence_edges.csv' lists these gene pairs and their co-occurrence counts."
            .format(MIN_COOC_SAMPLE_COUNT, n_pairs, len(genes_in_pairs))
        )
        lines.append(
            "The network visualization 'stage3_cooccurrence_network.png' focuses on the most highly connected "
            "determinants, while 'stage3_cooccurrence_heatmap_subgraph.png' shows the corresponding "
            "co-occurrence matrix. GraphML files ('stage3_cooccurrence_network_full.graphml' and "
            "'stage3_cooccurrence_network_subgraph.graphml') enable interactive exploration in Cytoscape, where "
            "node attributes such as AMR class, number of isolates, degree, and k-core membership can be inspected."
        )
        lines.append(
            "Per-group co-occurrence networks captured in 'stage3_group_edges_long.csv' and "
            "'stage3_group_edges_wide.csv' allow comparison of which AMR gene pairs tend to co-occur within each "
            "group and which co-occurrences may be group-specific."
        )
    except Exception:
        lines.append(
            "Stage 3 network outputs could not be summarized from 'stage3_cooccurrence_edges.csv'. "
            "It is possible that no gene pairs reached the co-occurrence threshold."
        )
    lines.append("")

    # Stage 4
    lines.append("Stage 4 - Relationship between virulence status and AMR burden")
    try:
        stats = pd.read_csv(outdir / "stage4_amr_burden_summary_stats.csv", index_col=0)
        means = stats.loc["mean"]
        if "Virulence+" in means.index and "Virulence-" in means.index:
            mean_pos = means["Virulence+"]
            mean_neg = means["Virulence-"]
            lines.append(
                "Stage 4 compared the number of distinct AMR genes per isolate between virulence-positive and "
                "virulence-negative genomes. The boxplot 'stage4_amr_burden_boxplot.png' visualizes this comparison."
            )
            lines.append(
                "On average, virulence-positive isolates carried {:.2f} AMR genes, whereas virulence-negative isolates "
                "carried {:.2f} AMR genes."
                .format(mean_pos, mean_neg)
            )
        else:
            lines.append(
                "Stage 4 summary statistics are available, but virulence-positive and virulence-negative groups "
                "could not be clearly distinguished from 'stage4_amr_burden_summary_stats.csv'."
            )
    except Exception:
        lines.append(
            "Stage 4 analysis was skipped or the required files were not generated (for example, if no VIRULENCE "
            "entries were present in the AMRFinder 'type' field)."
        )
    lines.append("")

    # Stage 5
    lines.append("Stage 5 - Candidate novel AMR variants")
    try:
        cand_hits = pd.read_csv(outdir / "stage5_candidate_novel_amr_hits.csv")
        cand_counts = pd.read_csv(outdir / "stage5_candidate_counts_by_gene.csv", index_col=0)["candidate_count"]
        n_cand_hits = cand_hits.shape[0]
        n_cand_genes = cand_counts.shape[0]
        top_cand = cand_counts.sort_values(ascending=False).head(5)

        lines.append(
            "Stage 5 flagged potential novel or divergent AMR variants using simple thresholds on AMRFinder percent "
            "identity (< 98 percent) and/or coverage (< 95 percent). A total of {} candidate hits were identified "
            "across {} AMR genes. These are listed in 'stage5_candidate_novel_amr_hits.csv', and gene-level "
            "frequencies are summarized in 'stage5_candidate_counts_by_gene.csv'."
            .format(n_cand_hits, n_cand_genes)
        )
        lines.append(
            "The AMR genes with the highest number of candidate divergent hits include: {}."
            .format(", ".join(["{} (n={})".format(g, int(c)) for g, c in top_cand.items()]))
        )
    except Exception:
        lines.append(
            "Stage 5 did not identify any hits meeting the candidate variant thresholds, or the corresponding "
            "CSV files are not available."
        )
    lines.append("")

    # Stage 6
    lines.append("Stage 6 - Machine-learning analysis of AMR profiles")
    try:
        rep_df = pd.read_csv(outdir / "stage6_rf_classification_report.csv", index_col=0)
        accuracy = None
        if "accuracy" in rep_df.index:
            try:
                accuracy = float(rep_df.loc["accuracy", "precision"])
            except Exception:
                accuracy = None

        lines.append(
            "Stage 6 used a random forest classifier to predict the group label from binary AMR gene presence across "
            "isolates. The confusion matrix is displayed in 'stage6_confusion_matrix_heatmap.png', while "
            "'stage6_feature_importances_barplot.png' and 'stage6_feature_importances_all_genes.csv' summarize which "
            "AMR determinants contributed most strongly to group discrimination. The PCA plot "
            "'stage6_pca_scatter.png' provides a low-dimensional visualization of AMR profiles coloured by group."
        )

        if accuracy is not None:
            lines.append(
                "On the held-out test set, the overall classification accuracy was {:.3f}."
                .format(accuracy)
            )

        try:
            top_feats = pd.read_csv(
                outdir / "stage6_top_feature_importances.csv",
                index_col=0
            )["importance"]
            top_names = list(top_feats.index[:5])
            lines.append(
                "Based on feature importance, a small set of AMR genes drives most of the discriminatory power of "
                "the classifier. Examples include: {}."
                .format(", ".join(top_names))
            )
        except Exception:
            lines.append(
                "Top feature importance values for Stage 6 could not be read from 'stage6_top_feature_importances.csv'."
            )

    except Exception:
        lines.append(
            "Stage 6 random forest outputs could not be summarized. Please check that the classification report, "
            "confusion matrix, and feature importance files were generated successfully."
        )
    lines.append("")

    with path_out.open("w") as f:
        for line in lines:
            f.write(line + "\n")

    print("[INFO] Wrote all-stages interpretation to:", path_out)

# ----------------------------------------------------------------------
# 12. Cytoscape instructions
# ----------------------------------------------------------------------
def write_cytoscape_instructions(outdir):
    """
    Write a plain-text guide explaining how to use the Stage 3 outputs in Cytoscape.
    """
    path_out = outdir / "instructions_for_cytoscape.txt"
    lines = []
    lines.append("Instructions: Visualizing AMR gene co-occurrence in Cytoscape")
    lines.append("==============================================================")
    lines.append("")
    lines.append("Files produced by Stage 3 that are relevant for Cytoscape:")
    lines.append("  - stage3_cooccurrence_network_full.graphml")
    lines.append("      Full AMR gene co-occurrence network; nodes are AMR genes, edges weighted by")
    lines.append("      the number of isolates where the two genes co-occur.")
    lines.append("  - stage3_cooccurrence_network_subgraph.graphml")
    lines.append("      Subnetwork of the most highly connected genes (top nodes by weighted degree).")
    lines.append("  - stage3_group_<GROUP>_network.graphml")
    lines.append("      Per-group co-occurrence networks (GROUP is the group label with non-alphanumeric")
    lines.append("      characters replaced by underscores).")
    lines.append("  - stage3_group_edges_long.csv and stage3_group_edges_wide.csv")
    lines.append("      Long and wide tables that summarize edge weights across groups.")
    lines.append("")
    lines.append("A. Importing the global network into Cytoscape")
    lines.append("---------------------------------------------")
    lines.append("1. Open Cytoscape.")
    lines.append("2. Go to: File -> Import -> Network -> File...")
    lines.append("3. Select 'stage3_cooccurrence_network_full.graphml' and click 'Open'.")
    lines.append("4. Cytoscape will load the network with node and edge attributes:")
    lines.append("      Node attributes:")
    lines.append("        - amr_class: AMR class most frequently associated with the gene.")
    lines.append("        - n_isolates: number of isolates carrying this gene.")
    lines.append("        - degree: unweighted degree in the co-occurrence graph.")
    lines.append("        - weighted_degree: sum of edge weights connected to the node.")
    lines.append("        - core_number: k-core index for core-periphery structure.")
    lines.append("        - is_core: 1 for nodes in the highest k-core, 0 otherwise.")
    lines.append("      Edge attributes:")
    lines.append("        - weight: cooccurrence_count (number of isolates with both genes).")
    lines.append("")
    lines.append("B. Basic styling suggestions")
    lines.append("----------------------------")
    lines.append("1. Use the 'Style' panel in Cytoscape to map attributes to visual properties:")
    lines.append("   - Map 'weighted_degree' to node size (larger nodes for highly connected genes).")
    lines.append("   - Map 'amr_class' to node fill color (different colors for different classes).")
    lines.append("   - Map 'weight' to edge width (thicker edges for stronger co-occurrence).")
    lines.append("   - Map 'is_core' to node border width or color (highlight core nodes).")
    lines.append("2. Apply a force-directed layout, for example:")
    lines.append("   - Layout -> Prefuse Force Directed")
    lines.append("   - or Layout -> yFiles Organic (if available).")
    lines.append("")
    lines.append("C. Exploring the top subnetwork")
    lines.append("-------------------------------")
    lines.append("1. Import 'stage3_cooccurrence_network_subgraph.graphml' in the same way as above.")
    lines.append("2. This smaller network focuses on the most connected AMR genes and is useful for")
    lines.append("   generating clear publication-quality figures.")
    lines.append("")
    lines.append("D. Per-group network comparison")
    lines.append("-------------------------------")
    lines.append("1. For each group, a file named 'stage3_group_<GROUP>_network.graphml' is created.")
    lines.append("   Here, <GROUP> is your group label with non-alphanumeric characters replaced by")
    lines.append("   underscores (for example, 'HEALTH_TRACHEAL' or 'CLINICAL_HOSPITAL').")
    lines.append("2. Import each per-group graph with:")
    lines.append("      File -> Import -> Network -> File...")
    lines.append("3. You can compare:")
    lines.append("   - which nodes (AMR genes) appear only in certain groups;")
    lines.append("   - which edges (gene pairs) are present or absent in different groups;")
    lines.append("   - how node degree, weighted degree, or core_number differ across groups.")
    lines.append("4. The tables 'stage3_group_edges_long.csv' and 'stage3_group_edges_wide.csv' can")
    lines.append("   be opened in Excel or R/Python to quantify differences between groups.")
    lines.append("")
    lines.append("E. Exporting figures from Cytoscape")
    lines.append("-------------------------------")
    lines.append("1. After adjusting layout and style, export figures using:")
    lines.append("      File -> Export -> Network View as Graphics...")
    lines.append("2. Choose formats such as PNG, PDF, or SVG for inclusion in publications.")
    lines.append("")
    lines.append("This guide is generic; you can adapt node color schemes, layouts, and thresholds")
    lines.append("to match the focus of your analysis (for example, highlighting specific AMR classes")
    lines.append("or genes identified as important in Stage 6).")

    with path_out.open("w") as f:
        for line in lines:
            f.write(line + "\n")

    print("[INFO] Wrote Cytoscape instructions to:", path_out)

# ----------------------------------------------------------------------
# 13. Run log and directory tree (ASCII-only tree)
# ----------------------------------------------------------------------
def write_run_log(outdir):
    log_path = outdir / "run_log_and_file_structure.txt"

    def _tree(root, prefix=""):
        lines = []
        entries = sorted(list(root.iterdir()), key=lambda p: (p.is_file(), p.name))
        for i, p in enumerate(entries):
            connector = "`-- " if i == len(entries) - 1 else "+-- "
            lines.append(prefix + connector + p.name)
            if p.is_dir():
                extension = "    " if i == len(entries) - 1 else "|   "
                lines.extend(_tree(p, prefix + extension))
        return lines

    with log_path.open("w") as f:
        f.write("Command run:\n")
        f.write("  " + " ".join(sys.argv) + "\n\n")
        f.write("Output directory structure:\n")
        f.write(str(outdir.resolve()) + "\n")
        for line in _tree(outdir):
            f.write(line + "\n")

    print("[INFO] Wrote run log and file structure to:", log_path)

# ----------------------------------------------------------------------
# 14. Main
# ----------------------------------------------------------------------
def main():
    sample_group_override = load_sample_groups_from_tsv(SAMPLE_TSV)
    df_combined = load_combined_amr(AMR_SUMMARY_XLSX, COMBINED_SHEET, sample_group_override)
    gene_matrix, sample_groups = make_gene_presence_matrix(df_combined)

    stage1_core_summary(df_combined, OUTDIR)
    stage2_group_comparison(df_combined, OUTDIR)
    stage3_cooccurrence_network(df_combined, gene_matrix, sample_groups, OUTDIR)
    stage4_virulence_amr_link(df_combined, gene_matrix, OUTDIR)
    stage5_candidate_novel_mutations(df_combined, OUTDIR)
    stage6_ml_without_phenotype(gene_matrix, sample_groups, OUTDIR)
    write_figure_captions_and_methods(df_combined, OUTDIR)
    write_all_stages_interpretation(df_combined, OUTDIR)
    write_cytoscape_instructions(OUTDIR)
    write_run_log(OUTDIR)

    print("[DONE] All results written to:", OUTDIR.resolve())


if __name__ == "__main__":
    main()
