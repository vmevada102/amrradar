#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
amr_network_mge_enrichment.py

Auto-discovers AMRFinder TSVs and Prokka GFFs under:
  results/tools/*_AMRFinder/*/<sample>/*.tsv
  results/tools/*_prokka/*/<sample>/*.gff

Combines them, then performs:
  - co-occurrence networks (binary & abundance via Spearman)
  - node centralities, Louvain communities
  - interactive Plotly network (colored by dominant group)
  - per-group subgraphs (network_<group>.gexf + interactive HTML)
  - simple MGE association mapping from AMR TSVs (contig/start/end heuristics)
  - mechanism/drug-class enrichment by group (Fisher + BH)

Outputs are written into a new numbered directory:
  results/tools/<N>_AMRNetwork
where N = max numeric prefix in results/tools/* + 1 (so numbering follows other tool folders).

Usage:
    python3 amr_network_mge_enrichment.py results/OneHealthAMR_AMRFinder_summary.xlsx

Optional flags:
    --min-rho FLOAT        (default 0.6) threshold for edges
    --min-samples INT      (default 3) min samples for per-group subgraphs
    --no-combine           don't auto-combine found AMR/GFF files
    --amr-tsv PATH         explicit combined AMR TSV (overrides auto-discovery)
    --gff PATH             explicit combined GFF (overrides auto-discovery)
    --keep-temp            keep combined temp files
    --force-outdir PATH    use exact output dir instead of auto-numbered results/tools/*_AMRNetwork
"""

from pathlib import Path
import argparse
import sys
import os
import glob
import csv
import math
import traceback
import json

import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as smm

import networkx as nx
try:
    import community as community_louvain
except Exception:
    community_louvain = None

# interactive plotting
try:
    import plotly.graph_objects as go
    import plotly.express as px
    PLOTLY = True
except Exception:
    PLOTLY = False

# -------------------------
# Utilities
# -------------------------
def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def get_next_numbered_outdir(base_root: Path, suffix: str):
    """
    Scan base_root for directories with a leading integer prefix (e.g. '1_qc', '9_AMRStatistics').
    Return base_root / '<max+1>_<suffix>'.
    """
    ensure_dir(base_root)
    max_num = 0
    for d in base_root.iterdir():
        if not d.is_dir():
            continue
        name = d.name
        parts = name.split('_', 1)
        if parts and parts[0].isdigit():
            try:
                num = int(parts[0])
                if num > max_num:
                    max_num = num
            except Exception:
                pass
    next_num = max_num + 1
    return base_root / f"{next_num}_{suffix}"

def find_amr_tsvs_auto():
    patterns = [
        "results/tools/*_AMRFinder/*/*amrfinder*.tsv",
        "results/tools/*_AMRFinder/*/*amrfinder*.txt",
        "results/tools/*_AMRFinder/*/*amr*.tsv",
        "results/tools/*_AMRFinder/*/*.tsv",
    ]
    found = []
    for p in patterns:
        found.extend(glob.glob(p))
    return sorted(set(found))

def find_prokka_gffs_auto():
    patterns = [
        "results/tools/*_prokka/*/*.gff",
        "results/tools/*_prokka/*/*.gff3",
    ]
    found = []
    for p in patterns:
        found.extend(glob.glob(p))
    return sorted(set(found))

def write_combined_amr(amr_files, outpath):
    if not amr_files:
        return None
    outpath.parent.mkdir(parents=True, exist_ok=True)
    header_written = False
    with open(outpath, 'w', newline='', encoding='utf-8') as fout:
        writer = None
        for f in sorted(amr_files):
            try:
                with open(f, 'r', encoding='utf-8', errors='replace') as fh:
                    lines = fh.readlines()
                if not lines:
                    continue
                header_line = None
                start_idx = 0
                for i, l in enumerate(lines):
                    if l.strip() == '':
                        continue
                    if l.startswith('#'):
                        continue
                    header_line = l.rstrip('\n')
                    start_idx = i
                    break
                if header_line is None:
                    continue
                cols = header_line.strip().split('\t')
                sample = Path(f).stem
                parent = Path(f).parent.name
                # prefer parent dir name if it looks like the sample folder
                if parent and parent not in ('.', 'amrfinder_out', 'AMRFinder'):
                    sample = parent
                if not header_written:
                    cols_out = cols + ['sample']
                    writer = csv.writer(fout, delimiter='\t')
                    writer.writerow(cols_out)
                    header_written = True
                for l in lines[start_idx+1:]:
                    if l.strip() == '':
                        continue
                    if l.startswith('#'):
                        continue
                    row = l.rstrip('\n').split('\t')
                    if len(row) < len(cols):
                        row += [''] * (len(cols) - len(row))
                    row_out = row + [sample]
                    writer.writerow(row_out)
            except Exception as e:
                print(f"Warning: failed to read/append {f}: {e}", file=sys.stderr)
    print(f"Combined AMR TSV written to: {outpath}")
    return outpath

def write_combined_gff(gff_files, outpath):
    if not gff_files:
        return None
    outpath.parent.mkdir(parents=True, exist_ok=True)
    header_written = False
    with open(outpath, 'w', encoding='utf-8') as outf:
        for f in sorted(gff_files):
            sample = Path(f).parent.name
            try:
                with open(f, 'r', encoding='utf-8', errors='replace') as fh:
                    for line in fh:
                        if line.startswith('##'):
                            if not header_written:
                                outf.write(line)
                            continue
                        if line.startswith('#'):
                            if not header_written:
                                outf.write(line)
                            else:
                                continue
                        else:
                            parts = line.rstrip('\n').split('\t')
                            if len(parts) < 9:
                                continue
                            seqid = parts[0]
                            seqid_pref = f"{sample}|{seqid}"
                            attrs = parts[8]
                            attrs = attrs.replace('ID=', f'ID={sample}|')
                            attrs = attrs.replace('Parent=', f'Parent={sample}|')
                            parts[0] = seqid_pref
                            parts[8] = attrs
                            outf.write('\t'.join(parts) + '\n')
                header_written = True
            except Exception as e:
                print(f"Warning: failed to read/append GFF {f}: {e}", file=sys.stderr)
    print(f"Combined (prefixed) GFF written to: {outpath}")
    return outpath

# -------------------------
# Correlation / network helpers
# -------------------------
def compute_spearman_network(mat_df, min_rho=0.6, max_p=0.01):
    if mat_df is None or mat_df.shape[1] < 2:
        return pd.DataFrame(columns=['class1','class2','rho','p']), pd.DataFrame()
    cols = list(mat_df.columns)
    rows = []
    n = len(cols)
    for i in range(n):
        for j in range(i+1, n):
            a = mat_df.iloc[:, i].values
            b = mat_df.iloc[:, j].values
            try:
                rho, p = stats.spearmanr(a, b, nan_policy='omit')
            except Exception:
                rho, p = np.nan, np.nan
            if pd.notna(rho) and abs(rho) >= min_rho and (math.isnan(p) or p <= max_p):
                rows.append((cols[i], cols[j], float(rho), float(p if not math.isnan(p) else np.nan)))
    edges_df = pd.DataFrame(rows, columns=['class1','class2','rho','p'])
    corr_df = mat_df.corr(method='spearman')
    return edges_df, corr_df

def build_network_from_edges(edges_df):
    G = nx.Graph()
    for _, r in edges_df.iterrows():
        u = r['class1']; v = r['class2']; w = float(r['rho'])
        G.add_node(u); G.add_node(v)
        G.add_edge(u, v, weight=w)
    return G

def compute_node_stats(G, meta_df=None, binary_df=None):
    """
    Compute node centrality stats and optionally annotate each node (AMR class/gene)
    with dominant group presence if metadata and binary matrices are available.
    """
    if G is None or G.number_of_nodes() == 0:
        return pd.DataFrame(columns=['node','degree','betweenness','closeness'])
    
    deg = dict(G.degree(weight='weight'))
    bet = nx.betweenness_centrality(G, weight='weight')
    clos = nx.closeness_centrality(G)
    rows = []

    group_info = {}
    if meta_df is not None and binary_df is not None and not binary_df.empty:
        shared_samples = binary_df.index.intersection(meta_df.index)
        for feat in G.nodes():
            if feat in binary_df.columns:
                subset = binary_df.loc[shared_samples, feat]
                if subset.empty:
                    continue
                grp_counts = subset.groupby(meta_df.loc[shared_samples, 'group']).sum()
                if not grp_counts.empty:
                    grp_counts = grp_counts[grp_counts > 0]
                    if not grp_counts.empty:
                        dominant_group = grp_counts.idxmax()
                        group_info[feat] = {
                            'dominant_group': dominant_group,
                            'group_counts': grp_counts.to_dict()
                        }

    for n in G.nodes():
        row = {
            'node': n,
            'degree': deg.get(n, 0),
            'betweenness': bet.get(n, 0),
            'closeness': clos.get(n, 0)
        }
        if n in group_info:
            row.update({
                'dominant_group': group_info[n]['dominant_group'],
                'group_counts': json.dumps(group_info[n]['group_counts'])
            })
        else:
            row.update({'dominant_group': None, 'group_counts': None})
        rows.append(row)

    df = pd.DataFrame(rows).sort_values('degree', ascending=False)
    return df

def run_louvain(G):
    if community_louvain is None or G.number_of_nodes() == 0:
        return {}
    partition = community_louvain.best_partition(G, weight='weight')
    return partition

def export_gexf(G, outpath, node_stats_df=None):
    try:
        # attach node attributes from node_stats_df if provided
        if node_stats_df is not None and not node_stats_df.empty:
            attr_map = node_stats_df.set_index('node').to_dict('index')
            for n in G.nodes():
                if n in attr_map:
                    for k,v in attr_map[n].items():
                        G.nodes[n][k] = v
        nx.write_gexf(G, str(outpath))
        print(f"Wrote GEXF: {outpath}")
    except Exception as e:
        print("Failed to write GEXF:", e, file=sys.stderr)

def interactive_network_plot(G, outpath, node_stats_df=None, title="AMR co-occurrence network (colored by dominant group)"):
    if not PLOTLY or G is None or G.number_of_nodes() == 0:
        return None
    pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)
    edge_x, edge_y = [], []
    for u, v, d in G.edges(data=True):
        x0, y0 = pos[u]; x1, y1 = pos[v]
        edge_x += [x0, x1, None]; edge_y += [y0, y1, None]
    edge_trace = go.Scatter(x=edge_x, y=edge_y, mode='lines', line=dict(width=0.5, color='#aaa'), hoverinfo='none')

    # build color palette for groups
    group_colors = {}
    palette = px.colors.qualitative.Plotly if PLOTLY else ['#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A','#19D3F3','#FF6692']
    if node_stats_df is not None and 'dominant_group' in node_stats_df.columns:
        unique_groups = node_stats_df['dominant_group'].dropna().unique()
        for i, g in enumerate(unique_groups):
            group_colors[g] = palette[i % len(palette)]

    node_x, node_y, node_text, node_size, node_color = [], [], [], [], []
    for n in G.nodes():
        x, y = pos[n]
        node_x.append(x); node_y.append(y)
        node_text.append(n)
        deg = G.degree(n, weight='weight')
        node_size.append(max(6, min(28, int(deg * 2))))
        color = '#888'
        if node_stats_df is not None and 'dominant_group' in node_stats_df.columns:
            dfn = node_stats_df[node_stats_df['node'] == n]
            if not dfn.empty and pd.notna(dfn.iloc[0].get('dominant_group')):
                color = group_colors.get(dfn.iloc[0]['dominant_group'], '#888')
        node_color.append(color)

    node_trace = go.Scatter(x=node_x, y=node_y, mode='markers+text', text=node_text, textposition='top center',
                            hoverinfo='text', marker=dict(size=node_size, color=node_color, line=dict(width=0.5)))
    layout = go.Layout(title=title, showlegend=False, margin=dict(l=10,r=10,b=10,t=50))
    fig = go.Figure(data=[edge_trace, node_trace], layout=layout)
    fig.write_html(str(outpath))
    print(f"Wrote interactive network HTML: {outpath}")
    return outpath

# -------------------------
# MGE mapping & mechanism enrichment
# -------------------------
def map_genes_to_mge_from_combined_amr(df_combined):
    if df_combined is None or df_combined.empty:
        return pd.DataFrame()
    cols = df_combined.columns.str.lower().tolist()
    seq_col = None; start_col=None; end_col=None; gene_col=None; sample_col=None; note_col=None
    for c in df_combined.columns:
        low = c.lower()
        if low in ('seqname','sequence_name','contig','target_name','seqid','refseq'):
            seq_col = c
        if low in ('start','startpos','location_start'):
            start_col = c
        if low in ('end','endpos','location_end'):
            end_col = c
        if 'gene' in low or 'symbol' in low or 'allele' in low or 'name' in low:
            if low not in ('gene_count','geneid'):
                gene_col = gene_col or c
        if low in ('sample','isolate','isolate_name'):
            sample_col = c
        if low in ('product','note','element','description','class'):
            note_col = note_col or c
    rows=[]
    for _, r in df_combined.iterrows():
        sample = r[sample_col] if sample_col and sample_col in r.index else None
        seq = r[seq_col] if seq_col and seq_col in r.index else None
        start = r[start_col] if start_col and start_col in r.index else None
        end = r[end_col] if end_col and end_col in r.index else None
        gene = r[gene_col] if gene_col and gene_col in r.index else None
        note = r[note_col] if note_col and note_col in r.index else None
        rows.append({'sample':sample, 'gene':gene, 'contig':seq, 'start':start, 'end':end, 'note':note})
    out = pd.DataFrame(rows)
    def find_mge_annotation(s):
        if not isinstance(s, str): return None
        s2 = s.lower()
        keys = [ 'plasmid','replicon','transpos','integron','insertion sequence','is element','integrase','transposase' ]
        hits = [k for k in keys if k in s2]
        return ';'.join(hits) if hits else None
    if not out.empty:
        out['nearby_MGE'] = out['note'].astype(str).apply(find_mge_annotation)
    return out

def mechanism_enrichment(binary_df, meta_df, mapping_df, group_col='group'):
    if meta_df is None or meta_df.empty or group_col not in meta_df.columns:
        print("No metadata groups available; skipping mechanism enrichment.")
        return pd.DataFrame()
    if mapping_df is None or mapping_df.empty:
        print("No mechanism mapping provided.")
        return pd.DataFrame()
    mech_map = {}
    if 'mechanism' in mapping_df.columns:
        mech_map = mapping_df['mechanism'].to_dict()
    else:
        mech_map = mapping_df.iloc[:,0].to_dict()
    features = [f for f in binary_df.columns if f in mech_map]
    if not features:
        keys = list(mech_map.keys())
        for f in binary_df.columns:
            for k in keys:
                if str(k) in str(f):
                    features.append(f)
                    break
    if not features:
        print("No features matched mapping for enrichment; skipping.")
        return pd.DataFrame()
    mech_list=[]
    for feat in features:
        mech = mech_map.get(feat, mech_map.get(str(feat), None))
        if mech is None: continue
        mech_list.append((feat, mech))
    mech_df = pd.DataFrame(mech_list, columns=['feature','mechanism']).dropna()
    mech_groups = mech_df.groupby('mechanism')['feature'].apply(list).to_dict()
    results=[]
    groups = meta_df[group_col].unique()
    for mech, feats in mech_groups.items():
        present = (binary_df[feats].sum(axis=1) > 0).astype(int)
        for g in groups:
            a = int(((present==1) & (meta_df[group_col]==g)).sum())
            b = int(((present==0) & (meta_df[group_col]==g)).sum())
            c = int(((present==1) & (meta_df[group_col]!=g)).sum())
            d = int(((present==0) & (meta_df[group_col]!=g)).sum())
            try:
                odr, p = stats.fisher_exact([[a,b],[c,d]])
            except Exception:
                odr, p = (np.nan, np.nan)
            results.append({'mechanism':mech, 'group':g, 'a':a,'b':b,'c':c,'d':d,'oddsratio':odr, 'pvalue':p})
    if not results:
        return pd.DataFrame()
    df = pd.DataFrame(results)
    df['qvalue'] = smm.multipletests(df['pvalue'].fillna(1.0), method='fdr_bh')[1]
    return df.sort_values('qvalue')

# -------------------------
# Main
# -------------------------
def main():
    p = argparse.ArgumentParser(description="AMR co-occurrence/network + MGE mapping + mechanism enrichment (auto outputs to results/tools/*_AMRNetwork)")
    p.add_argument('xlsx', help='OneHealthAMR_AMRFinder_summary.xlsx (path)')
    p.add_argument('--min-rho', type=float, default=0.6)
    p.add_argument('--min-samples', type=int, default=3)
    p.add_argument('--no-combine', action='store_true')
    p.add_argument('--amr-tsv', default=None)
    p.add_argument('--gff', default=None)
    p.add_argument('--keep-temp', action='store_true')
    p.add_argument('--force-outdir', default=None)
    args = p.parse_args()

    # determine output dir (numbered)
    if args.force_outdir:
        outdir = Path(args.force_outdir).expanduser().resolve()
        ensure_dir(outdir)
    else:
        base_root = Path('results/tools')
        outdir = get_next_numbered_outdir(base_root, 'AMRNetwork')
        ensure_dir(outdir)
    print("OUTPUT DIRECTORY:", outdir)

    xlsx = Path(args.xlsx)
    if not xlsx.exists():
        print("ERROR: Excel input not found:", xlsx, file=sys.stderr)
        sys.exit(2)

    combined_dir = Path('amrfinder_combined')
    ensure_dir(combined_dir)
    combined_amr = combined_dir / 'combined_amrfinder.tsv'
    combined_gff = combined_dir / 'combined_prokka_prefixed.gff'

    amr_tsv_path = Path(args.amr_tsv) if args.amr_tsv else None
    gff_path = Path(args.gff) if args.gff else None

    # Auto-discover and combine unless disabled
    if not args.no_combine:
        if amr_tsv_path is None:
            amr_files = find_amr_tsvs_auto()
            if not amr_files:
                print("No AMRFinder TSVs found under results/tools/*_AMRFinder/*/*.tsv", file=sys.stderr)
            else:
                print(f"Found {len(amr_files)} AMRFinder TSVs; combining into {combined_amr}")
                write_combined_amr(amr_files, combined_amr)
                amr_tsv_path = combined_amr
        else:
            print("Using provided AMR TSV:", amr_tsv_path)

        if gff_path is None:
            gff_files = find_prokka_gffs_auto()
            if not gff_files:
                print("No Prokka GFFs found under results/tools/*_prokka/*/*.gff", file=sys.stderr)
            else:
                print(f"Found {len(gff_files)} Prokka GFFs; combining into {combined_gff}")
                write_combined_gff(gff_files, combined_gff)
                gff_path = combined_gff
        else:
            print("Using provided GFF:", gff_path)
    else:
        print("Skipping auto-combine (--no-combine). Ensure you provided --amr-tsv/--gff if needed.")

    # load combined AMR TSV if available
    df_combined = None
    if amr_tsv_path and Path(amr_tsv_path).exists():
        try:
            df_combined = pd.read_csv(amr_tsv_path, sep='\t', dtype=str, low_memory=False)
            print("Loaded combined AMR TSV rows:", df_combined.shape[0])
        except Exception as e:
            print("Failed to read combined AMR TSV:", e, file=sys.stderr)

    # attempt to load counts/binary matrices from any existing AMRStatistics outputs
    counts = None; binary = None
    stats_candidates = sorted(Path('results/tools').glob('*_AMRStatistics'))
    for sc in stats_candidates:
        cf = sc / 'counts_class_by_sample.csv'
        bf = sc / 'binary_class_by_sample.csv'
        if cf.exists() and counts is None:
            try:
                counts = pd.read_csv(cf, index_col=0)
                print("Loaded counts matrix from", sc)
            except Exception:
                counts = None
        if bf.exists() and binary is None:
            try:
                binary = pd.read_csv(bf, index_col=0)
                print("Loaded binary matrix from", sc)
            except Exception:
                binary = None

    # fallback: try Excel sheets
    try:
        xls = pd.ExcelFile(str(xlsx))
        if 'counts_class_by_sample' in xls.sheet_names and (counts is None or counts.empty):
            counts = pd.read_excel(xls, sheet_name='counts_class_by_sample', index_col=0)
            print("Loaded counts_class_by_sample from Excel")
        if 'binary_class_by_sample' in xls.sheet_names and (binary is None or binary.empty):
            binary = pd.read_excel(xls, sheet_name='binary_class_by_sample', index_col=0)
            print("Loaded binary_class_by_sample from Excel")
    except Exception:
        pass

    # pivot df_combined -> counts if needed
    if (counts is None or counts.empty) and df_combined is not None:
        sample_col = None; class_col = None
        for c in df_combined.columns:
            if c.lower() == 'sample':
                sample_col = c
            if 'class' in c.lower() or 'gene' in c.lower() or 'subclass' in c.lower():
                class_col = class_col or c
        if sample_col and class_col:
            try:
                pivot = df_combined.groupby([sample_col, class_col]).size().unstack(fill_value=0)
                counts = pivot
                print("Built counts matrix from combined AMR TSV by pivoting sample x class.")
            except Exception:
                pass

    if (binary is None or binary.empty) and (counts is not None and not counts.empty):
        binary = (counts.fillna(0) > 0).astype(int)
        print("Derived binary matrix from counts matrix.")

    # load metadata (combined_raw sheet or diversity_per_sample)
    meta = pd.DataFrame()
    try:
        xls = pd.ExcelFile(str(xlsx))
        if 'combined_raw' in xls.sheet_names:
            combined_raw = pd.read_excel(xls, sheet_name='combined_raw')
            if {'sample','group'}.issubset(set(combined_raw.columns)):
                meta = combined_raw[['sample','group']].drop_duplicates().set_index('sample')
                print("Loaded metadata (combined_raw) with groups.")
    except Exception:
        pass
    if meta.empty:
        # maybe diversity_per_sample exists in previous AMRStatistics; try recent folder
        prev_stats = sorted(Path('results/tools').glob('*_AMRStatistics'))
        if prev_stats:
            latest = prev_stats[-1]
            divf = latest / 'diversity_per_sample.csv'
            if divf.exists():
                try:
                    meta = pd.read_csv(divf, index_col=0).filter(items=['group'])
                    print("Loaded metadata from diversity_per_sample.csv in", latest)
                except Exception:
                    pass

    # --- Binary co-occurrence (global)
    if binary is None or binary.empty:
        print("Binary matrix not available; skipping global binary co-occurrence.")
    else:
        print("Running global binary Spearman co-occurrence (rho >= {})".format(args.min_rho))
        edges_bin, corr_bin = compute_spearman_network(binary, min_rho=args.min_rho, max_p=0.01)
        edges_bin.to_csv(outdir / 'network_edges_binary_spearman.csv', index=False)
        corr_bin.to_csv(outdir / 'network_corr_binary_spearman.csv')
        print("Wrote binary network edges and correlation matrix.")
        G_bin = build_network_from_edges(edges_bin)
        nodes_stats = compute_node_stats(G_bin, meta_df=meta if not meta.empty else None, binary_df=binary)
        nodes_stats.to_csv(outdir / 'network_nodes_binary_stats.csv', index=False)
        if G_bin.number_of_nodes() > 0:
            export_gexf(G_bin, outdir / 'network_binary.gexf', node_stats_df=nodes_stats)
            if PLOTLY:
                interactive_network_plot(G_bin, outdir / 'network_binary_interactive.html', node_stats_df=nodes_stats,
                                         title="AMR co-occurrence network (global, binary)")
        if community_louvain is not None and G_bin.number_of_nodes() > 0:
            partition = run_louvain(G_bin)
            pd.Series(partition, name='community').to_csv(outdir / 'network_binary_communities.csv')
            print("Saved binary communities.")

    # --- Abundance co-occurrence (global) ---
    if counts is None or counts.empty:
        print("Counts matrix not available; skipping global abundance-based co-occurrence.")
    else:
        print("Running global abundance Spearman co-occurrence")
        edges_ab, corr_ab = compute_spearman_network(counts, min_rho=args.min_rho, max_p=0.01)
        edges_ab.to_csv(outdir / 'network_edges_abundance_spearman.csv', index=False)
        corr_ab.to_csv(outdir / 'network_corr_abundance_spearman.csv')
        print("Wrote abundance-based network edges and correlation matrix.")
        G_ab = build_network_from_edges(edges_ab)
        nodes_stats_ab = compute_node_stats(G_ab, meta_df=meta if not meta.empty else None, binary_df=binary if binary is not None else None)
        nodes_stats_ab.to_csv(outdir / 'network_nodes_abundance_stats.csv', index=False)
        if G_ab.number_of_nodes() > 0:
            export_gexf(G_ab, outdir / 'network_abundance.gexf', node_stats_df=nodes_stats_ab)
            if PLOTLY:
                interactive_network_plot(G_ab, outdir / 'network_abundance_interactive.html', node_stats_df=nodes_stats_ab,
                                         title="AMR co-occurrence network (global, abundance)")
        if community_louvain is not None and G_ab.number_of_nodes() > 0:
            part_ab = run_louvain(G_ab)
            pd.Series(part_ab, name='community').to_csv(outdir / 'network_abundance_communities.csv')
            print("Saved abundance communities.")

    # --- Per-group subgraphs ---
    if meta is None or meta.empty:
        print("No group metadata available; skipping per-group subgraphs.")
    else:
        groups = sorted(meta['group'].dropna().unique())
        print("Generating per-group subgraphs for groups:", groups)
        for g in groups:
            samples_in_g = meta.index[meta['group'] == g].tolist()
            if len(samples_in_g) < args.min_samples:
                print(f"Skipping group {g}: only {len(samples_in_g)} samples (min {args.min_samples})")
                continue
            if binary is None or binary.empty:
                print(f"No binary matrix; cannot build subgraph for group {g}")
                continue
            # subset binary to samples in this group
            sub_bin = binary.reindex(index=samples_in_g).dropna(axis=1, how='all')
            if sub_bin.shape[0] < args.min_samples or sub_bin.shape[1] < 2:
                print(f"Insufficient data for group {g}; skipping")
                continue
            edges_g, corr_g = compute_spearman_network(sub_bin, min_rho=args.min_rho, max_p=0.01)
            edges_g.to_csv(outdir / f'network_edges_binary_{g}.csv', index=False)
            corr_g.to_csv(outdir / f'network_corr_binary_{g}.csv')
            Gg = build_network_from_edges(edges_g)
            nodes_g = compute_node_stats(Gg, meta_df=meta.loc[samples_in_g], binary_df=sub_bin)
            nodes_g.to_csv(outdir / f'network_nodes_binary_{g}.csv', index=False)
            if Gg.number_of_nodes() > 0:
                export_gexf(Gg, outdir / f'network_{g}.gexf', node_stats_df=nodes_g)
                if PLOTLY:
                    interactive_network_plot(Gg, outdir / f'network_{g}_interactive.html', node_stats_df=nodes_g,
                                             title=f"AMR network ({g})")
            if community_louvain is not None and Gg.number_of_nodes() > 0:
                partg = run_louvain(Gg)
                pd.Series(partg, name='community').to_csv(outdir / f'network_{g}_communities.csv')
            print(f"Wrote per-group network outputs for group: {g}")

    # --- MGE mapping from combined AMR TSV (naive) ---
    amr_coords = None
    if df_combined is not None:
        try:
            amr_coords = map_genes_to_mge_from_combined_amr(df_combined)
            if not amr_coords.empty:
                amr_coords.to_csv(outdir / 'amr_gene_locations.csv', index=False)
                amr_coords.groupby('contig').size().reset_index(name='n_amr_genes').to_csv(outdir / 'amr_genes_per_contig.csv', index=False)
                print("Saved AMR gene locations and contig counts.")
        except Exception as e:
            print("MGE mapping failed:", e, file=sys.stderr)
            traceback.print_exc()
    else:
        print("No combined AMR TSV loaded; skipped MGE mapping.")

    # --- Mechanism enrichment ---
    mech_map = None
    descf = outdir / 'descriptive_by_class.csv'
    if descf.exists():
        try:
            desc = pd.read_csv(descf, index_col=0)
            for c in desc.columns:
                if c.lower() in ('mechanism','drug_class','drugclass','mechanism/drug'):
                    mech_map = desc[[c]].rename(columns={c:'mechanism'})
                    break
        except Exception:
            mech_map = None
    if mech_map is None and df_combined is not None:
        map_candidates = [col for col in df_combined.columns if col.lower() in ('mechanism','drug_class','drug','amr_class','class')]
        if map_candidates and 'gene' in df_combined.columns:
            mech_map = df_combined[['gene'] + [map_candidates[0]]].drop_duplicates().set_index('gene').rename(columns={map_candidates[0]:'mechanism'})
    if mech_map is None:
        print("Mechanism mapping not found; skipping mechanism enrichment.")
    else:
        if binary is None or binary.empty or meta is None or meta.empty:
            print("Binary matrix or metadata (groups) missing; cannot perform mechanism enrichment.")
        else:
            enrich = mechanism_enrichment(binary, meta, mech_map, group_col='group')
            if not enrich.empty:
                enrich.to_csv(outdir / 'mechanism_enrichment_by_group.csv', index=False)
                print("Saved mechanism_enrichment_by_group.csv")

    # Summary print
    print("\nCompleted network/MGE/enrichment. Outputs in:", outdir)
    for fname in sorted(outdir.iterdir()):
        print("  -", fname.name)

    # cleanup combined files unless requested
    if not args.keep_temp:
        try:
            if combined_amr.exists():
                combined_amr.unlink()
            if combined_gff.exists():
                combined_gff.unlink()
            try:
                combined_amr.parent.rmdir()
            except Exception:
                pass
        except Exception:
            pass

if __name__ == '__main__':
    main()
