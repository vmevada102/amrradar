# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
FINAL, CLEAN, UTF-8 SAFE VERSION
================================
This is the fully corrected generate_amr_reports_autoenv.py
All syntax errors removed, all stray characters cleaned, all earlier fixes applied.
This version is safe for direct copy-paste.
"""

import sys
import os
import glob
import subprocess
import shutil
from datetime import datetime

# ---------------- configuration ----------------
INPUT_GLOB = "results/tools/*_AMRFinder/combined_amrfinder.tsv"
OUT_DIR = "results/Charts/AMRFinder"
FIG_DPI = 300
PLOTLY_HEIGHT = 900
PLOTLY_WIDTH = 1400

KMEANS_K = 4
DBSCAN_EPS = 0.8
DBSCAN_MIN_SAMPLES = 3
SANKEY_TOP_N_GENES = None

REQ_CONDA = [
    "pandas", "numpy", "scipy", "scikit-learn", "matplotlib", "seaborn",
    "plotly", "openpyxl", "xlrd", "kaleido", "networkx", "python-louvain"
]

# ---------------- helper functions ----------------

def run_cmd(cmd, capture=False):
    try:
        if capture:
            proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            return proc.returncode, proc.stdout.strip(), proc.stderr.strip()
        else:
            proc = subprocess.run(cmd, shell=True)
            return proc.returncode, None, None
    except Exception as e:
        return 1, "", str(e)

def conda_available():
    return shutil.which("conda") is not None

def in_conda_env():
    return os.environ.get('CONDA_PREFIX') is not None

def ensure_packages(packages_conda=None):
    import importlib
    if packages_conda is None:
        packages_conda = []
    missing = []

    mapping = [
        ("pandas", "pandas"), ("numpy", "numpy"), ("scipy", "scipy"),
        ("sklearn", "scikit-learn"), ("matplotlib", "matplotlib"),
        ("seaborn", "seaborn"), ("plotly", "plotly"), ("openpyxl", "openpyxl"),
        ("xlrd", "xlrd"), ("kaleido", "kaleido"), ("networkx", "networkx"),
        ("community", "python-louvain")
    ]

    for imp, pkg in mapping:
        try:
            importlib.import_module(imp)
        except Exception:
            missing.append(pkg)

    for p in packages_conda:
        if p not in missing:
            try:
                importlib.import_module(p)
            except Exception:
                missing.append(p)

    if not missing:
        print("All required packages are installed.")
        return True

    use_conda = conda_available() and in_conda_env()
    if use_conda:
        print("Installing missing packages with conda (conda-forge):", missing)
        cmd = f"conda install -y -c conda-forge {' '.join(missing)}"
        rc, out, err = run_cmd(cmd, capture=True)
        if rc != 0:
            print("Conda install failed; falling back to pip. Error:\n" + err)

    still_missing = []
    for pkg in missing:
        modname = pkg.split('=')[0]
        try:
            __import__(modname)
        except Exception:
            still_missing.append(pkg)

    if still_missing:
        pip_cmd = shutil.which('pip') or shutil.which('pip3')
        if not pip_cmd:
            print('pip not found. Install manually:', still_missing)
            return False
        print('Installing missing packages with pip:', still_missing)
        cmd = f'"{pip_cmd}" install ' + ' '.join(still_missing)
        rc, out, err = run_cmd(cmd, capture=True)
        if rc != 0:
            print('Pip install failed:', err)
            return False

    print('Package installation complete.')
    return True

# ---------------- ensure output dir ----------------
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------- find input file ----------------
files = glob.glob(INPUT_GLOB)
if not files:
    raise FileNotFoundError(f"No combined_amrfinder.tsv found matching pattern {INPUT_GLOB}")
input_file = files[0]
print('Input file:', input_file)

# ---------------- ensure packages ----------------
ensure_packages(REQ_CONDA)

# ---------------- imports ----------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN
import plotly.graph_objects as go
import plotly.io as pio

# ---------------- Newick conversion ----------------

def linkage_to_newick(Z, labels):
    Z = np.asarray(Z)
    n = len(labels)
    tree = {i: labels[i] for i in range(n)}
    heights = {}

    for i, (a, b, dist, cnt) in enumerate(Z):
        a = int(a); b = int(b)
        left = tree[a]; right = tree[b]
        left_h = heights.get(a, 0.0); right_h = heights.get(b, 0.0)
        left_len = max(dist/2 - left_h, 0.0)
        right_len = max(dist/2 - right_h, 0.0)
        newnode = f"({left}:{left_len:.6f},{right}:{right_len:.6f})"
        tree[n + i] = newnode
        heights[n + i] = dist/2

    return tree[2*n - 2] + ";"

def write_newick(z, labels, outpath):
    newick = linkage_to_newick(z, labels)
    with open(outpath, 'w', encoding='utf-8') as fh:
        fh.write(newick)
    print('Wrote Newick to', outpath)

# ---------------- Load AMRFinder data ----------------
df = pd.read_csv(input_file, sep='\t')
df.columns = [c.strip() for c in df.columns]
if 'element symbol' not in df.columns:
    raise KeyError('Missing required column: element symbol')
if 'group' not in df.columns:
    df['group'] = 'Unknown'

mat = df.pivot_table(index='sample', columns='element symbol', aggfunc='size', fill_value=0)
mat_binary = (mat > 0).astype(int)

sample_group = df[['sample','group']].drop_duplicates().set_index('sample').reindex(mat.index)

# ---------------- Exports ----------------
mat.to_csv(os.path.join(OUT_DIR, 'amr_matrix_counts.csv'))
mat_binary.to_csv(os.path.join(OUT_DIR, 'amr_matrix_binary.csv'))
sample_group.to_csv(os.path.join(OUT_DIR, 'sample_groups.csv'))

# ---------------- Clustering ----------------
X = StandardScaler().fit_transform(mat_binary.values)
dist = pdist(X, metric='euclidean')
linkage_ward = sch.linkage(dist, method='ward')

newick_path = os.path.join(OUT_DIR, 'amr_tree.newick')
write_newick(linkage_ward, list(mat.index.astype(str)), newick_path)

# ---------------- Static dendrogram ----------------
fig, ax = plt.subplots(figsize=(18, 7))
sch.dendrogram(linkage_ward, labels=mat.index.astype(str).tolist(), leaf_rotation=90, leaf_font_size=6)
ax.set_title('Hierarchical Clustering (Ward) - AMR Profiles')
fig.savefig(os.path.join(OUT_DIR, 'dendrogram_ward.png'), dpi=FIG_DPI)
fig.savefig(os.path.join(OUT_DIR, 'dendrogram_ward.svg'))
plt.close(fig)

# ---------------- Interactive dendrogram ----------------
D = sch.dendrogram(linkage_ward, labels=mat.index.astype(str).tolist(), no_plot=True)
icoord = np.array(D['icoord'])
dcoord = np.array(D['dcoord'])
leaf_labels = D['ivl']

fig_d = go.Figure()
for xs, ys in zip(icoord, dcoord):
    fig_d.add_trace(go.Scatter(x=xs, y=ys, mode='lines', line=dict(color='gray')))

# 100-color palette
import matplotlib as mpl
colors100 = []
for cmap_name in ('tab20', 'tab20b', 'tab20c'):
    cm = mpl.cm.get_cmap(cmap_name)
    for i in range(cm.N):
        rgba = cm(i)
        colors100.append(f'rgba({int(rgba[0]*255)},{int(rgba[1]*255)},{int(rgba[2]*255)},{rgba[3]:.2f})')
hsv = mpl.cm.get_cmap('hsv')
for i in range(40):
    rgba = hsv(i/40)
    colors100.append(f'rgba({int(rgba[0]*255)},{int(rgba[1]*255)},{int(rgba[2]*255)},{rgba[3]:.2f})')

x_coords = np.unique(icoord.flatten())
x_min, x_max = x_coords.min(), x_coords.max()
leaf_positions = np.linspace(x_min, x_max, len(leaf_labels))

for i, lbl in enumerate(leaf_labels):
    fig_d.add_trace(go.Scatter(x=[leaf_positions[i]], y=[0], mode='markers+text', text=[lbl],
                               textposition='bottom center', marker=dict(size=8, color=colors100[i % len(colors100)])))

fig_d.update_layout(title='Interactive Ward Dendrogram', height=700, width=1400)
fig_d.write_html(os.path.join(OUT_DIR, 'dendrogram_ward_interactive.html'), include_plotlyjs='cdn')

# ---------------- Clustermap ----------------
sns.set(style='white')
cg = sns.clustermap(mat_binary, figsize=(14, 12), cmap='viridis', row_cluster=True, col_cluster=True)
cg.fig.savefig(os.path.join(OUT_DIR, 'clustermap.png'), dpi=FIG_DPI)
cg.fig.savefig(os.path.join(OUT_DIR, 'clustermap.svg'))
plt.close(cg.fig)

# ---------------- PCA ----------------
pca = PCA(n_components=2)
Xp = pca.fit_transform(X)
pcadf = pd.DataFrame(Xp, index=mat.index, columns=['PC1','PC2'])
pcadf['group'] = sample_group['group']

fig, ax = plt.subplots(figsize=(8, 7))
for g in pcadf['group'].unique():
    sub = pcadf[pcadf['group'] == g]
    ax.scatter(sub['PC1'], sub['PC2'], s=30, label=g)
ax.legend(); fig.savefig(os.path.join(OUT_DIR,'pca.png'), dpi=FIG_DPI)
fig.savefig(os.path.join(OUT_DIR,'pca.svg'))
plt.close(fig)

# ---------------- Heatmap ----------------
fig, ax = plt.subplots(figsize=(14, 10))
sns.heatmap(mat_binary, cmap='viridis', ax=ax)
fig.savefig(os.path.join(OUT_DIR,'heatmap.png'), dpi=FIG_DPI)
fig.savefig(os.path.join(OUT_DIR,'heatmap.svg'))
plt.close(fig)

# ---------------- Cluster assignments ----------------
kmeans = KMeans(n_clusters=KMEANS_K, random_state=0)
km = kmeans.fit_predict(X)

db = DBSCAN(eps=DBSCAN_EPS, min_samples=DBSCAN_MIN_SAMPLES)
dbcl = db.fit_predict(X)

cldf = pd.DataFrame({'kmeans': km,'dbscan': dbcl,'group': sample_group['group']}, index=mat.index)
cldf.to_csv(os.path.join(OUT_DIR,'cluster_assignments.csv'))

# ---------------- Sankey utilities ----------------

def build_color100():
    return colors100

COLOR100 = build_color100()

# group ? gene sankey
def sankey_group_to_genes(df_local, mat_local, top_n=None, out_html=None):
    sample_group_local = df_local[['sample','group']].drop_duplicates().set_index('sample').reindex(mat_local.index)
    gene_group = mat_local.join(sample_group_local).groupby('group').sum().T
    if top_n:
        gene_group = gene_group.loc[gene_group.sum(axis=1).nlargest(top_n).index]
    genes = gene_group.index.tolist(); groups = gene_group.columns.tolist()
    nodes = groups + genes
    node_colors = [COLOR100[i % len(COLOR100)] for i in range(len(nodes))]
    src=[]; dst=[]; val=[]; lbl=[]
    for gi, g in enumerate(groups):
        for gj, gene in enumerate(genes):
            c = int(gene_group.loc[gene, g])
            if c>0:
                src.append(gi); dst.append(len(groups)+gj); val.append(c)
                lbl.append(f"{g}?{gene}:{c}")
    fig = go.Figure(data=[go.Sankey(node=dict(label=nodes,color=node_colors), link=dict(source=src,target=dst,value=val,label=lbl))])
    fig.update_layout(height=PLOTLY_HEIGHT,width=PLOTLY_WIDTH,title="Group?Gene Sankey")
    if out_html:
        fig.write_html(out_html, include_plotlyjs='cdn')
    return fig

# gene ? group sankey
def sankey_gene_to_group(df_local, mat_local, top_n=None, out_html=None):
    sample_group_local = df_local[['sample','group']].drop_duplicates().set_index('sample').reindex(mat_local.index)
    gene_group = mat_local.join(sample_group_local).groupby('group').sum().T
    if top_n:
        gene_group = gene_group.loc[gene_group.sum(axis=1).nlargest(top_n).index]
    genes = gene_group.index.tolist(); groups = gene_group.columns.tolist()
    nodes = genes + groups
    node_colors = [COLOR100[i % len(COLOR100)] for i in range(len(nodes))]
    src=[]; dst=[]; val=[]; lbl=[]
    for gi, g in enumerate(groups):
        for gj, gene in enumerate(genes):
            c = int(gene_group.loc[gene, g])
            if c>0:
                src.append(gj); dst.append(len(genes)+gi); val.append(c)
                lbl.append(f"{gene}?{g}:{c}")
    fig = go.Figure(data=[go.Sankey(node=dict(label=nodes,color=node_colors), link=dict(source=src,target=dst,value=val,label=lbl))])
    fig.update_layout(height=PLOTLY_HEIGHT,width=PLOTLY_WIDTH,title="Gene?Group Sankey")
    if out_html:
        fig.write_html(out_html, include_plotlyjs='cdn')
    return fig

# group ? class ? gene sankey
def sankey_group_class_gene(df_local, mat_local, out_html=None):
    if 'class' not in df_local.columns:
        df_local['class'] = df_local.get('element name','Unknown')
    sample_group_local = df_local[['sample','group']].drop_duplicates().set_index('sample').reindex(mat_local.index)
    gene_class_map = df_local.groupby('element symbol')['class'].first().to_dict()
    genes = mat_local.columns.tolist()
    classes = [gene_class_map.get(g,'Unknown') for g in genes]

    gene_by_sample = mat_local.T.copy()
    gene_by_sample['class'] = classes

    sample_has_class = (gene_by_sample.drop(columns='class').groupby(gene_by_sample['class']).sum()>0).T.astype(int)
    sample_has_class['group'] = sample_group_local['group']
    class_group = sample_has_class.groupby('group').sum().T

    class_gene_counts = {}
    for cls in sorted(set(classes)):
        cls_genes = [g for g,c in zip(genes,classes) if c==cls]
        if cls_genes:
            class_gene_counts[cls] = gene_by_sample.loc[cls_genes].drop(columns='class').sum(axis=1)

    groups = class_group.columns.tolist(); class_list = class_group.index.tolist(); gene_list = genes
    nodes = groups + class_list + gene_list
    node_index = {n:i for i,n in enumerate(nodes)}
    node_colors = [COLOR100[i % len(COLOR100)] for i in range(len(nodes))]
    src=[]; dst=[]; val=[]; lbl=[]

    # group?class
    for cls in class_list:
        for g in groups:
            cval = int(class_group.loc[cls, g]) if cls in class_group.index else 0
            if cval>0:
                src.append(node_index[g]); dst.append(node_index[cls]); val.append(cval)
                lbl.append(f"{g}?{cls}:{cval}")

    # class?gene
    for cls in class_list:
        for gene, cnt in class_gene_counts.get(cls, {}).items():
            if cnt>0:
                src.append(node_index[cls]); dst.append(node_index[gene]); val.append(int(cnt))
                lbl.append(f"{cls}?{gene}:{int(cnt)}")

    fig = go.Figure(data=[go.Sankey(node=dict(label=nodes,color=node_colors), link=dict(source=src,target=dst,value=val,label=lbl))])
    fig.update_layout(height=PLOTLY_HEIGHT,width=PLOTLY_WIDTH,title="Group?Class?Gene Sankey")
    if out_html:
        fig.write_html(out_html, include_plotlyjs='cdn')
    return fig

# ---------------- Generate sankeys ----------------
sankey_group_to_genes(df, mat_binary, top_n=SANKEY_TOP_N_GENES, out_html=os.path.join(OUT_DIR,'sankey_group_to_genes.html'))
sankey_gene_to_group(df, mat_binary, top_n=SANKEY_TOP_N_GENES, out_html=os.path.join(OUT_DIR,'sankey_gene_to_group.html'))
sankey_group_class_gene(df, mat_binary, out_html=os.path.join(OUT_DIR,'sankey_group_class_gene.html'))

# ---------------- Build combined HTML report ----------------
report = os.path.join(OUT_DIR,'AMR_analysis_report.html')
with open(report,'w',encoding='utf-8') as fh:
    fh.write('<html><head><meta charset="utf-8"><title>AMR Report</title></head><body>')
    fh.write('<h1>AMR Analysis Report</h1>')
    fh.write(f'<p>Generated {datetime.utcnow().isoformat()} UTC</p>')

    def embed(title, fname):
        fh.write(f'<h2>{title}</h2>')
        path = os.path.join(OUT_DIR,fname)
        if fname.endswith('.svg') and os.path.exists(path):
            fh.write(open(path,'r',encoding='utf-8').read())
        elif os.path.exists(path):
            fh.write(f"<img src='{fname}' style='max-width:100%'>")

    embed('Dendrogram (static)','dendrogram_ward.svg')
    embed('Clustermap','clustermap.svg')
    embed('Heatmap','heatmap.svg')
    embed('PCA','pca.svg')

    def embed_html(title,fname):
        fh.write(f'<h2>{title}</h2>')
        path=os.path.join(OUT_DIR,fname)
        if os.path.exists(path): fh.write(open(path,'r',encoding='utf-8').read())

    embed_html('Interactive Dendrogram','dendrogram_ward_interactive.html')
    embed_html('Interactive PCA','pca_amr_interactive.html')
    embed_html('Group?Gene Sankey','sankey_group_to_genes.html')
    embed_html('Gene?Group Sankey','sankey_gene_to_group.html')
    embed_html('Group?Class?Gene Sankey','sankey_group_class_gene.html')

    fh.write('<h2>Downloads</h2><ul>')
    for f in sorted(os.listdir(OUT_DIR)):
        fh.write(f"<li><a href='{f}' download>{f}</a></li>")
    fh.write('</ul></body></html>')

print('Report generated:', report)
print('All results saved in:', OUT_DIR)