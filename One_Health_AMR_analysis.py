# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
One_Health_AMR_analysis.py  (updated: group PCA + HTML report + interactive PCA)

Usage:
    python One_Health_AMR_analysis.py /path/to/OneHealthAMR_AMRFinder_summary.xlsx [output_dir]

Outputs:
 - CSVs & PNG/SVG (as before)
 - Pairwise PCA plots for each pair of sample groups (PCA_groupA_vs_groupB.png/.svg)
 - Interactive PCA HTML (PCA_interactive.html) if plotly is installed
 - report.html combining key images/tables and linking interactive PCA
"""
import os, sys, traceback
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

try:
    import pandas as pd
    import numpy as np
    from scipy import stats
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    import matplotlib.pyplot as plt
    import seaborn as sns
except Exception as e:
    print("ERROR: Missing required packages. Install pandas,numpy,scipy,scikit-learn,matplotlib,seaborn", file=sys.stderr)
    raise

# optional
try:
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except Exception:
    PLOTLY_AVAILABLE = False

# optional PERMANOVA
try:
    from skbio.stats.distance import permanova
    from skbio.diversity import beta_diversity
    SKBIO_AVAILABLE = True
except Exception:
    SKBIO_AVAILABLE = False

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def savefig(fig, outpath_stem: Path, dpi=300):
    png = outpath_stem.with_suffix('.png')
    svg = outpath_stem.with_suffix('.svg')
    fig.savefig(png, dpi=dpi, bbox_inches='tight')
    fig.savefig(svg, dpi=dpi, bbox_inches='tight')
    print(f"Saved: {png} and {svg}")

def shannon(counts):
    freqs = np.array(counts, dtype=float)
    if freqs.sum() == 0:
        return 0.0
    p = freqs / freqs.sum()
    p = p[p > 0]
    return -np.sum(p * np.log(p))

def simpson(counts):
    freqs = np.array(counts, dtype=float)
    if freqs.sum() == 0:
        return 0.0
    p = freqs / freqs.sum()
    return 1.0 - np.sum(p ** 2)

def write_error(outdir: Path, exc: Exception):
    try:
        ensure_dir(outdir)
        ef = outdir / "analysis_error.txt"
        with open(ef, "w") as f:
            f.write("Exception during One_Health_AMR_analysis.py run\n")
            f.write("".join(traceback.format_exception(type(exc), exc, exc.__traceback__)))
        print(f"Wrote exception trace to {ef}", file=sys.stderr)
    except Exception as e:
        print("Failed to write exception trace:", e, file=sys.stderr)

def make_interactive_pca(df_counts, meta, outpath):
    """Create an interactive PCA html file using plotly (if available)."""
    try:
        X = df_counts.fillna(0).astype(float)
        X_log = np.log1p(X)
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_log)
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(X_scaled)
        explained = pca.explained_variance_ratio_ * 100
        pca_df = pd.DataFrame(pcs, index=X.index, columns=['PC1','PC2'])
        if not meta.empty:
            pca_df = pca_df.join(meta)
        if not PLOTLY_AVAILABLE:
            print("plotly not available - skipping interactive PCA")
            return None
        fig = px.scatter(pca_df.reset_index(), x='PC1', y='PC2', color=meta.columns[0] if not meta.empty else None,
                         hover_data=['index'] + (meta.columns.tolist() if not meta.empty else []),
                         title=f'Interactive PCA (PC1 {explained[0]:.1f}% | PC2 {explained[1]:.1f}%)')
        html_path = outpath / 'PCA_interactive.html'
        fig.write_html(str(html_path), include_plotlyjs='cdn')
        print(f"Saved interactive PCA to {html_path}")
        return html_path
    except Exception as e:
        print("Interactive PCA creation failed:", e)
        return None

def build_html_report(outdir: Path, summary_txt: str, images, interactive_link=None):
    """
    Build a simple HTML report that embeds PNG images and includes the summary text.
    images: list of tuples (title, relative_path)
    interactive_link: relative path to interactive html (optional)
    """
    try:
        report = []
        report.append("<html><head><meta charset='utf-8'><title>AMR Analysis Report</title>")
        report.append("<style>body{font-family:Arial,Helvetica,sans-serif;margin:20px;} h1{color:#1f4f97} .img{max-width:100%;height:auto;border:1px solid #ddd;padding:4px;margin-bottom:10px}</style>")
        report.append("</head><body>")
        report.append("<h1>One Health AMR Analysis Report</h1>")
        report.append(f"<p>Generated: {outdir.name}</p>")
        if interactive_link:
            report.append(f"<h2>Interactive PCA</h2><p><a href='{interactive_link}'>Open interactive PCA (new tab)</a></p>")
        report.append("<h2>Summary</h2>")
        report.append("<pre style='background:#f8f8f8;padding:10px;border:1px solid #eee;'>")
        report.append(summary_txt)
        report.append("</pre>")
        report.append("<h2>Figures</h2>")
        for title, img in images:
            report.append(f"<h3>{title}</h3>")
            report.append(f"<img class='img' src='{img}' alt='{title}' />")
        report.append("<hr><p>Report generated by One_Health_AMR_analysis.py</p>")
        report.append("</body></html>")
        report_path = outdir / "report.html"
        with open(report_path, "w") as f:
            f.write("\n".join(report))
        print(f"Wrote HTML report to {report_path}")
        return report_path
    except Exception as e:
        print("Failed to write HTML report:", e)
        return None

def main():
    try:
        if len(sys.argv) < 2:
            print("Usage: python One_Health_AMR_analysis.py /path/to/OneHealthAMR_AMRFinder_summary.xlsx [output_dir]", file=sys.stderr)
            sys.exit(1)

        input_xlsx = Path(sys.argv[1]).expanduser().resolve()
        if len(sys.argv) >= 3:
            outdir = Path(sys.argv[2]).expanduser().resolve()
        else:
            outdir = Path('./results/OneHealth_AMR_analysis').resolve()

        if not input_xlsx.exists():
            print(f"ERROR: Input file not found: {input_xlsx}", file=sys.stderr)
            sys.exit(2)

        ensure_dir(outdir)
        print(f"Input: {input_xlsx}")
        print(f"Output directory: {outdir}")

        from matplotlib.colors import LinearSegmentedColormap
        cmap_blue_white = LinearSegmentedColormap.from_list('blue_white', ['white', '#1f4f97'])

        # read excel
        xls = pd.ExcelFile(str(input_xlsx))
        print("Sheets found:", xls.sheet_names)

        def safe_read(sheet, index_col=None):
            try:
                return pd.read_excel(xls, sheet_name=sheet, index_col=index_col)
            except Exception:
                return pd.DataFrame()

        counts_class_by_sample = safe_read('counts_class_by_sample', index_col=0)
        binary_class_by_sample = safe_read('binary_class_by_sample', index_col=0)
        counts_class_by_group = safe_read('counts_class_by_group', index_col=0)
        binary_class_by_group = safe_read('binary_class_by_group', index_col=0)
        combined_raw = safe_read('combined_raw', index_col=None)

        # fallback sheet names
        if counts_class_by_sample.empty:
            for alt in ['counts_type_by_sample','counts_subclass_by_sample','counts_type_by_sample']:
                if alt in xls.sheet_names:
                    counts_class_by_sample = safe_read(alt, index_col=0)
                    break

        # metadata
        if not combined_raw.empty and {'sample','group'}.issubset(set(combined_raw.columns)):
            meta = combined_raw[['sample','group']].drop_duplicates().set_index('sample')
        else:
            meta = pd.DataFrame(index=counts_class_by_sample.index)
            if 'group' in counts_class_by_sample.columns:
                meta = counts_class_by_sample[['group']].copy()
                counts_class_by_sample = counts_class_by_sample.drop(columns=['group'])

        if counts_class_by_sample.empty:
            print("ERROR: counts_class_by_sample missing; aborting", file=sys.stderr)
            sys.exit(3)

        class_cols = [c for c in counts_class_by_sample.columns if c not in ('group',)]
        counts = counts_class_by_sample.copy()

        if not meta.empty:
            common = [s for s in counts.index if s in meta.index]
            counts = counts.loc[common]
            meta = meta.loc[counts.index]
        else:
            meta = pd.DataFrame(index=counts.index)

        if binary_class_by_sample.empty:
            binary = (counts[class_cols] > 0).astype(int)
        else:
            binary = binary_class_by_sample.copy()
            if set(counts.index).issuperset(set(binary.index)):
                binary = binary.reindex(counts.index).fillna(0).astype(int)
            else:
                common = [s for s in counts.index if s in binary.index]
                binary = binary.loc[common]
                counts = counts.loc[common]
                meta = meta.loc[common]

        counts['total_AMR_genes'] = counts[class_cols].sum(axis=1)

        # descriptive
        desc = counts[class_cols].describe().T
        desc.to_csv(outdir / 'descriptive_by_class.csv')
        print("Saved descriptive_by_class.csv")

        if 'group' in meta.columns:
            counts_with_group = counts.join(meta)
            group_summary = counts_with_group.groupby('group')[class_cols].agg(['mean','median','std','count'])
            group_summary.to_csv(outdir / 'descriptive_by_class_by_group.csv')
            print("Saved descriptive_by_class_by_group.csv")
        else:
            group_summary = None

        # diversity
        diversity = pd.DataFrame(index=counts.index)
        diversity['richness'] = (counts[class_cols] > 0).sum(axis=1)
        diversity['shannon'] = counts[class_cols].apply(shannon, axis=1)
        diversity['simpson'] = counts[class_cols].apply(simpson, axis=1)
        diversity['total_AMR_genes'] = counts['total_AMR_genes']
        diversity = diversity.join(meta)
        diversity.to_csv(outdir / 'diversity_per_sample.csv')
        print("Saved diversity_per_sample.csv")

        sns.set(context='notebook', style='whitegrid', font_scale=1.2)
        for metric in ['richness','shannon','simpson','total_AMR_genes']:
            fig, ax = plt.subplots(figsize=(7,5))
            if 'group' in diversity.columns:
                sns.boxplot(x='group', y=metric, data=diversity.reset_index(), ax=ax)
                sns.stripplot(x='group', y=metric, data=diversity.reset_index(), color='k', size=3, alpha=0.5, ax=ax)
                ax.set_title(f'{metric} by group')
                ax.set_xlabel('Group')
                ax.set_ylabel(metric)
                plt.xticks(rotation=45)
            else:
                sns.histplot(diversity[metric], kde=True, ax=ax)
                ax.set_title(metric)
            savefig(fig, outdir / f'diversity_{metric}')
            plt.close(fig)

        if 'group' in diversity.columns:
            tests = []
            groups = diversity['group'].unique()
            for metric in ['richness','shannon','simpson','total_AMR_genes']:
                samples_by_group = [diversity.loc[diversity['group']==g, metric].values for g in groups]
                try:
                    stat,p = stats.kruskal(*samples_by_group)
                except Exception:
                    stat,p = np.nan, np.nan
                tests.append({'metric':metric,'test':'Kruskal-Wallis','stat':stat,'pvalue':p})
            pd.DataFrame(tests).to_csv(outdir / 'diversity_group_tests.csv', index=False)
            print("Saved diversity_group_tests.csv")

        # PCA full
        X = counts[class_cols].fillna(0).astype(float)
        X_log = np.log1p(X)
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_log)
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(X_scaled)
        explained = pca.explained_variance_ratio_ * 100
        pca_df = pd.DataFrame(pcs, index=X.index, columns=['PC1','PC2']).join(meta)
        pca_df.to_csv(outdir / 'PCA_coordinates.csv')

        fig, ax = plt.subplots(figsize=(8,6))
        if 'group' in pca_df.columns:
            sns.scatterplot(x='PC1', y='PC2', hue='group', data=pca_df, s=70, ax=ax)
            ax.legend(bbox_to_anchor=(1.05,1), loc='upper left')
        else:
            sns.scatterplot(x='PC1', y='PC2', data=pca_df, s=70, ax=ax)
        ax.set_xlabel(f'PC1 ({explained[0]:.1f}% var)')
        ax.set_ylabel(f'PC2 ({explained[1]:.1f}% var)')
        ax.set_title('PCA of AMR class counts (log1p, scaled)')
        savefig(fig, outdir / 'PCA_AMR_classes')
        plt.close(fig)
        print("Saved PCA plot and coordinates")

        # PERMANOVA
        permanova_res = None
        if SKBIO_AVAILABLE and 'group' in meta.columns:
            try:
                dist = beta_diversity('braycurtis', X.values, ids=X.index)
                perma_res = permanova(dist, meta['group'].loc[X.index], permutations=999)
                with open(outdir / 'permanova.txt','w') as f:
                    f.write(str(perma_res))
                permanova_res = str(perma_res)
                print("Saved PERMANOVA result to permanova.txt")
            except Exception as e:
                print("PERMANOVA failed:", e)

        # clustermap
        heat_data = X_log.copy()
        row_colors = None
        if 'group' in meta.columns:
            uniq = list(meta['group'].unique())
            palette = dict(zip(uniq, sns.color_palette('tab10', n_colors=len(uniq))))
            row_colors = meta['group'].map(palette)
        try:
            cg = sns.clustermap(heat_data, metric='euclidean', method='average', cmap=cmap_blue_white, figsize=(12,10), row_colors=row_colors)
            cg.fig.suptitle('Clustermap of AMR class counts (log1p)')
            cg.fig.subplots_adjust(top=0.93)
            cg.fig.savefig(outdir / 'clustermap_AMR_classes.png', dpi=300, bbox_inches='tight')
            cg.fig.savefig(outdir / 'clustermap_AMR_classes.svg', dpi=300, bbox_inches='tight')
            plt.close(cg.fig)
            print("Saved clustermap")
        except Exception as e:
            print("Clustermap failed:", e)

        # group heatmap
        if not counts_class_by_group.empty:
            try:
                group_norm = counts_class_by_group.div(counts_class_by_group.sum(axis=1)+1e-9, axis=0)
                fig, ax = plt.subplots(figsize=(8, max(4, len(group_norm)/2)))
                sns.heatmap(group_norm.T, cmap=cmap_blue_white, cbar_kws={'label':'Proportion of genes'}, ax=ax)
                ax.set_title('AMR classes by group (proportions)')
                plt.yticks(rotation=0)
                savefig(fig, outdir / 'heatmap_group_proportions')
                plt.close(fig)
                print("Saved heatmap_group_proportions")
            except Exception as e:
                print("Failed to create group heatmap:", e)

        # correlation
        try:
            corr = binary[class_cols].corr(method='spearman')
            corr.to_csv(outdir / 'spearman_corr_binary_classes.csv')
            strong = []
            for i in corr.index:
                for j in corr.columns:
                    if i < j:
                        val = corr.loc[i,j]
                        if pd.notna(val) and abs(val) >= 0.6:
                            strong.append((i, j, float(val)))
            pd.DataFrame(strong, columns=['class1','class2','rho']).to_csv(outdir / 'strong_cooccurrences.csv', index=False)
            fig, ax = plt.subplots(figsize=(10,8))
            sns.heatmap(corr, cmap='vlag', center=0, ax=ax)
            ax.set_title('Spearman correlation (binary presence/absence)')
            savefig(fig, outdir / 'corr_heatmap_binary')
            plt.close(fig)
            print("Saved correlation heatmap and strong co-occurrences")
        except Exception as e:
            print("Correlation analysis failed:", e)

        # chi-square
        if 'group' in meta.columns:
            chi_results = []
            for cls in binary.columns:
                try:
                    contingency = pd.crosstab(meta['group'].loc[binary.index], binary[cls].loc[binary.index])
                    if contingency.shape == (2,2) or contingency.values.sum() > 0:
                        chi2, p, dof, ex = stats.chi2_contingency(contingency)
                    else:
                        chi2, p = np.nan, np.nan
                except Exception:
                    chi2, p = np.nan, np.nan
                chi_results.append({'class':cls, 'chi2':chi2, 'pvalue':p})
            pd.DataFrame(chi_results).to_csv(outdir / 'chi2_presence_by_group.csv', index=False)
            print("Saved chi2_presence_by_group.csv")

        # --- NEW: Pairwise PCA (group vs group) ---
        pairwise_images = []
        min_group_samples = 3
        if 'group' in meta.columns:
            groups = sorted(meta['group'].unique())
            from itertools import combinations
            for a,b in combinations(groups, 2):
                idx = meta.index[meta['group'].isin([a,b])]
                if len(idx) >= min_group_samples:
                    X_sub = counts.loc[idx, class_cols].fillna(0).astype(float)
                    X_log_sub = np.log1p(X_sub)
                    try:
                        scaler_sub = StandardScaler()
                        Xs = scaler_sub.fit_transform(X_log_sub)
                        pca_sub = PCA(n_components=2)
                        pcs_sub = pca_sub.fit_transform(Xs)
                        explained_sub = pca_sub.explained_variance_ratio_*100
                        psub = pd.DataFrame(pcs_sub, index=X_sub.index, columns=['PC1','PC2']).join(meta)
                        fname = f"PCA_{a}_vs_{b}"
                        fig, ax = plt.subplots(figsize=(7,6))
                        sns.scatterplot(x='PC1', y='PC2', hue='group', data=psub, s=80, ax=ax)
                        ax.set_xlabel(f'PC1 ({explained_sub[0]:.1f}% var)')
                        ax.set_ylabel(f'PC2 ({explained_sub[1]:.1f}% var)')
                        ax.set_title(f'PCA: {a} vs {b}')
                        savefig(fig, outdir / fname)
                        plt.close(fig)
                        pairwise_images.append((f"PCA {a} vs {b}", f"{fname}.png"))
                    except Exception as e:
                        print(f"Pairwise PCA {a} vs {b} failed:", e)
                else:
                    print(f"Skipping pairwise PCA for {a} vs {b} - only {len(idx)} samples (min {min_group_samples})")

        # interactive PCA (plotly) if possible
        interactive_path = None
        if PLOTLY_AVAILABLE:
            try:
                interactive_path = make_interactive_pca(counts[class_cols], meta, outdir)
            except Exception as e:
                print("Interactive PCA generation failed:", e)

        # produce run_summary.txt (same as before)
        try:
            summary_lines = []
            summary_lines.append("=== One_Health_AMR_analysis summary ===")
            # sample counts by group
            try:
                d = pd.read_csv(outdir / 'diversity_per_sample.csv', index_col=0)
                if 'group' in d.columns:
                    summary_lines.append("\nSample counts by group:")
                    vc = d['group'].value_counts()
                    for g,v in vc.items():
                        summary_lines.append(f"  {g}: {v}")
                else:
                    summary_lines.append("No group metadata found in diversity_per_sample.csv")
            except Exception:
                summary_lines.append("Could not read diversity_per_sample.csv")

            # top classes overall
            try:
                desc = pd.read_csv(outdir / 'descriptive_by_class.csv', index_col=0)
                if 'mean' in desc.columns:
                    summary_lines.append("\nTop 10 AMR classes by mean abundance:")
                    for cls,val in desc['mean'].sort_values(ascending=False).head(10).items():
                        summary_lines.append(f"  {cls}: {val:.2f}")
                else:
                    summary_lines.append("descriptive_by_class.csv missing mean column")
            except Exception:
                summary_lines.append("Could not read descriptive_by_class.csv")

            # strong co-occurrences
            try:
                co = pd.read_csv(outdir / 'strong_cooccurrences.csv')
                if not co.empty:
                    summary_lines.append("\nTop strong co-occurrences (rho):")
                    for _, r in co.sort_values('rho', ascending=False).head(10).iterrows():
                        summary_lines.append(f"  {r['class1']} - {r['class2']}: {r['rho']:.2f}")
                else:
                    summary_lines.append("\nNo strong co-occurrences found (threshold rho>=0.6)")
            except Exception:
                summary_lines.append("Could not read strong_cooccurrences.csv")

            # PERMANOVA
            if permanova_res is not None:
                summary_lines.append("\nPERMANOVA result:")
                summary_lines.append(permanova_res)
            else:
                try:
                    with open(outdir / 'permanova.txt') as f:
                        txt = ''.join(f.readlines()[:20])
                        summary_lines.append("\nPERMANOVA (file):")
                        summary_lines.append(txt)
                except Exception:
                    summary_lines.append("\nPERMANOVA not run or result not available")

            # write
            summary_path = outdir / 'run_summary.txt'
            with open(summary_path, 'w') as sf:
                sf.write("\n".join(summary_lines) + "\n")
            print("Saved run_summary.txt")
        except Exception as e:
            print("Failed to produce run summary:", e)

        # HTML report: embed some key PNGs and include pairwise PCA images and interactive link if present
        try:
            # collect images: prefer PCA, clustermap, heatmaps, diversity
            images = [
                ("PCA (full)", "PCA_AMR_classes.png"),
                ("Clustermap", "clustermap_AMR_classes.png"),
                ("Heatmap (group proportions)", "heatmap_group_proportions.png"),
                ("Correlation heatmap", "corr_heatmap_binary.png"),
                ("Diversity (richness)", "diversity_richness.png"),
            ]
            # append pairwise images
            images.extend(pairwise_images)
            # read summary text
            try:
                with open(outdir / 'run_summary.txt') as f:
                    summary_txt = f.read()
            except Exception:
                summary_txt = "No summary available."
            interlink = None
            if interactive_path is not None:
                interlink = interactive_path.name
            report_path = build_html_report(outdir, summary_txt, images, interactive_link=interlink)
        except Exception as e:
            print("Failed to build HTML report:", e)

        # DONE & generated files
        try:
            import glob
            files = sorted(glob.glob(str(outdir / '*')))
            with open(outdir / 'DONE.txt', 'w') as f:
                f.write('One_Health_AMR_analysis completed\n')
                f.write(f'Generated file count: {len(files)}\n')
            with open(outdir / 'generated_files.csv', 'w') as gf:
                gf.write('filename\n')
                for p in files:
                    gf.write(os.path.basename(p) + '\n')
            print('Wrote DONE.txt and generated_files.csv')
        except Exception as e:
            print('Failed to write completion markers:', e)

        print("Analysis complete. Files written to:", outdir)

    except Exception as exc:
        # ensure any unexpected exception is written to outdir
        try:
            write_error(outdir, exc)
        except Exception:
            print("Fatal error and failed to write error file. Traceback below:", file=sys.stderr)
            traceback.print_exc()
        raise

if __name__ == "__main__":
    main()
