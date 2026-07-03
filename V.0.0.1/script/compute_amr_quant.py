#!/usr/bin/env python3
# compute_amr_quant.py
# Input: combined counts TSV with header sample,gene,length,mapped_reads
# Produces per-sample RPK,RPKM,TPM and writes combined TSV
import sys, pandas as pd
if len(sys.argv) < 3:
    print("Usage: compute_amr_quant.py <combined_counts.tsv> <out_tpm.tsv>"); sys.exit(1)
inf, outf = sys.argv[1], sys.argv[2]
df = pd.read_csv(inf, sep='\t')
df['length_kb'] = df['length'].replace(0,1)/1000.0
df['RPK'] = df['mapped_reads'] / df['length_kb']
out_rows = []
for s, g in df.groupby('sample'):
    g = g.copy()
    sum_rpk = g['RPK'].sum()
    if sum_rpk <= 0:
        g['TPM'] = 0.0
    else:
        g['TPM'] = (g['RPK'] / sum_rpk) * 1e6
    total_mapped = g['mapped_reads'].sum()
    if total_mapped <= 0:
        g['RPKM'] = 0.0
    else:
        g['RPKM'] = g['RPK'] / (total_mapped/1e6)
    out_rows.append(g)
res = pd.concat(out_rows, ignore_index=True)
res.to_csv(outf, sep='\t', index=False)
print("Wrote", outf)
