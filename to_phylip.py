#!/usr/bin/env python3
"""
to_phylip.py
Convert Mash pairwise distances (two path columns + distance) into PHYLIP lower-triangular.
Input: TSV with at least 3 columns: fileA, fileB, distance
Also requires a list of taxa (one per line) to define order and ensure completeness.
"""
import sys, math
from pathlib import Path

if len(sys.argv) < 4:
    print("Usage: to_phylip.py <pairs.tsv> <taxa.list> <out.phylip>", file=sys.stderr)
    sys.exit(1)

pairs = Path(sys.argv[1]).read_text().strip().splitlines()
taxa = [Path(t).stem for t in Path(sys.argv[2]).read_text().strip().splitlines()]
# Build map from original paths to taxa names by stem
def stem(p): return Path(p).stem

idx = {t:i for i,t in enumerate(taxa)}
n = len(taxa)
# Initialize matrix with zeros on diagonal, NaN elsewhere
M = [[0.0 if i==j else math.nan for j in range(n)] for i in range(n)]

for line in pairs:
    if not line.strip(): continue
    parts = line.strip().split()
    if len(parts) < 3: continue
    a,b,d = parts[0], parts[1], parts[2]
    ta, tb = stem(a), stem(b)
    if ta not in idx or tb not in idx: continue
    i, j = idx[ta], idx[tb]
    try:
        dv = float(d)
    except Exception:
        continue
    M[i][j] = dv
    M[j][i] = dv

# Fill any remaining NaNs with a small positive distance
for i in range(n):
    for j in range(n):
        if math.isnan(M[i][j]):
            M[i][j] = 0.001 if i!=j else 0.0

# Write lower-triangular PHYLIP with names
with open(sys.argv[3], "w") as out:
    out.write(f"{n}\n")
    for i in range(n):
        name = taxa[i]
        row = [f"{M[i][k]:.6f}" for k in range(i+1)]  # include diagonal
        out.write(f"{name}\t" + "\t".join(row) + "\n")
