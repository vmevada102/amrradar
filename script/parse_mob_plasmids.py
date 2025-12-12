#!/usr/bin/env python3
# parse_mob_plasmids.py <mob_recon_dir> -> prints plasmid contig names (one per line)
import sys, os
d = sys.argv[1] if len(sys.argv)>1 else "."
out=set()
# mob_recon produces plasmids.fna and mob_recon_report.txt depending on version
pfn = os.path.join(d, "plasmids.fna")
if os.path.isfile(pfn):
    with open(pfn) as fh:
        for line in fh:
            if line.startswith(">"):
                out.add(line[1:].split()[0])
rpt = os.path.join(d, "mob_recon_report.txt")
if os.path.isfile(rpt):
    with open(rpt) as fh:
        for line in fh:
            parts=line.rstrip("\n").split("\t")
            # try to capture contig ids in report lines
            for p in parts:
                if p.startswith("NODE_") or p.startswith("contig_") or p.count(":")>=1:
                    out.add(p.strip())
for x in sorted(out):
    print(x)
