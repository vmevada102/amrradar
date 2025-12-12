#!/usr/bin/env python3
"""
extract_gene_and_class_amr.py

Automatically:
 - Reads summary from: results/OneHealthAMR_AMRFinder_summary.xlsx
 - Finds assemblies in: results/tools/*_assembly/<sample>/<sample>.assembly.fasta
 - Chooses output directory under results/tools as:
     - NEW RUN (default):
         Look at all directories in results/tools named "N_*"
         (e.g. 1_QC, 2_trim, 3_assembly) and create:
             (max(N) + 1)_amr_extract
         e.g. 4_amr_extract
     - RESUME (--resume):
         If any N_amr_extract exists, use the one with highest N.
         Otherwise behave like a new run (create next N_amr_extract).

Inside the chosen output dir (N_amr_extract) it creates:
 - per_gene/<CLEAN_GENE_NAME>.fa      (all sample hits for each gene)
 - per_class/<CLEAN_CLASS_NAME>.fa    (all sample hits for each antibiotic class)
 - amr_gene_extraction_mapping.csv    (optional mapping CSV)
 - command.txt                        (the command used to run this script)

Features:
 - Auto-detect common column names in summary (xlsx/csv/tsv)
   * sample
   * contig / contig id
   * start / stop
   * gene from element symbol / element name / gene symbol ...
   * class from class / type / etc.
 - Options:
    --dedup      : exact nucleotide dedup per gene/class file
    --flank N    : add N bp flanks around hits
    --translate  : also write protein FASTAs
    --mapping    : write mapping CSV
    --cdhit X    : optional cd-hit-est clustering at identity X (e.g. 0.99)
    --resume     : resume in latest N_amr_extract, skipping already-done genes/classes

Dependencies:
  pip install biopython pandas openpyxl
  (Optional) cd-hit-est on PATH for --cdhit

Usage (basic):
  python extract_gene_and_class_amr.py

Show help:
  python extract_gene_and_class_amr.py -h
  python extract_gene_and_class_amr.py -help
  python extract_gene_and_class_amr.py --help
"""

import os
import sys
import argparse
import csv
import hashlib
import subprocess
import re
import glob
import shlex
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# Candidate column names (common variants)
SAMPLE_CANDS = ["sample", "Sample", "SAMPLE"]

# Your file: contig id
CONTIG_CANDS = [
    "contig", "Contig", "contig id", "contig_id",
    "sequence-id", "seqid", "seq-id", "sseqid", "subject"
]

START_CANDS = ["start", "Start", "sstart", "subject_start", "hit_start"]

# Your file: stop
END_CANDS = ["end", "End", "stop", "Stop", "send", "subject_end", "hit_end"]

STRAND_CANDS = ["strand", "orientation", "frame"]

# Your file: element symbol / element name
GENE_CANDS = [
    "element symbol", "element name",
    "gene", "gene symbol", "gene_symbol", "Gene",
    "name", "closest reference name", "symbol", "gene_name"
]

CLASS_CANDS = ["class", "Class", "type", "category", "drug_class",
               "antibiotic_class", "name"]


def choose_column(cols, candidates, default=None):
    # exact match
    for c in candidates:
        if c in cols:
            return c
    # case-insensitive
    lower2orig = {c.lower(): c for c in cols}
    for c in candidates:
        if c.lower() in lower2orig:
            return lower2orig[c.lower()]
    return default


def clean_name(s):
    return "".join([c if c.isalnum() or c in "_-." else "_" for c in str(s)]).strip("_")


def load_assembly_index(fa_path):
    return SeqIO.to_dict(SeqIO.parse(fa_path, "fasta"))


def sha256_seq(s):
    return hashlib.sha256(s.encode()).hexdigest()


def run_cdhit_est(in_fasta, out_fasta, identity):
    """Run cd-hit-est to cluster and write representative sequences to out_fasta.
       Returns True if succeeded."""
    cmd = [
        "cd-hit-est",
        "-i", in_fasta,
        "-o", out_fasta,
        "-c", str(identity),
        "-n", "10",
        "-M", "0",
        "-T", "2",
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except Exception as e:
        print("cd-hit-est failed or not found:", e)
        return False


def write_fasta_records(path, records):
    # records: list of (header, seqstr)
    if not records:
        return False
    with open(path, "w") as fh:
        for h, s in records:
            fh.write(f">{h}\n")
            fh.write(s + "\n")
    return True


def find_numbered_tool_dirs(tools_dir):
    """
    Return list of (N, name, full_path) for directories named 'N_*'
    e.g. 1_QC, 2_trim, 3_assembly, 4_amr_extract
    """
    out = []
    if not os.path.isdir(tools_dir):
        return out
    for name in os.listdir(tools_dir):
        full = os.path.join(tools_dir, name)
        if not os.path.isdir(full):
            continue
        m = re.match(r"(\d+)_.*$", name)
        if m:
            n = int(m.group(1))
            out.append((n, name, full))
    return out


def determine_outdir(tools_dir, resume):
    """
    If resume:
      - use highest N_amr_extract directory if any exist
      - otherwise create (max N + 1)_amr_extract based on all N_* dirs
    Else:
      - create new directory with N = max existing N among all N_* dirs + 1
        name = f"{N}_amr_extract"
    """
    os.makedirs(tools_dir, exist_ok=True)
    existing = find_numbered_tool_dirs(tools_dir)

    if resume:
        # choose highest N among N_* where name endswith '_amr_extract'
        amr_dirs = [(n, name, path) for (n, name, path) in existing if name.endswith("_amr_extract")]
        if amr_dirs:
            max_n = max(n for n, _, _ in amr_dirs)
            outdir = [p for n, _, p in amr_dirs if n == max_n][0]
            print(f"Resuming in existing directory: {outdir}")
            return outdir
        else:
            # no existing amr_extract dir; behave as new run
            max_n_all = max([n for n, _, _ in existing], default=0)
            new_n = max_n_all + 1
            new_name = f"{new_n}_amr_extract"
            outdir = os.path.join(tools_dir, new_name)
            os.makedirs(outdir, exist_ok=False)
            print(f"No existing *_amr_extract dirs, created: {outdir}")
            return outdir
    else:
        max_n_all = max([n for n, _, _ in existing], default=0)
        new_n = max_n_all + 1
        new_name = f"{new_n}_amr_extract"
        outdir = os.path.join(tools_dir, new_name)
        os.makedirs(outdir, exist_ok=False)
        print(f"Output directory created: {outdir}")
        return outdir


def main():
    # custom help to support -h, -help, --help
    parser = argparse.ArgumentParser(
        description="Extract per-gene and per-class AMR sequences using fixed results/ directory structure",
        add_help=False,
    )
    parser.add_argument(
        "-h", "-help", "--help",
        action="help",
        help="show this help message and exit",
    )
    parser.add_argument(
        "--summary",
        default=os.path.join("results", "OneHealthAMR_AMRFinder_summary.xlsx"),
        help="Path to AMR summary file (default: results/OneHealthAMR_AMRFinder_summary.xlsx)",
    )
    parser.add_argument(
        "--tools-dir",
        default=os.path.join("results", "tools"),
        help="Path to results/tools directory (default: results/tools)",
    )
    parser.add_argument(
        "--dedup",
        action="store_true",
        help="Deduplicate exact identical sequences within each gene/class file",
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=0,
        help="Add N bp upstream/downstream around hits (default 0)",
    )
    parser.add_argument(
        "--translate",
        action="store_true",
        help="Also create protein FASTA per gene and per class (simple translation)",
    )
    parser.add_argument(
        "--mapping",
        action="store_true",
        help="Write mapping CSV (sample,gene,header,seq_len)",
    )
    parser.add_argument(
        "--cdhit",
        type=float,
        default=None,
        help="Optional: run cd-hit-est on each output file with given identity (e.g. 0.99)",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume in latest N_amr_extract directory and skip genes/classes that already have output files",
    )

    args = parser.parse_args()

    summary_path = args.summary
    tools_dir = args.tools_dir

    if not os.path.exists(summary_path):
        print(f"ERROR: Summary file not found at: {summary_path}", file=sys.stderr)
        sys.exit(1)

    # Read summary
    sfx = summary_path.lower()
    if sfx.endswith((".xls", ".xlsx")):
        df = pd.read_excel(summary_path)
    elif sfx.endswith(".csv"):
        df = pd.read_csv(summary_path)
    else:
        # try TSV/CSV autodetect
        df = pd.read_csv(summary_path, sep=None, engine="python")

    cols = list(df.columns)

    col_sample = choose_column(cols, SAMPLE_CANDS)
    col_contig = choose_column(cols, CONTIG_CANDS)
    col_start = choose_column(cols, START_CANDS)
    col_end = choose_column(cols, END_CANDS)
    col_strand = choose_column(cols, STRAND_CANDS, default=None)
    col_gene = choose_column(cols, GENE_CANDS)
    col_class = choose_column(cols, CLASS_CANDS)

    # Mandatory
    if not (col_sample and col_contig and col_start and col_end):
        print("ERROR: required columns not found in summary.", file=sys.stderr)
        print("Columns present:", cols, file=sys.stderr)
        print("Need at least: sample, contig, start, end (start/stop OK)", file=sys.stderr)
        sys.exit(2)

    print("Detected columns mapping:")
    print("  sample:", col_sample)
    print("  contig:", col_contig)
    print("  start :", col_start)
    print("  end   :", col_end)
    print("  strand:", col_strand)
    print("  gene  :", col_gene)
    print("  class :", col_class)

    # Create or reuse output dir: results/tools/N_amr_extract
    outdir = determine_outdir(tools_dir, resume=args.resume)
    per_gene_dir = os.path.join(outdir, "per_gene")
    per_class_dir = os.path.join(outdir, "per_class")
    os.makedirs(per_gene_dir, exist_ok=True)
    os.makedirs(per_class_dir, exist_ok=True)

    # Save command used into command.txt
    cmdline = " ".join(shlex.quote(a) for a in sys.argv)
    cmdfile = os.path.join(outdir, "command.txt")
    cmd_mode = "a" if os.path.exists(cmdfile) else "w"
    with open(cmdfile, cmd_mode) as cf:
        cf.write(cmdline + "\n")
    print(f"Saved command to {cmdfile}")

    # Build groups
    gene_groups = defaultdict(list)
    class_groups = defaultdict(list)
    mapping_rows = []

    for idx, row in df.iterrows():
        sample = row.get(col_sample)
        contig = row.get(col_contig)
        start = row.get(col_start)
        end = row.get(col_end)
        strand = row.get(col_strand) if col_strand else "+"
        gene = row.get(col_gene) if col_gene is not None else None
        cls = row.get(col_class) if col_class is not None else None

        if pd.isna(sample) or pd.isna(contig) or pd.isna(start) or pd.isna(end):
            print(f"Skipping incomplete row {idx} (missing sample/contig/start/end)")
            continue

        try:
            start = int(start)
            end = int(end)
        except Exception:
            print(f"Skipping row {idx} due to non-integer start/end: start={start}, end={end}")
            continue

        gene_key = str(gene) if gene is not None and not pd.isna(gene) else f"{contig}:{start}-{end}"
        class_key = str(cls) if cls is not None and not pd.isna(cls) else None

        hit = {
            "sample": str(sample),
            "contig": str(contig),
            "start": start,
            "end": end,
            "strand": str(strand),
            "gene": gene_key,
            "class": class_key,
        }
        gene_groups[gene_key].append(hit)
        if class_key:
            class_groups[class_key].append(hit)

    # cache assemblies by sample
    assembly_cache = {}

    def get_assembly_dict(sample_name):
        if sample_name in assembly_cache:
            return assembly_cache[sample_name]

        # Assemblies under: results/tools/*_assembly/<sample>/<sample>.assembly.fasta
        pattern = os.path.join(tools_dir, "*_assembly", sample_name, f"{sample_name}.assembly.fasta")
        matches = glob.glob(pattern)

        if not matches:
            print(f"Assembly not found for sample '{sample_name}' with pattern: {pattern}")
            assembly_cache[sample_name] = None
            return None

        if len(matches) > 1:
            print(f"Warning: multiple assemblies found for sample '{sample_name}'. Using first: {matches[0]}")
        path = matches[0]

        try:
            assembly_cache[sample_name] = load_assembly_index(path)
            return assembly_cache[sample_name]
        except Exception as e:
            print(f"Failed to load assembly {path} for sample {sample_name}: {e}")
            assembly_cache[sample_name] = None
            return None

    def extract_hits_to_records(hits):
        records = []
        seen_hashes = set()
        for h in hits:
            sample_name = h["sample"]
            seqdict = get_assembly_dict(sample_name)
            if seqdict is None:
                continue
            contig = h["contig"]
            if contig not in seqdict:
                print(f"Contig '{contig}' not found in sample '{sample_name}', skipping hit")
                continue
            seq = seqdict[contig].seq
            s = max(1, h["start"] - args.flank)
            e = min(len(seq), h["end"] + args.flank)
            subseq = seq[s - 1 : e]
            if str(h["strand"]).strip() in ("-", "-1", "minus"):
                subseq = subseq.reverse_complement()
            seqstr = str(subseq)
            header = f"{sample_name}|{h['gene']}|{contig}:{s}-{e}|strand={h['strand']}"
            if args.dedup:
                hsh = sha256_seq(seqstr)
                if hsh in seen_hashes:
                    continue
                seen_hashes.add(hsh)
            records.append((header, seqstr))
            if args.mapping:
                mapping_rows.append(
                    {"gene": h["gene"], "sample": sample_name, "header": header, "seq_len": len(seqstr)}
                )
        return records

    # Per-gene FASTAs
    print("\nExtracting per-gene FASTAs...")
    for gene, hits in gene_groups.items():
        safe_gene = clean_name(gene)
        out_fa = os.path.join(per_gene_dir, f"{safe_gene}.fa")

        # resume: skip if file already exists and non-empty
        if args.resume and os.path.exists(out_fa) and os.path.getsize(out_fa) > 0:
            print(f"[RESUME] Skipping existing gene file: {out_fa}")
            continue

        recs = extract_hits_to_records(hits)
        if not recs:
            continue
        write_fasta_records(out_fa, recs)
        print(f"Wrote {len(recs)} records to {out_fa}")

        # optional translation
        if args.translate:
            prot_fa = os.path.join(per_gene_dir, f"{safe_gene}_protein.fa")
            if args.resume and os.path.exists(prot_fa) and os.path.getsize(prot_fa) > 0:
                print(f"[RESUME] Skipping existing protein gene file: {prot_fa}")
            else:
                with open(prot_fa, "w") as pf:
                    for header, seqstr in recs:
                        prot = str(Seq(seqstr).translate(to_stop=False))
                        pf.write(f">{header}\n{prot}\n")
                print(f"Wrote protein translations to {prot_fa}")

        # optional cd-hit clustering
        if args.cdhit:
            clustered_path = os.path.join(per_gene_dir, f"{safe_gene}.cdhit.fa")
            if run_cdhit_est(out_fa, clustered_path, args.cdhit):
                os.replace(clustered_path, out_fa)
                print(f"cd-hit-est clustered {out_fa} at {args.cdhit} identity")

    # Per-class FASTAs
    print("\nExtracting per-class FASTAs...")
    for cls, hits in class_groups.items():
        safe_cls = clean_name(cls)
        out_fa = os.path.join(per_class_dir, f"{safe_cls}.fa")

        if args.resume and os.path.exists(out_fa) and os.path.getsize(out_fa) > 0:
            print(f"[RESUME] Skipping existing class file: {out_fa}")
            continue

        recs = extract_hits_to_records(hits)
        if not recs:
            continue
        write_fasta_records(out_fa, recs)
        print(f"Wrote {len(recs)} records to {out_fa}")

        if args.translate:
            prot_fa = os.path.join(per_class_dir, f"{safe_cls}_protein.fa")
            if args.resume and os.path.exists(prot_fa) and os.path.getsize(prot_fa) > 0:
                print(f"[RESUME] Skipping existing protein class file: {prot_fa}")
            else:
                with open(prot_fa, "w") as pf:
                    for header, seqstr in recs:
                        prot = str(Seq(seqstr).translate(to_stop=False))
                        pf.write(f">{header}\n{prot}\n")
                print(f"Wrote protein translations to {prot_fa}")

        if args.cdhit:
            clustered_path = os.path.join(per_class_dir, f"{safe_cls}.cdhit.fa")
            if run_cdhit_est(out_fa, clustered_path, args.cdhit):
                os.replace(clustered_path, out_fa)
                print(f"cd-hit-est clustered {out_fa} at {args.cdhit} identity")

    # Mapping CSV
    if args.mapping:
        map_csv = os.path.join(outdir, "amr_gene_extraction_mapping.csv")
        mode = "a" if args.resume and os.path.exists(map_csv) else "w"
        write_header = not (mode == "a" and os.path.exists(map_csv) and os.path.getsize(map_csv) > 0)

        with open(map_csv, mode, newline="") as csvf:
            w = csv.DictWriter(csvf, fieldnames=["gene", "sample", "header", "seq_len"])
            if write_header:
                w.writeheader()
            for r in mapping_rows:
                w.writerow(r)
        print("Updated mapping CSV:", map_csv)

    print("\nDone.")
    print(" Output root:", outdir)
    print(" Per-gene FASTAs:", per_gene_dir)
    print(" Per-class FASTAs:", per_class_dir)


if __name__ == "__main__":
    main()
