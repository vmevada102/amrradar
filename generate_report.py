#!/usr/bin/env python3
"""
generate_report.py

Generates an Excel summary from One Health-AMR-Radar outputs (tools layout).
Kraken2 processing removed.
Usage:
  python3 generate_report.py -i sample.tsv -o OUT.xlsx -r OUTDIR
"""
import argparse
from pathlib import Path
import gzip
import sys
import subprocess
import pandas as pd

def open_maybe_gzip(path):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(str(p))
    if str(p).endswith('.gz'):
        return gzip.open(str(p), 'rt', errors='ignore')
    else:
        return open(str(p), 'r', errors='ignore')

def count_fastq_seqs(path):
    p = Path(path)
    if not p.exists():
        return None
    # try seqkit if available
    try:
        subprocess.run(["seqkit","--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out = subprocess.check_output(["seqkit","stats","-a","-T","-j","1", str(p)], universal_newlines=True)
        # seqkit stats columns: file, format, type, num_seqs, ...
        vals = out.strip().split("\n")
        total = 0
        for line in vals:
            cols = line.split('\t')
            if len(cols) >= 4:
                try:
                    total += int(cols[3])
                except:
                    pass
        return total
    except Exception:
        # fallback: count lines and divide by 4
        try:
            if str(p).endswith('.gz'):
                cmd = ["gzip","-cd", str(p)]
            else:
                cmd = ["cat", str(p)]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            wc_out = subprocess.check_output(["wc","-l"], stdin=proc.stdout)
            proc.wait()
            lines = int(wc_out.strip().split()[0])
            return lines // 4
        except Exception:
            # slow Python fallback
            cnt = 0
            with open_maybe_gzip(p) as fh:
                for _ in fh:
                    cnt += 1
            return cnt // 4

def count_fasta_records(path):
    p = Path(path)
    if not p.exists():
        return None
    # try seqkit
    try:
        subprocess.run(["seqkit","--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out = subprocess.check_output(["seqkit","stats","-a","-T","-j","1", str(p)], universal_newlines=True)
        vals = out.strip().split("\n")
        total = 0
        for line in vals:
            cols = line.split('\t')
            if len(cols) >= 4:
                try:
                    total += int(cols[3])
                except:
                    pass
        return total
    except Exception:
        # fallback: count '>' lines
        c = 0
        with open_maybe_gzip(p) as fh:
            for line in fh:
                if line.startswith('>'):
                    c += 1
        return c

def read_checkm_parsed(checkm_path):
    """
    Read parsed_checkm.tsv robustly:
    - If file has header with column names (bin, completeness, contamination), use column name.
    - If file is headerless, fall back to positional indexing (first row, second column).
    Returns completeness (float) or None.
    """
    p = Path(checkm_path)
    if not p.exists():
        return None
    try:
        # peek first line to see if header present
        with p.open('r', encoding='utf-8') as fh:
            first = fh.readline().strip()
            if not first:
                return None
            parts = first.split('\t')
            # if second column is non-numeric, assume header present
            second = parts[1] if len(parts) > 1 else ""
            header_present = False
            try:
                float(second)
                header_present = False
            except Exception:
                # not numeric -> header likely present
                header_present = True

        if header_present:
            df = pd.read_csv(p, sep='\t', header=0)
            # try common column names
            for colname in ['completeness','Completeness','COMPLETE','checkm_completeness']:
                if colname in df.columns:
                    val = df[colname].iloc[0]
                    try:
                        return float(val)
                    except Exception:
                        return None
            # fallback: use second column by position
            if df.shape[1] >= 2:
                try:
                    return float(df.iloc[0,1])
                except Exception:
                    return None
            return None
        else:
            # headerless - read as no-header and fetch second column of first row
            df = pd.read_csv(p, sep='\t', header=None)
            if df.shape[0] >= 1 and df.shape[1] >= 2:
                try:
                    return float(df.iloc[0,1])
                except Exception:
                    # try to coerce numeric
                    try:
                        return float(str(df.iloc[0,1]).strip())
                    except:
                        return None
            return None
    except Exception:
        return None

def main():
    parser = argparse.ArgumentParser(description="Generate Excel summary from One Health-AMR-Radar outputs.")
    parser.add_argument('-i','--samples', required=True, help='sample TSV (tab-delimited) with columns: sample r1 r2 group')
    parser.add_argument('-o','--out', required=True, help='output Excel path (.xlsx)')
    parser.add_argument('-r','--root', required=True, help='pipeline OUTDIR root (contains tools/...)')
    args = parser.parse_args()

    samples_tsv = Path(args.samples)
    out_root = Path(args.root)
    tools_root = out_root / "tools"
    out_xlsx = Path(args.out)

    if not samples_tsv.exists():
        print(f"ERROR: samples TSV not found: {samples_tsv}", file=sys.stderr)
        sys.exit(1)

    df_samples = pd.read_csv(samples_tsv, sep='\t', dtype=str).fillna("")
    for col in ['sample','r1','r2','group']:
        if col not in df_samples.columns:
            print(f"ERROR: missing column {col} in sample TSV", file=sys.stderr)
            sys.exit(1)

    rows = []
    for _, row in df_samples.iterrows():
        sample = row['sample']
        r1 = row['r1']
        r2 = row['r2'] if 'r2' in row and row['r2'] not in ("", "NA", "na") else ""
        group = row.get('group','')

        paired = "Yes" if r2 else "No"
        # raw counts
        raw_r1 = None; raw_r2 = None
        try:
            if r1 and Path(r1).exists():
                raw_r1 = count_fastq_seqs(r1)
        except Exception:
            raw_r1 = None
        try:
            if r2 and Path(r2).exists():
                raw_r2 = count_fastq_seqs(r2)
        except Exception:
            raw_r2 = None

        # trimmed counts - try read_counts.tsv first, else trimmed file counts
        trim_dir = tools_root / "2_trim" / sample
        trim_r1_path = trim_dir / f"{sample}_R1.trim.fastq.gz"
        trim_r2_path = trim_dir / f"{sample}_R2.trim.fastq.gz"
        rc_file = trim_dir / "read_counts.tsv"
        trimmed_r1 = None; trimmed_r2 = None
        if rc_file.exists():
            try:
                rc_df = pd.read_csv(rc_file, sep='\t')
                if 'trimmed_R1' in rc_df.columns:
                    v = rc_df['trimmed_R1'].iloc[-1]; trimmed_r1 = int(v) if pd.notnull(v) else None
                if 'trimmed_R2' in rc_df.columns:
                    v = rc_df['trimmed_R2'].iloc[-1]; trimmed_r2 = int(v) if pd.notnull(v) else None
                if 'trimmed_reads' in rc_df.columns and (trimmed_r1 is None and trimmed_r2 is None):
                    tot = int(rc_df['trimmed_reads'].iloc[-1]) if pd.notnull(rc_df['trimmed_reads'].iloc[-1]) else None
                    if tot is not None:
                        if r2:
                            trimmed_r1 = trimmed_r2 = tot // 2
                        else:
                            trimmed_r1 = tot
            except Exception:
                pass
        if trimmed_r1 is None and trim_r1_path.exists():
            try: trimmed_r1 = count_fastq_seqs(trim_r1_path)
            except: trimmed_r1 = None
        if r2 and trimmed_r2 is None and trim_r2_path.exists():
            try: trimmed_r2 = count_fastq_seqs(trim_r2_path)
            except: trimmed_r2 = None

        # assembly stats
        stats_file = tools_root / "4_stats" / sample / "assembly_stats.tsv"
        assembly_fasta = tools_root / "3_assembly" / sample / f"{sample}.assembly.fasta"
        totalbp = None; contigs = None; n50 = None; l50 = None; gc = None; largest = None
        if stats_file.exists():
            try:
                st = pd.read_csv(stats_file, sep='\t')
                if st.shape[0] >= 1:
                    info = st.iloc[-1]
                    totalbp = int(info.get('TotalBP')) if pd.notnull(info.get('TotalBP')) else None
                    contigs = int(info.get('Contigs')) if pd.notnull(info.get('Contigs')) else None
                    n50 = int(info.get('N50')) if pd.notnull(info.get('N50')) else None
                    l50 = int(info.get('L50')) if pd.notnull(info.get('L50')) else None
                    gc = float(info.get('GC')) if pd.notnull(info.get('GC')) else None
                    largest = int(info.get('Largest')) if pd.notnull(info.get('Largest')) else None
            except Exception:
                pass
        if contigs is None and assembly_fasta.exists():
            try:
                contigs = count_fasta_records(assembly_fasta)
            except Exception:
                contigs = None

        assembled_seq_count = contigs

        # CheckM completeness
        checkm_file = tools_root / "5_checkm" / sample / "parsed_checkm.tsv"
        checkm_compl = None
        if checkm_file.exists():
            try:
                checkm_compl = read_checkm_parsed(checkm_file)
            except Exception:
                checkm_compl = None

        rows.append({
            'Name of group': group,
            'Name of sample': sample,
            'Paired (Yes/No)': "Yes" if r2 else "No",
            'Total No of raw sequences in R1': raw_r1,
            'Total No of raw sequences in R2': raw_r2,
            'Total No of sequences in Assembled genome': assembled_seq_count,
            'Total Assembly Size in Base pairs': totalbp,
            'Number of Contigs/Scaffolds': contigs,
            'N50': n50,
            'L50': l50,
            'GC Content': gc,
            'Largest Contig': largest,
            'CheckM completeness (%)': checkm_compl
        })

    out_df = pd.DataFrame(rows)
    cols = [
        'Name of group','Name of sample','Paired (Yes/No)',
        'Total No of raw sequences in R1','Total No of raw sequences in R2',
        'Total No of sequences in Assembled genome',
        'Total Assembly Size in Base pairs','Number of Contigs/Scaffolds','N50','L50',
        'GC Content','Largest Contig','CheckM completeness (%)'
    ]
    out_df = out_df[cols]
    try:
        out_df.to_excel(out_xlsx, index=False)
        print(f"Wrote summary to {out_xlsx}")
    except Exception as e:
        print(f"ERROR writing Excel: {e}", file=sys.stderr)
        csv_out = out_xlsx.with_suffix('.csv')
        out_df.to_csv(csv_out, index=False)
        print(f"Wrote CSV fallback to {csv_out}")

if __name__ == '__main__':
    main()
