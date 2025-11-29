#!/usr/bin/env python3
# parse_checkm_dict.py
# Usage: ./parse_checkm_dict.py path/to/bin_stats_ext.tsv > parsed_checkm.tsv

import sys, ast, io

if len(sys.argv) != 2:
    print("Usage: parse_checkm_dict.py bin_stats_ext.tsv", file=sys.stderr)
    sys.exit(2)

infile = sys.argv[1]

out = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
print("bin\tcompleteness\tcontamination", file=out)

with open(infile, 'r', encoding='utf-8') as fh:
    for raw in fh:
        line = raw.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        # split only on the first tab (some lines may have large dict in 2nd column)
        parts = line.split("\t", 1)
        if len(parts) == 1:
            # line contains only one column — skip or try to parse as JSON-like string
            continue
        bin_name, dict_str = parts[0].strip(), parts[1].strip()
        try:
            # bin_stats_ext uses Python-ish dicts with single quotes — use ast.literal_eval
            d = ast.literal_eval(dict_str)
            # sometimes checkm prints nested dict under 'marker lineage' etc.
            # completeness key may be 'Completeness' (case as shown)
            comp = d.get('Completeness') if isinstance(d, dict) else None
            cont = d.get('Contamination') if isinstance(d, dict) else None
            # fallback: try lowercase keys
            if comp is None:
                comp = d.get('completeness') if isinstance(d, dict) else None
            if cont is None:
                cont = d.get('contamination') if isinstance(d, dict) else None
            # If still None, try to search inside string for numbers (best-effort)
            if comp is None or cont is None:
                s = dict_str
                # very naive regex fallback
                import re
                if comp is None:
                    m = re.search(r"['\"]?Completeness['\"]?\s*:\s*([0-9]+(?:\.[0-9]+)?)", s, re.IGNORECASE)
                    comp = float(m.group(1)) if m else ""
                if cont is None:
                    m = re.search(r"['\"]?Contamination['\"]?\s*:\s*([0-9]+(?:\.[0-9]+)?)", s, re.IGNORECASE)
                    cont = float(m.group(1)) if m else ""
        except Exception as e:
            # if parsing fails, try to extract numbers with regex as last resort
            import re
            m1 = re.search(r"Completeness'\s*:\s*([0-9]+(?:\.[0-9]+)?)", dict_str)
            m2 = re.search(r"Contamination'\s*:\s*([0-9]+(?:\.[0-9]+)?)", dict_str)
            comp = float(m1.group(1)) if m1 else ""
            cont = float(m2.group(1)) if m2 else ""
        print(f"{bin_name}\t{comp}\t{cont}", file=out)
