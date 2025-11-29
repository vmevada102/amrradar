#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# 0. Basic paths
# -----------------------------
TOOLS_ROOT="results/tools"
SAMPLE_TSV="sample.tsv"   # path to your sample metadata file

if [[ ! -d "$TOOLS_ROOT" ]]; then
  echo "ERROR: $TOOLS_ROOT directory not found. Run this script from the project root."
  exit 1
fi

# -----------------------------
# 1. Find next N for N_biomarker
# -----------------------------
max_n=$(find "$TOOLS_ROOT" -maxdepth 1 -type d -regex '.*/[0-9]+_.*' \
  | sed 's#.*/\([0-9]\+\)_.*#\1#' \
  | sort -n | tail -n1 || true)

if [[ -z "${max_n:-}" ]]; then
  N=1
else
  N=$((max_n + 1))
fi

BIOMARKER_DIR="${TOOLS_ROOT}/${N}_biomarker"
mkdir -p "$BIOMARKER_DIR"
echo "Using biomarker directory: $BIOMARKER_DIR"

mkdir -p "$BIOMARKER_DIR/prokka_gff"
mkdir -p "$BIOMARKER_DIR/roary"
mkdir -p "$BIOMARKER_DIR/scripts"

# -----------------------------
# 2. Check tools in current conda env
# -----------------------------
required_tools=(roary panaroo snp-sites FastTree Rscript pyseer)

missing=()
for t in "${required_tools[@]}"; do
  if ! command -v "$t" >/dev/null 2>&1; then
    missing+=("$t")
  fi
done

if ((${#missing[@]} > 0)); then
  echo "The following tools are MISSING in your current conda environment:"
  for t in "${missing[@]}"; do
    case "$t" in
      roary)
        echo "  - roary (install:  conda install -c bioconda roary)"
        ;;
      panaroo)
        echo "  - panaroo (install:  conda install -c bioconda panaroo)"
        ;;
      snp-sites)
        echo "  - snp-sites (install:  conda install -c bioconda snp-sites)"
        ;;
      FastTree)
        echo "  - FastTree (install:  conda install -c bioconda fasttree)"
        ;;
      Rscript)
        echo "  - Rscript (install:  conda install -c conda-forge r-base)"
        ;;
      pyseer)
        echo "  - pyseer (install:  conda install -c bioconda pyseer)"
        ;;
      *)
        echo "  - $t (please install from bioconda/conda-forge as appropriate)"
        ;;
    esac
  done
  echo
  echo "Please install the missing tools in your active conda environment and rerun this script."
  exit 1
fi

echo "All required tools are available in the current conda environment."

# -----------------------------
# 3. Collect Prokka GFF files
# -----------------------------
echo "Linking Prokka GFF files into $BIOMARKER_DIR/prokka_gff ..."

prokka_found=0
for d in "$TOOLS_ROOT"/*_prokka/*; do
  [[ -d "$d" ]] || continue
  sample=$(basename "$d")
  gff=$(ls "$d"/*.gff 2>/dev/null | head -n1 || true)
  if [[ -z "$gff" ]]; then
    echo "  WARNING: No .gff found for sample $sample in $d, skipping."
    continue
  fi
  ln -sf "$(realpath "$gff")" "$BIOMARKER_DIR/prokka_gff/${sample}.gff"
  prokka_found=$((prokka_found+1))
done

if [[ $prokka_found -eq 0 ]]; then
  echo "ERROR: No Prokka GFF files were found under ${TOOLS_ROOT}/*_prokka/*."
  echo "Make sure Prokka has been run and try again."
  exit 1
fi

echo "Linked $prokka_found Prokka GFF files."

# -----------------------------
# 4. Run Roary pangenome
# -----------------------------
echo "Running Roary pangenome analysis..."

roary -e -n -p 8 \
  -f "$BIOMARKER_DIR/roary" \
  "$BIOMARKER_DIR"/prokka_gff/*.gff

echo "Roary finished. Outputs in: $BIOMARKER_DIR/roary"

# -----------------------------
# 5. Core-genome SNPs and tree
# -----------------------------
CORE_ALN="$BIOMARKER_DIR/roary/core_gene_alignment.aln"

if [[ -f "$CORE_ALN" ]]; then
  echo "Calling core-genome SNPs with snp-sites..."
  snp-sites -v -o "$BIOMARKER_DIR/roary/core_snps.vcf" "$CORE_ALN"

  echo "Building FastTree tree from core alignment..."
  FastTree -nt "$CORE_ALN" > "$BIOMARKER_DIR/roary/core_tree.nwk"
else
  echo "WARNING: core_gene_alignment.aln not found. Roary may have failed or there are too few core genes."
fi

# -----------------------------
# 6. Create phenotypes.txt from sample.tsv (col1: sample, col4: phenotype)
# -----------------------------
if [[ -f "$SAMPLE_TSV" ]]; then
  echo "Creating phenotypes.txt from $SAMPLE_TSV (1st col = sample, 4th col = phenotype)..."
  awk -F'\t' '
    NR==1 {
      pheno=$4;
      print "sample\t" pheno;
      next
    }
    {
      if ($1 != "" && $4 != "")
        print $1 "\t" $4;
    }
  ' "$SAMPLE_TSV" > "$BIOMARKER_DIR/phenotypes.txt"
  echo "Phenotype file written to: $BIOMARKER_DIR/phenotypes.txt"
else
  echo "WARNING: $SAMPLE_TSV not found."
  echo "No phenotypes.txt created. Please create $BIOMARKER_DIR/phenotypes.txt manually."
  echo "Expected format (tab-delimited):"
  echo "  sample<TAB>phenotype"
fi

# -----------------------------
# 7. Create R script: test_gene_biomarkers.R
# -----------------------------
cat <<'EOF' > "$BIOMARKER_DIR/scripts/test_gene_biomarkers.R"
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript test_gene_biomarkers.R gene_presence_absence.csv phenotypes.txt [output.tsv]\n", call. = FALSE)
}

gpa_file  <- args[1]
pheno_file <- args[2]
out_file  <- ifelse(length(args) >= 3, args[3], "gene_biomarker_results.tsv")

cat("Reading Roary gene_presence_absence from:", gpa_file, "\n")
gpa <- read.csv(gpa_file, check.names = FALSE, stringsAsFactors = FALSE)

cat("Reading phenotype data from:", pheno_file, "\n")
pheno <- read.table(pheno_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

if (ncol(pheno) < 2) {
  stop("Phenotype file must have at least two columns: sample <TAB> phenotype", call. = FALSE)
}

sample_col <- colnames(pheno)[1]
pheno_col  <- colnames(pheno)[2]

# Sample names from phenotype file
pheno_samples <- pheno[[sample_col]]

# Identify which columns in gene_presence_absence correspond to samples
gpa_cols <- colnames(gpa)
sample_cols <- intersect(gpa_cols, pheno_samples)

if (length(sample_cols) == 0) {
  stop("No overlapping sample names between gene_presence_absence.csv and phenotypes.txt", call. = FALSE)
}

cat("Number of samples with both genotype and phenotype:", length(sample_cols), "\n")

# Subset phenotype to the samples we actually have in gpa
pheno_sub <- pheno[match(sample_cols, pheno[[sample_col]]), , drop = FALSE]
y <- pheno_sub[[pheno_col]]

if (!all(y %in% c(0, 1))) {
  warning("Phenotype column is not strictly 0/1, attempting to coerce to numeric.")
  y <- as.numeric(as.character(y))
}

# Roary usually has a column "Gene" with gene IDs
if (!"Gene" %in% colnames(gpa)) {
  stop("Expected a column named 'Gene' in gene_presence_absence.csv, but it was not found.", call. = FALSE)
}

genes <- gpa[["Gene"]]
n_genes <- length(genes)
cat("Number of genes in Roary table:", n_genes, "\n")

results <- data.frame(
  Gene = character(0),
  OR   = numeric(0),
  pval = numeric(0),
  stringsAsFactors = FALSE
)

for (i in seq_len(n_genes)) {
  gene <- genes[i]
  
  # presence = non-empty string in that sample's cell
  row <- gpa[i, sample_cols, drop = TRUE]
  pres <- as.numeric(row != "")

  # skip genes that are all present or all absent (no variation)
  if (all(pres == 0) || all(pres == 1)) {
    next
  }
  
  tab <- table(pres, y)
  if (any(dim(tab) != c(2, 2))) {
    # if for some reason we don't get a 2x2 table, skip
    next
  }
  
  ft <- fisher.test(tab)
  or <- as.numeric(ft$estimate)
  p  <- ft$p.value
  
  results <- rbind(
    results,
    data.frame(Gene = gene, OR = or, pval = p, stringsAsFactors = FALSE)
  )
}

if (nrow(results) == 0) {
  stop("No informative genes (variable presence and valid 2x2 tables) found.", call. = FALSE)
}

results$FDR <- p.adjust(results$pval, method = "BH")

# Sort by p-value
results <- results[order(results$pval), ]

cat("Writing results to:", out_file, "\n")
write.table(results, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done.\n")
EOF

chmod +x "$BIOMARKER_DIR/scripts/test_gene_biomarkers.R"

# -----------------------------
# 8. Create Python script: roary_csv_to_pyseer_matrix.py
# -----------------------------
cat <<'EOF' > "$BIOMARKER_DIR/scripts/roary_csv_to_pyseer_matrix.py"
#!/usr/bin/env python3
import sys
import csv

if len(sys.argv) < 3:
    sys.stderr.write(
        "Usage: roary_csv_to_pyseer_matrix.py gene_presence_absence.csv output_matrix.txt\n"
    )
    sys.exit(1)

gpa_file = sys.argv[1]
out_file = sys.argv[2]

# Roary metadata columns (may vary slightly between versions, but these are standard)
metadata_cols = {
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory order with Fragment",
    "QC",
    "Min sequence length",
    "Max sequence length",
    "Avg sequence length",
}

with open(gpa_file, newline="") as fh:
    reader = csv.reader(fh)
    header = next(reader)

    # Identify sample columns (those not in metadata_cols)
    sample_cols_idx = []
    sample_names = []
    for idx, col in enumerate(header):
        if col not in metadata_cols:
            sample_cols_idx.append(idx)
            sample_names.append(col)

    if not sample_cols_idx:
        sys.stderr.write("ERROR: Could not identify any sample columns in Roary CSV.\n")
        sys.exit(1)

    sys.stderr.write(
        f"Identified {len(sample_cols_idx)} sample columns: {', '.join(sample_names)}\n"
    )

    genes = []
    # presence[ sample ][ gene ] = 0/1
    presence = {sample: [] for sample in sample_names}

    for row in reader:
        if not row:
            continue
        gene_name = row[0]  # "Gene" column
        if gene_name == "":
            continue
        genes.append(gene_name)

        # For each sample column, mark 1 if non-empty, else 0
        for sample, idx in zip(sample_names, sample_cols_idx):
            cell = row[idx].strip() if idx < len(row) else ""
            val = 1 if cell != "" else 0
            presence[sample].append(val)

# Write pyseer-style matrix: samples in rows, genes in columns
with open(out_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    header_row = ["strain"] + genes
    writer.writerow(header_row)
    for sample in sample_names:
        row = [sample] + presence[sample]
        writer.writerow(row)

sys.stderr.write(f"Wrote pyseer matrix to: {out_file}\n")
EOF

chmod +x "$BIOMARKER_DIR/scripts/roary_csv_to_pyseer_matrix.py"

# -----------------------------
# 9. Final message
# -----------------------------
cat <<EOF

==================================================================
PANGENOME + CORE SNP PIPELINE COMPLETED
Biomarker analysis folder:  $BIOMARKER_DIR
------------------------------------------------------------------
Generated files:
- Prokka GFF links:     $BIOMARKER_DIR/prokka_gff/
- Roary outputs:        $BIOMARKER_DIR/roary/
  - gene_presence_absence.csv
  - core_gene_alignment.aln
  - core_snps.vcf (if core alignment present)
  - core_tree.nwk (if core alignment present)
- Phenotypes (if sample.tsv found):
  - $BIOMARKER_DIR/phenotypes.txt (sample name from col1, phenotype from col4)
- Helper scripts:
  - $BIOMARKER_DIR/scripts/test_gene_biomarkers.R
  - $BIOMARKER_DIR/scripts/roary_csv_to_pyseer_matrix.py

Next examples:

1) Test gene presence vs phenotype (Fisher test, R):
   cd $BIOMARKER_DIR
   Rscript scripts/test_gene_biomarkers.R \\
       roary/gene_presence_absence.csv \\
       phenotypes.txt \\
       gene_biomarker_results.tsv

2) Create pyseer gene matrix and run pyseer:
   cd $BIOMARKER_DIR
   python scripts/roary_csv_to_pyseer_matrix.py \\
       roary/gene_presence_absence.csv \\
       gene_matrix.txt

   pyseer --phenotypes phenotypes.txt \\
          --phenotype-col \$(head -1 phenotypes.txt | cut -f2) \\
          --pres gene_matrix.txt \\
          > pyseer_gene_results.txt

==================================================================
EOF
