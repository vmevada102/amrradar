#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Configuration
###############################################################################

SAMPLES_TSV="sample.tsv"        # sample<tab>r1<tab>r2<tab>group
WINDOW_SIZE=1000
TOOLS_DIR="results/tools"

THREADS_PER_JOB=1
JOBS=1

###############################################################################
# Parse arguments (support -resume)
###############################################################################

RESUME_FORCE=0
if [[ "${1-}" == "-resume" || "${1-}" == "--resume" ]]; then
  RESUME_FORCE=1
  shift || true
fi

###############################################################################
# Reference FASTA (default: results/sequence.fasta)
###############################################################################

REFERENCE_DEFAULT="results/sequence.fasta"
REFERENCE="${REFERENCE:-$REFERENCE_DEFAULT}"

if [[ ! -f "$REFERENCE" ]]; then
  echo "Default reference not found at: $REFERENCE"
  read -rp "Enter path to reference genome FASTA: " REFERENCE
fi

if [[ ! -f "$REFERENCE" ]]; then
  echo "ERROR: reference FASTA not found at: $REFERENCE"
  exit 1
fi

echo "Using reference FASTA: $REFERENCE"

###############################################################################
# Auto-detect cores and set threads / parallel samples
###############################################################################

if command -v nproc &>/dev/null; then
  CORES=$(nproc)
elif [[ "$(uname)" == "Darwin" ]] && command -v sysctl &>/dev/null; then
  CORES=$(sysctl -n hw.ncpu)
else
  CORES=4
fi

USABLE=$((CORES - 2))
(( USABLE < 1 )) && USABLE=1

DEFAULT_THREADS_PER_JOB=4
if (( USABLE <= DEFAULT_THREADS_PER_JOB )); then
  THREADS_PER_JOB=$USABLE
  JOBS=1
else
  THREADS_PER_JOB=$DEFAULT_THREADS_PER_JOB
  JOBS=$(( USABLE / THREADS_PER_JOB ))
  (( JOBS < 1 )) && JOBS=1
fi

echo "Total cores        : $CORES"
echo "Usable (cores - 2) : $USABLE"
echo "Threads per sample : $THREADS_PER_JOB"
echo "Parallel samples   : $JOBS"

###############################################################################
# Check required CLI tools
###############################################################################

REQUIRED_TOOLS=(minimap2 samtools bcftools bedtools awk snp-sites snp-dists FastTree)

MISSING=()
for t in "${REQUIRED_TOOLS[@]}"; do
  if ! command -v "$t" &>/dev/null; then
    MISSING+=("$t")
  fi
done

if (( ${#MISSING[@]} > 0 )); then
  echo "ERROR: Missing required tools: ${MISSING[*]}"
  if command -v mamba &>/dev/null || command -v conda &>/dev/null; then
    echo "Current conda env: ${CONDA_DEFAULT_ENV:-<unknown>}"
    read -rp "Attempt to install missing tools into THIS env and then exit? [y/N] " ans
    if [[ "$ans" =~ ^[Yy]$ ]]; then
      if command -v mamba &>/dev/null; then
        INSTALLER=mamba
      else
        INSTALLER=conda
      fi
      echo "Installing with $INSTALLER..."
      "$INSTALLER" install -y "${MISSING[@]}"
      echo "Installation complete. Please re-run this script; it will resume."
      exit 0
    else
      echo "Please install: ${MISSING[*]} and re-run."
      exit 1
    fi
  else
    echo "Conda/mamba not found. Install tools manually and re-run."
    exit 1
  fi
fi

###############################################################################
# Check Rscript and required R packages, auto-install if missing (non-interactive)
###############################################################################

REQUIRED_R_PACKAGES=(
  data.table
  ggplot2
  viridis
  ggpubr
  FSA
  vegan
  factoextra
  openxlsx
  patchwork
  gridExtra
  reshape2
  ape
)

if command -v Rscript &>/dev/null; then
  missing_r_pkgs=$(Rscript - <<'RSCRIPT' || echo ""
    pkgs <- c(
      "data.table","ggplot2","viridis","ggpubr","FSA",
      "vegan","factoextra","openxlsx","patchwork",
      "gridExtra","reshape2","ape"
    )
    inst <- rownames(installed.packages())
    miss <- pkgs[!(pkgs %in% inst)]
    cat(paste(miss, collapse = " "))
RSCRIPT
  )

  if [[ -n "$missing_r_pkgs" ]]; then
    echo "The following R packages are missing:"
    echo "  $missing_r_pkgs"
    echo "Attempting automatic installation via R (non-interactive)..."
    Rscript - <<RSCRIPT || echo "WARNING: R package installation step encountered an error; will re-check remaining packages."
      pkgs <- strsplit("$missing_r_pkgs", " +")[[1]]
      pkgs <- pkgs[pkgs != ""]
      if (length(pkgs) > 0) {
        repos <- "https://cloud.r-project.org"
        install.packages(pkgs, repos = repos)
      }
RSCRIPT

    remaining=$(Rscript - <<'RSCRIPT' || echo ""
      pkgs <- c(
        "data.table","ggplot2","viridis","ggpubr","FSA",
        "vegan","factoextra","openxlsx","patchwork",
        "gridExtra","reshape2","ape"
      )
      inst <- rownames(installed.packages())
      miss <- pkgs[!(pkgs %in% inst)]
      cat(paste(miss, collapse = " "))
RSCRIPT
    )

    if [[ -n "$remaining" ]]; then
      echo "WARNING: Some R packages are still missing after attempted installation:"
      echo "  $remaining"
      echo "R-based plots and statistics that require these packages may fail."
    else
      echo "All required R packages are now installed."
    fi
  else
    echo "All required R packages appear to be installed."
  fi
else
  echo "Rscript not found; all R-based plots and statistics will be skipped."
fi

###############################################################################
# Choose mapping directory with resume facility
# - Scan results/tools/<N>_* for max N (any suffix)
# - If -resume: force latest *_mapping
# - Else: resume unfinished *_mapping, or create <N+1>_mapping
###############################################################################

shopt -s nullglob

max_any_n=0
for d in "$TOOLS_DIR"/*/; do
  base=$(basename "$d")
  if [[ "$base" =~ ^([0-9]+)_ ]]; then
    num=${BASH_REMATCH[1]}
    (( num > max_any_n )) && max_any_n=$num
  fi
done

max_mapping_n=0
for d in "$TOOLS_DIR"/*_mapping/; do
  base=$(basename "$d")
  if [[ "$base" =~ ^([0-9]+)_mapping$ ]]; then
    num=${BASH_REMATCH[1]}
    (( num > max_mapping_n )) && max_mapping_n=$num
  fi
done

shopt -u nullglob

RESUME=0
MAP_DIR=""

if (( RESUME_FORCE == 1 && max_mapping_n > 0 )); then
  MAP_DIR="$TOOLS_DIR/${max_mapping_n}_mapping"
  RESUME=1
  echo "Forced resume (-resume) in: $MAP_DIR"
elif (( max_mapping_n > 0 )); then
  last_mapping_dir="$TOOLS_DIR/${max_mapping_n}_mapping"
  if [[ -d "$last_mapping_dir" && ! -f "$last_mapping_dir/RUN_COMPLETE" ]]; then
    MAP_DIR="$last_mapping_dir"
    RESUME=1
    echo "Resuming existing mapping run in: $MAP_DIR (no RUN_COMPLETE marker found)"
  fi
fi

if (( RESUME == 0 )); then
  NEW_N=$((max_any_n + 1))
  MAP_DIR="$TOOLS_DIR/${NEW_N}_mapping"
  echo "Starting new mapping run in: $MAP_DIR"
fi

mkdir -p "$MAP_DIR"/{bam,vcf,density,plots,logs,phylo}

# Store the command used to run this script
{
  printf '# Command used on %s\n' "$(date)"
  printf '%q ' "$0" "$@"
  echo
} > "$MAP_DIR/command.txt"

###############################################################################
# Index reference
###############################################################################

if [[ ! -f "${REFERENCE}.mmi" ]]; then
  echo "Building minimap2 index for reference..."
  minimap2 -d "${REFERENCE}.mmi" "$REFERENCE"
else
  echo "Found existing minimap2 index: ${REFERENCE}.mmi"
fi

if [[ ! -f "${REFERENCE}.fai" ]]; then
  echo "Indexing reference with samtools faidx..."
  samtools faidx "$REFERENCE"
fi

###############################################################################
# Scan sample.tsv for groups (expected: sample, r1, r2, group)
###############################################################################

SAMPLES_FOR_MAPPING="$MAP_DIR/samples_for_mapping.tsv"

if [[ -f "$SAMPLES_TSV" ]]; then
  echo "Reading sample groups from: $SAMPLES_TSV"
  cp "$SAMPLES_TSV" "$MAP_DIR/original_sample.tsv"
  echo -e "sample\tgroup" > "$SAMPLES_FOR_MAPPING"

  while IFS=$'\t' read -r sample r1 r2 group rest; do
    [[ -z "$sample" ]] && continue
    [[ "$sample" == sample ]] && continue
    [[ "$sample" =~ ^# ]] && continue
    [[ -z "$group" ]] && group="UNGROUPED"
    echo -e "${sample}\t${group}" >> "$SAMPLES_FOR_MAPPING"
  done < "$SAMPLES_TSV"

  echo "Simplified sample-group file written to: $SAMPLES_FOR_MAPPING"
else
  echo "WARNING: $SAMPLES_TSV not found. Samples will be treated as ungrouped."
  SAMPLES_FOR_MAPPING=""
fi

###############################################################################
# Find assemblies
###############################################################################

echo "Scanning for assemblies under: $TOOLS_DIR/*_assembly/<sample>/<sample>.assembly.fasta"

mapfile -t ASSEMBLIES < <(find "$TOOLS_DIR" -maxdepth 3 -type f -name "*.assembly.fasta" | sort)

if (( ${#ASSEMBLIES[@]} == 0 )); then
  echo "ERROR: No *.assembly.fasta files found."
  exit 1
fi

echo "Found ${#ASSEMBLIES[@]} assemblies."

###############################################################################
# Per-sample processing function
###############################################################################

process_one() {
  local fasta="$1"

  local basename_f
  basename_f=$(basename "$fasta")
  local sample="${basename_f%.assembly.fasta}"

  local group="UNGROUPED"
  if [[ -n "${SAMPLES_FOR_MAPPING:-}" && -f "$SAMPLES_FOR_MAPPING" ]]; then
    group=$(awk -v S="$sample" 'NR>1 && $1==S {print $2; exit}' "$SAMPLES_FOR_MAPPING")
    [[ -z "$group" ]] && group="UNGROUPED"
  fi

  local out_bam="$MAP_DIR/bam/${sample}.sorted.bam"
  local out_vcf="$MAP_DIR/vcf/${sample}.vcf.gz"
  local out_density="$MAP_DIR/density/${sample}.density.${WINDOW_SIZE}bp.bed"
  local log="$MAP_DIR/logs/${sample}.log"

  {
    echo "=============================="
    echo "Sample: $sample"
    echo "Group : $group"
    echo "FASTA : $fasta"
    echo "Start : $(date)"
    echo "=============================="

    if [[ ! -f "$out_bam" ]]; then
      echo "[$sample] Mapping with minimap2..."
      minimap2 -t "$THREADS_PER_JOB" -ax asm5 "${REFERENCE}.mmi" "$fasta" \
        | samtools view -bS - \
        | samtools sort -o "$out_bam"
      samtools index "$out_bam"
    else
      echo "[$sample] BAM exists, skipping mapping."
    fi

    if [[ ! -f "$out_vcf" ]]; then
      echo "[$sample] Calling variants with bcftools..."
      bcftools mpileup -Ou -f "$REFERENCE" "$out_bam" \
        | bcftools call -mv -Oz -o "$out_vcf"
      bcftools index "$out_vcf"
    else
      echo "[$sample] VCF exists, skipping variant calling."
    fi

    if [[ ! -f "$out_density" ]]; then
      echo "[$sample] Computing variant density with window size ${WINDOW_SIZE} bp..."
      local tmp_bed="$MAP_DIR/density/${sample}.variants.tmp.bed"
      local windows_bed="$MAP_DIR/density/windows.${WINDOW_SIZE}bp.bed"

      bcftools view -H "$out_vcf" | awk '{print $1"\t"$2-1"\t"$2}' > "$tmp_bed"

      if [[ ! -f "$windows_bed" ]]; then
        bedtools makewindows -g "${REFERENCE}.fai" -w "$WINDOW_SIZE" > "$windows_bed"
      fi

      bedtools intersect -a "$windows_bed" -b "$tmp_bed" -c > "$out_density"
      rm -f "$tmp_bed"
    else
      echo "[$sample] Density file exists, skipping."
    fi

    echo "[$sample] Done at $(date)."
  } &>> "$log"
}

export -f process_one
export MAP_DIR REFERENCE WINDOW_SIZE THREADS_PER_JOB SAMPLES_FOR_MAPPING

###############################################################################
# Run samples in parallel (sample-level) with THREADS_PER_JOB each
###############################################################################

echo "Starting processing with up to $JOBS samples in parallel (THREADS_PER_SAMPLE=$THREADS_PER_JOB)..."

if command -v parallel &>/dev/null; then
  parallel -j "$JOBS" process_one ::: "${ASSEMBLIES[@]}"
else
  running=0
  for fasta in "${ASSEMBLIES[@]}"; do
    process_one "$fasta" &
    ((running++))
    if (( running >= JOBS )); then
      wait -n
      ((running--))
    fi
  done
  wait
fi

echo "All per-sample jobs finished."

###############################################################################
# Merge VCFs across all samples (combined)
###############################################################################

MERGED_VCF="$MAP_DIR/vcf/all_samples.merged.vcf.gz"
if [[ ! -f "$MERGED_VCF" ]]; then
  echo "Merging all sample VCFs into multi-sample VCF..."
  bcftools merge "$MAP_DIR"/vcf/*.vcf.gz -Oz -o "$MERGED_VCF"
  bcftools index "$MERGED_VCF"
else
  echo "Merged VCF already exists, skipping."
fi

###############################################################################
# Group-wise VCF merges
###############################################################################

if [[ -f "$SAMPLES_TSV" ]]; then
  echo "Building group-wise VCF merges from $SAMPLES_TSV"
  rm -f "$MAP_DIR"/vcf/group_*.vcf.list || true

  while IFS=$'\t' read -r sample r1 r2 group rest; do
    [[ -z "$sample" ]] && continue
    [[ "$sample" == sample ]] && continue
    [[ "$sample" =~ ^# ]] && continue
    [[ -z "$group" ]] && continue
    vcf="$MAP_DIR/vcf/${sample}.vcf.gz"
    [[ -f "$vcf" ]] && echo "$vcf" >> "$MAP_DIR/vcf/group_${group}.vcf.list"
  done < "$SAMPLES_TSV"

  for list in "$MAP_DIR"/vcf/group_*.vcf.list; do
    [[ -e "$list" ]] || continue
    grp=$(basename "$list")
    grp=${grp#group_}
    grp=${grp%.vcf.list}
    out_gvcf="$MAP_DIR/vcf/${grp}.group.merged.vcf.gz"
    if [[ ! -f "$out_gvcf" ]]; then
      echo "Merging VCFs for group $grp..."
      bcftools merge -l "$list" -Oz -o "$out_gvcf"
      bcftools index "$out_gvcf"
    else
      echo "Group VCF for $grp exists, skipping."
    fi
  done
else
  echo "WARNING: $SAMPLES_TSV not found; skipping group-wise VCF merging."
fi

###############################################################################
# Per-sample and per-group variant counts
###############################################################################

VARCOUNT_SAMPLE="$MAP_DIR/variant_counts_per_sample.tsv"
VARCOUNT_GROUP="$MAP_DIR/variant_counts_per_group.tsv"

if [[ -f "$MERGED_VCF" ]]; then
  echo "Computing per-sample and per-group non-ref variant counts..."

  echo -e "sample\tgroup\tnonref_variants" > "$VARCOUNT_SAMPLE"

  if [[ -n "${SAMPLES_FOR_MAPPING:-}" && -f "$SAMPLES_FOR_MAPPING" ]]; then
    while IFS=$'\t' read -r sample group; do
      [[ -z "$sample" ]] && continue
      [[ "$sample" == "sample" ]] && continue
      [[ "$sample" =~ ^# ]] && continue
      [[ -z "$group" ]] && group="UNGROUPED"
      count=$(bcftools view -s "$sample" -g ^ref -H "$MERGED_VCF" 2>/dev/null | wc -l || true)
      echo -e "${sample}\t${group}\t${count}" >> "$VARCOUNT_SAMPLE"
    done < "$SAMPLES_FOR_MAPPING"
  else
    for sample in $(bcftools query -l "$MERGED_VCF"); do
      count=$(bcftools view -s "$sample" -g ^ref -H "$MERGED_VCF" 2>/dev/null | wc -l || true)
      echo -e "${sample}\tUNGROUPED\t${count}" >> "$VARCOUNT_SAMPLE"
    done
  fi

  if [[ -n "${SAMPLES_FOR_MAPPING:-}" && -f "$SAMPLES_FOR_MAPPING" ]]; then
    echo -e "group\tsamples_in_group\tnonref_variants" > "$VARCOUNT_GROUP"
    while read -r grp; do
      [[ -z "$grp" ]] && continue
      sample_list=$(awk -v G="$grp" 'NR>1 && $2==G {print $1}' "$SAMPLES_FOR_MAPPING" | paste -sd, -)
      n_samples=$(awk -v G="$grp" 'NR>1 && $2==G {c++} END{print c+0}' "$SAMPLES_FOR_MAPPING")
      if [[ -n "$sample_list" ]]; then
        gcount=$(bcftools view -s "$sample_list" -g ^ref -H "$MERGED_VCF" 2>/dev/null \
          | awk '{print $1"\t"$2}' | sort -u | wc -l || true)
      else
        gcount=0
      fi
      echo -e "${grp}\t${n_samples}\t${gcount}" >> "$VARCOUNT_GROUP"
    done < <(awk 'NR>1 {print $2}' "$SAMPLES_FOR_MAPPING" | sort -u)
  else
    echo "No group information available; skipping group-wise variant counts."
  fi

  echo "Variant count tables written:"
  echo "  - Per-sample : $VARCOUNT_SAMPLE"
  echo "  - Per-group  : $VARCOUNT_GROUP (if groups defined)"
else
  echo "Merged VCF not found; cannot compute variant counts."
fi

###############################################################################
# Diversity plots (sample & group level)
###############################################################################

if command -v Rscript &>/dev/null; then
  VARCOUNT_PLOT_R="$MAP_DIR/plots/plot_variant_counts.R"
  cat > "$VARCOUNT_PLOT_R" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: plot_variant_counts.R <variant_counts_per_sample.tsv> <out_prefix>")
infile <- args[1]
out_prefix <- args[2]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(viridis)
})

dt <- fread(infile)
if (!all(c("sample", "group", "nonref_variants") %in% names(dt))) {
  stop("variant_counts_per_sample.tsv must have columns: sample, group, nonref_variants")
}

dt[, group := as.factor(group)]
dt[, sample := as.factor(sample)]
dt[, sample := reorder(sample, nonref_variants)]

p_sample_bar <- ggplot(dt, aes(x = sample, y = nonref_variants, fill = group)) +
  geom_col() +
  scale_fill_viridis(discrete = TRUE, option = "viridis") +
  labs(
    x = "Sample",
    y = "Non-reference variant count",
    fill = "Group",
    title = "Per-sample variant diversity"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )

ggsave(paste0(out_prefix, "_sample_diversity_bar.pdf"), p_sample_bar, width = 9, height = 5)
ggsave(paste0(out_prefix, "_sample_diversity_bar.png"), p_sample_bar, width = 9, height = 5, dpi = 300)

p_group_box <- ggplot(dt, aes(x = group, y = nonref_variants, fill = group)) +
  geom_boxplot(outlier.alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE, option = "viridis") +
  labs(
    x = "Group",
    y = "Non-reference variants (per sample)",
    fill = "Group",
    title = "Per-group distribution of sample diversity"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggsave(paste0(out_prefix, "_group_diversity_boxplot.pdf"), p_group_box, width = 6, height = 4)
ggsave(paste0(out_prefix, "_group_diversity_boxplot.png"), p_group_box, width = 6, height = 4, dpi = 300)

dt_group_sum <- dt[, .(total_nonref_variants = sum(nonref_variants)), by = group]

p_group_bar <- ggplot(dt_group_sum, aes(x = group, y = total_nonref_variants, fill = group)) +
  geom_col() +
  scale_fill_viridis(discrete = TRUE, option = "viridis") +
  labs(
    x = "Group",
    y = "Total non-reference variants (sum across samples)",
    fill = "Group",
    title = "Total diversity per group"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggsave(paste0(out_prefix, "_group_diversity_bar.pdf"), p_group_bar, width = 6, height = 4)
ggsave(paste0(out_prefix, "_group_diversity_bar.png"), p_group_bar, width = 6, height = 4, dpi = 300)
EOF

  if [[ -f "$VARCOUNT_SAMPLE" ]]; then
    echo "Creating diversity plots from $VARCOUNT_SAMPLE..."
    Rscript "$VARCOUNT_PLOT_R" "$VARCOUNT_SAMPLE" "$MAP_DIR/plots/diversity" || echo "WARNING: diversity plotting failed; continuing."
  else
    echo "variant_counts_per_sample.tsv not found; skipping diversity plots."
  fi
else
  echo "Rscript not found; skipping diversity plots from variant counts."
fi

###############################################################################
# Statistical analysis – Kruskal-Wallis, Dunn (skip if not possible)
###############################################################################

if command -v Rscript &>/dev/null; then
  STATS_R="$MAP_DIR/plots/variant_diversity_stats.R"
  cat > "$STATS_R" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: variant_diversity_stats.R <variant_counts_per_sample.tsv>")
infile <- args[1]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(FSA)
})

dt <- fread(infile)
if (!all(c("sample", "group", "nonref_variants") %in% names(dt))) {
  stop("File must contain: sample, group, nonref_variants")
}

sink(file = paste0(infile,"_normality_tests.txt"))
cat("Normality Tests (Shapiro-Wilk):\n")
for (g in unique(dt$group)) {
  cat("\nGroup: ", g, "\n")
  x <- dt[group == g]$nonref_variants
  if (length(x) < 3) {
    cat("Too few samples for normality test (n < 3)\n")
  } else if (sd(x) == 0) {
    cat("All values identical in this group; Shapiro-Wilk not applicable\n")
  } else {
    print(shapiro.test(x))
  }
}
sink(NULL)

kruskal_res <- kruskal.test(nonref_variants ~ group, data = dt)
sink(file = paste0(infile,"_kruskal_test.txt"))
cat("Kruskal-Wallis rank sum test on nonref_variants ~ group\n\n")
print(kruskal_res)
sink(NULL)

dunn_out <- tryCatch(
  {
    dunnTest(nonref_variants ~ group, data = dt, method = "bh")
  },
  error = function(e) {
    msg_file <- paste0(infile, "_dunn_posthoc_test_SKIPPED.txt")
    sink(msg_file)
    cat("Dunn post-hoc test could not be computed:\n")
    cat(e$message, "\n")
    sink(NULL)
    NULL
  }
)

if (!is.null(dunn_out)) {
  write.csv(dunn_out$res,
            file = paste0(infile,"_dunn_posthoc_test.csv"),
            row.names = FALSE)
}

p <- ggboxplot(dt, x = "group", y = "nonref_variants", add = "jitter") +
  stat_compare_means(method = "kruskal.test") +
  labs(title = "Group-wise Diversity with Significance") +
  theme_bw(base_size=11)

ggsave(paste0(infile,"_stats_boxplot.pdf"), p, width=7, height=5)
ggsave(paste0(infile,"_stats_boxplot.png"), p, width=7, height=5, dpi=300)
EOF

  if [[ -f "$VARCOUNT_SAMPLE" ]]; then
    echo "Running statistical tests on diversity..."
    Rscript "$STATS_R" "$VARCOUNT_SAMPLE" || echo "WARNING: statistical tests failed; continuing."
  else
    echo "variant_counts_per_sample.tsv missing — skipping statistics."
  fi
fi

###############################################################################
# Automatic comprehensive statistical report (robust to zero-variance data)
###############################################################################

if command -v Rscript &>/dev/null; then
  REPORT_R="$MAP_DIR/plots/generate_statistical_report.R"
  cat > "$REPORT_R" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: generate_statistical_report.R <variant_counts.tsv> <output_prefix>")
}
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(viridis)
  library(FSA)
  library(vegan)
  library(factoextra)
  library(openxlsx)
  library(patchwork)
  library(gridExtra)
})

infile <- args[1]
out_prefix <- args[2]
dt <- fread(infile)

if (!all(c("sample","group","nonref_variants") %in% colnames(dt))) {
  stop("variant_counts file must have: sample, group, nonref_variants")
}
dt[, group := factor(group)]

p_box <- ggboxplot(dt, x="group", y="nonref_variants",
                   color="group", fill="group", palette="viridis",
                   add="jitter", title="Variant Diversity per Group") +
         stat_compare_means(method="kruskal.test") +
         theme_bw(base_size=11)

permanova_res <- tryCatch(
  {
    dist_obj <- dist(dt$nonref_variants)
    adonis2(dist_obj ~ group, data=dt, permutations=999)
  },
  error = function(e) {
    msg_file <- paste0(out_prefix,"_PERMANOVA_SKIPPED.txt")
    sink(msg_file)
    cat("PERMANOVA could not be computed:\n")
    cat(e$message, "\n")
    sink(NULL)
    NULL
  }
)

if (!is.null(permanova_res)) {
  write.table(permanova_res, file=paste0(out_prefix,"_PERMANOVA.txt"),
              sep="\t", quote=FALSE)
}

pca_ok <- (length(unique(dt$nonref_variants)) > 1) &&
          (sd(dt$nonref_variants, na.rm = TRUE) > 0) &&
          (nrow(dt) >= 3)

if (pca_ok) {
  vals <- dt$nonref_variants
  mat  <- matrix(vals, ncol = 1)
  pca  <- prcomp(mat, scale. = TRUE, center = TRUE)
  df_pca <- data.frame(
    PC1   = pca$x[,1],
    group = dt$group,
    sample = dt$sample
  )
  p_pca <- ggplot(df_pca, aes(PC1, 0, color=group, label=sample)) +
    geom_point(size=3) +
    geom_text(hjust=1.2, vjust=1) +
    theme_bw(base_size=11) +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    labs(title="PCA of Variant Diversity (PC1)")
} else {
  df_pca <- data.frame(dummy=1)
  p_pca <- ggplot(df_pca) +
    annotate("text", x=1, y=1,
             label="PCA not computed:\nnonref_variants has zero or insufficient variance",
             size=4) +
    theme_void() +
    labs(title="PCA of Variant Diversity")
}

dunn_out <- tryCatch(
  {
    dunnTest(nonref_variants ~ group, data=dt, method="bh")
  },
  error = function(e) {
    msg_file <- paste0(out_prefix, "_DunnTest_SKIPPED.txt")
    sink(msg_file)
    cat("Dunn post-hoc test in report could not be computed:\n")
    cat(e$message, "\n")
    sink(NULL)
    NULL
  }
)

if (!is.null(dunn_out)) {
  write.csv(dunn_out$res, paste0(out_prefix,"_DunnTest.csv"), row.names=FALSE)
}

xlsx_file <- paste0(out_prefix,"_all_results.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "Per_sample")
writeData(wb, "Per_sample", dt)

if (!is.null(permanova_res)) {
  addWorksheet(wb, "PERMANOVA")
  writeData(wb, "PERMANOVA", permanova_res)
}

if (!is.null(dunn_out)) {
  addWorksheet(wb, "Dunn_test")
  writeData(wb, "Dunn_test", dunn_out$res)
}

saveWorkbook(wb, xlsx_file, overwrite = TRUE)

combined <- p_box | p_pca
ggsave(paste0(out_prefix,"_combined_figure.pdf"), combined, width=9, height=5)
ggsave(paste0(out_prefix,"_combined_figure.png"), combined, width=9, height=5, dpi=300)

pdf_report <- paste0(out_prefix,"_statistical_report.pdf")
pdf(pdf_report, width=8, height=10)
grid.arrange(p_box, p_pca, ncol=1)
dev.off()

cat("\nSTATISTICAL ANALYSIS COMPLETED (with PCA/PERMANOVA/Dunn skipped if not applicable).\n")
cat("Excel workbook: ", xlsx_file, "\n")
cat("PDF report:     ", pdf_report, "\n")
cat("Combined fig:   ", paste0(out_prefix,"_combined_figure.pdf"), "\n")
EOF

  if [[ -f "$VARCOUNT_SAMPLE" ]] ; then
    echo "Generating full statistical report..."
    Rscript "$REPORT_R" "$VARCOUNT_SAMPLE" "$MAP_DIR/plots/stats" || echo "WARNING: report generation failed; continuing."
  else
    echo "variant_counts_per_sample.tsv missing — skipping report."
  fi
fi

###############################################################################
# Build combined density table
###############################################################################

COMBINED_DENSITY="$MAP_DIR/density/all_samples.combined.${WINDOW_SIZE}bp.tsv"
if [[ ! -f "$COMBINED_DENSITY" ]]; then
  echo "Building combined density table: $COMBINED_DENSITY"
  echo -e "sample\tchrom\tstart\tend\tcount" > "$COMBINED_DENSITY"
  for f in "$MAP_DIR"/density/*.density.${WINDOW_SIZE}bp.bed; do
    [[ "$f" == *"windows."* ]] && continue
    s=$(basename "$f")
    s=${s%.density.*}
    awk -v S="$s" '{print S"\t"$1"\t"$2"\t"$3"\t"$4}' "$f" >> "$COMBINED_DENSITY"
  done
else
  echo "Combined density table already exists, skipping."
fi

###############################################################################
# Variant density plots (sample-wise, group-wise, heatmap)
###############################################################################

if command -v Rscript &>/dev/null; then
  PLOT_R="$MAP_DIR/plots/plot_variant_density.R"
  cat > "$PLOT_R" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: plot_variant_density.R <combined_density.tsv> <out_prefix> [sample_group_tsv]")
infile <- args[1]
out_prefix <- args[2]
sample_tsv <- if (length(args) >= 3) args[3] else NA

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(viridis)
})

dt <- fread(infile)
dt[, mid := (start + end)/2 / 1e6]

if (!is.na(sample_tsv) && file.exists(sample_tsv)) {
  st <- fread(sample_tsv)
  if (!all(c("sample", "group") %in% names(st))) {
    stop("sample_group_tsv must have columns: sample and group")
  }
  dt <- merge(dt, st[, .(sample, group)], by = "sample", all.x = TRUE)
} else {
  dt[, group := "UNGROUPED"]
}
dt[is.na(group), group := "UNGROUPED"]

p_sample <- ggplot(dt, aes(x = mid, y = count, colour = group, group = sample)) +
  geom_line(linewidth = 0.4, alpha = 0.9) +
  facet_wrap(~ chrom, scales = "free_x") +
  labs(
    x = "Genomic position (Mb)",
    y = "Variant count per window",
    colour = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(colour = "black", fill = "grey95"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(paste0(out_prefix, "_sample_variant_density.pdf"), p_sample, width = 8.5, height = 6)
ggsave(paste0(out_prefix, "_sample_variant_density.png"), p_sample, width = 8.5, height = 6, dpi = 300)

dt_group <- dt[, .(count = sum(count)), by = .(group, chrom, start, end)]
dt_group[, mid := (start + end)/2 / 1e6]
fwrite(dt_group, paste0(out_prefix, "_group_density.tsv"))

p_group <- ggplot(dt_group, aes(x = mid, y = count, colour = group, group = group)) +
  geom_line(linewidth = 0.6, alpha = 0.9) +
  facet_wrap(~ chrom, scales = "free_x") +
  labs(
    x = "Genomic position (Mb)",
    y = "Total variants per window (group)",
    colour = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(colour = "black", fill = "grey95"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(paste0(out_prefix, "_group_variant_density.pdf"), p_group, width = 8.5, height = 6)
ggsave(paste0(out_prefix, "_group_variant_density.png"), p_group, width = 8.5, height = 6, dpi = 300)

p_heat <- ggplot(dt_group, aes(x = mid, y = group, fill = count)) +
  geom_tile() +
  facet_wrap(~ chrom, scales = "free_x") +
  scale_fill_viridis(option = "magma") +
  labs(
    x = "Genomic position (Mb)",
    y = "Group",
    fill = "Variants / window"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(colour = "black", fill = "grey95"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(paste0(out_prefix, "_group_heatmap.pdf"), p_heat, width = 8.5, height = 6)
ggsave(paste0(out_prefix, "_group_heatmap.png"), p_heat, width = 8.5, height = 6, dpi = 300)
EOF

  echo "Creating publication-style variant density plots..."
  if [[ -n "$SAMPLES_FOR_MAPPING" && -f "$SAMPLES_FOR_MAPPING" ]]; then
    Rscript "$PLOT_R" "$COMBINED_DENSITY" "$MAP_DIR/plots/bacteria" "$SAMPLES_FOR_MAPPING" || echo "WARNING: density plotting failed; continuing."
  else
    Rscript "$PLOT_R" "$COMBINED_DENSITY" "$MAP_DIR/plots/bacteria" || echo "WARNING: density plotting failed; continuing."
  fi
else
  echo "Rscript not found; skipping variant density plotting."
fi

###############################################################################
# Phylogeny + SNP distance matrix
###############################################################################

PHYLO_DIR="$MAP_DIR/phylo"
SNPS_FASTA="$PHYLO_DIR/all_samples.snps.fa"
DIST_TSV="$PHYLO_DIR/snp_distance.tsv"
TREE_NWK="$PHYLO_DIR/snp_tree.nwk"

if [[ ! -f "$SNPS_FASTA" ]]; then
  echo "Building consensus sequences for each sample from merged VCF..."
  SAMPLES_IN_VCF=$(bcftools query -l "$MERGED_VCF")
  mkdir -p "$PHYLO_DIR"
  for S_RAW in $SAMPLES_IN_VCF; do
    S_BASENAME=$(basename "$S_RAW")
    S_LABEL=${S_BASENAME%.sorted.bam}
    S_LABEL=${S_LABEL%.bam}
    S_LABEL=${S_LABEL%.vcf.gz}
    outfa="$PHYLO_DIR/${S_LABEL}.fa"
    if [[ ! -f "$outfa" ]]; then
      echo "  - $S_RAW  (label: $S_LABEL)"
      bcftools consensus -s "$S_RAW" -f "$REFERENCE" "$MERGED_VCF" \
        | sed "s/^>.*/>${S_LABEL}/" > "$outfa"
    else
      echo "  - $S_LABEL consensus exists, skipping."
    fi
  done

  echo "Concatenating all consensus FASTAs..."
  cat "$PHYLO_DIR"/*.fa > "$PHYLO_DIR/all_samples.fa"

  echo "Extracting SNP-only alignment with snp-sites..."
  snp-sites -o "$SNPS_FASTA" "$PHYLO_DIR/all_samples.fa"
else
  echo "SNP alignment already exists: $SNPS_FASTA"
fi

if [[ ! -f "$DIST_TSV" ]]; then
  echo "Computing pairwise SNP distance matrix..."
  snp-dists "$SNPS_FASTA" > "$DIST_TSV"
else
  echo "Distance matrix already exists: $DIST_TSV"
fi

if [[ ! -f "$TREE_NWK" ]]; then
  echo "Building phylogenetic tree with FastTree..."
  FastTree -nt "$SNPS_FASTA" > "$TREE_NWK"
else
  echo "Phylogenetic tree already exists: $TREE_NWK"
fi

echo "Phylogeny + distance outputs:"
echo "  - SNP alignment : $SNPS_FASTA"
echo "  - SNP distances : $DIST_TSV"
echo "  - Tree (Newick) : $TREE_NWK"

###############################################################################
# SNP distance heatmap + tree plot (R)
###############################################################################

if command -v Rscript &>/dev/null; then
  PHYLO_PLOT_R="$MAP_DIR/plots/plot_phylo_heatmap.R"
  cat > "$PHYLO_PLOT_R" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: plot_phylo_heatmap.R <snp_distance.tsv> <tree.nwk>")
dist_file <- args[1]
tree_file <- args[2]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(viridis)
  library(ape)
})

out_prefix <- sub("\\.tsv$", "", dist_file)
dmat <- as.matrix(read.table(dist_file, header = TRUE, row.names = 1, check.names = FALSE))
dd <- melt(dmat)
colnames(dd) <- c("sample1", "sample2", "distance")

hc <- hclust(as.dist(dmat))
ord <- hc$labels[hc$order]
dd$sample1 <- factor(dd$sample1, levels = ord)
dd$sample2 <- factor(dd$sample2, levels = ord)

p_heat <- ggplot(dd, aes(x = sample1, y = sample2, fill = distance)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  labs(x = "", y = "", fill = "SNP distance") +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(paste0(out_prefix, "_heatmap.pdf"), p_heat, width = 6, height = 5)
ggsave(paste0(out_prefix, "_heatmap.png"), p_heat, width = 6, height = 5, dpi = 300)

tr <- read.tree(tree_file)

pdf(paste0(out_prefix, "_tree.pdf"), width = 7, height = 7)
plot(tr, cex = 0.5)
axis(1)
dev.off()

png(paste0(out_prefix, "_tree.png"), width = 7, height = 7, units = "in", res = 300)
plot(tr, cex = 0.5)
axis(1)
dev.off()
EOF

  echo "Creating SNP distance heatmap and tree plots..."
  Rscript "$PHYLO_PLOT_R" "$DIST_TSV" "$TREE_NWK" || echo "WARNING: phylo plotting failed; continuing."
else
  echo "Rscript not found; skipping SNP heatmap/tree plotting."
fi

###############################################################################
# Transition / transversion analysis (overall, per group, per sample)
# -> use bcftools stats on each per-sample VCF
###############################################################################

TITV_SAMPLE_TSV="$MAP_DIR/titv_per_sample.tsv"

echo "Computing transition/transversion statistics per sample using bcftools stats..."

# Header
echo -e "sample\tgroup\ttransitions\ttransversions\ttotal\tTiTv_ratio" > "$TITV_SAMPLE_TSV"

shopt -s nullglob
for vcf in "$MAP_DIR"/vcf/*.vcf.gz; do
  base=$(basename "$vcf")
  sample="${base%.vcf.gz}"

  # Skip combined / group VCFs
  if [[ "$sample" == "all_samples.merged" ]] || [[ "$sample" == *.group.merged ]]; then
    continue
  fi

  # Determine group
  group="UNGROUPED"
  if [[ -n "${SAMPLES_FOR_MAPPING:-}" && -f "$SAMPLES_FOR_MAPPING" ]]; then
    g=$(awk -v S="$sample" 'NR>1 && $1==S {print $2; exit}' "$SAMPLES_FOR_MAPPING")
    [[ -n "$g" ]] && group="$g"
  fi

  # Use bcftools stats TSTV lines:
  # Format is roughly:
  # TSTV  <chrom>  <type>  <ts>  <tv>  <ts/tv> ...
  # We want the line where type=="all" (or ALL) and read ts, tv.
  counts=$(
    bcftools stats "$vcf" 2>/dev/null \
      | awk '
          $1=="TSTV" && ($3=="all" || $3=="ALL") {
            ts=$4; tv=$5;
          }
          END {
            if (ts == "" || ts == "NA") ts = 0;
            if (tv == "" || tv == "NA") tv = 0;
            print ts "\t" tv;
          }
        ' \
      || echo -e "0\t0"
  )

  ti=$(echo "$counts" | cut -f1)
  tv=$(echo "$counts" | cut -f2)
  [[ -z "$ti" ]] && ti=0
  [[ -z "$tv" ]] && tv=0

  total=$((ti + tv))
  if [[ "$tv" -gt 0 ]]; then
    titv_ratio=$(awk -v ti="$ti" -v tv="$tv" 'BEGIN{print ti/tv}')
  else
    titv_ratio="NA"
  fi

  echo -e "${sample}\t${group}\t${ti}\t${tv}\t${total}\t${titv_ratio}" >> "$TITV_SAMPLE_TSV"
done
shopt -u nullglob

nrows=$(( $(wc -l < "$TITV_SAMPLE_TSV") - 1 ))
echo "Ti/Tv per-sample table written to: $TITV_SAMPLE_TSV  (samples: $nrows)"


###############################################################################
# Ti/Tv summary and plots (overall, per group, per sample)
###############################################################################

if command -v Rscript &>/dev/null && [[ -f "$TITV_SAMPLE_TSV" ]]; then
  TITV_R="$MAP_DIR/plots/titv_plots.R"

  cat > "$TITV_R" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: titv_plots.R <titv_per_sample.tsv> <out_prefix>")
}
infile <- args[1]
out_prefix <- args[2]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(viridis)
  library(reshape2)
})

dt <- fread(infile)
if (!all(c("sample","group","transitions","transversions","total","TiTv_ratio") %in% names(dt))) {
  stop("titv_per_sample.tsv must contain: sample, group, transitions, transversions, total, TiTv_ratio")
}

num_cols <- c("transitions","transversions","total","TiTv_ratio")
for (nc in num_cols) {
  dt[[nc]] <- suppressWarnings(as.numeric(dt[[nc]]))
}

dt[, group := factor(group)]

dt_overall <- dt[, .(
  transitions    = sum(transitions, na.rm = TRUE),
  transversions  = sum(transversions, na.rm = TRUE)
)]
dt_overall[, total := transitions + transversions]
dt_overall[, TiTv_ratio := ifelse(transversions > 0, transitions / transversions, NA_real_)]

dt_group <- dt[, .(
  transitions   = sum(transitions, na.rm = TRUE),
  transversions = sum(transversions, na.rm = TRUE)
), by = group]

dt_group[, total := transitions + transversions]
dt_group[, TiTv_ratio := ifelse(transversions > 0, transitions / transversions, NA_real_)]

fwrite(dt_overall, paste0(out_prefix, "_overall.tsv"), sep = "\t")
fwrite(dt_group,   paste0(out_prefix, "_per_group.tsv"), sep = "\t")
fwrite(dt,         paste0(out_prefix, "_per_sample.tsv"), sep = "\t")

dt_long_sample <- melt(
  dt,
  id.vars = c("sample","group"),
  measure.vars = c("transitions","transversions"),
  variable.name = "type",
  value.name = "count"
)

p_sample <- ggplot(dt_long_sample, aes(x = sample, y = count, fill = type)) +
  geom_col() +
  facet_wrap(~ group, scales = "free_x") +
  scale_fill_viridis(discrete = TRUE, option = "viridis") +
  labs(
    title = "Transitions / Transversions per sample",
    x = "Sample",
    y = "Count",
    fill = "Mutation type"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )

ggsave(paste0(out_prefix, "_per_sample_titv_counts.pdf"), p_sample, width = 10, height = 6)
ggsave(paste0(out_prefix, "_per_sample_titv_counts.png"), p_sample, width = 10, height = 6, dpi = 300)

dt_long_group <- melt(
  dt_group,
  id.vars = c("group"),
  measure.vars = c("transitions","transversions"),
  variable.name = "type",
  value.name = "count"
)

p_group_counts <- ggplot(dt_long_group, aes(x = group, y = count, fill = type)) +
  geom_col() +
  scale_fill_viridis(discrete = TRUE, option = "viridis") +
  labs(
    title = "Transitions / Transversions per group",
    x = "Group",
    y = "Count",
    fill = "Mutation type"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank()
  )

ggsave(paste0(out_prefix, "_per_group_titv_counts.pdf"), p_group_counts, width = 7, height = 4)
ggsave(paste0(out_prefix, "_per_group_titv_counts.png"), p_group_counts, width = 7, height = 4, dpi = 300)

p_group_ratio <- ggplot(dt_group, aes(x = group, y = TiTv_ratio)) +
  geom_col() +
  labs(
    title = "Ti/Tv ratio per group",
    x = "Group",
    y = "Ti/Tv ratio"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank()
  )

ggsave(paste0(out_prefix, "_per_group_titv_ratio.pdf"), p_group_ratio, width = 7, height = 4)
ggsave(paste0(out_prefix, "_per_group_titv_ratio.png"), p_group_ratio, width = 7, height = 4, dpi = 300)
EOF

  echo "Creating Ti/Tv plots and summary tables..."
  Rscript "$TITV_R" "$TITV_SAMPLE_TSV" "$MAP_DIR/plots/titv" \
    || echo "WARNING: Ti/Tv plotting failed; continuing."
else
  echo "Rscript not found or Ti/Tv table missing; skipping Ti/Tv plots."
fi

###############################################################################
# Mark run complete and summary
###############################################################################

touch "$MAP_DIR/RUN_COMPLETE"

echo "All done. Results in: $MAP_DIR"
echo "Key outputs:"
echo "  - BAMs           : $MAP_DIR/bam/"
echo "  - VCFs           : $MAP_DIR/vcf/"
echo "  - Merged VCF     : $MERGED_VCF"
echo "  - Variant counts : $VARCOUNT_SAMPLE, $VARCOUNT_GROUP"
echo "  - Ti/Tv table    : $TITV_SAMPLE_TSV"
echo "  - Plots          : $MAP_DIR/plots/"
echo "  - Phylo          : $PHYLO_DIR/"
echo "  - Command used   : $MAP_DIR/command.txt"
echo "  - Resume marker  : $MAP_DIR/RUN_COMPLETE"
