# amrradar
One Health AMR Radar Pipeline
Version: 1.0.0

📋 Overview
The One Health AMR Radar Pipeline is a fully automated, end-to-end bioinformatics workflow designed for the comprehensive surveillance and analysis of Antimicrobial Resistance (AMR) from high-throughput sequencing data.

This pipeline integrates identification, extraction, assembly, annotation, pangenome analysis, and advanced statistical reporting into a single, seamless execution. It is optimized for multi-core systems, automatically detecting available resources to maximize throughput.

🚀 Key Features
Automated Resource Management: Auto-detects CPU cores and optimizes threading.

Taxonomic Filtering: Extracts specific organism reads (e.g., Pseudomonas aeruginosa) from metagenomic data using Kraken2.

Assembly & QC: Performs de novo assembly (Spades/Shovill) and quality assessment (CheckM, QUAST).

AMR Detection: Identifies resistance genes using AMRFinderPlus.

Advanced Analytics: Generates co-occurrence networks, heatmaps, PCA plots, and statistical association tests.

Comparative Genomics: Includes FastANI, Mashtree, and Pangenome analysis (Panaroo/Roary).

Variant Calling: Mapping-based SNP detection and phylogeny.

MLST: Automated Sequence Typing.

🛠️ Prerequisites & Installation
1. System Requirements
OS: Linux (Ubuntu/CentOS recommended) or macOS.

Shell: Bash (v4+).

Resources: Minimum 16GB RAM recommended (varies by dataset size).

2. Dependencies (Conda/Mamba)
It is highly recommended to run this pipeline within a Conda environment containing the following tools:

Core Tools:

kraken2

fastp / trimmomatic

spades / shovill

checkm-genome

quast

amrfinder

prokka

panaroo / roary

iqtree

fasttree

mash / mashtree

fastani

minimap2

samtools

bcftools

bedtools

snp-sites / snp-dists

mlst

seqkit

Python Libraries:

pandas, numpy, scipy, matplotlib, seaborn, plotly, scikit-learn, networkx, biopython

R Libraries:

ggplot2, dplyr, pheatmap (if using R-based plotting modules)

3. Installation
Clone this repository and ensure all scripts are executable:

Bash

git clone https://github.com/yourusername/one-health-amr-radar.git
cd one-health-amr-radar
chmod +x *.sh *.py
🏃 Usage
The pipeline is controlled by a single master script: run_complete_pipeline.sh.

1. Prepare Input Data
Raw Reads: Place your paired-end FASTQ files (e.g., Sample_R1.fastq.gz, Sample_R2.fastq.gz) in a directory (e.g., raw_data/).

Mapping File (CSV): Prepare a CSV file linking samples to metadata groups (Format: SampleName,Group).

2. Run the Pipeline
Execute the master script:

Bash

./run_complete_pipeline.sh
3. Configuration Prompts
The script will interactively ask for:

Input Directory: Path to your raw FASTQ folder.

Mapping File: Path to your sample_group.csv.

Kraken2 DB: Path to your Kraken2 database.

Target Organism: Scientific name to extract (e.g., Pseudomonas aeruginosa).

Parallel Samples: How many samples to process simultaneously (default: 4).

Reference Accession: NCBI Accession number for the reference genome (e.g., NC_002516.2) for mapping.

🧬 Pipeline Workflow
The pipeline is divided into three major phases:

🔹 Phase 1: Identification & Extraction
Kraken2 Classification: Taxonomically classifies raw reads.

Extraction: Extracts reads belonging to the Target Organism using extract_kraken2_by_name_auto_fastq.sh.

Staging: Moves specific reads to a staging directory (input/) for downstream analysis.

Visualization: Generates interactive HTML charts of organism abundance.

🔹 Phase 2: Assembly & AMR Analysis
Assembly (AMR Radar):

QC: Fastp/Trimmomatic.

Assembly: Spades or Shovill.

Genome QC: CheckM (completeness/contamination).

AMR Detection: Runs AMRFinderPlus on assembled contigs.

Statistics & Reporting:

Generates heatmaps of gene presence/absence.

Performs statistical tests (Chi-square/T-tests) on AMR burden between groups.

Builds Co-occurrence Networks (Gephi/Cytoscape compatible).

Analyzes MGE (Mobile Genetic Elements) enrichment.

🔹 Phase 3: Comparative Genomics & Pangenome
Annotation: Annotates assemblies using Prokka.

Comparative Metrics:

FastANI: Computes Average Nucleotide Identity (ANI) heatmap.

Mashtree: Constructs distance-based phylogenetic trees.

QUAST: Generates assembly quality metrics.

Pangenome Analysis:

Runs Panaroo (or Roary) to determine core vs. accessory genome.

Generates Core Gene Phylogeny (IQ-TREE/FastTree).

Mapping & SNP Calling:

Downloads reference genome (NCBI).

Maps reads using Minimap2.

Calls variants using Bcftools.

Generates SNP distance matrix and phylogeny.

Typing: Performs MLST (Multi-Locus Sequence Typing).

📂 Output Directory Structure
Results are stored in the results/ directory, organized by tool:

Plaintext

results/
├── sample.tsv               # Final sample sheet used
├── Charts/                  # Final HTML/PNG visualizations (AMR, Taxonomy)
├── summary_kraken2.xls      # Kraken2 summary
├── OneHealthAMR_summary.xlsx # Comprehensive AMR Excel Report
└── tools/
    ├── 0_kraken2/           # Kraken2 reports and outputs
    ├── 1_fastqc/            # FastQC reports
    ├── 2_trim/              # Trimmed reads
    ├── 3_assembly/          # Assembled Contigs (FASTA)
    ├── 5_checkm/            # Genome Quality Stats
    ├── *_AMRFinder/         # AMRFinderPlus results
    ├── *_AMRStatistics/     # Statistical analysis tables
    ├── *_AMRNetwork/        # Co-occurrence networks (GraphML)
    ├── *_comparativegenomics/ # QUAST, FastANI, Mashtree results
    ├── *_prokka/            # GFF/GBK Annotation files
    ├── *_pangenome/         # Panaroo/Roary output (gene_presence_absence.csv)
    ├── *_mapping/           # BAM files, VCFs, and SNP trees
    └── *_mlst/              # MLST typing results
📜 License
This project is licensed under the MIT License - see the LICENSE file for details.

🤝 Acknowledgments
Developed by Dr. Vishal Mevada for One Health AMR surveillance initiatives.

For issues or contributions, please open an issue in the GitHub repository.
