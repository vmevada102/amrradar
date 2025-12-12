# One Health AMR Radar Pipeline

**Author:** Dr. Vishal Mevada  
**Repository:** https://github.com/vmevada102/amrradar.git  
**Version:** 1.0.0

## 📋 Overview

The **One Health AMR Radar Pipeline** is a fully automated, end-to-end bioinformatics workflow designed for the comprehensive surveillance and analysis of Antimicrobial Resistance (AMR) from high-throughput sequencing data.

This pipeline integrates identification, extraction, assembly, annotation, pangenome analysis, and advanced statistical reporting into a single, seamless execution. It is optimized for multi-core systems, automatically detecting available resources to maximize throughput.

### 🚀 Key Features
* **Automated Resource Management:** Auto-detects CPU cores and optimizes threading.
* **Taxonomic Filtering:** Extracts specific organism reads (e.g., *Pseudomonas aeruginosa*) from metagenomic data using **Kraken2**.
* **Assembly & QC:** Performs de novo assembly (**Spades/Shovill**) and quality assessment (**CheckM, QUAST**).
* **AMR Detection:** Identifies resistance genes using **AMRFinderPlus**.
* **Advanced Analytics:** Generates co-occurrence networks, heatmaps, PCA plots, and statistical association tests.
* **Comparative Genomics:** Includes **FastANI**, **Mashtree**, and **Pangenome analysis** (Panaroo/Roary).
* **Variant Calling:** Mapping-based SNP detection and phylogeny.
* **MLST:** Automated Sequence Typing.

---

## 🛠️ Prerequisites & Installation

### 1. System Requirements
* **OS:** Linux (Ubuntu/CentOS recommended) or macOS.
* **Shell:** Bash (v4+).
* **Resources:** Minimum 16GB RAM recommended (varies by dataset size).

### 2. Dependencies (Conda/Mamba)
It is highly recommended to run this pipeline within the provided Conda environment to ensure all tool versions are compatible.

**Setup Instructions:**

1.  Make sure you have [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/en/latest/) installed.
2.  Clone the repository:
    ```bash
    git clone https://github.com/vmevada102/amrradar.git
    cd amrradar
    ```
3.  Create the environment using the included `environment.yml` file:
    ```bash
    conda env create -f environment.yml
    ```
4.  Activate the environment:
    ```bash
    conda activate amr_radar
    ```
5.  Ensure all scripts are executable:
    ```bash
    chmod +x *.sh *.py
    ```

---

## 🏃 Usage

The pipeline is controlled by a single master script: `run_complete_pipeline.sh`.

### 1. Prepare Input Data
* **Raw Reads:** Place your paired-end FASTQ files (e.g., `Sample_R1.fastq.gz`, `Sample_R2.fastq.gz`) in a directory (e.g., `raw_data/`).
* **Mapping File (CSV):** Prepare a CSV file **(keep the name of file as sample_group.mapping.csv)** linking samples to metadata groups (Format: `SampleName,Group`).

### 2. Run the Pipeline
Ensure your conda environment is active, then execute the master script:

```bash
./run_complete_pipeline.sh

### 3. Configuration Prompts
The script will interactively ask for the following inputs:
1.  **Input Directory:** Path to your raw FASTQ folder.
2.  **Mapping File:** Path to your `sample_group_mapping.csv`.
3.  **Kraken2 DB:** Path to your Kraken2 database.
4.  **Target Organism:** Scientific name to extract (e.g., *Pseudomonas aeruginosa*).
5.  **Parallel Samples:** How many samples to process simultaneously (default: 4).
6.  **Reference Accession:** NCBI Accession number for the reference genome (e.g., `NC_002516.2`) used for SNP mapping.

---

## 🧬 Pipeline Workflow

The pipeline is divided into three major phases:

### 🔹 Phase 1: Identification & Extraction
* **Kraken2 Classification:** Taxonomically classifies raw reads.
* **Extraction:** Extracts reads belonging to the **Target Organism** using `extract_kraken2_by_name_auto_fastq.sh`.
* **Staging:** Moves specific reads to a staging directory (`input/`) for downstream analysis.
* **Visualization:** Generates interactive HTML charts of organism abundance.

### 🔹 Phase 2: Assembly & AMR Analysis
* **Assembly (AMR Radar):**
    * **QC:** Fastp/Trimmomatic.
    * **Assembly:** Spades or Shovill.
    * **Genome QC:** CheckM (completeness/contamination).
* **AMR Detection:** Runs **AMRFinderPlus** on assembled contigs.
* **Statistics & Reporting:**
    * Generates heatmaps of gene presence/absence.
    * Performs statistical tests (Chi-square/T-tests) on AMR burden between groups.
    * Builds Co-occurrence Networks (Gephi/Cytoscape compatible).
    * Analyzes MGE (Mobile Genetic Elements) enrichment.

### 🔹 Phase 3: Comparative Genomics & Pangenome
* **Annotation:** Annotates assemblies using **Prokka**.
* **Comparative Metrics:**
    * **FastANI:** Computes Average Nucleotide Identity (ANI) heatmap.
    * **Mashtree:** Constructs distance-based phylogenetic trees.
    * **QUAST:** Generates assembly quality metrics.
* **Pangenome Analysis:**
    * Runs **Panaroo** (or Roary) to determine core vs. accessory genome.
    * Generates **Core Gene Phylogeny** (IQ-TREE/FastTree).
* **Mapping & SNP Calling:**
    * Downloads reference genome (NCBI).
    * Maps reads using **Minimap2**.
    * Calls variants using **Bcftools**.
    * Generates SNP distance matrix and phylogeny.
* **Typing:** Performs **MLST** (Multi-Locus Sequence Typing).

---

## 📂 Output Directory Structure

Results are stored in the `results/` directory, organized by tool:

```text
results/
├── sample.tsv               # Final sample sheet used by the pipeline
├── Charts/                  # Final HTML/PNG visualizations (AMR, Taxonomy, Heatmaps)
├── summary_kraken2.xls      # Kraken2 summary report
├── OneHealthAMR_summary.xlsx # Comprehensive AMR Excel Report
└── tools/
    ├── 0_kraken2/           # Kraken2 reports and outputs
    ├── 1_fastqc/            # FastQC reports
    ├── 2_trim/              # Trimmed reads
    ├── 3_assembly/          # Assembled Contigs (FASTA)
    ├── 5_checkm/            # Genome Quality Stats
    ├── *_AMRFinder/         # AMRFinderPlus results
    ├── *_AMRStatistics/     # Statistical analysis tables and plots
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

