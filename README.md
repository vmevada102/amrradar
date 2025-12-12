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
    git clone [https://github.com/vmevada102/amrradar.git](https://github.com/vmevada102/amrradar.git)
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
* **Mapping File (CSV):** Prepare a CSV file linking samples to metadata groups (Format: `SampleName,Group`).

### 2. Run the Pipeline
Ensure your conda environment is active, then execute the master script:

```bash
./run_complete_pipeline.sh
