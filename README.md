# amrradar

# 🧬 One Health AMR-Radar – Full Genome AMR Analysis Pipeline

This repository contains **AMR-Radar full automation** pipeline for:
✔ Kraken2 Taxonomy  
✔ AMR Detection & Statistics  
✔ Resistome & Comparative Genomics  
✔ Pangenome & Biomarkers Discovery  

---

## 🚀 Install & Download from GitHub

```bash
git clone https://github.com/vmevada102/amrradar.git


**## Create enviornment**

conda env create -f environment.yml
conda activate amr-radar

cd AMR-Radar
chmod +x run_full_amr_radar_pipeline.sh


🧬 AMR-Radar: Automated Whole-Genome AMR & Resistome Analysis Pipeline

One Health Approach | WGS-based AMR Surveillance | Comparative Genomics

🚀 Overview

AMR-Radar is a fully automated pipeline for bacterial whole-genome sequencing (WGS) analysis, integrating multiple tools to perform:

Functionality	Tools Used
Taxonomic Classification	Kraken2, Bracken
AMR Detection & Gene Extraction	AMRFinderPlus, Custom Scripts
Resistome Profiling	Python + R
Comparative Genomics	fastANI, mashtree, QUAST
Annotation	Prokka
Pangenome Analysis	Panaroo, Roary
Biomarker Discovery	Machine Learning (Python)
📂 Repository Structure
AMR-Radar/
 ├─ raw/                      → Input FASTQ files (user-provided)
 ├─ input/                    → Auto-filled after Kraken2 extraction
 ├─ results/                  → Output results and charts
 ├─ sample.tsv                → Sample metadata file (required)
 ├─ run_full_amr_radar_pipeline.sh   → MAIN pipeline script (interactive)
 ├─ environment.yml           → Conda environment setup
 ├─ Dockerfile                → Docker support
 └─ README.md                 → Main documentation

📌 Pipeline Flowchart
flowchart TD
    A[raw FASTQ files] --> B[Kraken2]
    B --> C[Extract Organism Reads]
    C --> D[Generate sample.tsv]
    D --> E[AMR Detection: AMRFinderPlus]
    E --> F[Resistome & Statistics]
    F --> G[Comparative Genomics]
    G --> H[Prokka Annotation]
    H --> I[Pangenome Analysis]
    I --> J[Biomarker Discovery (ML)]
    J --> K[Final Reports & Charts]

⚙️ Installation
🅰️ Conda (Recommended)
conda env create -f environment.yml
conda activate amr-radar


Check:

which kraken2
which prokka

🅱️ Docker (Fully Reproducible)
docker build -t amr-radar:latest .
docker run --rm -it -v /path/to/data:/data amr-radar:latest


Inside the container:

./run_full_amr_radar_pipeline.sh

🧾 Required Input File – sample.tsv
sample	forward_read	reverse_read	group

Generate using:

./make_samples_tsv.sh

▶️ Run the FULL Pipeline
chmod +x run_full_amr_radar_pipeline.sh
./run_full_amr_radar_pipeline.sh


During execution, interactive questions will appear:

Variable	Example Input
KRAKEN_DB	/data/db/kraken2/standard
THREADS	40
READS_DIR	raw
ORGANISM	Pseudomonas aeruginosa
PERCENT_OF	total

Results are saved in: /results/

📊 Output Overview
Output	Location
QC reports	results/summary_input_count.xls
AMRFinder gene results	results/OneHealthAMR_AMRFinder_summary.xlsx
HTML charts	results/charts/
Comparative Genomics	results/comparative/
Prokka Annotation	results/prokka/
Pangenome Analysis	results/pangenome/
Biomarker Analysis	results/biomarker/
