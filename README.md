# 16S rRNA Amplicon Analysis Pipeline

A reproducible pipeline for 16S rRNA amplicon sequencing analysis using QIIME 2 and DADA2, built for microbiome community profiling.

Connects to dissertation work on microbial community analysis at IISER Pune (16S rRNA/ITS amplification, primer design, qPCR).

## What This Pipeline Does
1. Import and demultiplex raw FASTQ reads
2. Quality filtering and denoising with DADA2
3. Taxonomic classification against SILVA database
4. Alpha and beta diversity analysis
5. Differential abundance testing

## Tools Used
- QIIME 2 — core amplicon analysis framework
- DADA2 — denoising and ASV generation
- SILVA — 16S taxonomic reference database
- R (phyloseq, vegan) — diversity analysis and visualization

## Repository Structure# 16S-rRNA-amplicon-pipeline
16S rRNA amplicon analysis pipeline using QIIME 2 and DADA2 for microbiome community profiling

16S-rRNA-amplicon-pipeline/

├── data/

│   ├── raw/          # Raw FASTQ files (not tracked by Git)

│   └── processed/    # Intermediate files

├── scripts/          # Analysis scripts

├── notebooks/        # Jupyter/R notebooks with results

├── results/

│   ├── taxonomy/     # Taxonomic classification outputs

│   └── diversity/    # Alpha and beta diversity outputs

└── docs/             # Methods documentation

## Status
🔧 In development — pipeline scripts being added progressively

## Author
Taruna Gupta | MSc Bioinformatics | Mumbai University
