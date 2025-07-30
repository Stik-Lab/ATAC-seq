# ATAC-seq Pipeline

This repository contains a complete ATAC-seq analysis pipeline designed for SLURM environments using Bash scripts.

## Contents
- Quality Control with FastQC
- Adapter Trimming with Trim Galore
- Alignment with Bowtie2
- Filtering, Deduplication, Blacklist Removal
- Peak Calling with MACS2
- Signal Tracks with deepTools

## How to Run

1. Modify the variables inside `pipeline.sh` with your paths and sample names.
2. Submit to SLURM:
   ```bash
   sbatch pipeline.sh
