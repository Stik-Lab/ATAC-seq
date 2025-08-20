# ATAC-seq Pipeline

This repository contains a complete ATAC-seq analysis pipeline designed for SLURM environments using Bash scripts.

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#pre-processing-step-by-step-description">Pre-processing Step-by-Step Description</a>
      <ul>
        <li><a href="#1-quality-control-fastqc">1. Quality Control (FastQC)</a></li>
        <li><a href="#2-adapter-trimming-trim-galore">2. Adapter Trimming (Trim Galore)</a></li>
        <li><a href="#3-alignment-bowtie2">3. Alignment (Bowtie2)</a></li>
        <li><a href="#4-sam-to-bam-conversion-and-sorting-samtools">4. SAM to BAM Conversion and Sorting (SAMtools)</a></li>
        <li><a href="#5-remove-mitochondrial-reads-chrm">5. Remove Mitochondrial Reads (chrM)</a></li>
        <li><a href="#6-filter-low-quality-reads">6. Filter Low-Quality Reads</a></li>
        <li><a href="#7-remove-duplicates-picard">7. Remove Duplicates (Picard)</a></li>
        <li><a href="#8-remove-blacklist-regions-bedtools">8. Remove Blacklist Regions (bedtools)</a></li>
        <li><a href="#9-generate-signal-tracks-deeptools">9. Generate Signal Tracks (deepTools)</a></li>
        <li><a href="#10-peak-calling-macs2">10. Peak Calling (MACS2)</a></li>
      </ul>
    </li>
    <li><a href="#differential-accessibility-analysis-description">Differential Accessibility Analysis Description</a></li>
    <li><a href="#motif-analysis-description">Motif Analysis Description</a></li>
  </ol>
</details>
   
<!-- Pre-procesing Step-by-Step Description -->
## Pre-procesing Step-by-Step Description

### How to Run

1. Modify the variables inside `pipeline.sh` with your paths and sample names.
2. Submit to SLURM:
   ```bash
   sbatch ATAC_PIPELINE.sh
   ```
### 1. Quality Control (FastQC)
Runs a quality control check on the raw FASTQ files.

Command:

```bash
fastqc sample.fastq.gz -o output_directory

```
Arguments:

- sample.fastq.gz: Input FASTQ file
- -o: Output directory



### 2. Adapter Trimming (Trim Galore)
Removes adapters and low-quality bases from paired-end FASTQ files.

Command:

```bash
trim_galore --fastqc --output_dir output_dir --paired sample_1.fastq.gz sample_2.fastq.gz

```
Arguments:

- --fastqc: Runs FastQC after trimming
- --output_dir: Output directory
- --paired: Paired-end read mode
- sample_1.fastq.gz, sample_2.fastq.gz: Input files

### 3. Alignment (Bowtie2)
Aligns trimmed reads to the reference genome.

Command:

```bash
bowtie2 --very-sensitive -x reference_index --threads 8 -1 sample_1.fq.gz -2 sample_2.fq.gz -S output.sam
```

Arguments:
- --very-sensitive: High-sensitivity alignment
- -x: Bowtie2 reference index
- -1, -2: Paired-end FASTQ files
- -S: Output SAM file

### 4. SAM to BAM Conversion and Sorting (SAMtools)
Converts the SAM file to a sorted BAM file and creates an index.

Command:

```bash
samtools view -bS input.sam | samtools sort -T temp_prefix -o output.bam
samtools index output.bam
```
Arguments:

- view: Converts SAM to BAM
- sort: Sorts BAM by coordinates
- index: Creates index for BAM

### 5. Remove Mitochondrial Reads (chrM)
Filters out mitochondrial reads from the BAM file.

Command:

```bash
samtools view -h input.bam | grep -v chrM | samtools sort -O bam -o output.bam -T temp_prefix

```

### 6. Filter Low-Quality Reads
Removes low-quality and non-primary/duplicate reads.

Command:

```bash
samtools view -F 2304 -b -q 10 input.bam > output.bam

```
Arguments:

- -F 2304: Filters out non-primary and duplicate reads
- -q 10: Minimum mapping quality

### 7. Remove Duplicates (Picard)
Marks and removes PCR duplicates.

Command:

```bash
java -jar picard.jar MarkDuplicates \
  I=input.bam \
  O=dedup.bam \
  M=metrics.txt \
  REMOVE_DUPLICATES=true
```

### 8. Remove Blacklist Regions (bedtools)
Filters out reads mapping to blacklisted genomic regions.

Command:

```bash
bedtools intersect -nonamecheck -v -abam input.bam -b blacklist.bed > clean.bam
samtools index clean.bam
```
### 9. Generate Signal Tracks (deepTools)
Creates a normalized bigWig file for visualization.

Command:

```bash
bamCoverage \
  --bam clean.bam \
  --outFileName output.bw \
  --effectiveGenomeSize genome_size \
  --outFileFormat bigwig \
  --binSize 1 \
  --normalizeUsing RPGC
```
Arguments:

- --bam: Input BAM file (cleaned, filtered, and deduplicated)
- --outFileName: Output bigWig file name
- --effectiveGenomeSize: Effective genome size (e.g. 2913022398 for human hg38)
- --outFileFormat: Output format, typically bigwig
- --binSize: Bin size for signal aggregation (e.g. 1 bp resolution)
- --normalizeUsing: Normalization method (e.g. RPGC = Reads Per Genomic Content)


### 10. Peak Calling (MACS2)
Identifies open chromatin regions (peaks) from the cleaned BAM.

Command:

```bash
macs2 callpeak \
  --format BAMPE \
  -t clean.bam \
  -g hs \
  -n sample_name \
  -B \
  -q 0.05 \
  --outdir peak_output
```

Arguments:

- --format BAMPE: Specifies paired-end BAM format
- -t: Input treatment BAM file (clean.bam)
- -g hs: Genome size (hs for human, or use effective genome size like 2.7e9) 
- -n: Sample name prefix for output files
- -B: Generates signal pileup files (*.bdg) 
- -q 0.05: FDR cutoff for peak detection (q-value) 
- --outdir: Output directory for peak files

## Differential Accessibility Analysis Description



## Motif Analysis Descirption




Notes
- Make sure all required modules are loaded or installed in your environment.
- The pipeline assumes paired-end sequencing data and a pre-built Bowtie2 genome index.
