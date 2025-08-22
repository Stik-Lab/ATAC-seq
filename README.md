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
    <li>
      <a href="#motif-analysis-description">Motif Analysis Description</a>
      <ul>
        <li><a href="#how-to-run">How to Run</a></li>
        <li><a href="#step-by-step-description">Step-by-Step Description</a></li>
      </ul>
    </li>
    <li>
      <a href="#differential-accessibility-analysis-description">Differential Accessibility Analysis Description</a>
      <ul>
        <li><a href="#how-to-run">How to Run</a></li>
        <li><a href="#1-peak-quantification-multicovsh">1. Peak quantification (multicov.sh)</a></li>
        <li><a href="#2-differential-analysis-and-visualization-atac_diffanalysisrmd">2. Differential analysis and visualization (ATAC_diffanalysis.Rmd)</a></li>
      </ul>
    </li>
  </ol>
</details>


   
<!-- Pre-procesing Step-by-Step Description -->
## Pre-processing Step-by-Step Description

### Input files

- Raw paired-end FASTQ files (sample_1.fastq.gz, sample_2.fastq.gz)
- Reference genome index for Bowtie2.
- Blacklist BED file (genomic regions to exclude)

### How to Run

1. Modify the variables inside `ATAC_PIPELINE.sh` with your paths and sample names.
2. Submit to SLURM:
```bash
sbatch ATAC_PIPELINE.sh
```
### 1. Quality Control (FastQC)
Runs a quality control check on the raw FASTQ files.

#### Command line:

```bash
fastqc sample.fastq.gz -o output_directory
```
#### Arguments:

- sample.fastq.gz: Input FASTQ file
- -o: Output directory



### 2. Adapter Trimming (Trim Galore)
Removes adapters and low-quality bases from paired-end FASTQ files.

#### Command line

```bash
trim_galore --fastqc --output_dir output_dir --paired sample_1.fastq.gz sample_2.fastq.gz
```
#### Arguments:

- --fastqc: Runs FastQC after trimming
- --output_dir: Output directory
- --paired: Paired-end read mode
- sample_1.fastq.gz, sample_2.fastq.gz: Input files

### 3. Alignment (Bowtie2)
Aligns trimmed reads to the reference genome.

#### Command line:

```bash
bowtie2 --very-sensitive -x reference_index --threads 8 -1 sample_1.fq.gz -2 sample_2.fq.gz -S output.sam
```

#### Arguments:
- --very-sensitive: High-sensitivity alignment
- -x: Bowtie2 reference index
- -1, -2: Paired-end FASTQ files
- -S: Output SAM file

### 4. SAM to BAM Conversion and Sorting (SAMtools)
Converts the SAM file to a sorted BAM file and creates an index.

#### Command line:

```bash
samtools view -bS input.sam | samtools sort -T temp_prefix -o output.bam
samtools index output.bam
```
#### Arguments:

- view: Converts SAM to BAM
- sort: Sorts BAM by coordinates
- index: Creates index for BAM

### 5. Remove Mitochondrial Reads (chrM)
Filters out mitochondrial reads from the BAM file.

#### Command:

```bash
samtools view -h input.bam | grep -v chrM | samtools sort -O bam -o output.bam -T temp_prefix

```

### 6. Filter Low-Quality Reads
Removes low-quality and non-primary/duplicate reads.

#### Command line:

```bash
samtools view -F 2304 -b -q 10 input.bam > output.bam

```
#### Arguments:

- -F 2304: Filters out non-primary and duplicate reads
- -q 10: Minimum mapping quality

### 7. Remove Duplicates (Picard)
Marks and removes PCR duplicates.

#### Command line:

```bash
java -jar picard.jar MarkDuplicates \
  I=input.bam \
  O=dedup.bam \
  M=metrics.txt \
  REMOVE_DUPLICATES=true
```
#### Arguments:
- I → Input BAM file (input.bam) 
- O → Output BAM file after duplicate removal (dedup.bam)
- M → Metrics file reporting duplication statistics (metrics.txt)
- REMOVE_DUPLICATES=true → Removes duplicates instead of just marking them


### 8. Remove Blacklist Regions (bedtools)
Filters out reads mapping to blacklisted genomic regions.

#### Command line:

```bash
bedtools intersect -nonamecheck -v -abam input.bam -b blacklist.bed > clean.bam
samtools index clean.bam
```
#### Arguments (bedtools intersect):
- -nonamecheck → Ignores mismatches between sequence names in the BAM and BED files
- -v → Reports entries in input.bam that have no overlap with blacklist.bed
- -abam input.bam → Input BAM file
- -b blacklist.bed → BED file of regions to exclude


### 9. Generate Signal Tracks (deepTools)
Creates a normalized bigWig file for visualization.

#### Command line:

```bash
bamCoverage \
  --bam clean.bam \
  --outFileName output.bw \
  --effectiveGenomeSize genome_size \
  --outFileFormat bigwig \
  --binSize 1 \
  --normalizeUsing RPGC
```
#### Arguments:

- --bam: Input BAM file (cleaned, filtered, and deduplicated)
- --outFileName: Output bigWig file name
- --effectiveGenomeSize: Effective genome size (e.g. 2913022398 for human hg38)
- --outFileFormat: Output format, typically bigwig
- --binSize: Bin size for signal aggregation (e.g. 1 bp resolution)
- --normalizeUsing: Normalization method (e.g. RPGC = Reads Per Genomic Content)


### 10. Peak Calling (MACS2)
Identifies open chromatin regions (peaks) from the cleaned BAM.

#### Command line:

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

#### Arguments:

- --format BAMPE: Specifies paired-end BAM format
- -t: Input treatment BAM file (clean.bam)
- -g hs: Genome size (hs for human, or use effective genome size like 2.7e9) 
- -n: Sample name prefix for output files
- -B: Generates signal pileup files (*.bdg) 
- -q 0.05: FDR cutoff for peak detection (q-value) 
- --outdir: Output directory for peak files

### Output files

- FastQC quality reports (.html, .zip)
- Trimmed FASTQ files
- Alignment files (.sam, .bam, .bai)
- Filtered and deduplicated BAM files (clean.bam, dedup.bam)
- Duplicate metrics (metrics.txt)
- Normalized signal tracks (.bw)
- Peak calling results from MACS2 (.narrowPeak, .xls, .bdg, etc.)

## Motif Analysis Description

This pipeline identifies transcription factor motifs using that show differential binding footprints between two ATAC-seq conditions, providing insight into regulatory drivers of chromatin accessibility changes. This pipeline uses TOBIAS to performing this analysis, taking BAM files ans a merged set of accessible peaks.

### Input files

- Merged peaks file: mergedpks.bed
  - Contains the consensus accessible regions across all samples.
    
*Command line to generate a merged peak file.*

```bash
cat *_peaks.narrowPeak | bedtools merge | bedtools sort > ${describer}pk.merged.sort.bed
```

- BAM files: {cond1}_clean.bam and {cond2}_clean.bam
  - Cleaned, deduplicated, and properly indexed ATAC-seq alignments for each condition.

- Reference genome: genome.fa
  - FASTA file of the reference genome (must match BAM alignment).

- Motif database: JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt
  - Position Frequency Matrices (PFMs) used for motif scanning.



### How to Run

1. Modify the variables inside `tobias_motif_analysis.sh` with your paths and sample names.
2. Submit to SLURM:

```bash
sbatch tobias_motif_analysis.sh
```
### Step-by-Step Description

#### 1. ATAC correct
Corrects ATAC-seq signal for Tn5 insertion bias.

```bash
TOBIAS ATACorrect --bam condX_clean.bam \
                  --genome genome.fa \
                  --peaks mergedpks.bed \
                  --outdir path_tobias \
                  --cores 8
```

#### 2. FootprintScores
Computes footprint scores across the merged peak set.

```bash
TOBIAS FootprintScores --signal *_corrected.bw \
                       --regions mergedpks.bed \
                       --output *_footprints.bw \
                       --cores 8
```

#### 3. BINDetect

Performs motif binding site detection using footprint scores and tests for differential TF binding between **cond1** and **cond2**.

```bash
TOBIAS BINDetect --motifs JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
                 --signals cond1_footprints.bw cond2_footprints.bw \
                 --genome genome.fa \
                 --peaks mergedpks.bed \
                 --cond_names cond1 cond2 \
                 --cores 8 \
                 --outdir path_tobias
```

## Differential Accessibility Analysis Description

Differential accessibility analysis is the process of identifying genomic regions (peaks) that show significant differences in chromatin accessibility between two or more biological conditions. In ATAC-seq or other chromatin profiling experiments, “accessibility” reflects how open or closed a region of chromatin is, which can influence transcription factor binding and gene regulation.

By comparing read counts between conditions, differentially accessible peaks can be detected, showing regions that are more open in condition 1 (**cond1**) or in condition 2 (**cond2**).

This analysis is divided into two main steps: **(1) peak quantification across samples** and **(2) differential accessibility analysis with visualization and functional annotation**.

### Input files

- Cleaned BAM files.
- MACS2 peak files.

### 1. Peak quantification (multicov.sh)

#### How to run

**1.** Edit the variables inside multicov.sh in the differential accessibility analysis folder to specify your paths and sample names.
**2.** Submit to SLURM:
   
```bash
sbatch multicov.sh
```
#### Step-by-step description

**1.** All peak files (*_peaks.narrowPeak) are merged and sorted into a consensus peak set.

##### Command line
```bash
cat ${path_macs2}/*_peaks.narrowPeak | \
  bedtools merge | \
  bedtools sort > ${path_output}/${describer}pk.merged.sort.bed
```

**2.** For each BAM file (*_clean.bam), read counts are quantified across the merged peaks using bedtools multicov.
   
##### Command line
```bash
bedtools multicov \
  -bams ${path_bam}/*_clean.bam \
  -bed ${path_output}/${describer}pk.merged.sort.bed \
  > ${path_output}/${describer}_coverage.tmp
```
  
- The output file serves as the input for downstream differential accessibility analysis in R.

### 2. Differential analysis and visualization (ATAC_diffanalysis.Rmd)

#### How to Run

**1.** Open RStudio and load the script ATAC_diffanalysis.Rmd.
**2.** Edit the input file name (replace "file.bed" with your multicov count matrix).
**3.** Adjust the sample conditions (e.g., X vs. Y) in the metadata section of the script.
**4.** Run the script step by step (chunk by chunk) in RStudio.
**5.** At the end, all results and plots will be saved automatically in the working directory.


#### Step-by-step description

- **Data import and DESeq2 setup** : The multicov output is read into R, unique peak IDs are generated, and a DESeq2 dataset is constructed with condition metadata.
- **Normalization and PCA** :  Variance stabilizing transformation (VST) is applied, followed by PCA to assess clustering between conditions.
- **Differential accessibility analysis**: Using DESeq2, log2 fold changes, p-values, and adjusted FDR values are computed. Results are merged with raw and normalized counts.  
- **Visualization**:
  - Volcano plots are generated to highlight significantly up- or down-regulated peaks.
  - Boxplots compare accessibility changes between conditions.
- **Peak annotation**: Significant UP and DOWN peaks are annotated relative to nearby genes using ChIPseeker.
- **Functional analysis**: Annotated gene lists are tested for Gene Ontology (GO) term enrichment, highlighting biological processes associated with accessibility changes.

  
### Output files

- Merged peak set (*.bed)
- Normalized counts and DESeq2 results (resultsX.tab)
- PCA plot (pca_X.pdf)
- Volcano plot (volcano_X.pdf)
- Peak annotation files (UP_peaks.bed, DOWN_peaks.bed)
- GO enrichment plots


#### Notes
- Make sure all required modules are loaded or installed in your environment.
- The pipeline assumes paired-end sequencing data and a pre-built Bowtie2 genome index.
