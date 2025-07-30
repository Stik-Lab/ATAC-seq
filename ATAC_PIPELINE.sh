#!/bin/bash

#SBATCH --job-name=ATAC_PIPELINE
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=atac_pipeline_%A-%a.log
#SBATCH --array=1-X

# ========== VARIABLES ==========
describer_list=(X)  # Replace X with actual sample names
describer=${describer_list[$SLURM_ARRAY_TASK_ID-1]}

path_fq="path_to_fastq_files"
path_bam="path_to_output_files"
path_temp="path_to_temp_files"
path_bw="path_to_bigwig_files"
path_macs2="path_to_macs2"
refgenome="path_to_refgenome"
blacklist_file="path_to_blacklist"
effective_genome_size=2913022398

mkdir -p ${path_fq} ${path_bam} ${path_temp} ${path_bw} ${path_macs2}

# ========== MODULES ==========
module load fastqc-0.11.9-gcc-11.2.0-dd2vd2m
module load Trim_Galore/0.6.6-foss-2021b-Python-3.8.5
module load Bowtie2/2.4.4.1-GCC-11.2.0
module load SAMtools/1.13-foss-2021b
module load picard/2.26.3-Java-11
module load MACS2/2.2.5-foss-2021b-Python-3.8.5
module load deepTools/3.5.1-foss-2021b
module load bedtools2-2.30.0-gcc-11.2.0-q7z4zez

# ========== STEP 1: FASTQC ==========

echo "................................................................ 1. START_FASTQC ${describer} ................................................................"

fastqc ${path_fq}/${describer}.fastq.gz -o ${path_fq}

echo "................................................................ 1. END_FASTQC ${describer} ................................................................"

# ========== STEP 2: TRIM GALORE ==========

echo "................................................................ 2. START_TRIM_GALORE ${describer} ................................................................"

trim_galore --fastqc --output_dir ${path_fq} --paired \
    ${path_fq}/${describer}_1.fastq.gz ${path_fq}/${describer}_2.fastq.gz

echo "................................................................ 2. END_TRIM_GALORE ${describer} ................................................................"

# ========== STEP 3: ALIGNMENT ==========

echo "................................................................ 3. START_ALIGNMENT ${describer} ................................................................"

bowtie2 --very-sensitive -x ${refgenome} --threads 8 \
    -1 ${path_fq}/${describer}_1_val_1.fq.gz \
    -2 ${path_fq}/${describer}_2_val_2.fq.gz \
    -S ${path_bam}/${describer}.sam

echo "................................................................ 3. END_ALIGNMENT ${describer} ................................................................"

# ========== STEP 4: SAM TO BAM ==========

echo "................................................................ 4. START_SAM_TO_BAM ${describer} ................................................................"

samtools view -bS ${path_bam}/${describer}.sam | \
    samtools sort -T ${path_bam}/tmp_${describer} -o ${path_bam}/${describer}.bam -
samtools index ${path_bam}/${describer}.bam

echo "................................................................ 4. END_SAM_TO_BAM ${describer} ................................................................"

# ========== STEP 5: REMOVE chrM ==========

echo "................................................................ 5. START_REMOVE_chrM ${describer} ................................................................"

samtools view -h ${path_bam}/${describer}.bam | grep -v chrM | \
    samtools sort -O bam -o ${path_temp}/${describer}.rmChrM.bam -T ${path_temp}

echo "................................................................ 5. END_REMOVE_chrM ${describer} ................................................................"

# ========== STEP 6: FILTER LOW-QUALITY READS ==========

echo "................................................................ 6. START_FILTER_LOW_QUALITY ${describer} ................................................................"

samtools view -F 2304 -b -q 10 ${path_temp}/${describer}.rmChrM.bam \
    > ${path_temp}/${describer}.qual.bam

echo "................................................................ 6. END_FILTER_LOW_QUALITY ${describer} ................................................................"

# ========== STEP 7: REMOVE DUPLICATES ==========

echo "................................................................ 7. START_REMOVE_DUPLICATES ${describer} ................................................................"

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=${path_temp}/${describer}.qual.bam \
    O=${path_temp}/${describer}_removed_duplicates.bam \
    M=${path_temp}/${describer}_marked_dup_metrics.txt \
    REMOVE_DUPLICATES=true ASSUME_SORTED=true VERBOSITY=WARNING

echo "................................................................ 7. END_REMOVE_DUPLICATES ${describer} ................................................................"

# ========== STEP 8: REMOVE BLACKLIST REGIONS ==========

echo "................................................................ 8. START_REMOVE_BLACKLIST ${describer} ................................................................"

bedtools intersect -nonamecheck -v -abam ${path_temp}/${describer}_removed_duplicates.bam \
    -b ${blacklist_file} > ${path_temp}/${describer}_clean.bam
samtools index ${path_temp}/${describer}_clean.bam

echo "................................................................ 8. END_REMOVE_BLACKLIST ${describer} ................................................................"

# ========== STEP 9: BAM TO BIGWIG ==========

echo "................................................................ 9. START_BAMCOVERAGE ${describer} ................................................................"

bamCoverage --bam ${path_temp}/${describer}_clean.bam \
    --outFileName ${path_bw}/${describer}.bw \
    --effectiveGenomeSize ${effective_genome_size} \
    --outFileFormat bigwig --binSize 1 --normalizeUsing RPGC \
    > ${path_bw}/${describer}.log

echo "................................................................ 9. END_BAMCOVERAGE ${describer} ................................................................"

# ========== STEP 10: PEAK CALLING ==========

echo "................................................................ 10. START_PEAK_CALLING ${describer} ................................................................"

macs2 callpeak --format BAMPE -t ${path_temp}/${describer}_clean.bam \
    -g hs -n ${describer} -B -q 0.05 --outdir ${path_macs2}

echo "................................................................ 10. END_PEAK_CALLING ${describer} ................................................................"

