#!/bin/bash

#SBATCH --job-name=ATAC_PIPELINE
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=atac_pipeline_%A-%a.log
#SBATCH --array=1-X # change X for the actual number of samples before submiting

# ========== VARIABLES ==========
describer_list=(X)  # Replace X with actual sample names
describer=${describer_list[$SLURM_ARRAY_TASK_ID-1]}

#  Folder paths 
path_fq="path_to_fastq_files"
path_bam="path_to_bam_files"
path_temp="path_to_temp_files"
path_bw="path_to_bigwig_files"
path_macs2="path_to_macs2"

# Files and programs 
indexgenome='bowtie2/GRCh38_noalt_as/GRCh38_noalt_as' # Bowtie2 genome index
blacklist_file="path_to_blacklist"
effective_genome_size=2913022398

# Find R1 and R2
R1_RAW=$(ls ${path_fq}/${describer}*{_1,_*1}*.f*q.gz 2>/dev/null | head -n 1)
R2_RAW=$(ls ${path_fq}/${describer}*{_2,_*2}*.f*q.gz 2>/dev/null | head -n 1)

SAM=${path_bam}/${describer}.sam
BAM_CLEAN=${path_bam}/${describer}_clean.bam
BW=${path_bw}/${describer}.bw


for dir in "${path_bam}" "${path_temp}" "${path_bw}" "${path_macs2}"; do
  if [ ! -d "${dir}" ]; then
    mkdir -p "${dir}"
  fi
done


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

FQC_DONE=$(ls ${path_fq}/${describer}*_fastqc.zip 2>/dev/null | head -n 1)
if [ -n "${FQC_DONE}" ]; then
    echo "SKIP FastQC: reports already exist for ${describer}"
else
  fastqc ${path_fq}/${describer}*.fastq.gz -o ${path_fq}
fi

echo "................................................................ 1. END_FASTQC ${describer} ................................................................"

# ========== STEP 2: TRIM GALORE ==========

echo "................................................................ 2. START_TRIM_GALORE ${describer} ................................................................"
R1_TRIM=$(ls ${path_fq}/${describer}*val_1.f*q.gz 2>/dev/null | head -n 1)
R2_TRIM=$(ls ${path_fq}/${describer}*val_2.f*q.gz 2>/dev/null | head -n 1)
if [ -s "${R1_TRIM}" ] && [ -s "${R2_TRIM}" ]; then
    echo "SKIP Trim Galore: trimmed files already exist for ${describer}"
else
  trim_galore --fastqc --output_dir ${path_fq} --paired \
    "${R1_RAW}" "${R2_RAW}"
fi

echo "................................................................ 2. END_TRIM_GALORE ${describer} ................................................................"

# ========== STEP 3: ALIGNMENT ==========

echo "................................................................ 3. START_ALIGNMENT ${describer} ................................................................"

if [ -s "${SAM}" ]; then
    echo "SKIP Bowtie2: SAM already exists for ${describer}"
else
  bowtie2 --very-sensitive -x ${indexgenome} --threads 8 \
      -1 ${path_fq}/${describer}_*1_val_1.fq.gz \
      -2 ${path_fq}/${describer}_*2_val_2.fq.gz \
      -S ${path_bam}/${describer}.sam
fi

echo "................................................................ 3. END_ALIGNMENT ${describer} ................................................................"

# ========== STEP 4: SAM TO BAM ==========

echo "................................................................ 4. START_SAM_TO_BAM ${describer} ................................................................"

if [ -s "${path_bam}/${describer}.bam" ] && [ -s "${path_bam}/${describer}.bam.bai" ]; then
    echo "SKIP SAM to BAM: BAM already exists for ${describer}"
else
  samtools view -bS ${path_bam}/${describer}.sam | \
      samtools sort -T ${path_bam}/tmp_${describer} -o ${path_bam}/${describer}.bam -
  samtools index ${path_bam}/${describer}.bam
fi

echo "................................................................ 4. END_SAM_TO_BAM ${describer} ................................................................"

# ========== STEP 5: REMOVE chrM ==========

echo "................................................................ 5. START_REMOVE_chrM ${describer} ................................................................"

if [ -s "${path_temp}/${describer}.rmChrM.bam" ]; then
    echo "SKIP chrM removal: file already exists for ${describer}"
else
  samtools view -h ${path_bam}/${describer}.bam | grep -v chrM | \
      samtools sort -O bam -o ${path_temp}/${describer}.rmChrM.bam -T ${path_temp}
fi

echo "................................................................ 5. END_REMOVE_chrM ${describer} ................................................................"

# ========== STEP 6: FILTER LOW-QUALITY READS ==========

echo "................................................................ 6. START_FILTER_LOW_QUALITY ${describer} ................................................................"

if [ -s "${path_temp}/${describer}.qual.bam" ]; then
    echo "SKIP quality filter: file already exists for ${describer}"
else
  samtools view -F 2304 -b -q 10 ${path_temp}/${describer}.rmChrM.bam \
      > ${path_temp}/${describer}.qual.bam
fi

echo "................................................................ 6. END_FILTER_LOW_QUALITY ${describer} ................................................................"

# ========== STEP 7: REMOVE DUPLICATES ==========

echo "................................................................ 7. START_REMOVE_DUPLICATES ${describer} ................................................................"
if [ -s "${path_temp}/${describer}_removed_duplicates.bam" ]; then
    echo "SKIP MarkDuplicates: file already exists for ${describer}"
else
  java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=${path_temp}/${describer}.qual.bam \
      O=${path_temp}/${describer}_removed_duplicates.bam \
      M=${path_temp}/${describer}_marked_dup_metrics.txt \
      REMOVE_DUPLICATES=true ASSUME_SORTED=true VERBOSITY=WARNING
fi
echo "................................................................ 7. END_REMOVE_DUPLICATES ${describer} ................................................................"

# ========== STEP 8: REMOVE BLACKLIST REGIONS ==========

echo "................................................................ 8. START_REMOVE_BLACKLIST ${describer} ................................................................"

if [ -s "${BAM_CLEAN}" ] && [ -s "${BAM_CLEAN}.bai" ]; then
    echo "SKIP blacklist removal: file already exists for ${describer}"
else
  bedtools intersect -nonamecheck -v -abam ${path_temp}/${describer}_removed_duplicates.bam \
      -b ${blacklist_file} > ${path_bam}/${describer}_clean.bam
  samtools index ${path_bam}/${describer}_clean.bam
fi

echo "................................................................ 8. END_REMOVE_BLACKLIST ${describer} ................................................................"

# ========== STEP 9: BAM TO BIGWIG ==========

echo "................................................................ 9. START_BAMCOVERAGE ${describer} ................................................................"

if [ -s "${BW}" ]; then
    echo "SKIP bamCoverage: bigwig already exists for ${describer}"
else
  bamCoverage --bam ${path_bam}/${describer}_clean.bam \
      --outFileName ${path_bw}/${describer}.bw \
      --effectiveGenomeSize ${effective_genome_size} \
      --outFileFormat bigwig --binSize 10 --normalizeUsing RPGC \
      > ${path_bw}/${describer}.log
fi

echo "................................................................ 9. END_BAMCOVERAGE ${describer} ................................................................"

# ========== STEP 10: PEAK CALLING ==========

echo "................................................................ 10. START_PEAK_CALLING ${describer} ................................................................"
PEAK_DONE=$(ls ${path_macs2}/${describer}_peaks.narrowPeak 2>/dev/null | head -n 1)
if [ -n "${PEAK_DONE}" ]; then
    echo "SKIP MACS2: peaks already exist for ${describer}"
else
  macs2 callpeak --format BAMPE -t ${path_bam}/${describer}_clean.bam \
      -g ${effective_genome_size} -n ${describer} -B -q 0.05  --keep-dup all --outdir ${path_macs2}
fi
echo "................................................................ 10. END_PEAK_CALLING ${describer} ................................................................"

