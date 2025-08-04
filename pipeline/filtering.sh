#!/bin/bash

#SBATCH --job-name=ATAC_analysis
#SBATCH --mem=50gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=atac_analysis_%A-%a.log
#SBATCH --array=1-X

names=( X )
describer=${names[$SLURM_ARRAY_TASK_ID-1]}


module load SAMtools/1.13-foss-2021b
module load picard/2.26.3-Java-11
module load bedtools2-2.30.0-gcc-11.2.0-q7z4zez

path_bam='path_to_output_files'
path_temp='path_to_temp_files'
blacklist_file='path_to_BLregions'
path_bw='path_to_bigwig_files'
path_macs2='path_macs2'

echo " ................................................................ 1. STAR_SAM -> BAM and index ${describer} ................................................................"

samtools view -bS ${path_bam}/${describer}.sam | samtools sort -T ${path_bam}/tmp_${describer} -o ${path_bam}/${describer}.bam -
samtools index ${path_bam}/${describer}.bam

echo " ................................................................ 1. END_SAM -> BAM and index ${describer} ................................................................"


echo "................................................................2. START_remove_mtDNA_READS ${describer} ................................................................"

samtools view -h ${path_bam}/${describer}.bam | grep -v chrM | samtools sort -O bam -o ${path_temp}/${describer}.rmChrM.bam -T ${path_temp}

echo "................................................................ 2. END_remove_mtDNA_READS ${describer} ................................................................"


echo "................................................................ 3. START_filtering_protperly_reads ${describer} ................................................................"

samtools view -F 2304 -b -q 10  ${path_temp}/${describer}.rmChrM.bam > ${path_temp}/${describer}.qual.bam

echo "................................................................ 3. END_filtering_protperly_reads ${describer} ................................................................"



echo "................................................................ 4. START_mark_duplicates ${describer} ................................................................"

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
         I=${path_temp}/${describer}.qual.bam \
         O=${path_temp}/${describer}_removed_duplicates.bam \
         M=${path_temp}/${describer}_marked_dup_metrics.txt \
         REMOVE_DUPLICATES=true ASSUME_SORTED=true VERBOSITY=WARNING

echo "................................................................ 4. END_mark_duplicates ${describer} ................................................................"

# remove ENCODE blacklist regions

echo "................................................................ 5. START_blacklist regions ${describer} ................................................................"

bedtools intersect -nonamecheck -v -abam ${path_temp}/${describer}_removed_duplicates.bam -b ${blacklist_file} > ${path_bam}/${describer}_clean.bam
 
echo "................................................................ 5. END_blacklist ${describer} ................................................................"

echo "................................................................ 6. START_index_duplicates ${describer} ................................................................"

# index the results

samtools index ${path_bam}/${describer}_clean.bam

echo "................................................................ 6. END_index_duplicates ${describer} ................................................................"


