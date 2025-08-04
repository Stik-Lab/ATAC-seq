#!/bin/bash

#SBATCH --job-name=ATAC_alligment
#SBATCH --mem=50gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=atac_alignment_%A-%a.log
#SBATCH --array=1-X

names=( X )
describer=${names[$SLURM_ARRAY_TASK_ID-1]}

#MODULES REQUIRED

module load Bowtie2/2.4.4.1-GCC-11.2.0
module load SAMtools/1.13-foss-2021b

refgenome='path_to_refgenome'
path_fq='path_to_fastq_files'
path_bam='path_to_output_files'

echo " ................................................................ START_BOWTIE2 alligment ${describer} ................................................................"

bowtie2 --very-sensitive -x ${refgenome} --threads 8 -1 ${path_fq}/${describer}_1_val_1.fq.gz -2 ${path_fq}/${describer}_2_val_2.fq.gz -S ${path_bam}/${describer}.sam

echo " ................................................................ END_BOWTIE2 alligment ${describer} ................................................................"
