#!/bin/bash

#SBATCH --job-name=trimming_atac
#SBATCH --mem=50gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=?
#SBATCH --output=atac_QC_%A-%a.log
#SBATCH --array=1-X

names=( X )
describer=${names[$SLURM_ARRAY_TASK_ID-1]}

module load Trim_Galore/0.6.6-foss-2021b-Python-3.8.5

path_fq='path_to_fastq_files'


echo "................................................................ START TRIM GALORE ${describer}  ................................................................"

trim_galore --fastqc --output_dir ${path_fq} --paired ${path_fq}/${describer}_1.fastq.gz ${path_fq}/${describer}_2.fastq.gz 

echo "................................................................ END TRIM GALORE ${describer}  ................................................................"

