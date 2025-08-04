#!/bin/bash

#SBATCH --job-name=ATAC_QC
#SBATCH --mem=50gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=atac_QC_%A-%a.log
#SBATCH --array=1-X

names=( X  )
describer=${names[$SLURM_ARRAY_TASK_ID-1]}

path_fq='path_to_fastq_files'

module load fastqc-0.11.9-gcc-11.2.0-dd2vd2m


echo    " ................................................................   START QC ${describer}     ................................................................"

fastqc ${path_fq}/${describer}.fastq.gz  -o ${path_fq}

echo    " ................................................................  END QC ${describer}     ................................................................" 

