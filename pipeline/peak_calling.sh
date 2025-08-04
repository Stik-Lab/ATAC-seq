#!/bin/bash

#SBATCH --job-name=ATAC_QC
#SBATCH --mem=50gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=atac_pkcalling_%A-%a.log
#SBATCH --array=1-X

names=( X  )
describer=${names[$SLURM_ARRAY_TASK_ID-1]}

path_macs2='path_macs2'
path_bam='path_to_output_files'

module load MACS2/2.2.5-foss-2021b-Python-3.8.5

echo "................................................................ 9. START_peak_calling ${describer} ................................................................"

macs2 callpeak --format BAMPE -t ${path_bam}/${describer}_clean.bam -g hs -n ${describer} -B -q 0.05 --outdir ${path_macs2}

echo "................................................................ 9. END_peak_calling ${describer} ................................................................"


