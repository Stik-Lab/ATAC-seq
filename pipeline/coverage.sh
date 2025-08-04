#!/bin/bash

#SBATCH --job-name=ATAC_coverage
#SBATCH --mem=50gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=atac_cov_%A-%a.log
#SBATCH --array=1-X

names=( X  )
describer=${names[$SLURM_ARRAY_TASK_ID-1]}

path_bam='path_to_output_files'
path_bw='path_to_bigwig_files'

module load deepTools/3.5.1-foss-2021b

echo "................................................................ 8. START_bamcoverage ${describer} ................................................................"

bamCoverage --bam ${path_bam}/${describer}_clean.bam --outFileName ${path_bw}/${describer}.bw --effectiveGenomeSize 2913022398 --outFileFormat ${path_bw} --binSize 1 --normalizeUsing RPGC > ${path_bw}/${describer}.log

echo "................................................................ 8. END_bamcoverage ${describer} ................................................................"
