#!/bin/bash

#SBATCH --job-name=motif
#SBATCH --mem=60gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=motif_%A-%a.log
#SBATCH --array=1-3

names=( X Y Z )
describer=${names[$SLURM_ARRAY_TASK_ID-1]}

module load homer/0.1

path_bed='macs2'
path_motifs='motifs'

if [ ! -d "${path_motifs}" ]; then
    mkdir -p "${path_motifs}"
fi

# ========== RUN HOMER ==========
echo "................................................................ START_HOMER ${describer} ................................................................"

findMotifsGenome.pl ${path_bed}/${describer}.bed hg38 \
    ${path_motifs}/${describer} \
    -size 200 -p 4

echo "................................................................ END_HOMER ${describer} ................................................................"
