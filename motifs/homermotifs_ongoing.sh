#!/bin/bash

#SBATCH --job-name=MOTIF_analysis
#SBATCH --mem=50gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=homer_%A.log
#SBATCH --array=1-1

# ========== VARIABLES ==========
names=( X )
describer=${names[$SLURM_ARRAY_TASK_ID-1]}

input_bed='bam_atac/${describer}.bed'
genome='hg38'
output_dir='motifs_loops/${describer}'

# ========== MODULES ==========
module load homer/0.1

# ========== CHECKS & DIRECTORIES ==========
if [ ! -f "${input_bed}" ]; then
    echo "Error: Input BED file '${input_bed}' not found."
    exit 1
fi

if [ ! -d "${output_dir}" ]; then
    mkdir -p "${output_dir}"
fi

# ========== RUN HOMER ==========
echo ">> Running HOMER on $input_bed..."
findMotifsGenome.pl "${input_bed}" "${genome}" "${output_dir}" -size 200 -nomotif

echo ">> HOMER analysis completed for ${describer}."
