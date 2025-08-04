#!/bin/sh
#SBATCH --job-name=bedcov
#SBATCH --mem=40gb
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=ATAC_multucov_%A-%a.log
#SBATCH --array=1-1

names=( X )
describer=${names[$SLURM_ARRAY_TASK_ID-1]}

# ========== MODULES ==========
module load bedtools2-2.30.0-gcc-11.2.0-q7z4zez

path_bam='path_to_bam_files'
path_macs2='path_to_macs2'
path_output='path_to_outputdir'

# ========== Create all output directories ==========
for dir in "${path_output}" ; do
  if [ ! -d "${dir}" ]; then
    mkdir -p "${dir}"
  fi
done

# list all bam files 
bam_files=$(ls ${path_bam}/*_clean.bam | tr '\n' ' ')

# fist prepare bedfile with all merged peaks

cat *_peaks.narrowPeak | bedtools merge | bedtools sort > ${path_output}/${describer}pk.merged.sort.bed


bedtools multicov -bams  ${bam_files} \
     -bed ${path_macs2}/${describer}pk.merged.sort.bed > ${path_output}/${describer}_coverage.bed


# Create header
header="chr\tstart\tend"
for bam in ${bam_files}; do
    sample=$(basename $bam _clean.bam)
    header="${header}\t${sample}"
done

# add header to final output
echo -e "${header}" > ${path_output}/${describer}_coverage.bed
cat ${path_output}/${describer}_coverage.tmp >> ${path_output}/${describer}_coverage.bed
rm ${path_output}/${describer}_coverage.tmp