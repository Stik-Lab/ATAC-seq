#!/bin/sh
#SBATCH --job-name=TOBIAS
#SBATCH --mem=50gb
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=tobias_%A-%a.log
#SBATCH --array=1-1

# ========== VARIABLES ==========
merged_peaks_file="mergedpks.bed"
path_bam='path_to_bam_files'
path_tobias='path_to_tobias_files'
jaspar_motifs='JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt'
refgenome='genome.fa'
cond1='cond1'
cond2='cond2'

# ========== MODULES ==========
module load TOBIAS/0.14.0-foss-2021b

# ========== CHECK FILES ==========
if [ ! -f "$merged_peaks_file" ]; then
    echo "Error: merged peaks file '$merged_peaks_file' not found."
    exit 1
fi

for cond in "$cond1" "$cond2"; do
    if [ ! -f "${path_bam}/${cond}_clean.bam" ]; then
        echo "Error: BAM file for '$cond' not found."
        exit 1
    fi
done

# ========== Create all output directories ==========
if [ ! -d "${path_tobias}" ]; then
  mkdir -p "${path_output}"
fi

# ========== ATACorrect + FootprintScores ==========
for cond in ${cond1} ${cond2}
do

	echo ">> Processing ${cond}..."

TOBIAS ATACorrect --bam ${{path_bam}/${cond}_clean.bam --genome ${refgenome} \
	--peaks ${merged_peaks_file} \
	--outdir ${path_tobias} --cores 8

TOBIAS FootprintScores --signal ${path_tobias}/${cond}_clean_corrected.bw \
	--regions ${merged_peaks_file} \
	--output ${path_tobias}/${cond}_footprints.bw \
	--cores 8

done

echo ">> Running BINDetect..."
TOBIAS BINDetect --motifs ${jaspar_motifs} \
	--signals ${path_tobias}/${cond1}_footprints.bw ${path_tobias}/${cond2}_footprints.bw \
	--genome ${refgenome} \
	--peaks ${merged_peaks_file}/ \
	--cond_names ${cond1} ${cond2} --cores 8 \
	--outdir ${path_tobias}


echo ">> Analysis complete."