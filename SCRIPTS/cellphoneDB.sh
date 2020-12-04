#!/bin/bash

#SBATCH --job-name=cellphoneDB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=30GB
#SBATCH --time=5:00:00
#SBATCH --error=logs/db_error.txt
#SBATCH --output=logs/db_stdout.txt

module load anaconda3
conda activate /gpfs/data/ruggleslab/home/devlij03/sc_env
module unload anaconda3

cellphonedb method statistical_analysis OBJECTS/Interactions/all_meta.txt \
	OBJECTS/Interactions/all_counts.txt \
	--threads=36 \
	--counts-data hgnc_symbol \
	--output-path OBJECTS/Interactions

cellphonedb plot dot_plot --means-path OBJECTS/Interactions/means.txt \
	--pvalues-path OBJECTS/Interactions/pvalues.txt \
	--output-path OBJECTS/Interactions

cellphonedb plot heatmap_plot OBJECTS/Interactions/all_meta.txt \
	--pvalues-path OBJECTS/Interactions/pvalues.txt \
	--output-path OBJECTS/Interactions
