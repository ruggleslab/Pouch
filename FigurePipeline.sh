#!/bin/bash

#SBATCH --job-name=r_runner
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3GB
#SBATCH --partition=fn_medium
#SBATCH --time=24:00:00
#SBATCH --error=logs/error_e.txt
#SBATCH --output=logs/stdout_e.txt

module load anaconda3
conda activate /gpfs/data/ruggleslab/home/devlij03/sc_env
module unload anaconda3

#### Pre-processing

#for f in inputs/*; do
#   BN=`basename $f`
#   echo "Working on..."$BN
#   Rscript SCRIPTS/scrna-10x-seurat-3CD.R create OBJECTS/Samples/$BN $BN $f;
#done

#removed OBJECTS/Samples/EG \ OBJECTS/Samples/RR \

#Rscript scrna-10x-seurat-3CD.R integrate outputs/Pouch_UCI3 20 \
#    OBJECTS/Samples/mc \
#	OBJECTS/Samples/AM \
#	OBJECTS/Samples/AG \
#	OBJECTS/Samples/rk \
#	OBJECTS/Samples/rm \
#	OBJECTS/Samples/sh \
#	OBJECTS/Samples/am \
#	OBJECTS/Samples/lm \
#	OBJECTS/Samples/sb \
#	OBJECTS/Samples/kk \
#	OBJECTS/Samples/lb \
#	OBJECTS/Samples/CS-PB \
#	OBJECTS/Samples/KK \
#	OBJECTS/Samples/HC \
#	OBJECTS/Samples/GS-PB \
#	OBJECTS/Samples/NM-PB \
#	OBJECTS/Samples/AC-PI \
#	OBJECTS/Samples/CM-PI \
#	OBJECTS/Samples/KE \
#	OBJECTS/Samples/BK \
#	OBJECTS/Samples/JD \
#	OBJECTS/Samples/LB \
#	OBJECTS/Samples/KK2 \
#	OBJECTS/Samples/JD2 \
#	OBJECTS/Samples/SR-PB \
#	OBJECTS/Samples/TP-PB


##### Figure 1

#Rscript SCRIPTS/scrna-10x-seurat-3CD.R cluster OBJECTS 20
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R identify OBJECTS MajorPopulations
#Rscript SCRIPTS/Figure1b.R
#Rscript SCRIPTS/Figure1c.R
#Rscript SCRIPTS/Figure1d.R

#####

##### Subset T, B and Myeloid cells

#Rscript SCRIPTS/MajorSubset.R

#####

##### Figure 2

#Rscript SCRIPTS/scrna-10x-seurat-3CD.R cluster OBJECTS/Myeloid_cells 20
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R identify OBJECTS/Myeloid_cells MinorPopulations
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R de OBJECTS/Myeloid_cells MinorPopulations
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R identify OBJECTS/Myeloid_cells MinorPopulations2
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R de OBJECTS/Myeloid_cells MinorPopulations2

#Rscript SCRIPTS/Figure2ab.R
#Rscript SCRIPTS/Figure2c.R
#Rscript SCRIPTS/Figure2d.R
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R de OBJECTS/Myeloid_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/SOX4 Facs_pop
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R de OBJECTS/Myeloid_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/IL1B Facs_pop
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R de OBJECTS/Myeloid_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/APOE Facs_pop
#####

##### Figure 3

#Rscript SCRIPTS/scrna-10x-seurat-3CD.R cluster OBJECTS/T_cells 20
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R identify OBJECTS/T_cells MinorPopulations
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R de OBJECTS/T_cells MinorPopulations
#Rscript SCRIPTS/Figure3ab.R
#Rscript SCRIPTS/Figure3c_h.R
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R de OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/FOXP3 Facs_pop
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R de OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/CD8A Facs_pop

#####

##### Supplementary Figure ?

#Rscript SCRIPTS/scrna-10x-seurat-3CD.R cluster OBJECTS/B_cells 20
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R identify OBJECTS/B_cells MinorPopulations
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R de OBJECTS/B_cells MinorPopulations
#Rscript SCRIPTS/FigureSXa.R

#####

##### additional data

#Rscript SCRIPTS/atlas_integration.R
#Rscript SCRIPTS/epi_create.R

#Rscript SCRIPTS/scrna-10x-seurat-3CD_Smillie.R create OBJECTS/Smillie Smillie /gpfs/data/ruggleslab/home/devlij03/IBD/MANUSCRIPT_FILES/additional_data/Smillie_et_al/
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R cluster OBJECTS/Smillie 20

#####

#Rscript SCRIPTS/scrna-10x-seurat-3CD_Mitsialis.R create OBJECTS/Mitsialis Mitsialis /gpfs/data/ruggleslab/home/devlij03/IBD/MANUSCRIPT_FILES/additional_data/Mitsialis_et_al/
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R cluster OBJECTS/Mitsialis 20

#####

##### Figure 4
#Rscript SCRIPTS/Figure4a.R
#Rscript	SCRIPTS/Figure4b_d.R
#Rscript	SCRIPTS/Figure4e_g.R

#####

##### Figure 5
#Rscript SCRIPTS/Figure5a.R
#Rscript SCRIPTS/Figure5b.R

#Run cellPhoneDB single cell analysis (optional)
#Results stored in OBJECTS/Interactions
Rscript SCRIPTS/count_table_grab.R
sbatch SCRIPTS/cellphoneDB.sh


####atlas

#Rscript SCRIPTS/scrna-10x-seurat-3CD.R identify OBJECTS/Atlas/tcells Atlas_populations
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R identify OBJECTS/Atlas/bcells1 Atlas_populations
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R identify OBJECTS/Atlas/bcells2 Atlas_populations
#Rscript SCRIPTS/scrna-10x-seurat-3CD.R identify OBJECTS/Atlas/mcells Atlas_populations

