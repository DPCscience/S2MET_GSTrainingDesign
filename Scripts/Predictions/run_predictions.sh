#!/bin/bash

#### #PBS -l walltime=24:00:00,mem=62gb,nodes=1:ppn=24
#PBS -l walltime=08:00:00,mem=24gb,nodes=1:ppn=16
#### #PBS -N S2_MET_environmental_distance_pred
#PBS -N S2_MET_environmental_distance_realistic
#### #PBS -N S2_MET_environmental_distance_pred_window
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/Scripts/Predictions/

module load R/3.2.0_intel_mkl

# Predictions by environmental rank
# Rscript distance_rank_predictions.R

# Heritability by environmental rank
# Rscript distance_rank_heritability.R

# Predictions by environmental rank in a realistic breeding scenario
Rscript realistic_predictions.R
