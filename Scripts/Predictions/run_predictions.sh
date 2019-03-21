#!/bin/bash

#PBS -l walltime=24:00:00,mem=64gb,nodes=1:ppn=4
# #PBS -l walltime=24:00:00,mem=24gb,nodes=1:ppn=24
# #PBS -l walltime=24:00:00,mem=62gb,nodes=1:ppn=4
#PBS -N environmental_distance_pred
# #PBS -N environmental_covariance_predictions
# #PBS -N cluster_preds
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions/Scripts/Predictions/

module load R/3.5.0
# module load R/3.5.2_mkl

# Predictions by environmental rank
Rscript distance_rank_predictions.R

# # Predictions 
# Rscript environment_covariance_matrix_predictions.R

# Cluster predictions
# Rscript cluster_predictions.R

