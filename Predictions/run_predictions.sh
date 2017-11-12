#!/bin/bash

#PBS -l walltime=24:00:00,mem=24gb,nodes=1:ppn=16
#PBS -N S2_MET_prediction_random_env_add
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET

module load R/3.4.0

# Additive GEBLUP predictions
#Rscript Predictions/MSI/aGEBLUP_predictions.R

# Interation GEBLUP predictions
#Rscript Predictions/MSI/iGEBLUP_predictions.R

# Clustering and prediction
Rscript Predictions/S2MET_cluster_prediction.R
