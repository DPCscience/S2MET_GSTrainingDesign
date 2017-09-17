#!/bin/bash

#PBS -l walltime=12:00:00,mem=16gb,nodes=1:ppn=1
#PBS -N S2_MET_iGEBLUP_predictios
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
Rscript Predictions/MSI/S2_MET_cluster_prediction.R
