#!/bin/bash

#PBS -l walltime=24:00:00,mem=64gb,nodes=1:ppn=8
# #PBS -l walltime=24:00:00,mem=24gb,nodes=1:ppn=24
# #PBS -l walltime=24:00:00,mem=62gb,nodes=1:ppn=4
#PBS -N environmental_distance_pred
# # #PBS -N environmental_covariance_predictions
# #PBS -N cluster_preds
#PBS -M username@emaildomain.com
#PBS -m abe
#PBS -r n

# Change the working directory
cd /path/to/directory/containing/folder/Scripts

module load R/3.5.0

# Predictions by environmental rank
Rscript distance_rank_predictions.R


# Cluster predictions
# Rscript all_data_cluster_predictions.R

