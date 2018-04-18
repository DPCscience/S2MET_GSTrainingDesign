#!/bin/bash

#PBS -l walltime=24:00:00,mem=24gb,nodes=1:ppn=16
#PBS -N S2_MET_model_based_clustering
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/Scripts/Clustering

module load R/3.2.0_intel_mkl

# # Cluster heritability calculations
# Rscript cluster_heritability.R

# Cluster lrt
Rscript model_based_clustering.R
