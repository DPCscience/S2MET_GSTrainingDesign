#!/bin/bash

# #PBS -l walltime=36:00:00,mem=62gb,nodes=1:ppn=8
#PBS -l walltime=24:00:00,mem=22gb,nodes=1:ppn=8
#PBS -N distance_heritability
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/Scripts/Clustering

module load R/3.5.0

# # Cluster heritability calculations
Rscript distance_rank_heritability.R
