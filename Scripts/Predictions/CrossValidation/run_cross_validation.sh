#!/bin/bash

# #PBS -l walltime=24:00:00,mem=62gb,nodes=1:ppn=24
#PBS -l walltime=24:00:00,mem=64gb,nodes=1:ppn=24
# #PBS -N cross-validation-CV12
#PBS -N cross-validation-CVzeroLOEO
# #PBS -N cross-validation-CVzeroFuture
# #PBS -N cross-validation-POV
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions/Scripts/Predictions/CrossValidation/

module load R/3.5.0
# module load R/3.5.2_mkl

# Cross validation
# Rscript cross_validation_CV12.R

Rscript cross_validation_CVzero_loeo.R

# Rscript cross_validation_CVzero_future.R

## Parent-offspring validation
# Rscript cross_validation_POV.R
