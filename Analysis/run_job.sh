#!/bin/bash

#PBS -l walltime=12:00:00,mem=24gb,nodes=1:ppn=1
#PBS -N S2_MET_heritability
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# This is a generic job-launching script.


# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET

module load R/3.4.0

# Heritability
Rscript Analysis/S2_MET_heritability_MSI.R
