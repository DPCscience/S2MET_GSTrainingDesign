## S2MET source
## 
## A script that automatically loads the data relevant for the S2MET project

library(tidyverse)
library(readxl)

## Directories
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET/"

# Geno, pheno, and enviro data
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/"
pheno_dir <- file.path(proj_dir, "Phenotype_Data/")
env_var_dir <- file.path(proj_dir, "Environmental_Variables")


# Other directories
fig_dir <- file.path(proj_dir, "Figures/")
pred_dir <- file.path(proj_dir, "Predictions")
entry_dir <- file.path(proj_dir, "Plant_Materials")
analysis_dir <- file.path(proj_dir, "Analysis")
result_dir <- file.path(proj_dir, "Results")



# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))
# Load environmental data
load(file.path(env_var_dir, "environmental_data_compiled.RData"))

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))

# Load the trial metadata
trial_info <- read_csv(file.path(pheno_dir, "trial_metadata.csv"))


# Grab the entry names that are not checks
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)

vp <- entry_list %>% 
  filter(Class == "S2C1R") %>% 
  pull(Line)

# Find the tp and vp that are genotypes
tp_geno <- intersect(tp, row.names(s2_imputed_mat))
vp_geno <- intersect(vp, row.names(s2_imputed_mat))

# Define the checks
checks <- entry_list %>% 
  filter(Class == "Check") %>% 
  pull(Line)

entries <- entry_list %>% 
  pull(Line)

# Extract the tp and vp from the G matrix
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno, vp_geno),]

# Filter the S2 MET BLUEs for non-irrigated trials
S2_MET_BLUEs <- S2_MET_BLUEs %>% 
  filter(!grepl(pattern = "BZI|HTM", x = environment))

