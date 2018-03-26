## S2MET source
## 
## A script that automatically loads the data relevant for the S2MET project

library(tidyverse)
library(readxl)
library(rrBLUP)

## Directories
proj_dir <- repo_dir

# Geno, pheno, and enviro data
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/"

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")



# Load the phenotypic data
load(file.path(data_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))
# Load environmental data
load(file.path(data_dir, "environmental_data_compiled.RData"))

# Load an entry file
entry_list <- read_excel(file.path(data_dir, "project_entries.xlsx"))

# Load the trial metadata
trial_info <- read_csv(file.path(data_dir, "trial_metadata.csv"))


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
  filter(Class != "Check") %>% 
  pull(Line)

# Extract the tp and vp from the G matrix
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno, vp_geno),]

# Calculate the K matrix
K <- A.mat(X = s2_imputed_mat_use, min.MAF = 0, max.missing = 1)

# Filter the S2 MET BLUEs for non-irrigated trials
S2_MET_BLUEs <- S2_MET_BLUEs %>% 
  filter(!grepl(pattern = "BZI|HTM|AID", x = environment))


# Create a replacement vector for each trait that adds a space between words and
# also included the unit
trait_replace <- unique(S2_MET_BLUEs$trait) %>%
  set_names(x = c("Grain Yield", "Heading Date", "Plant Height"), nm = .)
trait_replace_unit <- unique(S2_MET_BLUEs$trait) %>%
  set_names(x = c("Grain Yield\n(kg ha^-1)", "Heading Date\n(days)", "Plant Height\n(cm)"), nm = .)
