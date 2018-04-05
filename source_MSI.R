## S2MET source for MSI
## 
## A script that automatically loads the data relevant for the S2MET project


# Load packages
packages <- c("dplyr", "tidyr", "tibble", "stringr", "readxl", "readr", "parallel",
              "rrBLUP", "purrr")
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-unknown-linux-gnu-library/3.2/"
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"

geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Genos"
pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")


# Source the project functions
source(file.path(proj_dir, "source_functions.R"))


# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))

# Load an entry file
entry_list <- read_excel(file.path(data_dir, "project_entries.xlsx"))

# Load the trial metadata
trial_info <- read_csv(file.path(data_dir, "trial_metadata.csv"))

# Vector of relevant traits
traits <- c("GrainYield", "HeadingDate", "PlantHeight")


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

# Create a replacement vector for each trait that adds a space between words and
# also included the unit
trait_replace <- unique(S2_MET_BLUEs$trait) %>%
  setNames(object = c("Grain Yield", "Heading Date", "Plant Height"), nm = .)
trait_replace_unit <- unique(S2_MET_BLUEs$trait) %>%
  setNames(object = c("Grain Yield\n(kg ha^-1)", "Heading Date\n(days)", "Plant Height\n(cm)"), nm = .)

# Pull out the trials with irrigation
irrig_env <- subset(trial_info, irrigated == "yes", environment)$environment
# Filter the S2 MET BLUEs for non-irrigated trials
S2_MET_BLUEs <- S2_MET_BLUEs %>% 
  filter(!environment %in% irrig_env) %>%
  # Remove environments deemed failures (i.e. HNY16 for grain yield)
  filter(!(environment == "HNY16" & trait == "GrainYield"),
         !(environment == "BCW16" & trait == "PlantHeight"))



# Find environments in which just data on both the TP and VP is available
tp_vp_env <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) > 1, sum(line_name %in% vp_geno) > 1) %>% 
  distinct(environment) %>% 
  pull()

# Find environments in which just data on the TP is available
tp_only_env <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) > 1, sum(line_name %in% vp_geno) == 0) %>% 
  distinct(environment) %>% 
  pull()

# Find environments in which just data on the VP is available
vp_only_env <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) == 0, sum(line_name %in% vp_geno) > 1) %>% 
  distinct(environment) %>% 
  pull()

