## S2MET source
## 
## A script that automatically loads the data relevant for the S2MET project

library(tidyverse)
library(readxl)
library(rrBLUP)
library(neyhart)
library(boot)
library(pbr)
library(lmerTest)

## Directories
proj_dir <- repo_dir

# Geno, pheno, and enviro data
geno_dir <-  "path/to/directory/containing/genotype/data"
pheno_dir <-  "path/to/directory/containing/phenotypic/data"

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- pheno_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")



######
# MSI Source starts here
######


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

## Traits with units
traits_unit <- setNames(c("kg~ha^-1", "days", "cm"), traits)
traits_unit1 <- setNames(c("kg~ha^-1", "days", "AGDD", "cm"), unique(S2_MET_BLUEs$trait))


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


## Rank the environments according to heritability
env_herit_rank <- stage_one_data %>% 
  ungroup() %>% 
  filter(!str_detect(trial, "S2C1")) %>% 
  select(trait, environment, heritability) %>% 
  split(.$trait) %>%
  map(~arrange(., desc(heritability)) %>% 
        mutate(environment = factor(environment, levels = .$environment)))


# Pull out the trials with irrigation
irrig_env <- subset(trial_info, irrigated == "yes", environment)[[1]]

# Filter the S2 MET BLUEs for non-irrigated trials
S2_MET_BLUEs_temp <- S2_MET_BLUEs %>% 
  filter(!environment %in% irrig_env) %>%
  # Remove environments deemed failures (i.e. HNY16 for grain yield)
  filter(!(environment == "HNY16" & trait == "GrainYield"),
         !(environment == "BCW16" & trait == "PlantHeight"))

## Remove environments with low heritability
high_herit_environments <- env_herit_rank %>%
  map(mutate, environment = as.character(environment)) %>% 
  map(~filter(., heritability >= 0.10))

S2_MET_BLUEs_temp1 <- S2_MET_BLUEs_temp %>% 
  split(.$trait) %>%
  map2_df(.x = ., .y = high_herit_environments, ~filter(.x, environment %in% .y$environment))


## Reassign BLUEs
S2_MET_BLUEs <- S2_MET_BLUEs_temp1


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

## Split these vectors based on traits
tp_vp_env_trait <- S2_MET_BLUEs %>% 
  group_by(environment, trait) %>%
  filter(sum(line_name %in% tp_geno) > 1, sum(line_name %in% vp_geno) > 1) %>% 
  split(.$trait) %>% 
  map(~unique(.$environment))
  
tp_only_env_trait <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) > 1, sum(line_name %in% vp_geno) == 0) %>%
  split(.$trait) %>% 
  map(~unique(.$environment))
  
vp_only_env_trait <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) == 0, sum(line_name %in% vp_geno) > 1) %>% 
  split(.$trait) %>% 
  map(~unique(.$environment))





## Create two list of environments - for complete or realistic training and validation
complete_train_env <- tp_vp_env_trait

## Create a df that contains trait, set, train environmnents, and validation environments
train_val_environment_df <- tibble(trait = names(complete_train_env), trait_env = complete_train_env) %>%
  filter(trait %in% traits) %>%
  crossing(., set = c("complete", paste0("realistic", sort(unique(S2_MET_BLUEs$year))))) %>%
  # Filter out testing envs
  mutate(test_env = trait_env,
         pattern = str_extract(set, "[0-9]{2}$") %>% ifelse(is.na(.), ".*", .),
         test_env = map2(test_env, pattern, ~str_subset(string = .x, pattern = .y)),
         train_env = map2(trait_env, test_env, setdiff),
         train_env = map2(trait_env, train_env, ~if (length(.y) == 0) .x else .y)) %>%
  select(trait, set, train_env, test_env)
  

## Final filter of BLUEs
S2_MET_BLUEs <- filter(S2_MET_BLUEs, environment %in% tp_vp_env)


## Remove
rm(S2_MET_BLUEs_temp, S2_MET_BLUEs_temp1, s2_discrete_mat, s2_imputed_mat, s2_imputed_mat_use)


### Colors, abbreviations, etc.

## Significant level
alpha <- 0.05

## Create a factor for the distance methods
dist_method_replace <- c("pheno_dist" = "Phenotypic Distance", "pheno_loc_dist" = "Location Phenotypic Distance",  "great_circle_dist" = "Great Circle Distance", 
                         "OYEC_All" = "One Year All ECs", "OYEC_Mean" = "One Year Mean Cor EC", "OYEC_IPCA" = "One Year IPCA Cor EC", 
                         "MYEC_All" = "Multi Year All ECs", "MYEC_Mean" = "Multi Year Mean Cor EC", "MYEC_IPCA" = "Multi Year IPCA Cor EC",
                         "AMMI" = "AMMI", "sample" = "Random")
dist_method_abbr <- abbreviate(dist_method_replace)

## Alternative abbreviations
dist_method_abbr <- setNames(c("PD", "LocPD", "GCD", "1Yr-All-EC", "1Yr-Mean-EC", "1Yr-IPCA-EC", "All-EC", 
                               "Mean-EC", "IPCA-EC", "AMMI", "Random"), names(dist_method_replace))

## Color scheme for distance methods
## Warm colors for those that don't allow for new environments (AMMI and PD)
## Intermediate for LocPD (allows new years, but not new locations)
## Cool colors for those that do allow new environments (GCD, All-EC, Mean-EC, IPCA-EC)
## Grey for random
## Black for all data

## Create a cool colors palette
cool_colors <- rev(colorRampPalette(colors = c(umn_palette(4)[3], umn_palette(3)[3]))(4))

# colors <- wesanderson::wes_palette(name = "Darjeeling1", n = 9, type = "continuous")
# colors2 <- wesanderson::wes_palette(name = "Darjeeling2", n = 10, type = "continuous")
# colors_use <- c(colors[c(1:3, 5:7)], colors2[2:4], colors[4], "grey75")

dist_colors <- c("AMMI" = "#FFB71E", "PD" = umn_palette(3)[2], "LocPD" = umn_palette(3)[5], "GCD" = cool_colors[1],
                 "All-EC" = cool_colors[2], "Mean-EC" = cool_colors[3], "IPCA-EC" = cool_colors[4], 
                 "Random" = "grey75")
  
## Subset for use
dist_method_abbr_use <- dist_method_abbr[c("AMMI", "pheno_dist", "pheno_loc_dist", "great_circle_dist", 
                                           "MYEC_All", "MYEC_IPCA", "sample")]
dist_colors_use <- dist_colors[dist_method_abbr_use]



# ## Replacement vector for CV
cv_replace <- c("cv1", "pov1", "pocv1",  "cv2" , "pocv2", "cv0", "pov0", "pocv0", "cv00", "pov00", "pocv00") %>%
  setNames(object = toupper(.), .)

## Set replacement vector and function
set_replace <- c("complete" = "Leave-one-environment-out", "realistic" = "Leave-one-year-out")

f_set_replace <- function(x, abbr = FALSE) {
  x1 <- str_replace_all(string = x, set_replace)
  # Abbreviate?
  if (abbr) {
    x2 <- str_replace_all(string = x1, pattern = "-", replacement = " ") %>% abbreviate(4) %>% toupper()
  } else {
    x2 <- x1
  }
  
  str_replace(string = x2, pattern = "([A-Za-z-]*)([0-9]{4})", replacement = "\\1 (\\2)")

}
  
  
  
  
  
  
  
  
