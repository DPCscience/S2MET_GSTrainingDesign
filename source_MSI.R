## S2MET source for MSI
## 
## A script that automatically loads the data relevant for the S2MET project


## Change the library paths
if (version$minor == "5.2") .libPaths(gsub(pattern = "lab_library", replacement = "mesabi_library", x = .libPaths()))

# Load packages
packages <- c("sommer", "dplyr", "tidyr", "tibble", "stringr", "readxl", "readr", "parallel",
              "rrBLUP", "purrr", "boot", "pbr", "lme4", "modelr", "neyhart")


invisible(lapply(packages, library, character.only = TRUE))

## Directories
proj_dir <- "path/to/project/directory/on/supercomputer/"

geno_dir <-  "path/to/directory/on/supercomputer/containing/genotype/data"
pheno_dir <- "path/to/directory/on/supercomputer/containing/phenotype/data"

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")


## Source the 'source.R' script from a define starting point
source_lines <- readLines(file.path(repo_dir, "source.R"))
source_lines_discard <- seq(which(grepl(pattern = "^# MSI", x = source_lines)))
source_lines_run <- source(textConnection(source_lines[-source_lines_discard]))


