## S2MET source for MSI
## 
## A script that automatically loads the data relevant for the S2MET project


# Load packages
packages <- c("dplyr", "tidyr", "tibble", "stringr", "readxl", "readr", "parallel",
              "rrBLUP", "purrr", "boot")
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


## Source the 'source.R' script from a define starting point
source_lines <- readLines(file.path(repo_dir, "source.R"))
source_lines_discard <- seq(which(grepl(pattern = "^# MSI", x = source_lines)))
source_lines_run <- source(textConnection(source_lines[-source_lines_discard]))

