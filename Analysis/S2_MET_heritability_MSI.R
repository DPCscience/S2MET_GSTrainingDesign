## Script to calculate heritability across all environments
## 
## 
## Script was designed to run on MSI

# List of packages
packages <- c("tidyverse", "stringr", "readxl", "pbr", "regress", "rrBLUP")

# Set the directory of the R packages
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/" 

fig_dir <- file.path(proj_dir, "Figures/")
pred_dir <- file.path(proj_dir, "Predictions")
geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/GBS_Genos"
env_var_dir <- file.path(proj_dir, "Environmental_Variables/")
pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"
entry_dir <- file.path(proj_dir, "Plant_Materials")

# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))


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



### Calculate the heritability across all environments
# Group by trait
stage_two <- S2_MET_BLUEs %>%
  mutate(gt = interaction(line_name, environment)) %>%
  mutate_if(is.character, as.factor) %>%
  group_by(trait) %>%
  do(mod = {
    
    # Create an object for the data.frame
    df <- .
    
    # Extract the standard errors
    R <- solve(diag(df$std_error^2))
    R <- structure(R, dimnames = replicate(2, str_c(df$line_name, df$environment, sep = ":"), 
                                           simplify = FALSE))
    
    # Table of lines by environments (i.e. plots)
    plot_table <- xtabs(formula = ~ line_name + environment, data = df)
    
    # Find the number of environments
    n_e <- plot_table %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Now replicates
    n_r <- plot_table %>% 
      rowSums() %>%
      harm_mean()
    
    # Use regress function
    fit <- regress(formula = value ~ environment,
                   Vformula = ~ line_name + gt + R, 
                   data = df, identity = FALSE, pos = rep(TRUE, 3))
    
    # Reorganize the BLUPs
    blup_out <- BLUP(model = fit)
    
    BLUPs <- blup_out$Mean %>% 
      data_frame(term = names(.), value = .) %>% 
      separate(col = term, into = c("term", "line_name", "environment"), 
               sep = "\\.", fill = "right") %>%
      select(-term)
    
    data_frame(mod = list(fit), BLUPs = list(BLUPs), BLUP_Covar = list(blup_out$Covariance), 
               n_e = n_e, n_r = n_r)  })


## Now fit models without the GxE term in order to conduct a LRT

reduced_mod <- S2_MET_BLUEs %>%
  mutate_if(is.character, as.factor) %>%
  group_by(trait) %>%
  do(mod = {
    
    # Create an object for the data.frame
    df <- .
    
    # Extract the standard errors
    R <- solve(diag(df$std_error^2))
    R <- structure(R, dimnames = replicate(2, str_c(df$line_name, df$environment, sep = ":"), 
                                           simplify = FALSE))
    
    # Use regress function
    fit <- regress(formula = value ~ environment,
                   Vformula = ~ line_name + R, 
                   data = df, identity = FALSE, pos = rep(TRUE, 2))
    
    data_frame(mod = list(fit))  })














