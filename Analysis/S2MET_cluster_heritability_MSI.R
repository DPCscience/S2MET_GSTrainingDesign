## S2MET Cluster Heritability
## 
## This script will perform different clustering procedures on the S2MET data
## and calculate the within and across-cluster heritabilities.

# List of packages
packages <- c("tidyverse", "stringr", "readxl", "pbr", "regress", "rrBLUP", "Matrix", 
              "psych", "parallel")

# Set the directory of the R packages
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"
# package_dir <- NULL

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/" 
# proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET/"

fig_dir <- file.path(proj_dir, "Figures/")
pred_dir <- file.path(proj_dir, "Predictions")
env_var_dir <- file.path(proj_dir, "Environmental_Variables/")
entry_dir <- file.path(proj_dir, "Plant_Materials")
analysis_dir <- file.path(proj_dir, "Analysis")

geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/GBS_Genos"
pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"

# # Alternative directories
# geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/Genotypic Data/GBS Genotype Data/"
# pheno_dir <- file.path(proj_dir, "Phenotype_Data/")

# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))
# Load the heritability results (i.e. BLUPs)
load(file.path(pheno_dir, "S2_MET_heritability_results.RData"))

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

# Detect cores
n_cores <- detectCores()




## Define some functions for calculating hertiability
# Function to calculate across and within cluster heritabilities
clust_herit <- function(mod, g, e, cl) {
  # Model frame
  mf <- model.frame(mod)
  # Unique environments
  n_e <- n_distinct(mf[[e]])
  # Reps in environments
  n_r <- pbr::harm_mean(table(subset(mf, , c(g, e))))
  # Clusters
  n_s <- n_distinct(mf[[cl]])
  
  varcomp <- as.data.frame(VarCorr(mod))
  varG <- subset(varcomp, grp == g, vcov, drop = TRUE)
  varGE <- subset(varcomp, grp == paste(g, e, sep = ":"), vcov, drop = TRUE)
  varR <- subset(varcomp, grp == "Residual", vcov, drop = TRUE)
  varGS <- subset(varcomp, grp == paste(g, cl, sep = ":"), vcov, drop = TRUE)
  
  # Heritability - across
  H_a <- varG / (varG + (varGS / n_s) + (varGE / n_e) + (varR / (n_e * n_r)))
  # Heritability - within
  H_w <- (varG + varGS) / ((varG + varGS) + (n_s * ( (varGE / n_e) + (varR / (n_e * n_r)) ) ))
  data.frame(selection = c("across", "within"), heritability = c(H_a, H_w))
}


# Create a new data.frame to hold the different datasets
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  group_by(trait) %>% 
  nest()

# For each dataset, fit the full model
S2_MET_BLUPs <- S2_MET_BLUEs_use %>%
  group_by(trait) %>%
  mutate(data = list(mutate(data[[1]], gt = interaction(line_name, environment, drop = TRUE)) %>%
                       mutate_if(is.character, as.factor)))  %>%
  do(mod = {
    
    # Create an object for the data.frame
    df <- unnest(.)
    
    # Extract the standard errors
    R <- solve(diag(df$std_error^2))
    R <- structure(R, dimnames = replicate(2, str_c(df$line_name, df$environment, sep = ":"), 
                                           simplify = FALSE))
    
    # Use regress function
    fit <- regress(formula = value ~ environment,
                   Vformula = ~ line_name + gt + R, 
                   data = df, identity = FALSE, pos = rep(TRUE, 3))
    
    # Reorganize the BLUPs
    blup_out <- BLUP(model = fit)
    
    # Return the BLUPs
    blup_out$Mean %>% 
      data_frame(term = names(.), value = .) %>% 
      separate(col = term, into = c("term", "line_name", "environment"), 
               sep = "\\.", fill = "right") %>%
      select(-term) })


# Save these results
save_file <- file.path(pheno_dir, "S2_MET_BLUPs.RData")
save("S2_MET_BLUPs", file = save_file)



# Perform the distance calculation on the line means in each environment
ge_mean_D <- S2_MET_BLUEs_use %>% 
  unnest() %>%
  group_by(trait) %>%
  do(D = {
    dist1 <- dist_env(x = ., gen.col = "line_name", 
                      env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })

# Factor analysis
## Fit a diagonal model to get the BLUP of each genotype in each environment
S2_MET_BLUPs_diag <- S2_MET_BLUEs_use %>%
  group_by(trait) %>%
  mutate(data = list(mutate(data[[1]], gt = interaction(line_name, environment, drop = TRUE)) %>%
                       mutate_if(is.character, as.factor)))  %>%
  do(mod = {
    
    # Create an object for the data.frame
    df <- unnest(.)
    
    # Extract the standard errors
    R <- solve(Diagonal(x = df$std_error^2))
    dimnames(R) <- replicate(2, str_c(df$line_name, df$environment, sep = ":"), 
                             simplify = FALSE)
    
    # Model frame
    mf <- model.frame(value ~ line_name + environment + gt, data = df)
    y <- model.response(data = mf)
    X <- model.matrix(~ environment, data = mf)
    Z <- gws::ranef_model_matrix(~ at(environment):line_name, data = mf)
    
    # Fit
    fit <- sommer::mmer(Y = y, X = X, Z = Z, R = list(res = R), constraint = TRUE)
    
    # Return the BLUPs
    fit$u.hat %>% 
      map(as.data.frame) %>% 
      map(rownames_to_column, "term") %>% 
      map(rename, value = T1) %>% 
      map(~separate(., col = "term", into = c("environment", "line_name"), sep = ":")) %>% 
      map(mutate, environment = str_replace_all(string = environment, 
                                                pattern = "at\\(|\\)", replacement = "")) %>% 
      bind_rows()
    
  })

# Save this data
save_file <- file.path(pheno_dir, "S2_MET_BLUPs_diag.RData")
save("S2_MET_BLUPs_diag", file = save_file)


























