## S2_MET Prediction Execution

# This R script will perform the predictions for the S2MET project.
# See the 'S2mET_prediction_analysis.Rmd' for a document of the analysis of these results


## Setup
# First load libraries and set directories
library(tidyverse)
library(stringr)
library(rrBLUP)
library(readxl)
library(neyhart)
library(sommer)
library(gws)
library(lme4qtl)

# The head directory
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET/"

# Prediction directory
pred_dir <- file.path(proj_dir, "Predictions")

# Other directories
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/Genotypic Data/GBS Genotype Data/"
env_var_dir <- file.path(proj_dir, "Environmental_Variables/")
pheno_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET/Phenotype_Data/"


## Load Data
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))

# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUE.RData"))
load(file.path(pheno_dir, "S2_MET_tidy.RData"))

# Load environmental data
load(file.path(env_var_dir, "environmental_data_compiled.RData"))

# Load an entry file
entry_file <- file.path(proj_dir, "Plant_Materials/S2MET_project_entries.xlsx")
entry_list <- read_excel(entry_file)

# Get rid of checks
entries <- entry_list %>% 
  filter(Class != "Check") %>%
  pull(Line)

# Pull out the tp lines that intersect with the genotype data
tp <- entries %>% 
  str_subset("^[0-9]{2}") %>%
  intersect(., row.names(s2_imputed_mat))

vp <- entries %>% 
  str_subset("^2MS") %>%
  intersect(., row.names(s2_imputed_mat))

# Keep only the relevant lines
# Other phenotype data manipulations

## We need to combine some environments and prepare others for prediction, including:
# The environments "CNY15" and "HNY15" are the same

pheno <- S2_MET_BLUE %>%
  filter(line_name %in% c(tp, vp)) %>%
  mutate(environment = ifelse(environment == "CNY15", "HNY15", as.character(environment))) %>%
  # Change line_name and environment to factors
  mutate_at(vars(line_name, environment), funs(as.factor)) %>%
  droplevels()



## Matrix construction
# Genomic relationship matrix
# Convert the matrix name
M <- s2_imputed_mat[c(tp, vp), ]

# Calculate the genomic relationship matrix
G_1 <- A.mat(X = M, min.MAF = 0, max.missing = 1)


# Environmental relationship matrix

# The one_year_mat uses environmental covariates from the year of the trial
# The multi_year_mat uses environmental covariates from the past 10 years
N <- one_year_mat
# N <- multi_year_mat

# Create the relationship matrix
# E_1 is formed using the same procedure as the genomic relationship matrix
# E_2 is formed using the correlation of PCs
E_1 <- tcrossprod(N) / ncol(N)
E_2 <- svd(N)$u %>% 
  structure(dimnames = list(row.names(N), paste("PC", seq(nrow(.)), sep = ""))) %>%
  t() %>%
  cor()
  
  # PCA
  one_year_svd <- svd(one_year_mat)

# Use the PCs as the matrix - n_PC = n_EC
N <- one_year_svd$u %>%
  structure(dimnames = list(row.names(one_year_mat), paste("PC", seq(nrow(.)), sep = "")))

# Create the relationship matrix
E_2 <- cor(t(N))


# Create the GEI relationship matrix as the kronecker product of the G and E matrices`
O_1 <- kronecker(E_1, G_1, make.dimnames = TRUE)
O_2 <- kronecker(E_2, G_1, make.dimnames = TRUE)



#### Predictions

# Extract data for predictions

env_train <- pheno %>% 
  filter(line_name %in% tp) %>% 
  distinct(environment) %>% 
  pull() %>% 
  parse_character()

# For testing methods
env_train <- c("AID16", "BZD15", "HNY15", "CRM15", "STP16")

env_pred <- pheno %>% 
  filter(line_name %in% vp) %>% 
  distinct(environment) %>% 
  pull() %>% 
  parse_character()

env_pred <- c("CRM15", "STP16")

# Trait vector
traits <- pheno %>% 
  distinct(trait) %>% 
  pull() %>% 
  parse_character()

# # Save data for MSI
# save_file <- file.path(pred_dir, "MSI/S2MET_MSI_prediction_material.RData")
# # List of objects to save
# object_list <- c("E_1", "E_2", "O_1", "O_2", "env_pred", "G", "pheno", "tp", "traits", "vp")
# save(list = object_list, file = save_file)



### Homogenous Variance

## Normal GBLUP
## i.e. only genotypic main effect

# Extract data on only the environments to test
pheno_to_use <- pheno %>%
  filter(environment %in% c(env_pred, env_train)) %>%
  droplevels()

# Iterate over all traits
gblup2_unst_results <- expand.grid(environment = env_pred, trait = traits, stringsAsFactors = FALSE) %>%
  group_by(environment, trait) %>%
  do(pred_out = {
    
    env <- .$environment
    tr <- .$trait
    
    # Training and prediction phenotypes
    pheno_train <- filter(pheno_to_use, environment != env, trait == tr) %>%
      droplevels() %>%
      filter(line_name %in% tp)
    pheno_pred <- filter(pheno_to_use, line_name %in% vp, environment == env, trait == tr)
    
    # Create a model frame
    model_frame_train <- model.frame(value ~ line_name * environment, data = pheno_train)
    
    # Create matrices
    Y <- model.response(data = model_frame_train) %>%
      as.matrix() %>%
      structure(dimnames = list(NULL, trait))
      
      
    
    X <- fixef_model_matrix(fixed = value ~ environment, data = pheno_train)
    
    ZETA <- ranef_model_matrix(random = ~ g(line_name), data = pheno_train, 
                            vcov = list(line_name = G_1))
    
    R <- resid_model_matrix(resid = ~ units, data = pheno_train)
    
    # solve the model
    solve_out <- mmer(Y = Y, X = X, Z = ZETA, R = R)
    
    # Extract the random effects and return
    ranefs <- randef(solve_out)
    
    u_hat_df <- ranefs %>% 
      reduce(., cbind) %>% 
      as.data.frame() %>% 
      structure(names = names(ranefs)) %>% 
      rownames_to_column("id") %>%
      # rename_at(vars(-id), str_extract, "[A-Z]{3}[0-9]{2}") %>%
      separate(data = ., col = id, into = c("environment", "line_name"), sep = ":", 
               fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "right", "left")) %>%
      mutate_if(., .predicate = is.character, .funs = str_replace, 
                pattern = "line_name|environment", replacement = "") 
    
    # Combine predictions with observations
    pheno_pred %>% 
      mutate_at(vars(environment, line_name), parse_character) %>% 
      left_join(., dplyr::select(u_hat_df, -environment), by = "line_name") %>% 
      rename(pred = `g(line_name)`) })






































### Appendix

## Prediction Testing

# Assign training and prediction environments
env_train <- sort(c("STP16", "CRM16"))
env_pred <- c("FND16")

# Only use height data
pheno_use <- pheno %>%
  filter(trait == "GrainYield", environment %in% c(env_train, env_pred)) %>%
  droplevels()

# Pull out the phenotypic data
pheno_train <- pheno_use %>% 
  filter(environment %in% env_train) %>%
  droplevels() %>% 
  filter(line_name %in% tp)
pheno_pred <- pheno_use %>% 
  filter(environment %in% env_pred, line_name %in% vp)
  
# Filter the E matrix
E_1_use <- E_1[c(env_train), c(env_train)]




# Solve
solve_out <- mmer2(fixed = value ~ 1, 
                   random = ~ at(environment):g(line_name), 
                   G = list(line_name = G_1),
                   data = pheno_train)
                   
solve_out <- mmer2_fixed(fixed = value ~ 1 + environment, 
                   random = ~ at(environment):g(line_name), 
                   G = list(line_name = G_1), 
                   rcov = ~ at(environment):units,
                   data = pheno_train)
                   





                   
# Extract the random effects
ranefs <- randef(solve_out)

u_hat_df <- ranefs %>% 
  reduce(., cbind) %>% 
  as.data.frame() %>% 
  structure(names = names(ranefs)) %>% 
  rownames_to_column("id") %>%
  rename_at(vars(-id), str_extract, "[A-Z]{3}[0-9]{2}") %>%
  separate(data = ., col = id, into = c("environment", "line_name"), sep = ":", 
           fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "right", "left")) %>%
  mutate_if(., .predicate = is.character, .funs = str_replace, pattern = "line_name|environment", replacement = "")
    
    
# Convert to data.frames
u_hat_df <- solve_out$u %>% 
  data.frame(id = names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE) %>%
  separate(data = ., col = id, into = c("line_name", "environment"), sep = ":", 
           fill = ifelse(str_detect(.$id[1], "environment"), "left", "right")) %>%
  mutate_if(., .predicate = is.character, .funs = str_replace, pattern = "line_name|environment", replacement = "") %>%
  filter(line_name %in% vp) 






## Examples provided by package
data(cornHybrid)
hybrid2 <- cornHybrid$hybrid # extract cross data
A <- cornHybrid$K
K1 <- A[levels(hybrid2$GCA1), levels(hybrid2$GCA1)]; dim(K1)
K2 <- A[levels(hybrid2$GCA2), levels(hybrid2$GCA2)]; dim(K2)
S <- kronecker(K1, K2, make.dimnames=TRUE) ; dim(S)

ans <- mmer2(Yield ~ 1,
             random = ~ at(Location):g(GCA2),
             rcov = ~ units,
             data=hybrid2,
             G = list(GCA2 = A))

Z <- ans$Zforvec
             
# Try fitting using design matrices
Y <- model.frame(Yield ~ 1, data = hybrid2)

X <- model.matrix(Yield ~ 1, data = hybrid2, na.action = "na.pass")

# Try creating Z matrices on own
V <- model.frame(~ GCA2 + at(Location):GCA2, data = hybrid2, na.action = "na.pass")



Z_me <- 


Z_use <- list(GCA2_1 = list(Z = t(Z[[1]]), K = K2),
              GCA2_2 = list(Z = t(Z[[2]]), K = K2),
              GCA2_3 = list(Z = t(Z[[3]]), K = K2),
              GCA2_4 = list(Z = t(Z[[4]]), K = K2))


Z_use <- list(GCA2_1 = list(Z = Z[[1]], K = K2))


ans2 <- mmer(Y = Y, X = X, Z = Z_use)


    