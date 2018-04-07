## S2_MET Prediction Execution

# This R script will perform the predictions for the S2MET project.
# See the 'S2mET_prediction_analysis.Rmd' for a document of the analysis of these results

## Note this is the testing script with limited predictions


## Setup
# First load libraries and set directories
library(tidyverse)
library(stringr)
library(readxl)
library(neyhart)
library(sommer)
library(gws)

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


# Levels for the different models
model_levels <- c("GBLUP2", "aGEBLUP2", "iGEBLUP2", "GBLUP2_Str", "aGEBLUP2_Str", 
                  "iGEBLUP2_Str", "aGEBLUP2_E2", "iGEBLUP2_E2", "aGEBLUP2_Str_E2", 
                  "iGEBLUP2_Str_E2")



###################
#### Use E_1 ######
###################




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
      structure(dimnames = list(NULL, tr))
      
      
    
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


## Analyze
# Perform a bootstrap correlation and plot
gblup2_unst_analysis <- gblup2_unst_results %>% 
  ungroup() %>% 
  unnest() %>% 
  group_by(environment, trait) %>% 
  do({pred_acc = boot_cor(x = .$value, y = .$pred, boot.reps = 1000)}) %>%
  mutate(model = "GBLUP2") %>% 
  dplyr::select(model, environment, trait, names(.))

gblup2_unst_analysis %>%
  ggplot(aes(x = environment, y = r_hat, fill = model, group = model, ymin = CI_lower, ymax = CI_upper)) +
  geom_col(position = "dodge") +
  geom_errorbar(position = "dodge", width = 0.5) +
  facet_grid(facets = . ~ trait) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Prediction Accuracy") +
  xlab("Environment") +
  scale_fill_manual(values = umn_palette("Secondary_Tier1")) +
  labs(title = "G-BLUP Prediction Accuracy")



## Additive G + E GBLUP
## i.e. only genotypic and environmental main effects
# Extract data on only the environments to test
pheno_to_use <- pheno %>%
  filter(environment %in% c(env_pred, env_train)) %>%
  droplevels()

# Iterate over all traits
a_geblup2_unst_results <- expand.grid(environment = env_pred, trait = traits, stringsAsFactors = FALSE) %>%
  group_by(environment, trait) %>%
  do(pred_out = {
    
    env <- .$environment
    tr <- .$trait
    
    # Training and prediction phenotypes
    pheno_train <- pheno_to_use %>% 
      filter(environment != env, trait == tr, line_name %in% tp)
    pheno_pred <- filter(pheno_to_use, line_name %in% vp, environment == env, trait == tr)
    
    # Create a model frame
    model_frame_train <- model.frame(value ~ line_name * environment, data = pheno_train)
    
    # Edit the covariance matrices
    E_1_use <- E_1[levels(pheno_train$environment), levels(pheno_train$environment)]
    
    # Create matrices
    Y <- model.response(data = model_frame_train) %>%
      as.matrix() %>%
      structure(dimnames = list(NULL, tr))
    
    X <- fixef_model_matrix(fixed = value ~ 1, data = pheno_train)
    
    ZETA <- ranef_model_matrix(random = ~ g(line_name) + g(environment), data = pheno_train, 
                               vcov = list(line_name = G_1, environment = E_1_use))
    
    R <- resid_model_matrix(resid = ~ units, data = pheno_train)
    
    # solve the model
    solve_out <- mmer(Y = Y, X = X, Z = ZETA, R = R)
    
    # Extract the random effects and return
    ranefs <- randef(solve_out)
    
    u_hat_df <- ranefs %>%
      map(~data.frame(id = row.names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE)) %>% 
      map(~rename_(.data = ., u_hat = tr)) %>%
      map(~separate(data = ., col = id, into = c("line_name", "environment"), sep = ":", 
                    fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "left", "right"))) %>%
      map(~mutate_if(., .predicate = is.character, .funs = str_replace, 
                     pattern = "line_name|environment", replacement = ""))
    
    # Combine with phenos of the predicted lines
    pheno_pred %>% 
      mutate_at(vars(environment, line_name), parse_character) %>% 
      left_join(., dplyr::select(u_hat_df$`g(line_name)`, -environment), by = "line_name") %>% 
      left_join(., dplyr::select(u_hat_df$`g(environment)`, -line_name), by = "environment") %>% 
      mutate(pred = u_hat.x + u_hat.y) %>%
      dplyr::select(-u_hat.x, -u_hat.y) })


## Analyze
# Perform a bootstrap correlation and plot
a_geblup2_unst_analysis <- a_geblup2_unst_results %>% 
  ungroup() %>% 
  unnest() %>% 
  group_by(environment, trait) %>% 
  do({pred_acc = boot_cor(x = .$value, y = .$pred, boot.reps = 1000)}) %>%
  mutate(model = "aGEBLUP2") %>% 
  dplyr::select(model, environment, trait, names(.))

# Combine graphs
# Combine
bind_rows(a_geblup2_unst_analysis, gblup2_unst_analysis) %>%
  mutate(model = parse_factor(model, levels = c("GBLUP2", "aGEBLUP2", "iGEBLUP2"))) %>%
  ggplot(aes(x = environment, y = r_hat, fill = model, group = model, ymin = CI_lower, ymax = CI_upper)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(position = position_dodge(1), width = 0.5) +
  facet_grid(facets = . ~ trait) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Prediction Accuracy") +
  xlab("Environment") +
  scale_fill_manual(values = umn_palette("Secondary_Tier1")) +
  labs(title = "G-BLUP Prediction Accuracy")


## Interaction G + E + GE GBLUP
# Extract data on only the environments to test
pheno_to_use <- pheno %>%
  filter(environment %in% c(env_pred, env_train)) %>%
  droplevels()

# Iterate over all traits
i_geblup2_unst_results <- expand.grid(environment = env_pred, trait = traits, stringsAsFactors = FALSE) %>%
  group_by(environment, trait) %>%
  do(pred_out = {
    
    env <- .$environment
    tr <- .$trait
    
    # Training and prediction phenotypes
    pheno_train <- pheno_to_use %>% 
      filter(environment != env, trait == tr, line_name %in% tp)
    pheno_pred <- filter(pheno_to_use, line_name %in% vp, environment == env, trait == tr)
    
    # Create a model frame
    model_frame_train <- model.frame(value ~ line_name * environment, data = pheno_train)
    
    # Edit the covariance matrices
    E_1_use <- E_1[levels(pheno_train$environment), levels(pheno_train$environment)]
    # Create the interaction matrix
    O_1_use <- kronecker(G_1, E_1_use, make.dimnames = TRUE)
    
    # Create matrices
    Y <- model.response(data = model_frame_train) %>%
      as.matrix() %>%
      structure(dimnames = list(NULL, tr))
    
    X <- fixef_model_matrix(fixed = value ~ 1, data = pheno_train)
    
    ZETA <- ranef_model_matrix(random = ~ g(line_name) + g(environment) + g(line_name:environment), 
                               data = pheno_train, 
                               vcov = list(line_name = G_1, environment = E_1_use,
                                           `line_name:environment` = O_1_use))
    
    R <- resid_model_matrix(resid = ~ units, data = pheno_train)
    
    # solve the model
    solve_out <- mmer(Y = Y, X = X, Z = ZETA, R = R)
    
    # Extract the random effects and return
    ranefs <- randef(solve_out)
    
    u_hat_df <- ranefs %>%
      map(~data.frame(id = row.names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE)) %>% 
      map(~rename_(.data = ., u_hat = tr)) %>%
      map(~separate(data = ., col = id, into = c("line_name", "environment"), sep = ":", 
                    fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "left", "right"))) %>%
      map(~mutate_if(., .predicate = is.character, .funs = str_replace, 
                     pattern = "line_name|environment", replacement = ""))
    
    # Combine with phenos of the predicted lines
    pheno_pred %>% 
      mutate_at(vars(environment, line_name), parse_character) %>% 
      left_join(., dplyr::select(u_hat_df$`g(line_name)`, -environment), by = "line_name") %>% 
      left_join(., dplyr::select(u_hat_df$`g(environment)`, -line_name), by = "environment") %>% 
      left_join(., u_hat_df$`g(line_name:environment)`, by = c("line_name", "environment")) %>%
      mutate(pred = u_hat.x + u_hat.y + u_hat) %>%
      dplyr::select(-u_hat.x, -u_hat.y, -u_hat) })


## Analyze
# Perform a bootstrap correlation and plot
i_geblup2_unst_analysis <- i_geblup2_unst_results %>% 
  ungroup() %>% 
  unnest() %>% 
  group_by(environment, trait) %>% 
  do({pred_acc = boot_cor(x = .$value, y = .$pred, boot.reps = 1000)}) %>%
  mutate(model = "iGEBLUP2") %>% 
  dplyr::select(model, environment, trait, names(.))

# Combine
bind_rows(i_geblup2_unst_analysis, a_geblup2_unst_analysis, gblup2_unst_analysis) %>%
  mutate(model = parse_factor(model, levels = c("GBLUP2", "aGEBLUP2", "iGEBLUP2"))) %>%
  ggplot(aes(x = environment, y = r_hat, fill = model, group = model, ymin = CI_lower, ymax = CI_upper)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(position = position_dodge(1), width = 0.5) +
  facet_grid(facets = . ~ trait) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Prediction Accuracy") +
  xlab("Environment") +
  scale_fill_manual(values = umn_palette("Secondary_Tier1")) +
  labs(title = "G-BLUP Prediction Accuracy")






### Heterogenous Variance

## Compound Symmetry and heterogenous variance
## We only use the base prediction of the genotypic value for prediction?

## This does not allow for the prediction of new environments.

# Extract data on only the environments to test
pheno_to_use <- pheno %>%
  filter(environment %in% c(env_pred, env_train)) %>%
  droplevels()

# Iterate over all traits
gblup2_str_results <- expand.grid(environment = env_pred, trait = traits, stringsAsFactors = FALSE) %>%
  group_by(environment, trait) %>%
  do(pred_out = {
    
    env <- .$environment
    tr <- .$trait
    
    # Training and prediction phenotypes
    pheno_train <- filter(pheno_to_use, environment != env, trait == tr) %>%
      filter(line_name %in% tp)
    pheno_pred <- filter(pheno_to_use, line_name %in% vp, environment == env, trait == tr)
    
    # Create a model frame
    model_frame_train <- model.frame(value ~ line_name * environment, data = pheno_train)
    
    
    # Create matrices
    Y <- model.response(data = model_frame_train) %>%
      as.matrix() %>%
      structure(dimnames = list(NULL, tr))
    
    X <- fixef_model_matrix(fixed = value ~ 1, data = pheno_train)
    
    ZETA <- ranef_model_matrix(random = ~ g(line_name) + at(environment):g(line_name), 
                               data = pheno_train, 
                               vcov = list(line_name = G_1))
    
    R <- resid_model_matrix(resid = ~ at(environment):units, data = pheno_train)
    
    
    
    # solve the model
    solve_out <- mmer(Y = Y, X = X, Z = ZETA, R = R)
    
    # Random effects to extract
    ranefs_tokeep <- c("g(line_name)")
    
    # Extract the random effects and return
    ranefs <- randef(solve_out)[ranefs_tokeep]
    
    # Only extract the base predictions of genotypic values (i.e. compound symmetry)
    u_hat_df <- ranefs %>% 
      map(~data.frame(id = row.names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE)) %>% 
      map(~rename_(.data = ., u_hat = tr)) %>%
      map(~separate(data = ., col = id, into = c("line_name", "environment"), sep = ":", 
                    fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "left", "right"))) %>%
      map(~mutate_if(., .predicate = is.character, .funs = str_replace, 
                     pattern = "line_name|environment", replacement = ""))
    
    # Combine predictions with observations
    pheno_pred %>% 
      mutate_at(vars(environment, line_name), parse_character) %>% 
      left_join(., dplyr::select(u_hat_df$`g(line_name)`, -environment), by = "line_name") %>% 
      rename(pred = u_hat) })


## Analyze
# Perform a bootstrap correlation and plot
gblup2_str_analysis <- gblup2_str_results %>% 
  ungroup() %>% 
  unnest() %>% 
  group_by(environment, trait) %>% 
  do({pred_acc = boot_cor(x = .$value, y = .$pred, boot.reps = 1000)}) %>%
  mutate(model = "GBLUP2_Str") %>% 
  dplyr::select(model, environment, trait, names(.))

# Compare with normal GBLUP
bind_rows(gblup2_str_analysis, gblup2_unst_analysis) %>%
  mutate(model = parse_factor(model, levels = c("GBLUP2", "aGEBLUP2", "iGEBLUP2", "GBLUP2_Str"))) %>%
  ggplot(aes(x = environment, y = r_hat, fill = model, group = model, ymin = CI_lower, ymax = CI_upper)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(position = position_dodge(1), width = 0.5) +
  facet_grid(facets = . ~ trait) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Prediction Accuracy") +
  xlab("Environment") +
  scale_fill_manual(values = umn_palette("Secondary_Tier1")) +
  labs(title = "G-BLUP Prediction Accuracy")






## Additive G + E GBLUP
## i.e. only genotypic and environmental main effects
# Extract data on only the environments to test
pheno_to_use <- pheno %>%
  filter(environment %in% c(env_pred, env_train)) %>%
  droplevels()

# Iterate over all traits
a_geblup2_str_results <- expand.grid(environment = env_pred, trait = traits, stringsAsFactors = FALSE) %>%
  group_by(environment, trait) %>%
  do(pred_out = {
    
    env <- .$environment
    tr <- .$trait
    
    # Training and prediction phenotypes
    pheno_train <- pheno_to_use %>% 
      filter(environment != env, trait == tr, line_name %in% tp)
    pheno_pred <- filter(pheno_to_use, line_name %in% vp, environment == env, trait == tr)
    
    # Create a model frame
    model_frame_train <- model.frame(value ~ line_name * environment, data = pheno_train)
    
    # Edit the covariance matrices
    E_1_use <- E_1[levels(pheno_train$environment), levels(pheno_train$environment)]
    
    # Create matrices
    Y <- model.response(data = model_frame_train) %>%
      as.matrix() %>%
      structure(dimnames = list(NULL, tr))
    
    X <- fixef_model_matrix(fixed = value ~ 1, data = pheno_train)
    
    ZETA <- ranef_model_matrix(random = ~ g(line_name) + at(environment):g(line_name) + g(environment), 
                               data = pheno_train, 
                               vcov = list(line_name = G_1, environment = E_1_use))
    
    R <- resid_model_matrix(resid = ~ at(environment):units, data = pheno_train)
    
    # solve the model
    solve_out <- mmer(Y = Y, X = X, Z = ZETA, R = R)
    
    # Random effects to extract
    ranefs_tokeep <- c("g(line_name)", "g(environment)")
    
    # Extract the random effects and return
    ranefs <- randef(solve_out)[ranefs_tokeep]
    
    u_hat_df <- ranefs %>%
      map(~data.frame(id = row.names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE)) %>% 
      map(~rename_(.data = ., u_hat = tr)) %>%
      map(~separate(data = ., col = id, into = c("line_name", "environment"), sep = ":", 
                    fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "left", "right"))) %>%
      map(~mutate_if(., .predicate = is.character, .funs = str_replace, 
                     pattern = "line_name|environment", replacement = ""))
    
    # Combine with phenos of the predicted lines
    pheno_pred %>% 
      mutate_at(vars(environment, line_name), parse_character) %>% 
      left_join(., dplyr::select(u_hat_df$`g(line_name)`, -environment), by = "line_name") %>% 
      left_join(., dplyr::select(u_hat_df$`g(environment)`, -line_name), by = "environment") %>% 
      mutate(pred = u_hat.x + u_hat.y) %>%
      dplyr::select(-u_hat.x, -u_hat.y) })


## Analyze
# Perform a bootstrap correlation and plot
a_geblup2_str_analysis <- a_geblup2_str_results %>% 
  ungroup() %>% 
  unnest() %>% 
  group_by(environment, trait) %>% 
  do({pred_acc = boot_cor(x = .$value, y = .$pred, boot.reps = 1000)}) %>%
  mutate(model = "aGEBLUP2_Str") %>% 
  dplyr::select(model, environment, trait, names(.))

# Combine graphs
# Combine
bind_rows(a_geblup2_str_analysis, a_geblup2_unst_analysis, gblup2_str_analysis, gblup2_unst_analysis) %>%
  mutate(model = parse_factor(model, levels = c("GBLUP2", "aGEBLUP2", "iGEBLUP2", "GBLUP2_Str", 
                                                "aGEBLUP2_Str"))) %>%
  ggplot(aes(x = environment, y = r_hat, fill = model, group = model, ymin = CI_lower, ymax = CI_upper)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(position = position_dodge(1), width = 0.5) +
  facet_grid(facets = . ~ trait) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Prediction Accuracy") +
  xlab("Environment") +
  scale_fill_manual(values = umn_palette("Secondary_Tier1")) +
  labs(title = "G-BLUP Prediction Accuracy")


## Interaction G + E + GE GBLUP
# Extract data on only the environments to test
pheno_to_use <- pheno %>%
  filter(environment %in% c(env_pred, env_train)) %>%
  droplevels()


# Iterate over all traits
i_geblup2_str_results <- expand.grid(environment = env_pred, trait = traits, stringsAsFactors = FALSE) %>%
  group_by(environment, trait) %>%
  do(pred_out = {
    
    env <- .$environment
    tr <- .$trait
    
    # Training and prediction phenotypes
    pheno_train <- pheno_to_use %>% 
      filter(environment != env, trait == tr, line_name %in% tp)
    pheno_pred <- filter(pheno_to_use, line_name %in% vp, environment == env, trait == tr)
    
    # Create a model frame
    model_frame_train <- model.frame(value ~ line_name * environment, data = pheno_train)
    
    # Edit the covariance matrices
    E_1_use <- E_1[levels(pheno_train$environment), levels(pheno_train$environment)]
    # Create the interaction matrix
    O_1_use <- kronecker(G_1, E_1_use, make.dimnames = TRUE)
    
    # Create matrices
    Y <- model.response(data = model_frame_train) %>%
      as.matrix() %>%
      structure(dimnames = list(NULL, tr))
    
    X <- fixef_model_matrix(fixed = value ~ 1, data = pheno_train)
    
    ZETA <- ranef_model_matrix(random = ~ g(line_name) + at(environment):g(line_name) + 
                                 g(environment) + g(line_name:environment), 
                               data = pheno_train, 
                               vcov = list(line_name = G_1, environment = E_1_use,
                                           `line_name:environment` = O_1_use))
    
    R <- resid_model_matrix(resid = ~ at(environment):units, data = pheno_train)
    
    # solve the model
    solve_out <- mmer(Y = Y, X = X, Z = ZETA, R = R)
    
    # Random effects to extract
    ranefs_tokeep <- c("g(line_name)", "g(environment)", "g(line_name:environment)")
    
    # Extract the random effects and return
    ranefs <- randef(solve_out)[ranefs_tokeep]
    
    u_hat_df <- ranefs %>%
      map(~data.frame(id = row.names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE)) %>% 
      map(~rename_(.data = ., u_hat = tr)) %>%
      # Reverse the order of the 'at' random effects
      map(~separate(data = ., col = id, into = c("line_name", "environment"), sep = ":", 
                    fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "left", "right"))) %>%
      map(~mutate_if(., .predicate = is.character, .funs = str_replace, 
                     pattern = "line_name|environment", replacement = ""))
    
    # Combine with phenos of the predicted lines
    pheno_pred %>% 
      mutate_at(vars(environment, line_name), parse_character) %>% 
      left_join(., dplyr::select(u_hat_df$`g(line_name)`, -environment), by = "line_name") %>% 
      left_join(., dplyr::select(u_hat_df$`g(environment)`, -line_name), by = "environment") %>% 
      left_join(., u_hat_df$`g(line_name:environment)`, by = c("line_name", "environment")) %>%
      mutate(pred = u_hat.x + u_hat.y + u_hat) %>%
      dplyr::select(-u_hat.x, -u_hat.y, -u_hat) })


## Analyze
# Perform a bootstrap correlation and plot
i_geblup2_str_analysis <- i_geblup2_str_results %>% 
  ungroup() %>% 
  unnest() %>% 
  group_by(environment, trait) %>% 
  do({pred_acc = boot_cor(x = .$value, y = .$pred, boot.reps = 1000)}) %>%
  mutate(model = "iGEBLUP2_Str") %>% 
  dplyr::select(model, environment, trait, names(.))

# Combine
bind_rows(i_geblup2_str_analysis, a_geblup2_str_analysis, gblup2_str_analysis,
          i_geblup2_unst_analysis, a_geblup2_unst_analysis, gblup2_unst_analysis) %>%
  mutate(model = parse_factor(model, levels = c("GBLUP2", "aGEBLUP2", "iGEBLUP2", "GBLUP2_Str", 
                                                "aGEBLUP2_Str", "iGEBLUP2_Str"))) %>%
  ggplot(aes(x = environment, y = r_hat, fill = model, group = model, ymin = CI_lower, ymax = CI_upper)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(position = position_dodge(1), width = 0.5) +
  facet_grid(facets = . ~ trait) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Prediction Accuracy") +
  xlab("Environment") +
  scale_fill_manual(values = umn_palette("Secondary_Tier1")) +
  labs(title = "G-BLUP Prediction Accuracy")














###################
#### Use E_2 ######
###################




### Homogenous Variance

## Additive G + E GBLUP
## i.e. only genotypic and environmental main effects
# Extract data on only the environments to test
pheno_to_use <- pheno %>%
  filter(environment %in% c(env_pred, env_train)) %>%
  droplevels()

# Iterate over all traits
a_geblup2_unst_results_E2 <- expand.grid(environment = env_pred, trait = traits, stringsAsFactors = FALSE) %>%
  group_by(environment, trait) %>%
  do(pred_out = {
    
    env <- .$environment
    tr <- .$trait
    
    # Training and prediction phenotypes
    pheno_train <- pheno_to_use %>% 
      filter(environment != env, trait == tr, line_name %in% tp)
    pheno_pred <- filter(pheno_to_use, line_name %in% vp, environment == env, trait == tr)
    
    # Create a model frame
    model_frame_train <- model.frame(value ~ line_name * environment, data = pheno_train)
    
    # Edit the covariance matrices
    E_2_use <- E_2[levels(pheno_train$environment), levels(pheno_train$environment)]
    
    # Create matrices
    Y <- model.response(data = model_frame_train) %>%
      as.matrix() %>%
      structure(dimnames = list(NULL, tr))
    
    X <- fixef_model_matrix(fixed = value ~ 1, data = pheno_train)
    
    ZETA <- ranef_model_matrix(random = ~ g(line_name) + g(environment), data = pheno_train, 
                               vcov = list(line_name = G_1, environment = E_2_use))
    
    R <- resid_model_matrix(resid = ~ units, data = pheno_train)
    
    # solve the model
    solve_out <- mmer(Y = Y, X = X, Z = ZETA, R = R)
    
    # Extract the random effects and return
    ranefs <- randef(solve_out)
    
    u_hat_df <- ranefs %>%
      map(~data.frame(id = row.names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE)) %>% 
      map(~rename_(.data = ., u_hat = tr)) %>%
      map(~separate(data = ., col = id, into = c("line_name", "environment"), sep = ":", 
                    fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "left", "right"))) %>%
      map(~mutate_if(., .predicate = is.character, .funs = str_replace, 
                     pattern = "line_name|environment", replacement = ""))
    
    # Combine with phenos of the predicted lines
    pheno_pred %>% 
      mutate_at(vars(environment, line_name), parse_character) %>% 
      left_join(., dplyr::select(u_hat_df$`g(line_name)`, -environment), by = "line_name") %>% 
      left_join(., dplyr::select(u_hat_df$`g(environment)`, -line_name), by = "environment") %>% 
      mutate(pred = u_hat.x + u_hat.y) %>%
      dplyr::select(-u_hat.x, -u_hat.y) })


## Analyze
# Perform a bootstrap correlation and plot
a_geblup2_unst_analysis_E2 <- a_geblup2_unst_results_E2 %>% 
  ungroup() %>% 
  unnest() %>% 
  group_by(environment, trait) %>% 
  do({pred_acc = boot_cor(x = .$value, y = .$pred, boot.reps = 1000)}) %>%
  mutate(model = "aGEBLUP2_E2") %>% 
  dplyr::select(model, environment, trait, names(.))

# Combine graphs
# Combine
bind_rows(a_geblup2_unst_analysis, a_geblup2_unst_analysis_E2) %>%
  mutate(model = parse_factor(model, levels = model_levels)) %>%
  ggplot(aes(x = environment, y = r_hat, fill = model, group = model, ymin = CI_lower, ymax = CI_upper)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(position = position_dodge(1), width = 0.5) +
  facet_grid(facets = . ~ trait) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Prediction Accuracy") +
  xlab("Environment") +
  scale_fill_manual(values = umn_palette("Secondary_Tier1")) +
  labs(title = "G-BLUP Prediction Accuracy")


## Interaction G + E + GE GBLUP
# Extract data on only the environments to test
pheno_to_use <- pheno %>%
  filter(environment %in% c(env_pred, env_train)) %>%
  droplevels()

# Iterate over all traits
i_geblup2_unst_results_E2 <- expand.grid(environment = env_pred, trait = traits, stringsAsFactors = FALSE) %>%
  group_by(environment, trait) %>%
  do(pred_out = {
    
    env <- .$environment
    tr <- .$trait
    
    # Training and prediction phenotypes
    pheno_train <- pheno_to_use %>% 
      filter(environment != env, trait == tr, line_name %in% tp)
    pheno_pred <- filter(pheno_to_use, line_name %in% vp, environment == env, trait == tr)
    
    # Create a model frame
    model_frame_train <- model.frame(value ~ line_name * environment, data = pheno_train)
    
    # Edit the covariance matrices
    E_2_use <- E_2[levels(pheno_train$environment), levels(pheno_train$environment)]
    # Create the interaction matrix
    O_2_use <- kronecker(G_1, E_2_use, make.dimnames = TRUE)
    
    # Create matrices
    Y <- model.response(data = model_frame_train) %>%
      as.matrix() %>%
      structure(dimnames = list(NULL, tr))
    
    X <- fixef_model_matrix(fixed = value ~ 1, data = pheno_train)
    
    ZETA <- ranef_model_matrix(random = ~ g(line_name) + g(environment) + g(line_name:environment), 
                               data = pheno_train, 
                               vcov = list(line_name = G_1, environment = E_2_use,
                                           `line_name:environment` = O_2_use))
    
    R <- resid_model_matrix(resid = ~ units, data = pheno_train)
    
    # solve the model
    solve_out <- mmer(Y = Y, X = X, Z = ZETA, R = R, silent = TRUE)
    
    # Extract the random effects and return
    ranefs <- randef(solve_out)
    
    u_hat_df <- ranefs %>%
      map(~data.frame(id = row.names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE)) %>% 
      map(~rename_(.data = ., u_hat = tr)) %>%
      map(~separate(data = ., col = id, into = c("line_name", "environment"), sep = ":", 
                    fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "left", "right"))) %>%
      map(~mutate_if(., .predicate = is.character, .funs = str_replace, 
                     pattern = "line_name|environment", replacement = ""))
    
    # Combine with phenos of the predicted lines
    pheno_pred %>% 
      mutate_at(vars(environment, line_name), parse_character) %>% 
      left_join(., dplyr::select(u_hat_df$`g(line_name)`, -environment), by = "line_name") %>% 
      left_join(., dplyr::select(u_hat_df$`g(environment)`, -line_name), by = "environment") %>% 
      left_join(., u_hat_df$`g(line_name:environment)`, by = c("line_name", "environment")) %>%
      mutate(pred = u_hat.x + u_hat.y + u_hat) %>%
      dplyr::select(-u_hat.x, -u_hat.y, -u_hat) })


## Analyze
# Perform a bootstrap correlation and plot
i_geblup2_unst_analysis_E2 <- i_geblup2_unst_results_E2 %>% 
  ungroup() %>% 
  unnest() %>% 
  group_by(environment, trait) %>% 
  do({pred_acc = boot_cor(x = .$value, y = .$pred, boot.reps = 1000)}) %>%
  mutate(model = "iGEBLUP2_E2") %>% 
  dplyr::select(model, environment, trait, names(.))

# Combine
bind_rows(i_geblup2_unst_analysis_E2, i_geblup2_unst_analysis, a_geblup2_unst_analysis, 
          a_geblup2_unst_analysis_E2, gblup2_unst_analysis) %>%
  mutate(model = parse_factor(model, levels = model_levels)) %>%
  ggplot(aes(x = environment, y = r_hat, fill = model, group = model, ymin = CI_lower, ymax = CI_upper)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(position = position_dodge(1), width = 0.5) +
  facet_grid(facets = . ~ trait) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Prediction Accuracy") +
  xlab("Environment") +
  scale_fill_manual(values = umn_palette("Secondary_Tier1")) +
  labs(title = "G-BLUP Prediction Accuracy")





### Heterogenous Variance

## Compound Symmetry and heterogenous variance
## We only use the base prediction of the genotypic value for prediction?


## Additive G + E GBLUP
## i.e. only genotypic and environmental main effects
# Extract data on only the environments to test
pheno_to_use <- pheno %>%
  filter(environment %in% c(env_pred, env_train)) %>%
  droplevels()

# Iterate over all traits
a_geblup2_str_results_E2 <- expand.grid(environment = env_pred, trait = traits, stringsAsFactors = FALSE) %>%
  group_by(environment, trait) %>%
  do(pred_out = {
    
    env <- .$environment
    tr <- .$trait
    
    # Training and prediction phenotypes
    pheno_train <- pheno_to_use %>% 
      filter(environment != env, trait == tr, line_name %in% tp)
    pheno_pred <- filter(pheno_to_use, line_name %in% vp, environment == env, trait == tr)
    
    # Create a model frame
    model_frame_train <- model.frame(value ~ line_name * environment, data = pheno_train)
    
    # Edit the covariance matrices
    E_2_use <- E_2[levels(pheno_train$environment), levels(pheno_train$environment)]
    
    # Create matrices
    Y <- model.response(data = model_frame_train) %>%
      as.matrix() %>%
      structure(dimnames = list(NULL, tr))
    
    X <- fixef_model_matrix(fixed = value ~ 1, data = pheno_train)
    
    ZETA <- ranef_model_matrix(random = ~ g(line_name) + at(environment):g(line_name) + g(environment), 
                               data = pheno_train, 
                               vcov = list(line_name = G_1, environment = E_2_use))
    
    R <- resid_model_matrix(resid = ~ at(environment):units, data = pheno_train)
    
    # solve the model
    solve_out <- mmer(Y = Y, X = X, Z = ZETA, R = R, silent = TRUE)
    
    # Random effects to extract
    ranefs_tokeep <- c("g(line_name)", "g(environment)")
    
    # Extract the random effects and return
    ranefs <- randef(solve_out)[ranefs_tokeep]
    
    u_hat_df <- ranefs %>%
      map(~data.frame(id = row.names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE)) %>% 
      map(~rename_(.data = ., u_hat = tr)) %>%
      map(~separate(data = ., col = id, into = c("line_name", "environment"), sep = ":", 
                    fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "left", "right"))) %>%
      map(~mutate_if(., .predicate = is.character, .funs = str_replace, 
                     pattern = "line_name|environment", replacement = ""))
    
    # Combine with phenos of the predicted lines
    pheno_pred %>% 
      mutate_at(vars(environment, line_name), parse_character) %>% 
      left_join(., dplyr::select(u_hat_df$`g(line_name)`, -environment), by = "line_name") %>% 
      left_join(., dplyr::select(u_hat_df$`g(environment)`, -line_name), by = "environment") %>% 
      mutate(pred = u_hat.x + u_hat.y) %>%
      dplyr::select(-u_hat.x, -u_hat.y) })


## Analyze
# Perform a bootstrap correlation and plot
a_geblup2_str_analysis_E2 <- a_geblup2_str_results_E2 %>% 
  ungroup() %>% 
  unnest() %>% 
  group_by(environment, trait) %>% 
  do({pred_acc = boot_cor(x = .$value, y = .$pred, boot.reps = 1000)}) %>%
  mutate(model = "aGEBLUP2_Str_E2") %>% 
  dplyr::select(model, environment, trait, names(.))

## Interaction G + E + GE GBLUP
# Extract data on only the environments to test
pheno_to_use <- pheno %>%
  filter(environment %in% c(env_pred, env_train)) %>%
  droplevels()


# Iterate over all traits
i_geblup2_str_results_E2 <- expand.grid(environment = env_pred, trait = traits, stringsAsFactors = FALSE) %>%
  group_by(environment, trait) %>%
  do(pred_out = {
    
    env <- .$environment
    tr <- .$trait
    
    # Training and prediction phenotypes
    pheno_train <- pheno_to_use %>% 
      filter(environment != env, trait == tr, line_name %in% tp)
    pheno_pred <- filter(pheno_to_use, line_name %in% vp, environment == env, trait == tr)
    
    # Create a model frame
    model_frame_train <- model.frame(value ~ line_name * environment, data = pheno_train)
    
    # Edit the covariance matrices
    E_1_use <- E_1[levels(pheno_train$environment), levels(pheno_train$environment)]
    # Create the interaction matrix
    O_1_use <- kronecker(G_1, E_1_use, make.dimnames = TRUE)
    
    # Create matrices
    Y <- model.response(data = model_frame_train) %>%
      as.matrix() %>%
      structure(dimnames = list(NULL, tr))
    
    X <- fixef_model_matrix(fixed = value ~ 1, data = pheno_train)
    
    ZETA <- ranef_model_matrix(random = ~ g(line_name) + at(environment):g(line_name) + 
                                 g(environment) + g(line_name:environment), 
                               data = pheno_train, 
                               vcov = list(line_name = G_1, environment = E_1_use,
                                           `line_name:environment` = O_1_use))
    
    R <- resid_model_matrix(resid = ~ at(environment):units, data = pheno_train)
    
    # solve the model
    solve_out <- mmer(Y = Y, X = X, Z = ZETA, R = R)
    
    # Random effects to extract
    ranefs_tokeep <- c("g(line_name)", "g(environment)", "g(line_name:environment)")
    
    # Extract the random effects and return
    ranefs <- randef(solve_out)[ranefs_tokeep]
    
    u_hat_df <- ranefs %>%
      map(~data.frame(id = row.names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE)) %>% 
      map(~rename_(.data = ., u_hat = tr)) %>%
      # Reverse the order of the 'at' random effects
      map(~separate(data = ., col = id, into = c("line_name", "environment"), sep = ":", 
                    fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "left", "right"))) %>%
      map(~mutate_if(., .predicate = is.character, .funs = str_replace, 
                     pattern = "line_name|environment", replacement = ""))
    
    # Combine with phenos of the predicted lines
    pheno_pred %>% 
      mutate_at(vars(environment, line_name), parse_character) %>% 
      left_join(., dplyr::select(u_hat_df$`g(line_name)`, -environment), by = "line_name") %>% 
      left_join(., dplyr::select(u_hat_df$`g(environment)`, -line_name), by = "environment") %>% 
      left_join(., u_hat_df$`g(line_name:environment)`, by = c("line_name", "environment")) %>%
      mutate(pred = u_hat.x + u_hat.y + u_hat) %>%
      dplyr::select(-u_hat.x, -u_hat.y, -u_hat) })


## Analyze
# Perform a bootstrap correlation and plot
i_geblup2_str_analysis_E2 <- i_geblup2_str_results_E2 %>% 
  ungroup() %>% 
  unnest() %>% 
  group_by(environment, trait) %>% 
  do({pred_acc = boot_cor(x = .$value, y = .$pred, boot.reps = 1000)}) %>%
  mutate(model = "iGEBLUP2_Str_E2") %>% 
  dplyr::select(model, environment, trait, names(.))

# Combine
map(cfind(what = "blup2.*analysis"), get) %>%
  bind_rows() %>%
  mutate(model = parse_factor(model, levels = model_levels)) %>%
  ggplot(aes(x = environment, y = r_hat, fill = model, group = model, ymin = CI_lower, ymax = CI_upper)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(position = position_dodge(1), width = 0.5) +
  facet_grid(facets = . ~ trait) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Prediction Accuracy") +
  xlab("Environment") +
  scale_fill_manual(values = umn_palette("Secondary_Tier1")) +
  labs(title = "G-BLUP Prediction Accuracy")



#####################


# Save everything
save(list = cfind(what = "blup2.*analysis"), 
     file = file.path(pred_dir, "Trial_Predictions/trial_prediction_results.RData"))




    