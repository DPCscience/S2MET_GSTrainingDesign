## Cluster environments in the S2MET project
## 
## Author: Jeff Neyhart
## Last modified: March 27, 2018
## 
## This script will create different distance matrices for use in clustering. The
## clustering algorithm will be constant
## 

# Load packages and the source script
library(tidyverse)
library(stringr)
library(readxl)
library(pbr)

# The head directory
repo_dir <- getwd()

source(file.path(repo_dir, "source.R"))


## We will test 6 different distance measurements:
## 
## 1. Great Circle Distance using lat/long
## 2. Genetic correlation of the TP
## 3. The environment distance metric (as described in Bernardo2010)
## 4. Factor analytic or PCA of GxE term
## 5. 1 year environmental covariables
## 6. 10 year environmental covariables
## 


# Create a new data.frame to hold the different datasets
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  group_by(trait) %>%
  nest() %>%
  mutate(data = map(data, droplevels))


## Great Circle Distance
# Use the geosphere package
# Subset the lat/long of the environment
trial_lat_long <- trial_info %>% 
  select(environment, longitude, latitude) %>% 
  distinct() %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>% # No NA's allowed
  as.data.frame() %>%
  column_to_rownames("environment") %>%
  as.matrix()


# Create pairwise combinations of environments and calculate the great circle distance
trial_great_circle <- row.names(trial_lat_long) %>%
  expand.grid(environment1 = ., environment2 = ., stringsAsFactors = FALSE) %>%
  mutate(dist = geosphere::distGeo(p1 = trial_lat_long[environment1,], p2 = trial_lat_long[environment2,]))
  

# Convert to a dist object
great_circle_dist <- trial_great_circle %>% 
  spread(environment2, dist) %>% 
  remove_rownames() %>% 
  column_to_rownames("environment1") %>%
  as.matrix() %>% 
  as.dist()

# Copy the list per trait (the distance is the same for each trait)
great_circle_dist_list <- rerun(length(traits), great_circle_dist) %>% set_names(traits)



## Genetic correlation of the TP or all lines
# Load the environmental genetic correlation matrices
load(file.path(result_dir, "environmental_genetic_correlations.RData"))

# Convert each correlation matrix to a distance matrix
env_cor_all_dist_list <- env_cor_all %>% map(as.dist)
env_cor_tp_dist_list <- env_cor_tp %>% map(as.dist)


## Environmental distance metric
# This is described by Bernardo2010

# First calculate the distance using all data
ge_mean_D_all <- S2_MET_BLUEs %>%
  # Split by trait
  split(.$trait) %>%
  map(function(df){
    dist1 <- dist_env(x = df, gen.col = "line_name", env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })

# Now calculate it based only on data from the TP
ge_mean_D_tp <- S2_MET_BLUEs %>%
  filter(line_name %in% tp_geno) %>%
  # Split by trait
  split(.$trait) %>%
  map(function(df){
    dist1 <- dist_env(x = df, gen.col = "line_name", env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })



## PCA of GxE BLUPs

# Fit a new model with compound symmetry and get the GxE BLUPs
S2_MET_BLUPs_ge_all <- S2_MET_BLUEs %>%
  split(.$trait) %>%
  map(function(df) {

    # Create an object for the data.frame
    ## Center and scale the data
    df <- mutate(df, value = scale(value))

    # Extract the standard errors
    wts <- df$std_error^2

    # lmer control
    lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore",
                                calc.derivs = FALSE)

    # Fit the model
    form <- value ~ (1|line_name) + (1|environment) + (1|line_name:environment)
    fit <- lmer(formula = form, data = df, control = lmer_control, weights = wts)

    # Extract the model.frame
    mf <- model.frame(fit)

    # Extract the BLUPs of each genotype in each environment
    blups <- ranef(fit)$`line_name:environment` %>%
      rownames_to_column("term") %>%
      separate(term, c("line_name", "environment"), sep = ":") %>%
      rename(value = "(Intercept)")

    # Extract the variance components
    var_comp <- as.data.frame(VarCorr(fit))

    # Return a data.frame
    data_frame(fit = list(fit), var_comp = list(var_comp), BLUP = list(blups)) })


# Fit a new model with compound symmetry and get the GxE BLUPs
# Use only the TP data
S2_MET_BLUPs_ge_tp <- S2_MET_BLUEs %>%
  filter(line_name %in% tp_geno) %>%
  split(.$trait) %>%
  map(function(df) {
    
    # Create an object for the data.frame
    ## Center and scale the data
    df <- mutate(df, value = scale(value))
    
    # Extract the standard errors
    wts <- df$std_error^2
    
    # lmer control
    lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore",
                                calc.derivs = FALSE)
    
    # Fit the model
    form <- value ~ (1|line_name) + (1|environment) + (1|line_name:environment)
    fit <- lmer(formula = form, data = df, control = lmer_control, weights = wts)
    
    # Extract the model.frame
    mf <- model.frame(fit)
    
    # Extract the BLUPs of each genotype in each environment
    blups <- ranef(fit)$`line_name:environment` %>%
      rownames_to_column("term") %>%
      separate(term, c("line_name", "environment"), sep = ":") %>%
      rename(value = "(Intercept)")
    
    # Extract the variance components
    var_comp <- as.data.frame(VarCorr(fit))
    
    # Return a data.frame
    data_frame(fit = list(fit), var_comp = list(var_comp), BLUP = list(blups)) })

# # Factor analysis using BLUEs
# ge_mean_FA <- S2_MET_BLUEs_use %>%
#   group_by(trait) %>%
#   mutate(BLUEs = list({
#     data[[1]] %>%
#       select(line_name, environment, value) %>%
#       complete(line_name, environment) %>%
#       group_by(environment) %>%
#       mutate(value = ifelse(is.na(value), mean(value, na.rm = TRUE), value)) %>%
#       spread(environment, value) %>%
#       as.data.frame() %>%
#       remove_rownames() %>%
#       column_to_rownames("line_name") %>%
#       as.matrix() })) %>%
#   do(fa_out = fa(r = .$BLUEs[[1]], nfactors = 2, rotate = "varimax")) %>%
#   # Grab the loadings
#   mutate(delta = list(structure(fa_out$loadings, class = "matrix"))) %>%
#   # Distance matrix
#   mutate(FA = list(dist(delta)))


# Number of PCs to choose for GxE BLUPs and ECs
n_PC <- 2

# Extract the GxE BLUPs and run PCA
ge_PCA_dist_all <- S2_MET_BLUPs_ge_all %>%
  map(function(out) {
    blup_mat <- out$BLUP[[1]] %>%
      spread(environment, value) %>%
      as.data.frame() %>%
      remove_rownames() %>%
      column_to_rownames("line_name") %>%
      as.matrix() 
    
    # Impute with mean
    blup_prcomp <- blup_mat %>%
      apply(MARGIN = 2, FUN = function(env) ifelse(is.na(env), mean(env, na.rm = T), env)) %>%
      t() %>%
      prcomp(center = TRUE, scale. = TRUE)
    
    # Extract all the PCs
    blup_pcs <- blup_prcomp$x
    # Subset the ones to use
    blup_pcs_use <- blup_pcs[,1:n_PC,drop = FALSE]
    # Calculate the distance between points using the PCs
    blup_pc_dist <- dist(blup_pcs_use)
    
    # Return a data.frame
    data_frame(blup_mat = list(blup_mat), blup_prcomp = list(blup_prcomp), 
               blup_pcs = list(blup_pcs), blup_pc_dist = list(blup_pc_dist))
    
  })
    
# Copy for the tp
ge_PCA_dist_tp <- S2_MET_BLUPs_ge_tp %>%
  map(function(out) {
    blup_mat <- out$BLUP[[1]] %>%
      spread(environment, value) %>%
      as.data.frame() %>%
      remove_rownames() %>%
      column_to_rownames("line_name") %>%
      as.matrix() 
    
    # Impute with mean
    blup_prcomp <- blup_mat %>%
      apply(MARGIN = 2, FUN = function(env) ifelse(is.na(env), mean(env, na.rm = T), env)) %>%
      t() %>%
      prcomp(center = TRUE, scale. = TRUE)
    
    # Extract all the PCs
    blup_pcs <- blup_prcomp$x
    # Subset the ones to use
    blup_pcs_use <- blup_pcs[,1:n_PC,drop = FALSE]
    # Calculate the distance between points using the PCs
    blup_pc_dist <- dist(blup_pcs_use)
    
    # Return a data.frame
    data_frame(blup_mat = list(blup_mat), blup_prcomp = list(blup_prcomp), 
               blup_pcs = list(blup_pcs), blup_pc_dist = list(blup_pc_dist))
    
  })



## Environmental covariates
# Use the PCs of the environmental covariates to form clusters
EC_one_PCA <- prcomp(x = one_year_env_mat, center = TRUE, scale. = TRUE)
EC_multi_PCA <- prcomp(x = multi_year_env_mat, center = TRUE, scale. = TRUE)

# Copy the list per trait (the distance is the same for each trait)
EC_one_PCA_dist_list <- rerun(length(traits), dist(EC_one_PCA$x[,1:n_PC, drop = FALSE])) %>% 
  set_names(traits)
EC_multi_PCA_dist_list <- rerun(length(traits), dist(EC_multi_PCA$x[,1:n_PC, drop = FALSE])) %>% 
  set_names(traits)



## Combine the lists into a data.frame
dist_method_df_all <- data_frame(
  trait = traits, 
  great_circle_dist = great_circle_dist_list, 
  env_cor_dist = env_cor_all_dist_list, 
  ge_mean_D = ge_mean_D_all, 
  ge_PCA_dist = ge_PCA_dist_all,
  ec_one_PCA_dist = EC_one_PCA_dist_list,
  ec_multi_PCA_dist = EC_multi_PCA_dist_list)


dist_method_df_tp <- data_frame(
  trait = traits, 
  great_circle_dist = great_circle_dist_list, 
  env_cor_dist = env_cor_tp_dist_list, 
  ge_mean_D = ge_mean_D_tp, 
  ge_PCA_dist = ge_PCA_dist_tp,
  ec_one_PCA_dist = EC_one_PCA_dist_list,
  ec_multi_PCA_dist = EC_multi_PCA_dist_list)


# Save this
save_file <- file.path(result_dir, "distance_methods_results.RData")
save("dist_method_df_all", "dist_method_df_tp", file = save_file)





