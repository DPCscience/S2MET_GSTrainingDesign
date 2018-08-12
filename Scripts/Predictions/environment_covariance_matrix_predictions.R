## Predictions using different environmental covariance matrices
## TESTING VERSION
## 
## Author: Jeff Neyhart
## Last Updated: August 12, 2018
## 
## This script will look at prediction accuracies from reaction norm models of GxE using different
## environmental covariance matrices. We will use LOEO prediction.
## 


### Run on MSI
# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"
source(file.path(repo_dir, "source_MSI.R"))

# # Load packages
# packages <- c("lme4")
# invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))

# # Source some scripts from pbr
# source("/panfs/roc/groups/6/smithkp/neyha001/R/my_packages/pbr/R/convenience_functions.R")
# source("/panfs/roc/groups/6/smithkp/neyha001/R/my_packages/pbr/R/herit.R")

### Run on a local machine

# # Run the source script
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))


## Number of cores
n_core <- ifelse(Sys.info()["sysname"] == "Windows", 1, detectCores())


# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))

# Modify the BLUEs for predictions
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate_at(vars(environment:line_name), as.factor)


## Construct different covariance matrices using distance, covariates, or phenotypic correlation
env_cov_mats <- clust_method_df %>% 
  filter(population == "all") %>% 
  mutate(env_cov_mat = map(dist, as.matrix)) %>% 
  select(trait, model, env_cov_mat)



## How accurate are predictions of each of the testing environments when using LOEO prediction and 
## different environmental covariance matrices?

# Load the testing environment sampling data
load(file.path(result_dir, "environmental_distance_predictions_TESTING.RData"))


# Create empty dfs for testing
loeo_predictions_df <- env_dist_predictions_out %>% 
  bind_rows() %>% 
  filter(!str_detect(model, "sample")) %>% 
  distinct(environment, trait, model) %>%
  left_join(., env_cov_mats) %>%
  mutate(predictions = list(NULL))

## Split by core
loeo_predictions_split <- loeo_predictions_df %>%
  assign_cores(n_core = n_core) %>% 
  split(.$core)

# Use mclapply to parallelize
loeo_predictions_out <- mclapply(X = loeo_predictions_split, FUN = function(core_df) {
  
  for (i in seq(nrow(core_df))) {
  
    pred_env <- core_df$environment[i]
    tr <- core_df$trait[i]
    mod <- core_df$model[i]
    
    # Extract the training set
    train_df <- S2_MET_BLUEs_use %>% 
      filter(trait == tr) %>%
      droplevels() %>%
      filter(environment != pred_env, line_name %in% tp_geno)
    
    # Create a model frame
    mf <- model.frame(value ~ line_name + environment, train_df)
    y <- model.response(mf)
    
    # Predict using the GxE as the only random effect
    Z_ge <- model.matrix(~ -1 + line_name:environment, mf)
    
    ## Create the GxE covariance matrix
    E <- core_df$env_cov_mat[[i]][levels(train_df$environment), levels(train_df$environment)]
    E <- cor(E)
    # E <- scale(E, center = rep(mean(E), ncol(E)), scale = rep(sd(E), ncol(E)))
    
    K_ge <- kronecker(E, K, make.dimnames = T)
    
    
    # Fit the model
    fit_ge <- mixed.solve(y = y, Z = Z_ge, K = K_ge)
    
    ## Extract the random effects
    preds <- fit_ge$u %>% data_frame(term = names(.), pred_value = .) %>% 
      separate(term, c("environment", "line_name"), sep = ":")
    
  
    # Get the testing set and combine the predictions
    validate_df <- S2_MET_BLUEs_use %>% 
      filter(trait == tr, environment == pred_env) %>%
      filter(line_name %in% vp_geno) %>%
      left_join(., preds)
    
    # Add this to a df
    core_df$predictions[[i]] <- validate_df
    
  }
  
  # Return the core_df
  select(core_df, -core)
  
}, mc.cores = n_core)







## Perform predictions using the genotype mean
## Split by core
loeo_g_predictions_split <- loeo_predictions_df %>%
  distinct(environment, trait) %>%
  assign_cores(n_core = n_core) %>% 
  split(.$core)



# Use mclapply to parallelize
loeo_g_predictions_out <- mclapply(X = loeo_g_predictions_split, FUN = function(core_df) {
  
  for (i in seq(nrow(core_df))) {
    
    pred_env <- core_df$environment[i]
    tr <- core_df$trait[i]
    mod <- core_df$model[i]
    
    # Extract the training set
    train_df <- S2_MET_BLUEs_use %>% 
      filter(trait == tr) %>%
      droplevels() %>%
      filter(environment != pred_env, line_name %in% tp_geno)
    
    # Create a model frame
    mf <- model.frame(value ~ line_name + environment, train_df)
    y <- model.response(mf)
    
    # Predict using the GxE as the only random effect
    # Z_ge <- model.matrix(~ -1 + line_name:environment, mf)
    Z_g <- model.matrix(~ -1 + line_name, mf)
  
    
    # Fit the model
    # fit_ge <- mixed.solve(y = y, Z = Z_ge, K = K_ge)
    fit_g <- mixed.solve(y = y, Z = Z_g, K = K)
    
    ## Extract the random effects
    ## Extract the random effects
    preds <- fit_g$u %>% 
      data_frame(line_name = names(.), pred_value = .)
    
    
    # Get the testing set and combine the predictions
    validate_df <- S2_MET_BLUEs_use %>% 
      filter(trait == tr, environment == pred_env) %>%
      filter(line_name %in% vp_geno) %>%
      left_join(., preds)
    
    # Add this to a df
    core_df$predictions[[i]] <- validate_df
    
  }
  
  # Return the core_df
  select(core_df, -core)
  
}, mc.cores = n_core)



## Combine the results
predictions_out <- list(
  ge = bind_rows(loeo_predictions_out),
  g = bind_rows(loeo_g_predictions_out)
)


# Save
save_file <- file.path(result_dir, "env_cov_mat_predictions.RData")
save("predictions_out", file = save_file)












