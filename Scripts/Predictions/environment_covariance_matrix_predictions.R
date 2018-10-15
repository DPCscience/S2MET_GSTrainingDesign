## Predictions using different environmental covariance matrices
## 
## Author: Jeff Neyhart
## Last Updated: October 10, 2018
## 
## This script will look at prediction accuracies from reaction norm models of GxE using different
## environmental covariance matrices. We will use LOEO and 60-40 environment predictions.
## 


### Run on MSI
# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"
source(file.path(repo_dir, "source_MSI.R"))


# # Run the source script
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# library(modelr)


# Number of cores
n_core <- 16
n_core <- detectCores()


# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))

# Modify the BLUEs for predictions
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate_at(vars(environment:line_name), as.factor)


## Construct different covariance matrices using distance, covariates, or phenotypic correlation
env_cov_mats <- clust_method_df %>% 
  filter(population == "tp", !map_lgl(cov, is.null)) %>% 
  filter(!str_detect(model, "Cor|Fstat")) %>% # For now, just focus on phenotypic distance and all covariates
  select(trait, model, env_cov_mat = cov) %>%
  arrange(trait, model)

## Subset the matrices for environments in which both the TP and VP were phenotyped
env_cov_mats_use <- env_cov_mats %>%
  mutate(env_cov_mat = map(env_cov_mat, ~.[row.names(.) %in% tp_vp_env, colnames(.) %in% tp_vp_env]),
         # Create the GxE covariance matrices
         env_cov_mat = map(env_cov_mat, ~kronecker(., K, make.dimnames = T)))



## Run different proportions of training environments
prop_train_env <- c(0.25, 0.5, 0.75)
n_iter <- 10

# Generate training and test sets
environment_mc_samples <- S2_MET_BLUEs_use %>%
  filter(environment %in% tp_vp_env) %>%
  group_by(trait) %>%
  do({
    df <- .
    # Number the rows
    df1 <- mutate(df, row = seq(nrow(df))) %>%
      droplevels()
    envs <- distinct(df1, environment)
    
    # Generate environment samples
    samples <- data_frame(pTrainEnv = prop_train_env) %>%
      mutate(envSample = map(pTrainEnv, ~rerun(.n = n_iter, sample_frac(tbl = envs, size = .)))) %>% 
      unnest() %>% 
      group_by(pTrainEnv) %>% 
      mutate(iter = seq(n_iter)) %>%
      ungroup()
    
    # Generate the samples
    samples %>% 
      mutate(train = map(envSample, ~left_join(., df1, by = "environment") %>% filter(line_name %in% tp_geno)), 
             test = map(train, ~setdiff(df1, .) %>% filter(line_name %in% vp_geno))) %>% 
      mutate_at(vars(train, test), ~map(., ~pull(., row) %>% resample(data = df1, .)))
    
  }) %>% ungroup()
    
  



## Run predictions
prediction_model_split <- environment_mc_samples %>%
  assign_cores(n_core) %>%
  split(.$core)



# Use mclapply to parallelize
environment_mc_predictions  <- mclapply(X = prediction_model_split, FUN = function(core_df) {
  
  # # For local machine
  # i <- 7
  # core_df <- prediction_model_split[[i]]
  # #
  
  results_out <- vector("list", nrow(core_df))
  
  
  # Iterate over rows
  for (i in seq_along(results_out)) {
    
    # Extract the covariance matrix
    K_mats <- env_cov_mats_use %>%
      filter(trait == core_df$trait[i]) %>%
      pull(env_cov_mat)
    
    # Create model matrices
    mf <- model.frame(value ~ line_name + environment, core_df$train[[i]])
    y <- model.response(mf)
    Z <- model.matrix(~ -1 + line_name:environment, mf)
    X <- model.matrix(~ 1 + environment, droplevels(mf))
    
    # Fit the model
    pgvs <- K_mats %>%
      map(~mixed.solve(y = y, Z = Z, X = X, K = .)) %>%
      map("u") %>%
      # Combine the PGVs with the observations
      map(~data.frame(term = names(.), pgv = ., row.names = NULL, stringsAsFactors = FALSE) %>%
            separate(term, c("environment", "line_name"), sep = ":") %>% 
            left_join(as.data.frame(core_df$test[[i]]), ., by = c("environment", "line_name")) %>%
            select(trait, environment, line_name, value, pgv))
    
    K_mat_pgv <- env_cov_mats_use %>% 
      filter(trait == core_df$trait[i]) %>% 
      mutate(predictions = pgvs) %>% 
      select(-env_cov_mat)
    
    
    
    ## Predict just using the mean
    Zg <- model.matrix(~ -1 + line_name, mf)
    pgv_mean <- mixed.solve(y = y, Z = Zg, X = X, K = K)$u %>% 
      data.frame(line_name = names(.), pgv = ., row.names = NULL, stringsAsFactors = FALSE) %>% 
      left_join(as.data.frame(core_df$test[[i]]), ., by = c("line_name")) %>%
      select(trait, environment, line_name, value, pgv)
    
    # Add to the K_mat df

    # Add to the K list
    results_out[[i]] <- K_mat_pgv %>% 
      add_row(trait = unique(.$trait), model = "mean", predictions = list(pgv_mean))
    
  } 
  
  core_df %>%
    mutate(out = results_out) %>%
    select(trait, pTrainEnv, iter, out)
  
})
    
environment_mc_predictions <- bind_rows(environment_mc_predictions)






### Leave-one-environment-out
# Generate training and test sets
environment_loeo_samples <- S2_MET_BLUEs_use %>%
  filter(environment %in% tp_vp_env) %>%
  group_by(trait) %>%
  do({
    df <- .
    # Number the rows
    df1 <- mutate(df, row = seq(nrow(df))) %>%
      droplevels()
    envs <- as.character(unique(df1$environment))
    
    # Generate environment samples
    samples <- data_frame(testEnv = envs) %>%
      mutate(train = map(testEnv, ~filter(df1, environment != ., line_name %in% tp_geno)),
             test = map(testEnv, ~filter(df1, environment == ., line_name %in% vp_geno))) %>%
      mutate_at(vars(train, test), ~map(., ~pull(., row) %>% resample(data = df1, .)))

    
  }) %>% ungroup()




## Run predictions
prediction_model_split <- environment_loeo_samples %>%
  assign_cores(n_core) %>%
  split(.$core)

# Use mclapply to parallelize
environment_loeo_predictions  <- mclapply(X = prediction_model_split, FUN = function(core_df) {
  
  # # For local machine
  # i <- 7
  # core_df <- prediction_model_split[[i]]
  # #
  
  results_out <- vector("list", nrow(core_df))
  
  
  # Iterate over rows
  for (i in seq_along(results_out)) {
    
    # Extract the covariance matrix
    K_mats <- env_cov_mats_use %>%
      filter(trait == core_df$trait[i]) %>%
      pull(env_cov_mat)
    
    # Create model matrices
    mf <- model.frame(value ~ line_name + environment, core_df$train[[i]])
    y <- model.response(mf)
    Z <- model.matrix(~ -1 + line_name:environment, mf)
    X <- model.matrix(~ 1 + environment, droplevels(mf))
    
    # Fit the model
    pgvs <- K_mats %>%
      map(~mixed.solve(y = y, Z = Z, X = X, K = .)) %>%
      map("u") %>%
      # Combine the PGVs with the observations
      map(~data.frame(term = names(.), pgv = ., row.names = NULL, stringsAsFactors = FALSE) %>%
            separate(term, c("environment", "line_name"), sep = ":") %>% 
            left_join(as.data.frame(core_df$test[[i]]), ., by = c("environment", "line_name")) %>%
            select(trait, environment, line_name, value, pgv))
    
    K_mat_pgv <- env_cov_mats_use %>% 
      filter(trait == core_df$trait[i]) %>% 
      mutate(predictions = pgvs) %>% 
      select(-env_cov_mat)
    
    
    
    ## Predict just using the mean
    Zg <- model.matrix(~ -1 + line_name, mf)
    pgv_mean <- mixed.solve(y = y, Z = Zg, X = X, K = K)$u %>% 
      data.frame(line_name = names(.), pgv = ., row.names = NULL, stringsAsFactors = FALSE) %>% 
      left_join(as.data.frame(core_df$test[[i]]), ., by = c("line_name")) %>%
      select(trait, environment, line_name, value, pgv)
    
    # Add to the K_mat df
    
    # Add to the K list
    results_out[[i]] <- K_mat_pgv %>% 
      add_row(trait = unique(.$trait), model = "mean", predictions = list(pgv_mean))
    
  } 
  
  core_df %>%
    mutate(out = results_out) %>%
    select(trait, pTrainEnv, iter, out)
  
})

environment_loeo_predictions <- bind_rows(environment_loeo_predictions)




# Save
save_file <- file.path(result_dir, "env_cov_mat_predictions.RData")
save("environment_loeo_predictions", "environment_mc_predictions", file = save_file)












