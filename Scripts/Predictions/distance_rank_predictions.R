## Predictions based on environmental distance
## 
## Author: Jeff Neyhart
## Last Updated: October 11, 2018
## 
## This script will look at prediction accuracies from adding environments after
## ranking based on specific distance metrics
## 


### Run for MSI

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"
source(file.path(repo_dir, "source_MSI.R"))



# ## Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))


## Number of cores
n_core <- 16
n_core <- detectCores()



# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))



### Testing Version ###


# Maximum number of training environments to test
max_env <- 15

## Prepare the BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate(line_name = as.factor(line_name))


## What is the heritability of each "training cluster" of environments as you add increasingly distant environments?

# First use the TP clusters
clusters_model <- pred_env_dist_rank$tp %>%
  filter(environment %in% sample_envs) %>%
  # Remove some covariate models
  filter(!str_detect(model, "Fstat|Top1")) %>%
  mutate(pred_environment = map(env_rank, names)) %>%
  select(-env_rank)

# Combine with the random clusters
clusters_rand <- pred_env_rank_random$tp %>%
  rename(pred_environment = env_rank) %>%
  mutate(pred_environment = map(pred_environment, names))

clusters_to_model <- bind_rows(clusters_model, clusters_rand) %>%
  mutate(pred_environment = map(pred_environment, ~head(., max_env))) %>%
  filter(!str_detect(model, "sample")) %>%
  tbl_df()


# Split by cores
clusters_to_model_split <- clusters_to_model %>%
  assign_cores(n_core = n_core) %>%
  split(.$core)


# Parallelize
cluster_pred_out <- mclapply(X = clusters_to_model_split, FUN = function(core_df) {
  
  # ##
  # i = 1
  # core_df <- clusters_to_model_split[[i]]
  # ##
  
  # Results list
  results_out <- vector("list", nrow(core_df))
  
  # Iterate over rows
  for (i in seq(nrow(core_df))) {
    
    # Vector of training environments, passed to an accumulation function
    envs <- core_df$pred_environment[[i]] %>%
      accumulate(., c)
    val_env <- core_df$environment[i]
    tr <- core_df$trait[i]
    
    # Create a list of training data
    train_data_list <- envs %>% 
      map(~filter(S2_MET_BLUEs_tomodel, trait == tr, environment %in% ., line_name %in% tp_geno)) %>%
      map(~mutate(., env = as.factor(environment)))
    
    test_data <- S2_MET_BLUEs_tomodel %>%
      filter(line_name %in% vp_geno, environment == val_env, trait == tr) %>%
      mutate(line_name = as.character(line_name))
    
    
    # Iterate over the list and predict
    pred_list <- train_data_list %>% map(~gblup(K = K, train = ., test = test_data))
    
    # Add to the results
    results_out[[i]] <- data.frame(n_e = map_dbl(envs, length), accuracy = map_dbl(pred_list, "accuracy"))

  }
  
  # Add the results list to the core_df and return
  core_df %>% 
    mutate(out = results_out) %>% 
    select(-pred_environment, -core)
  
})


cluster_pred_out <- bind_rows(cluster_pred_out)



#### Window-based calculation ####

# Window size
env_window <- 5 



# Parallelize
cluster_pred_out_window <- mclapply(X = clusters_to_model_split, FUN = function(core_df) {
  
  # ##
  # i = 1
  # core_df <- clusters_to_model_split[[i]]
  # ##
  
  # Results list
  results_out <- vector("list", nrow(core_df))
  
  # Iterate over rows
  for (i in seq(nrow(core_df))) {
    
    # Vector of training environments, passed to an accumulation function
    envs <- core_df$pred_environment[[i]]
    val_env <- core_df$environment[i]
    envs <- seq(env_window, length(envs)) %>% 
      seq_along() %>% 
      map(~. + seq(env_window) - 1) %>%
      map(~envs[.])
    
    
    tr <- core_df$trait[i]
    
    # Create a list of training data
    train_data_list <- envs %>% 
      map(~filter(S2_MET_BLUEs_tomodel, trait == tr, environment %in% ., line_name %in% tp_geno)) %>%
      map(~mutate(., env = as.factor(environment)))
    
    test_data <- S2_MET_BLUEs_tomodel %>%
      filter(line_name %in% vp_geno, environment == val_env, trait == tr) %>%
      mutate(line_name = as.character(line_name))
    
    # Iterate over the list and predict
    pred_list <- train_data_list %>%
      map(~gblup(K = K, train = ., test = test_data))
    
    # Add to the results
    results_out[[i]] <- data.frame(window = seq_along(envs), accuracy = map_dbl(pred_list, "accuracy"))
    
  }
  
  # Add the results list to the core_df and return
  core_df %>% 
    mutate(out = results_out) %>% 
    select(-pred_environment, -core)
  
})


cluster_pred_out_window <- bind_rows(cluster_pred_out_window)


# Save the results
save_file <- file.path(result_dir, "cluster_predictions_tp.RData")
save("cluster_pred_out_window", "cluster_pred_out", file = save_file)




