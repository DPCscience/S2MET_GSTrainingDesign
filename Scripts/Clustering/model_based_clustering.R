## Cluster environments in the S2MET project
## Model-based clustering strategies
## 
## Author: Jeff Neyhart
## Last modified: April 18, 2018
## 
## This script will explore model-based strategies to group environments based on
## (hopefully) difference distance metrics
## 


### Run for MSI

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"
source(file.path(repo_dir, "source_MSI.R"))

# Load packages
packages <- c("lme4")
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


# # The head directory
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# 
# library(lme4)

## Number of cores
n_core <- ifelse(Sys.info()["sysname"] == "Windows", 1, detectCores())
# Define a significance threshold for the LRT
alpha <- 0.05


## Load the environment covariable data
load(file.path(data_dir, "environmental_data_compiled.RData"))
load(file.path(result_dir, "distance_methods_results.RData"))


# Modify the BLUEs for predictions
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate_at(vars(environment:line_name), as.factor)

# Vectors of training and prediction environments
pred_envs <- c(tp_vp_env, vp_only_env)
train_envs <- c(tp_vp_env, tp_only_env)


# Summarize the traits available in those environments
pred_envs_traits <- S2_MET_BLUEs_use %>%
  filter(environment %in% pred_envs) %>% 
  group_by(environment) %>% 
  distinct(trait) %>%
  ungroup()

train_envs_traits <- S2_MET_BLUEs_use %>%
  filter(environment %in% train_envs) %>% 
  group_by(environment) %>% 
  distinct(trait) %>%
  ungroup()


## Use the distances calculated using 'all' data
pred_env_dist_rank <- pred_env_rank_random$all %>%
  # Use only the first 50 samples
  filter(!dist_method %in% str_c("sample", 51:100))


## Split the 'pred_env_rank_random' data.frame by core and then pipe to mclapply
pred_env_dist_rank_split <- pred_env_dist_rank %>% 
  assign_cores(n_core = n_core) %>%
  split(.$core)



# Run over multiple cores
env_dist_lrt_predictions_out <- mclapply(X = pred_env_dist_rank_split, FUN = function(core_df) {

  ## Create an empty results list
  results_out <- vector("list", nrow(core_df))

  # Iterate over the rows in the core_df
  for (i in seq_along(results_out)) {

    ## Given an ordered list of environments sorted by distance to the prediction
    ## environment, implement a function that calculates a cumulative mean distance
    ## and then use the information from those environments to run predictions
    pred_env <- core_df$environment[i]
    tr <- core_df$trait[i]
    dist_met <- core_df$dist_method[i]
    
    # Print a message
    print(paste(pred_env, tr, dist_met))

    # Subset the trainin environments for the trait
    tr_train_env <- train_envs_traits %>% filter(trait == tr) %>% pull(environment)

    # Remove the environments not in the training environments for that trait
    sorted_train_envs <- core_df$env_rank[[i]] %>% unlist() %>% .[names(.) %in% tr_train_env]
    
    # Subset the BLUEs for the traits and the training individuals
    train_data_full <- S2_MET_BLUEs_use %>%
      filter(environment %in% names(sorted_train_envs),
             trait == tr,
             line_name %in% tp_geno) %>%
      rename(., env = environment)

    
    # Model control
    control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
    
    ## Define the formula and fit the full model
    ## Use the last element in the above list, because it has all the data
    form_full <- value ~ 1 + env + (1|line_name) + (1|line_name:env)
    fit_full <- lmer(formula = form_full, data = train_data_full, control = control)
    

    ## Map over the data and fit the "fuller" models. Stop when the pvalue for the LRT
    ## reaches a local minimum
    
    # First create an empty list to store the output of the lrt
    lrt_list <- vector("list", length(sorted_train_envs))
    
    # Define the expanded model
    form_exp <- value ~ 1 + env + (1|line_name) + (1|line_name:env) + (1|line_name:cluster) + (1|line_name:env:cluster)
    
    # Setup the while loop conditions
    lrt_sig <- FALSE
    lrt_opt <- FALSE
    j <- 1
    
    # Iterate using a while loop
    # While lrt_sig or lrt_opt is false AND j <= length(sorted_train_envs)
    while (any(!lrt_sig, !lrt_opt) & j <= length(sorted_train_envs)) {
      # Get a vector of environments to be clustered
      clust_env <- names(sorted_train_envs[1:j])
      
      # Assign the environment(s) to a cluster
      train_data_full_clust <- train_data_full %>% 
        mutate(cluster = ifelse(env %in% clust_env, "train", "not_train"))
      
      ## Fit the model
      fit_exp <- lmer(formula = form_exp, data = train_data_full_clust, control = control)
      
      ## Run the likelihood ratio test
      lrt_out <- lr_test(model1 = fit_full, model2 = fit_exp)
      lrt_list[[j]] <- cbind(n_env = j, lrt_out)
      
      # Is the LRT significant?
      lrt_sig <- lrt_out$p_value <= alpha
      
      # If j is one, continue
      if (j != 1) {
        # Is the jth pvalue >= than the j - 1th pvalue
        lrt_opt <- lrt_out$p_value >= lrt_list[[j - 1]]$p_value
      }
      
      # Advance the counter
      j <- j + 1
    } # Close the while loop
    
    # Bind the rows of the list and return
    results_out[[i]] <- bind_rows(lrt_list)
    
  } # Close the for loop
  
  core_df1 <- core_df %>%
    select(-env_rank, -core) %>%
    mutate(results_out = results_out)
  
  # Return the results
  return(core_df1)
  
}, mc.cores = n_core)

# Save the results
save_file <- file.path(result_dir, "environmental_distance_lrt.RData")
save("env_dist_lrt_predictions_out", file = save_file)

