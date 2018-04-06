## Predictions based on environmental distance
## 
## Author: Jeff Neyhart
## Last Updated: April 5, 2018
## 
## This script will look at prediction accuracies from adding environments after
## ranking based on specific distance metrics
## 


### Run for MSI

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"
source(file.path(repo_dir, "source_MSI.R"))

# Load packages
packages <- c("lme4")
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))

# Source some scripts from pbr
source("/panfs/roc/groups/6/smithkp/neyha001/R/my_packages/pbr/R/convenience_functions.R")
source("/panfs/roc/groups/6/smithkp/neyha001/R/my_packages/pbr/R/herit.R")



### Run on a local machine

# # Run the source script
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))


## Number of cores
n_core <- ifelse(Sys.info()["sysname"] == "Windows", 1, detectCores())



# Load the clustering results
load(file.path(result_dir, "distance_methods_results.RData"))

# Modify the BLUEs for predictions
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate_at(vars(environment:line_name), as.factor)


### Question 1 - What is the trend in prediction accuracy as you add increasingly
### "distant" environment?
### 

## For each prediction environment (the tp+vp envs and just the vp envs), rank the 
## training environments by different distance metrics
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


## Rank the environments relative to each other
# Extract the "all" cluster metrics
dist_method_df_all <- clust_method_df %>%
  filter(population == "all") %>%
  select(-cluster) %>%
  mutate(dist = map(dist, as.matrix))

# Linearize the distance objects and rank
dist_method_df_rank <- dist_method_df_all %>% 
  mutate(dist_rank = map(dist, ~as.data.frame(.) %>% 
                           rownames_to_column("env") %>% 
                           filter(env %in% pred_envs) %>%
                           split(.$env) %>% 
                           map(~.[,!names(.) %in% .$env] %>% .[,-1] %>% sort() %>% 
                                 select(., which(names(.) %in% train_envs))) %>% 
                           data_frame(environment = names(.), env_rank = .))) %>%
  select(-dist) %>%
  unnest()


# Combine this data with the the trait information for the prediction environments
pred_env_dist_rank <- pred_envs_traits %>%
  left_join(., dist_method_df_rank,  by = c("environment", "trait"))



## Split the 'pred_env_dist_rank' data.frame by core and then pipe to mclapply
pred_env_dist_rank_split <- pred_env_dist_rank %>% 
  assign_cores(n_core = n_core) %>%
  split(.$core)

# Run over multiple cores
env_dist_predictions_out <- mclapply(X = pred_env_dist_rank_split, FUN = function(core_df) {

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
    
    # Subset the trainin environments for the trait
    tr_train_env <- train_envs_traits %>% filter(trait == tr) %>% pull(environment)

    # Subset the prediction environment data from the BLUEs
    pred_env_data <- S2_MET_BLUEs_use %>%
      filter(environment == pred_env, trait == tr, line_name %in% vp_geno)

    # Remove the environments not in the training environments for that trait
    sorted_train_envs <- core_df$env_rank[[i]] %>% unlist() %>% .[names(.) %in% tr_train_env]

    ## Use the accumulate function in dplyr to create lists of training environments
    ## and the cumulative mean distance from the prediction environment
    train_envs_accumulate <- sorted_train_envs %>%
      names() %>%
      accumulate(~c(.x, .y)) %>%
      map(~sorted_train_envs[.]) %>%
      map(~list(train_envs = names(.), cummean_dist = mean(.)))

    # Map over these environments and gather training data
    train_envs_data <- train_envs_accumulate %>%
      map(~filter(S2_MET_BLUEs_use, environment %in% .$train_envs, trait == tr, line_name %in% tp_geno))

    # Map over the data and predict
    predictions_out <- train_envs_data %>%
      map(~rename(., env = environment)) %>%
      map(~gblup(K = K, train = ., test = pred_env_data, bootreps = 100))

    # Return the bootstrap data.frame results
    predictions_boot <- predictions_out %>%
      map_df("boot") %>%
      rename(accuracy = cor)

    # Return a data_frame with the training environments, cumumative mean distance,
    # and the prediction accuracy results
    results_out[[i]] <- bind_cols(
      data_frame(train_envs = map(train_envs_accumulate, "train_envs"),
                 distance = map_dbl(train_envs_accumulate, "cummean_dist")),
      predictions_boot)

  } # Close the for loop
  
  # Add the results list to the original core DF
  core_df %>%
    select(-env_rank, -core) %>%
    mutate(results_out = results_out)

}, mc.cores = n_core)


# Save the results
save_file <- file.path(result_dir, "environmental_distance_predictions.RData")
save("env_dist_predictions_out", file = save_file)





### Use the distance objects to rank training environments relevative to prediction 
### environments and calculate the heritability within a "cluster"
### 
### 

# # Run over multiple cores
# env_dist_heritability_out <- mclapply(X = pred_env_dist_rank_split, FUN = function(core_df) {
#   
#   ## Create an empty results list
#   results_out <- vector("list", nrow(core_df))
#   
#   # Iterate over the rows in the core_df
#   for (i in seq_along(results_out)) {
#     
#     ## Given an ordered list of environments sorted by distance to the prediction
#     ## environment, implement a function that calculates a cumulative mean distance
#     ## and then use the information from those environments to run predictions
#     pred_env <- core_df$environment[i]
#     tr <- core_df$trait[i]
#     dist_met <- core_df$dist_method[i]
#     
#     sorted_train_envs <- core_df$env_rank[[i]] %>% unlist()
#     
#     ## Use the accumulate function in dplyr to create lists of training environments
#     ## and the cumulative mean distance from the prediction environment
#     train_envs_accumulate <- sorted_train_envs %>% 
#       names() %>% 
#       accumulate(~c(.x, .y)) %>% 
#       map(~sorted_train_envs[.]) %>% 
#       map(~list(train_envs = names(.), cummean_dist = mean(.)))
#     
#     # Map over these environments and gather training data
#     train_envs_data <- train_envs_accumulate %>%
#       map(~filter(S2_MET_BLUEs_use, environment %in% .$train_envs, trait == tr, line_name %in% tp_geno))
#     
#     # lmer control
#     lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore",
#                                 calc.derivs = FALSE)
#     
#     # Map over the data and calculate heritability
#     heritability_out <- train_envs_data %>%
#       map_dbl(function(pheno_tomodel) {
#         
#         # Extract the weights
#         wts <- pheno_tomodel$std_error^2
#         
#         # Determine the formula based on the number of environments
#         if (n_distinct(pheno_tomodel$environment) == 1) {
#           form <- value ~ (1|line_name) 
#         } else {
#           form <- value ~ (1|line_name) + (1|environment) + (1|line_name:environment)
#         }
#         
#         # Fit the model
#         fit <- lmer(formula = form, data = pheno_tomodel, weights = wts, control = lmer_control)
#         # Calculate heritability and return
#         training_heritability(object = fit)
#         
#       })
#     
#     # Return a data_frame with the training environments, cumumative mean distance,
#     # and the prediction accuracy results
#     results_out[[i]] <- data_frame(
#       train_envs = map(train_envs_accumulate, "train_envs"), 
#       distance = map_dbl(train_envs_accumulate, "cummean_dist"),
#       heritability_out)
#     
#   }
#   
#   # Combine the results with the original DF and return
#   core_df %>% 
#     select(-env_rank, -core) %>% 
#     mutate(results_out = results_out)
#   
# }, mc.cores = n_core)
# 
# 
# # Save the results
# save_file <- file.path(result_dir, "environmental_distance_heritability.RData")
# save("env_dist_heritability_out", file = save_file)













