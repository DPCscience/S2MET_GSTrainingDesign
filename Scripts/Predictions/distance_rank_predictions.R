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


# Maximum number of training environments to test
max_env <- 15
max_env <- Inf # Use all available environments

## Prepare the BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate(line_name = as.factor(line_name))

# Data_frame of possible environments
possible_envs <- bind_rows(
  complete_train_env %>% data_frame(trait = names(.), envs = .),
  realistic_train_env %>% data_frame(trait = names(.), envs = .)
) %>% mutate(set = c(rep("complete", length(complete_train_env)), rep("realistic", length(realistic_train_env))))
  
  
# Use the clusters based on the TP-only data
environment_rank_model <- pred_env_dist_rank %>%
  # Make sure the pool of training environments is consistent with the environments that are available
  left_join(., possible_envs) %>% mutate(rank = map2(.x = rank, .y = envs, ~intersect(.x, .y))) %>%
  filter(!mat_set %in% c("Jarquin", "MalosettiStand")) %>% # Filter out the relationship matrices
  # filter(environment %in% sample_envs) %>%
  # Remove some covariate models
  rename(training_environment = rank) %>%
  select(-envs)

# Combine with the random clusters
environment_rank_random <- pred_env_rank_random %>%
  # Make sure the pool of training environments is consistent with the environments that are available
  left_join(., possible_envs) %>% mutate(rank = map2(.x = rank, .y = envs, ~intersect(.x, .y))) %>%
  rename(training_environment = rank) %>%
  select(-envs)

environment_rank_tomodel <- bind_rows(environment_rank_model, environment_rank_random) %>%
  mutate(training_environment = map(training_environment, ~head(., max_env))) %>%  filter(!str_detect(model, "sample"))

# Split by cores
environment_rank_tomodel_split <- environment_rank_tomodel %>%
  assign_cores(n_core = n_core) %>%
  split(.$core)


## Should the fixed effect of environment be fitted?
fit_env <- FALSE


# ## Local machine
# ##
# environment_rank_pred_out <- environment_rank_tomodel %>%
#   group_by(mat_set, set, trait, validation_environment, model) %>%
#   do({
#     # df <- clusters_to_model %>% filter_at(vars(validation_environment, trait, model), all_vars(. == .[1]))
#     df <- .
    

### Comment below for local machine

# Parallelize
environment_rank_pred_out <- mclapply(X = environment_rank_tomodel_split, FUN = function(core_df) {
  # core_df <- environment_rank_tomodel_split[[1]]

  results_out <- vector("list", nrow(core_df))

  for (i in seq_along(results_out)) {
    df <- core_df[i,]

    ### Comment above here for local machine
    
    # Vector of training environments, passed to an accumulation function
    envs <- df$training_environment[[1]] %>% accumulate(., c)
    val_env <- df$validation_environment
    tr <- df$trait

    # Create a list of training data
    train_data_list <- envs %>%
      map(~filter(S2_MET_BLUEs_tomodel, trait == tr, environment %in% ., line_name %in% tp_geno)) %>%
      map(~mutate(., environment = as.factor(environment)))

    test_data <- S2_MET_BLUEs_tomodel %>%
      filter(line_name %in% vp_geno, environment == val_env, trait == tr) %>%
      mutate(line_name = as.character(line_name))
    
    if (!fit_env) {
    
      ## Calculate fixed effects, then predict
      # list of formulae
      form_list <- train_data_list %>%
        map(~{if (n_distinct(.$environment) == 1) formula(value ~ line_name) else formula(value ~ line_name + environment)})
  
      blues_list <- map2(.x = train_data_list, .y = form_list, ~lm(.y, data = .x)) %>%
        map(~effects::Effect("line_name", .) %>% as.data.frame()) %>%
        map(~select(., line_name, value = fit) %>% 
              mutate(std_error = 0, environment = "mean", line_name = factor(line_name, levels = c(tp_geno, vp_geno))))
      
      pred_list_fixed <- map(blues_list, ~gblup(K = K, train = ., test = test_data, fit.env = fit_env))
      

    } else {
    
      # Fit a mixed model that calculates the fixed effect of environment and random effect of genotype
      blues_list <- train_data_list
  
      # Iterate over the list and predict
      pred_list_fixed <- map(blues_list, ~gblup(K = K, train = ., test = test_data, fit.env = fit_env))

    }
    
  #   # Return results
  #   data.frame(n_e = map_dbl(envs, length), accuracy = map_dbl(pred_list_fixed, "accuracy"))
  # 
  # 
  # }) %>% ungroup()
    
    results_out[[i]] <- data.frame(n_e = map_dbl(envs, length), accuracy = map_dbl(pred_list_fixed, "accuracy"))

  }

  core_df %>%
    mutate(out = results_out) %>%
    select(-core)

}, mc.cores = n_core)




#### Window-based calculation ####

# Window size
env_window <- 5


# ## Local machine
# environment_window_pred_out <- environment_rank_tomodel %>%
#   group_by(set, trait, validation_environment, model) %>%
#   do({
#     # df <- clusters_to_model %>% filter_at(vars(validation_environment, trait, model), all_vars(. == .[1]))
#     df <- .

# Parallelize
environment_window_pred_out <- mclapply(X = environment_rank_tomodel_split, FUN = function(core_df) {
  # core_df <- environment_rank_tomodel_split[[1]]

  results_out <- vector("list", nrow(core_df))

  for (i in seq_along(results_out)) {

    df <- core_df[i,]
    
    # Vector of training environments, passed to an accumulation function
    envs <- df$training_environment[[1]]
    val_env <- df$validation_environment
    envs <- seq(env_window, length(envs)) %>%
      seq_along() %>%
      map(~. + seq(env_window) - 1) %>%
      map(~envs[.])
    
    tr <- df$trait
    
    # Create a list of training data
    train_data_list <- envs %>%
      map(~filter(S2_MET_BLUEs_tomodel, trait == tr, environment %in% ., line_name %in% tp_geno)) %>%
      map(~mutate(., environment = as.factor(environment)))
    
    test_data <- S2_MET_BLUEs_tomodel %>%
      filter(line_name %in% vp_geno, environment == val_env, trait == tr) %>%
      mutate(line_name = as.character(line_name))
    
    if (!fit_env) {
      
      ## Calculate fixed effects, then predict
      # list of formulae
      form_list <- train_data_list %>%
        map(~{if (n_distinct(.$environment) == 1) formula(value ~ line_name) else formula(value ~ line_name + environment)})
      
      blues_list <- map2(.x = train_data_list, .y = form_list, ~lm(.y, data = .x)) %>%
        map(~effects::Effect("line_name", .) %>% as.data.frame()) %>%
        map(~select(., line_name, value = fit) %>% 
              mutate(std_error = 0, environment = "mean", line_name = factor(line_name, levels = c(tp_geno, vp_geno))))
      
      pred_list_fixed <- map(blues_list, ~gblup(K = K, train = ., test = test_data, fit.env = fit_env))
      
      
    } else {
      
      # Fit a mixed model that calculates the fixed effect of environment and random effect of genotype
      blues_list <- train_data_list
      
      # Iterate over the list and predict
      pred_list_fixed <- map(blues_list, ~gblup(K = K, train = ., test = test_data, fit.env = fit_env))
      
    }
    
    # # Return results
    # data.frame(window = seq_along(envs), accuracy = map_dbl(pred_list_fixed, "accuracy"))
    # 
    # }) %>% ungroup()
    
    results_out[[i]] <- data.frame(n_e = map_dbl(envs, length), accuracy = map_dbl(pred_list_fixed, "accuracy"))
    
    
  }

  core_df %>%
    mutate(out = results_out) %>%
    select(-core)

}, mc.cores = n_core)
    
    

# Save the results
save_file <- file.path(result_dir, "distance_rank_predictions.RData")
save("environment_rank_pred_out", "environment_window_pred_out", file = save_file)






















# ### Appendix
# 
# 
# ## Test alternative prediction methods
# # 1. Test the normal method using REML (Fixed effects and prediction)
# # 2. Fixed effects, then prediction
# # 3. Scale, fixed effects, then prediction
# cluter_pred_testing <- clusters_to_model %>%
#   sample_n(20) %>%
#   group_by(environment, trait, model) %>%
#   do({
#     # df <- clusters_to_model %>% filter_at(vars(environment, trait, model), all_vars(. == .[1]))
#     df <- .
#     
#     # Vector of training environments, passed to an accumulation function
#     envs <- df$pred_environment[[1]] %>% accumulate(., c) %>% .[1:10]
#     val_env <- df$environment
#     tr <- df$trait
#     
#     # Create a list of training data
#     train_data_list <- envs %>% 
#       map(~filter(S2_MET_BLUEs_tomodel, trait == tr, environment %in% ., line_name %in% tp_geno)) %>%
#       map(~mutate(., env = as.factor(environment)))
#     
#     test_data <- S2_MET_BLUEs_tomodel %>%
#       filter(line_name %in% vp_geno, environment == val_env, trait == tr) %>%
#       mutate(line_name = as.character(line_name))
#     
#     ## First predict normally
#     pred_list_normal <- train_data_list %>% map(~gblup(K = K, train = ., test = test_data, fit.env = F))
#     
#     ## Calculate fixed effects, then predict
#     # list of formulae
#     form_list <- train_data_list %>% 
#       map(~if (n_distinct(.$environment) == 1) formula(value ~ -1 + line_name) else formula(value ~ -1 + line_name + environment))
#     
#     blues_list <- map2(.x = train_data_list, .y = form_list, ~lm(.y, data = .x, contrasts = list(environment = "contr.sum"))) %>% 
#       map(coef) %>% 
#       map(~data_frame(line_name = names(.), value = .) %>%
#             filter(str_detect(line_name, "line_name")) %>%
#             mutate(line_name = factor(str_remove_all(line_name, "line_name"), levels = c(tp_geno, vp_geno)),
#                    env = "mean", std_error = 0))
#     
#     # Iterate over the list and predict
#     pred_list_fixed <- blues_list %>% map(~gblup(K = K, train = ., test = test_data, fit.env = F))
#     
#     ## Scale, calculate fixed effects, then predict
#     blues_list <- train_data_list %>%
#       map(~group_by(., environment) %>% mutate(value = scale(value)) %>% ungroup()) %>% 
#       map(~lm(value ~ -1 + line_name, data = .)) %>% 
#       map(coef) %>% 
#       map(~data_frame(line_name = names(.), value = .) %>%
#             filter(str_detect(line_name, "line_name")) %>%
#             mutate(line_name = factor(str_remove_all(line_name, "line_name"), levels = c(tp_geno, vp_geno)),
#                    env = "mean", std_error = 0))
#     
#     pred_list_scaled <- blues_list %>% map(~gblup(K = K, train = ., test = test_data, fit.env = F))
#     
#     ## Compare accuracies
#     data_frame(normal = pred_list_normal, fixed = pred_list_fixed, scaled = pred_list_scaled) %>% 
#       mutate_all(~map_dbl(., "accuracy"))
#     
#   })
# 
# ## Compare
# cluter_pred_testing %>% 
#   summarize_at(vars(fixed, scaled), funs(cor(normal, .)))
# 
# 
# 
# 
# 
# 
# 
