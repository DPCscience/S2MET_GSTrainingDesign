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
max_env <- Inf # Use all available environments

## Prepare the BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate(line_name = as.factor(line_name))


# First use the TP clusters
clusters_model <- pred_env_dist_rank$tp %>%
  # filter(environment %in% sample_envs) %>%
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
  # filter(!str_detect(model, "sample")) %>%
  tbl_df()

# Split by cores
clusters_to_model_split <- clusters_to_model %>%
  assign_cores(n_core = n_core) %>%
  split(.$core)


# # Parallelize
# cluster_pred_out <- mclapply(X = clusters_to_model_split, FUN = function(core_df) {
#   
#   # ##
#   i = 1
#   core_df <- clusters_to_model_split[[i]]
#   # ##
#   
#   # Results list
#   results_out <- vector("list", nrow(core_df))
#   
#   # Iterate over rows
#   for (i in seq(nrow(core_df))) {
#     
#     # Vector of training environments, passed to an accumulation function
#     envs <- core_df$pred_environment[[i]] %>%
#       accumulate(., c)
#     val_env <- core_df$environment[i]
#     tr <- core_df$trait[i]
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
#     blues_list <- train_data_list %>% map(~lm(value ~ -1 + line_name, data = .)) %>% 
#       map(coef) %>% 
#       map(~data_frame(line_name = names(.), value = .) %>% 
#             mutate(line_name = factor(str_remove_all(line_name, "line_name"), levels = c(tp_geno, vp_geno)),
#                    env = "mean", std_error = 0))
#     
#     # Iterate over the list and predict
#     pred_list <- blues_list[1:10] %>% map(~gblup(K = K, train = ., test = test_data, fit.env = F))
#     # pred_list1 <- train_data_list[1:10] %>% map(~gblup(K = K, train = ., test = test_data, fit.env = F))
#     
#     
#     # Add to the results
#     results_out[[i]] <- data.frame(n_e = map_dbl(envs, length), accuracy = map_dbl(pred_list, "accuracy"))
# 
#   }
#   
#   # Add the results list to the core_df and return
#   core_df %>% 
#     mutate(out = results_out) %>% 
#     select(-pred_environment, -core)
#   
# })
# 
# 
# cluster_pred_out <- bind_rows(cluster_pred_out)




# ## Local machine
# cluster_pred_out <- clusters_to_model %>%
#   group_by(environment, trait, model) %>%
#   do({
#     # df <- clusters_to_model %>% filter_at(vars(environment, trait, model), all_vars(. == .[1]))
#     df <- .
    
# Parallelize
cluster_pred_out <- mclapply(X = clusters_to_model_split, FUN = function(core_df) {
  # core_df <- clusters_to_model_split[[1]]
  
  results_out <- vector("list", nrow(core_df))
  
  for (i in seq_along(results_out)) {
    
    df <- core_df[i,]

    # Vector of training environments, passed to an accumulation function
    envs <- df$pred_environment[[1]] %>% accumulate(., c)
    val_env <- df$environment
    tr <- df$trait

    # Create a list of training data
    train_data_list <- envs %>%
      map(~filter(S2_MET_BLUEs_tomodel, trait == tr, environment %in% ., line_name %in% tp_geno)) %>%
      map(~mutate(., env = as.factor(environment)))

    test_data <- S2_MET_BLUEs_tomodel %>%
      filter(line_name %in% vp_geno, environment == val_env, trait == tr) %>%
      mutate(line_name = as.character(line_name))
    
    ## Calculate fixed effects, then predict
    # list of formulae
    form_list <- train_data_list %>%
      map(~if (n_distinct(.$environment) == 1) formula(value ~ -1 + line_name) else formula(value ~ -1 + line_name + environment))

    blues_list <- map2(.x = train_data_list, .y = form_list, ~lm(.y, data = .x, contrasts = list(environment = "contr.sum"))) %>%
      map(coef) %>%
      map(~data_frame(line_name = names(.), value = .) %>%
            filter(str_detect(line_name, "line_name")) %>%
            mutate(line_name = factor(str_remove_all(line_name, "line_name"), levels = c(tp_geno, vp_geno)),
                   env = "mean", std_error = 0))

    # Iterate over the list and predict
    pred_list_fixed <- map(blues_list, ~gblup(K = K, train = ., test = test_data, fit.env = F))
    
    # Return results
    # data.frame(n_e = map_dbl(envs, length), accuracy = map_dbl(pred_list_fixed, "accuracy"))
    results_out[[i]] <- data.frame(n_e = map_dbl(envs, length), accuracy = map_dbl(pred_list_fixed, "accuracy"))

    
  # })
    
  }
  
  core_df %>%
    mutate(out = results_out) %>%
    select(-core)
  
}, mc.cores = n_core)


# cluster_pred_out <- ungroup(cluster_pred_out)
cluster_pred_out <- bind_rows(cluster_pred_out)





#### Window-based calculation ####

# Window size
env_window <- 5



# # Parallelize
# cluster_pred_out_window <- mclapply(X = clusters_to_model_split, FUN = function(core_df) {
# 
#   # ##
#   # i = 1
#   # core_df <- clusters_to_model_split[[i]]
#   # ##
# 
#   # Results list
#   results_out <- vector("list", nrow(core_df))
# 
#   # Iterate over rows
#   for (i in seq(nrow(core_df))) {
# 
#     # Vector of training environments, passed to an accumulation function
#     envs <- core_df$pred_environment[[i]]
#     val_env <- core_df$environment[i]
    # envs <- seq(env_window, length(envs)) %>%
    #   seq_along() %>%
    #   map(~. + seq(env_window) - 1) %>%
    #   map(~envs[.])
# 
# 
#     tr <- core_df$trait[i]
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
#     # Iterate over the list and predict
#     pred_list <- train_data_list %>%
#       map(~gblup(K = K, train = ., test = test_data))
# 
#     # Add to the results
#     results_out[[i]] <- data.frame(window = seq_along(envs), accuracy = map_dbl(pred_list, "accuracy"))
# 
#   }
# 
#   # Add the results list to the core_df and return
#   core_df %>%
#     mutate(out = results_out) %>%
#     select(-pred_environment, -core)
# 
# })
# 
# 
# cluster_pred_out_window <- bind_rows(cluster_pred_out_window)



# ## Local machine
# cluster_pred_out_window <- clusters_to_model %>%
#   group_by(environment, trait, model) %>%
#   do({
#     df <- .

# Parallelize
cluster_pred_out_window <- mclapply(X = clusters_to_model_split, FUN = function(core_df) {
  # core_df <- clusters_to_model_split[[1]]
  
  results_out <- vector("list", nrow(core_df))
  
  for (i in seq_along(results_out)) {
    
    df <- core_df[i,]
    
    # Vector of training environments, passed to an accumulation function
    envs <- df$pred_environment[[1]]
    val_env <- df$environment
    envs <- seq(env_window, length(envs)) %>%
      seq_along() %>%
      map(~. + seq(env_window) - 1) %>%
      map(~envs[.])
    
    tr <- df$trait
    
    # Create a list of training data
    train_data_list <- envs %>%
      map(~filter(S2_MET_BLUEs_tomodel, trait == tr, environment %in% ., line_name %in% tp_geno)) %>%
      map(~mutate(., env = as.factor(environment)))
    
    test_data <- S2_MET_BLUEs_tomodel %>%
      filter(line_name %in% vp_geno, environment == val_env, trait == tr) %>%
      mutate(line_name = as.character(line_name))
    
    ## Calculate fixed effects, then predict
    # list of formulae
    form_list <- train_data_list %>%
      map(~if (n_distinct(.$environment) == 1) formula(value ~ -1 + line_name) else formula(value ~ -1 + line_name + environment))
    
    blues_list <- map2(.x = train_data_list, .y = form_list, ~lm(.y, data = .x, contrasts = list(environment = "contr.sum"))) %>%
      map(coef) %>%
      map(~data_frame(line_name = names(.), value = .) %>%
            filter(str_detect(line_name, "line_name")) %>%
            mutate(line_name = factor(str_remove_all(line_name, "line_name"), levels = c(tp_geno, vp_geno)),
                   env = "mean", std_error = 0))
    
    # Iterate over the list and predict
    pred_list_fixed <- map(blues_list, ~gblup(K = K, train = ., test = test_data, fit.env = F))
    
    # Return results
    # data.frame(n_e = map_dbl(envs, length), accuracy = map_dbl(pred_list_fixed, "accuracy"))
    results_out[[i]] <- data.frame(n_e = map_dbl(envs, length), accuracy = map_dbl(pred_list_fixed, "accuracy"))
    
    
    # })
    
  }
  
  core_df %>%
    mutate(out = results_out) %>%
    select(-core)
  
}, mc.cores = n_core)
    
    
# cluster_pred_out_window <- ungroup(cluster_pred_out_window)
cluster_pred_out_window <- bind_rows(cluster_pred_out_window)


# Save the results
save_file <- file.path(result_dir, "cluster_predictions_tp.RData")
save("cluster_pred_out_window", "cluster_pred_out", file = save_file)






















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
