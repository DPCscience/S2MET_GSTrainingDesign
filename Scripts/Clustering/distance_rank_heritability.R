## Heritability based on environmental distance
## 
## Author: Jeff Neyhart
## Last Updated: October 9, 2018
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

## Subset only the TP
S2_MET_BLUEs_tp <- S2_MET_BLUEs %>%
  filter(line_name %in% tp) %>%
  droplevels()



## What is the heritability of each "training cluster" of environments as you add increasingly distant environments?

# First use the TP clusters
clusters_model <- pred_env_dist_rank$tp %>%
  mutate(environments = map(env_rank, names)) %>%
  select(-env_rank)

# Combine with the random clusters - only 25 of them
clusters_rand <- pred_env_rank_random$tp %>%
  filter(str_detect(model, paste0("sample", 1:25))) %>%
  mutate(environments = map(env_rank, names)) %>%
  select(-env_rank)

clusters_to_model <- bind_rows(clusters_model, clusters_rand)



# Model control
control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
# Heritability expression
herit_exp <- "line_name / (line_name + (line_name:environment / n_e) + (Residual / (n_e * n_r)))"


# Split by cores
clusters_to_model_split <- clusters_to_model %>%
  assign_cores(n_core = n_core) %>%
  split(.$core)

# Parallelize
cluster_herit_out <- mclapply(X = clusters_to_model_split, FUN = function(core_df) {
  
  # ##
  # i = 11
  # core_df <- clusters_to_model_split[[i]]
  # ##
  
  # Results list
  results_out <- vector("list", nrow(core_df))
  
  # Iterate over rows
  for (i in seq(nrow(core_df))) {
    
    # Vector of training environments, passed to an accumulation function
    envs <- core_df$environments[[i]] %>%
      accumulate(., c)
    tr <- core_df$trait[i]
    
    # Create a list of training data
    # Remove the first (only one environment)
    train_data_list <- envs[-1] %>% map(~filter(S2_MET_BLUEs_tp, trait == tr, environment %in% .))
    
    # Iterate over the list and fit a model
    fit_list <- train_data_list %>%
      map(~{
        wts <- .$std_error^2
        fit <- lmer(value ~ (1|line_name) + (1|environment) + (1|line_name:environment), data = ., control = control, weights = wts)
        
        plot_table <- xtabs(formula = ~ line_name + environment, data = .)
        
        harm_env <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
          ifelse(. > 1, 1, .) %>%
          rowSums() %>% 
          harm_mean()
        
        list(fit = fit, n_e = harm_env)
      })
    
    # Calculate heritability
    results_out[[i]] <- fit_list %>%
      map(~herit(object = .$fit, exp = herit_exp, n_e = .$n_e, n_r = 1))
    
  }
  
  # Add the results list to the core_df and return
  core_df %>% 
    mutate(out = results_out) %>% 
    select(-environments, -core)
  
})


cluster_herit_out <- bind_rows(cluster_herit_out)

# Save
save_file <- file.path(result_dir, "cluster_heritability_tp.RData")
save("cluster_herit_out", file = save_file)
  
  



























#### OLD CODE ####


# ## What is the change in the ratio of genetic variance to GxE variance as you add increasingly distant environments?
# 
# ## Load the distance prediction results
# ## Testing
# load(file.path(result_dir, "environmental_distance_predictions_TESTING.RData"))
# 
# ## Cumulative predictions
# # For each trait/method/environment, find the data.frame of training environments
# cumulative_predictions_df <- env_dist_predictions_out %>% 
#   bind_rows() %>% 
#   unnest()
# 
# 
# ## Look at GxE only in the training set
# training_set_var_comp <- cumulative_predictions_df %>%
#   filter(!str_detect(model, "sample")) %>%
#   mutate(n_train_env = map_dbl(train_envs, length)) %>%
#   mutate(var_comp_out = list(NULL))
# 
# # Map over the df
# for (i in seq(nrow(training_set_var_comp))) {
#   
#   # Get the training environment(s)
#   train_env <- training_set_var_comp$train_envs[[i]]
#   tr <- training_set_var_comp$trait[i]
#   
#   # Need at least 2 envs for GxE
#   if (length(train_env) < 2) {
#     training_set_var_comp$var_comp_out[[i]] <- as_data_frame(NULL)
#     
#   } else {
#     
#     # Get the data to model
#     pheno_to_model <- S2_MET_BLUEs_tp %>%
#       filter(environment %in% train_env, trait == tr)
#     
#     # Fit the model
#     control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
#     fit <- lmer(value ~ (1|line_name) + environment + (1|line_name:environment), data = pheno_to_model, control = control)
# 
#     # Heritability
#     h2 <- herit(object = fit, exp = "line_name / (line_name + (line_name:environment / ne) + (Residual / ne))", ne = training_set_var_comp$n_train_env[i])
#     # Find the proportion of G to total
#     propG <- subset(h2$var_comp, source == "line_name", variance, drop = T) / sum(subset(h2$var_comp, , variance, drop = T))
#     
#       
#     # add to the df
#     training_set_var_comp$var_comp_out[[i]] <- data_frame(h2 = h2$heritability, propG = propG)
#     
#   }
#   
# }
# 
# 
# training_set_var_comp_df <- training_set_var_comp %>% 
#   mutate(var_comp_out = map(var_comp_out, as_data_frame)) %>%
#   unnest(var_comp_out)
#   
# ## Plot propG versus accuracy
# training_set_var_comp_df %>%
#   filter(model == "pheno_dist") %>%
#   ggplot(aes(x = n_train_env)) + 
#   geom_line(aes(y = accuracy, color = "accurary")) + 
#   geom_line(aes(y = propG, color = "propG")) + 
#   facet_grid(trait ~ environment) +
#   theme_bw()
#   
# ## Plot h2 versus accuracy
# training_set_var_comp_df %>%
#   # filter(model == "pheno_dist") %>%
#   filter(trait == "GrainYield") %>%
#   ggplot(aes(x = n_train_env)) + 
#   geom_line(aes(y = accuracy, color = "accurary")) + 
#   geom_line(aes(y = h2, color = "h2")) + 
#   facet_grid(trait + model ~ environment) +
#   theme_bw()
# 
# 
# ## Compare these relationships - average over environments
# training_set_var_comp_df %>% 
#   group_by(trait, model, environment) %>% 
#   summarize_at(vars(h2, propG), funs(cor = cor(., accuracy))) %>%
#   ggplot(aes(x = model)) + 
#   geom_boxplot(aes(y = h2_cor, fill = "h2")) +
#   geom_boxplot(aes(y = propG_cor, fill = "propG")) +
#   facet_grid(~ trait)
# 
# 
# 
# 
# 
# ## For each prediction environment (the tp+vp envs and just the vp envs), rank the 
# ## training environments by different distance metrics
# pred_envs <- c(tp_vp_env, vp_only_env)
# train_envs <- c(tp_vp_env, tp_only_env)
# 
# # Summarize the traits available in those environments
# pred_envs_traits <- S2_MET_BLUEs_use %>%
#   filter(environment %in% pred_envs) %>% 
#   group_by(environment) %>% 
#   distinct(trait) %>%
#   ungroup()
# 
# train_envs_traits <- S2_MET_BLUEs_use %>%
#   filter(environment %in% train_envs) %>% 
#   group_by(environment) %>% 
#   distinct(trait) %>%
#   ungroup()
# 
# 
# ## Use the distances calculated using 'all' data
# pred_env_dist_rank <- pred_env_rank_random$all %>%
#   # Use only the first 50 samples
#   filter(!dist_method %in% str_c("sample", 51:100))
# 
# 
# ## Split the 'pred_env_dist_rank' data.frame by core and then pipe to mclapply
# pred_env_dist_rank_split <- pred_env_dist_rank %>% 
#   assign_cores(n_core = n_core) %>%
#   split(.$core)
# 
# 
# 
# 
# ### Use the distance objects to rank training environments relevative to prediction 
# ### environments and calculate the heritability within a "cluster"
# ### 
# ### 
# 
# # Run over multiple cores
# env_dist_heritability_out <- mclapply(X = pred_env_dist_rank_split, FUN = function(core_df) {
# 
# #core_df <- pred_env_dist_rank_split[[1]]
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
#     # Subset the trainin environments for the trait
#     tr_train_env <- train_envs_traits %>% filter(trait == tr) %>% pull(environment)
# 
# 
#     # Remove the environments not in the training environments for that trait
#     sorted_train_envs <- core_df$env_rank[[i]] %>% unlist() %>% .[names(.) %in% tr_train_env]
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
# 
# 
# ### Use the distance objects to rank training environments relevative to prediction 
# ### environments and calculate the heritability within a fixed sliding window of 
# ### environments
# ### 
# 
# # Number of environments in the sliding window
# window_size <- 5
# 
# 
# # Run over multiple cores
# env_dist_window_heritability_out <- mclapply(X = pred_env_dist_rank_split, FUN = function(core_df) {
#   
#   #core_df <- pred_env_dist_rank_split[[1]]
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
#     # Subset the trainin environments for the trait
#     tr_train_env <- train_envs_traits %>% filter(trait == tr) %>% pull(environment)
#     
#     
#     # Remove the environments not in the training environments for that trait
#     sorted_train_envs <- core_df$env_rank[[i]] %>% unlist() %>% .[names(.) %in% tr_train_env]
#     
#     # Create a list of length (n - k + 1) of the training environments in the sliding
#     # window
#     n <- length(sorted_train_envs) - window_size
#     
#     # Create the list
#     train_envs_accumulate <- seq(0, n) %>% 
#       map(~. + seq(window_size)) %>%
#       map(~sorted_train_envs[.]) %>%
#       map(~list(train_envs = names(.), mean_dist = mean(.)))
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
#         form <- value ~ (1|line_name) + (1|environment) + (1|line_name:environment)
#       
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
#       distance = map_dbl(train_envs_accumulate, "mean_dist"),
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
# save("env_dist_heritability_out", "env_dist_window_heritability_out", file = save_file)













