## Heritability based on environmental distance
## 
## Author: Jeff Neyhart
## Last Updated: October 9, 2018
## 
## This script will look at prediction accuracies from adding environments after
## ranking based on specific distance metrics
## 
## TESTING VERSION
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








### Testing Version ###


# Maximum number of training environments to test
max_env <- 15


## What is the heritability of each "training cluster" of environments as you add increasingly distant environments?

# First use the TP clusters
clusters_model <- pred_env_dist_rank$tp %>%
  filter(environment %in% sample_envs) %>%
  mutate(pred_environment = map(env_rank, names)) %>%
  select(-env_rank)

# Combine with the random clusters - only 25 of them
clusters_rand <- pred_env_rank_random$tp %>%
  filter(model %in% paste0("sample", 1:25),
         environment %in% sample_envs) %>%
  mutate(pred_environment = map(env_rank, names)) %>%
  select(-env_rank)

clusters_to_model <- bind_rows(clusters_model, clusters_rand) %>%
  mutate(pred_environment = map(pred_environment, ~head(., max_env)))



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
    envs <- core_df$pred_environment[[i]] %>%
      accumulate(., c)
    tr <- core_df$trait[i]
    
    # Create a list of training data
    # Remove the first (only one environment)
    train_data_list <- envs[-1] %>% 
      map(~filter(S2_MET_BLUEs_tp, trait == tr, environment %in% .))
    
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
    select(-pred_environment, -core)
  
})


cluster_herit_out <- bind_rows(cluster_herit_out)

# Save
save_file <- file.path(result_dir, "cluster_heritability_tp_TESTING.RData")
save("cluster_herit_out", file = save_file)

