## Analysis of predictions based on environmental distance
## 
## Author: Jeff Neyhart
## Last Updated: April 6, 2018
## 
## This script will look at prediction accuracies from adding environments after
## ranking based on specific distance metrics
## 


# # Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load the results
load(file.path(result_dir, "environmental_distance_predictions.RData"))
load(file.path(result_dir, "environmental_distance_heritability.RData"))

################################################################################
######### Remove this when the new results come it

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
  assign_cores(n_core = 16) %>%
  split(.$core)


################################################################################




### This code will need to be revised, as I will look at an initial version of
### the results that included errors.
# Which elements of the result list are useable?
tokeep <- map_lgl(env_dist_predictions_out, ~!inherits(., "try-error"))
env_dist_predictions_use <- env_dist_predictions_out[tokeep]

i <- which(tokeep)[2]

test_results <- pred_env_dist_rank_split[[i]] %>%
  select(environment:dist_method) %>%
  mutate(results_out = env_dist_predictions_out[[i]]) %>%
  unnest() %>%
  group_by(environment) %>%
  filter(n_distinct(dist_method) == 6) %>% 
  ungroup()

# Re-scale the distance measurements to unit variance
test_results_toplot <- test_results %>%
  group_by(environment, dist_method) %>% 
  mutate(scaled_distance = scale(distance)) %>%
  ungroup()
  

# Plot
test_results_toplot %>% 
  ggplot(aes(x = scaled_distance, y = accuracy, color = dist_method)) +
  geom_point() +
  geom_line() +
  facet_grid(environment + trait ~ ., scales = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank())





### Heritability across different ranks of environments

# Bind rows and unnest
env_dist_heritability <- env_dist_heritability_out %>%
  bind_rows() %>% 
  unnest()

# Scale the distance measurements by environment and distance method
# Then subtract the minimum to have the values start at 0
env_dist_heritability1 <- env_dist_heritability %>%
  group_by(environment, trait, dist_method) %>%
  mutate(scaled_distance = scale(distance),
         scaled_distance = scaled_distance - min(scaled_distance)) %>%
  ungroup()

# Plot 5 environments
g_distance_herit <- env_dist_heritability1 %>% 
  filter(environment %in% sample(unique(.$environment), 5), trait == "GrainYield") %>%
  ggplot(aes(x = scaled_distance, y = heritability_out, color = dist_method)) + 
  geom_point() +
  geom_line() + 
  facet_grid(environment + trait ~ .)












