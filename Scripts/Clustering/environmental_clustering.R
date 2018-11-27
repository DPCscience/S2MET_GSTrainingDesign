## Cluster environments in the S2MET project
## 
## Author: Jeff Neyhart
## Last modified: October 04, 2018
## 
## This script will create different distance matrices for use in clustering. The
## clustering algorithm will be constant
## 

# Load packages and the source script
library(pbr)
library(ggdendro)
library(cowplot)

# The head directory
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))



## Load the environment covariable distance matrices
load(file.path(result_dir, "environmental_covariable_distance_mat.RData"))
# Load the genetic correlation estimates
load(file.path(result_dir, "environmental_genetic_correlations.RData"))

# # Create a new data.frame to hold the different datasets
# S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
#   group_by(trait) %>%
#   nest() %>%
#   mutate(data = map(data, droplevels))


traits <- unique(S2_MET_BLUEs$trait)


## Great Circle Distance
# Use the geosphere package
# Subset the lat/long of the environment
trial_lat_long <- trial_info %>% 
  select(environment, longitude, latitude) %>% 
  distinct() %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>% # No NA's allowed
  as.data.frame() %>%
  column_to_rownames("environment") %>%
  as.matrix()


# Create pairwise combinations of environments and calculate the great circle distance
trial_great_circle <- row.names(trial_lat_long) %>%
  expand.grid(environment1 = ., environment2 = ., stringsAsFactors = FALSE) %>%
  mutate(dist = geosphere::distGeo(p1 = trial_lat_long[environment1,], p2 = trial_lat_long[environment2,]))
  

# Convert to a dist object
great_circle_dist <- trial_great_circle %>% 
  spread(environment2, dist) %>% 
  remove_rownames() %>% 
  column_to_rownames("environment1") %>%
  as.matrix() %>% 
  as.dist()

# Copy the list per trait (the distance is the same for each trait)
great_circle_dist_list <- rerun(length(traits), great_circle_dist) %>% set_names(traits)




## Phenotypic distance metric
# This is described by Bernardo2010

# First calculate the distance based on each environment
ge_mean_D <- S2_MET_BLUEs %>%
  filter(line_name %in% tp) %>%
  # Split by trait
  split(.$trait) %>%
  map(function(df){
    dist1 <- dist_env(x = df, gen.col = "line_name", env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })



control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")

## Now calculate the line means in each location
location_BLUEs <- S2_MET_BLUEs %>%
  group_by(trait, line_name, location) %>%
  mutate(value = mean(value)) %>%
  ungroup()

## Correlations among locations
env_loc_cor_complete <- location_BLUEs %>%
  split(.$trait) %>% 
  map(~select(., line_name, environment, value) %>% spread(environment, value) %>% as.data.frame() %>% 
        column_to_rownames("line_name") %>% as.matrix() %>% cor(., use = "pairwise.complete.obs"))

# Calculate the distance between locations using all data
gl_mean_D_complete <- location_BLUEs %>%
  filter(line_name %in% tp) %>%
  # Split by trait
  split(.$trait) %>%
  map(function(df){
    dist1 <- dist_env(x = df, gen.col = "line_name", env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })

## Now calculate the line means in each location, excluding 2017 data
location_BLUEs <- S2_MET_BLUEs %>%
  filter(year != 2017) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup() %>%
  left_join(distinct(S2_MET_BLUEs, trait, environment, location, year), .) %>%
  filter(!is.na(line_name))


## Correlations among locations
env_loc_cor_realistic <- location_BLUEs %>%
  split(.$trait) %>% 
  map(~select(., line_name, environment, value) %>% spread(environment, value) %>% as.data.frame() %>% 
        column_to_rownames("line_name") %>% as.matrix() %>% cor(., use = "pairwise.complete.obs"))

# Calculate the distance between locations using all data
gl_mean_D_realistic <- location_BLUEs %>%
  filter(line_name %in% tp) %>%
  # Split by trait
  split(.$trait) %>%
  map(function(df){
    dist1 <- dist_env(x = df, gen.col = "line_name", env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })



## Combine the one-year and multi-year ECs
# Create 3 different sets of matrices based on the distance calculations
sim_mat_group <- bind_rows(
  mutate(sim_mat_df, mat_set = "Malosetti"),
  mutate(sim_mat_df1, mat_set = "Jarquin"),
  mutate(sim_mat_df2, mat_set = "MalosettiStand")
) %>% select(-mat)


ec_sim_mat_df1 <- sim_mat_group %>%
  mutate(model = str_replace_all(ec_group, "_", " ") %>% str_to_title() %>% abbreviate(2),
         model = str_c(model, group),
         dist = map(sim_mat, dist)) %>%
  select(mat_set, set, model, trait, cov = sim_mat, dist)


## Combine - scale such that greater values mean "closer" or more similar
dist_method_df_complete <- ec_sim_mat_df1 %>% 
  filter(set == "complete") %>% 
  add_row(trait = traits, model = "great_circle_dist", dist = great_circle_dist_list, cov = map(great_circle_dist_list, ~1 - as.matrix(.))) %>% 
  add_row(trait = traits, model = "pheno_dist", dist = ge_mean_D, cov = env_cor_all) %>% 
  add_row(trait = traits, model = "pheno_loc_dist", dist = gl_mean_D_complete, cov = env_loc_cor_complete) %>%
  mutate(set = "complete")
  

dist_method_df_realistic <- ec_sim_mat_df1 %>% 
  filter(set == "realistic") %>% 
  add_row(trait = traits, model = "great_circle_dist", dist = great_circle_dist_list, cov = map(great_circle_dist_list, ~1 - as.matrix(.))) %>% 
  add_row(trait = traits, model = "pheno_loc_dist", dist = gl_mean_D_realistic, cov = env_loc_cor_realistic) %>%
  mutate(set = "realistic")


## Combine
dist_method_df <- bind_rows(dist_method_df_complete, dist_method_df_realistic)



# Tidy the the distance matrix data.frame, then combine the data.frames for
# TP and TP + VP

# For each trait, identify the most common set of environments
common_env_all <- dist_method_df %>% 
  group_by(set, trait) %>% 
  summarize(common_env = list(map(cov, ~row.names(as.matrix(.))) %>% reduce(intersect)))

# Combine with the distance metrics, convert to a matrix, subset the environments,
# then convert back to a dist object
# Then create cluster objects
clust_method_df <- full_join(dist_method_df, common_env_all, by = c("trait", "set")) %>% 
  mutate(cov = map2(.x = cov, .y = common_env, ~.x[.y,.y]),
         dist = map2(.x = dist, .y = common_env, ~subset_env(dist = .x, envs = .y)),
         cluster = map(dist, ~hclust(., method = "ward.D"))) %>%
  select(-common_env)


## For each model, plot the MDS of the distance matrix and the clustering
clust_method_df_toplot <- clust_method_df %>%
  group_by(trait, model, population) %>%
  nest(dist, cluster) %>%
  ## For each row, create a DF of the mds of the distance object
  mutate(data = map(data, ~{
    # Create the data.frame of MDS coordinates
    dist_df <- cmdscale(.$dist[[1]]) %>% 
      as.data.frame() %>% 
      rename_all(~c("x", "y")) %>% 
      rownames_to_column("environment") %>% 
      as_data_frame()
    
    # Return a list
    data_frame(dist_df = list(dist_df), cluster = .$cluster[1])
  }))

## Plot
clust_method_plot_list <- clust_method_df_toplot %>%
  mutate(model1 = model) %>%
  arrange(model, trait) %>%
  group_by(model, population) %>% 
  nest(trait, model1, data) %>%
  mutate(plot_obj = map(data, ~{
    # Unnest the df
    temp1 <- unnest(.)
    
    # Create the plot objects
    dist_plots <- temp1$dist_df %>% 
      map2(.x = ., .y = temp1$trait, ~mutate(.x, trait = .y)) %>%
      map(~ggplot(data = ., aes(x = x, y = y, label = environment)) + 
            geom_point() +
            ggrepel::geom_text_repel(size = 2) + 
            facet_wrap(~ trait, strip.position = "left") + 
            theme_acs() +
            theme(axis.title = element_blank(), axis.text = element_blank()))
    
    # Create the cluster plot objects
    clust_plots <- temp1$cluster %>% 
      map(~ggdendrogram(.) + 
            theme(axis.text.x = element_text(size = 6), axis.text.y = element_blank()))
    
    # Model name
    model <- unique(temp1$model1) %>% str_replace_all(pattern = "_", replacement = " ") %>% str_to_title() %>% str_remove_all(" ")
    
    # Create cowplots
    plots2 <- map2(.x = dist_plots, .y = clust_plots, plot_grid, nrow = 1, rel_widths = c(0.8, 1))
    # Bind these together and return
    plot_grid(plotlist = plots2, ncol = 1) %>% add_sub(label = model) %>% ggdraw()
    
    }))

## Iterate and save
for (i in seq(nrow(clust_method_plot_list))) {
  filename <- str_c("distance_cluster_plot_", clust_method_plot_list$model[i], "_pop", clust_method_plot_list$population[i], ".jpg")
  ggsave(filename = filename, plot = clust_method_plot_list$plot_obj[[i]], path = fig_dir, height = 6, width = 5, dpi = 1000)
}
    






### Rank environments according to a prediction environment
### This will be used for prediction and heritability calculations
### 

# Modify the BLUEs for predictions
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate_at(vars(environment:line_name), as.factor)


## For each prediction environment (the tp+vp envs and just the vp envs), rank the 
## training environments by different distance metrics
pred_envs <- c(tp_vp_env, vp_only_env)
train_envs <- c(tp_vp_env, tp_only_env)

# Summarize the traits available in those environments
val_envs_traits <- S2_MET_BLUEs_use %>%
  filter(environment %in% pred_envs) %>% 
  group_by(environment) %>% 
  distinct(trait) %>%
  ungroup() %>%
  mutate(environment = as.character(environment),
         set = "complete") %>%
  bind_rows(., filter(., str_detect(environment, "17")) %>% mutate(set = "realistic")) %>%
  nest(environment, .key = "val_environments")

train_envs_traits <- S2_MET_BLUEs_use %>%
  filter(environment %in% train_envs) %>% 
  group_by(environment) %>% 
  distinct(trait) %>%
  ungroup() %>%
  mutate(set = "complete") %>%
  bind_rows(., filter(., !str_detect(environment, "17")) %>% mutate(set = "realistic")) %>%
  nest(environment, .key = "train_environments")


## Rank the environments relative to each other
# Do this for both population groups
env_rank_df <- clust_method_df %>% 
  left_join(., val_envs_traits) %>%  ## Add validation and training environments
  left_join(., train_envs_traits) %>%
  mutate(env_rank = pmap(list(dist, val_environments, train_environments), ~{

    dmat <- ..1
    val_env_use <- ..2
    train_env_use <- ..3
    
    ddf <- as.matrix(dmat) %>% 
      broom::fix_data_frame(newcol = "environment") %>%
      filter(environment %in% val_env_use$environment) %>%
      select(environment, which(names(.) %in% train_env_use$environment))
    
    # Create a list of environment ranks
    ddf %>% 
      gather(validation_environment, distance, -environment) %>% 
      filter(environment != validation_environment) %>%
      split(.$environment) %>% 
      map(~.$validation_environment[order(.$distance, decreasing = F)]) %>% ## Decreasing = FALSE (add the most distant (largest) last)
      data_frame(environment = names(.), rank = .)

  })) %>% select(-contains("environments"))


# Combine this data with the the trait information for the prediction environments
# This will remove environment-trait combinations that were not observed.
pred_env_dist_rank <- unnest(env_rank_df, env_rank) %>%
  inner_join(unnest(val_envs_traits), .,  by = c("environment", "trait", "set")) %>%
  rename(validation_environment = environment)




## Randomly order non-prediction environments
# Number of random samples
n_sample <- 10

set.seed(153)
## For each prediction environment, take the environments in one of the distance
## matrices and randomly sample it
pred_env_rank_random <- pred_env_dist_rank %>%
  group_by(mat_set, set, trait, validation_environment) %>%
  do({
    df <- .
    smpls <- rerun(.n = n_sample, sample(df$rank[[1]]))
    data_frame(validation_environment = df$validation_environment[1], trait = df$trait[1],
               set = df$set[1], model = str_c("rank_sample", seq_along(smpls)), 
               rank = smpls)
    
  }) %>% ungroup()


# Now generate 100 samples using all environments
set.seed(1004)

n_sample <- 100
pred_env_random <- pred_env_dist_rank %>%
  group_by(mat_set, set, trait) %>%
  do({
    df <- .
    envs <- df$rank %>% reduce(union)
    smpls <- rerun(.n = n_sample, sample(envs))
    data_frame(validation_environment = NA, trait = df$trait[1],
               set = df$set[1], model = str_c("sample", seq_along(smpls)), 
               rank = smpls)
          
  }) %>% ungroup()





# Save this
save_file <- file.path(result_dir, "distance_method_results.RData")
save("env_rank_df", "pred_env_dist_rank","pred_env_rank_random", "pred_env_random", file = save_file)




