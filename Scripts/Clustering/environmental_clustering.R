## Cluster environments in the S2MET project
## 
## Author: Jeff Neyhart
## Last modified: October 04, 2018
## 
## This script will create different distance matrices for use in clustering. The
## clustering algorithm will be constant
## 

# Load packages and the source script
library(mclust)
library(pbr)
library(ggdendro)
library(cowplot)
library(modelr)

# The head directory
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))



## Load the environment covariable distance matrices
load(file.path(result_dir, "environmental_covariable_distance_mat.RData"))
# Load the genetic correlation estimates
load(file.path(result_dir, "environmental_genetic_correlations.RData"))
# Load the two-way geno/env tables
load(file.path(result_dir, "genotype_environment_ammi_analysis.RData"))

# # Create a new data.frame to hold the different datasets
# S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
#   group_by(trait) %>%
#   nest() %>%
#   mutate(data = map(data, droplevels))



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


## Create different lat/long matrices depending on the trait
gcd_mat_list <- S2_MET_BLUEs %>% 
  filter(trait %in% traits) %>% 
  split(.$trait) %>%
  map2(.x = ., .y = tp_vp_env_trait[names(.)], ~filter(.x, environment %in% .y)) %>%
  map(~unique(.$environment)) %>% 
  map(~trial_lat_long[.,])


# Create pairwise combinations of environments and calculate the great circle distance
trial_great_circle <- gcd_mat_list %>%
  map(~{
    mat <- .
    row.names(mat) %>% 
      expand.grid(environment1 = ., environment2 = ., stringsAsFactors = FALSE) %>%
      mutate(dist = geosphere::distGeo(p1 = mat[environment1,], p2 = mat[environment2,])) })
  

# Convert to a dist object
great_circle_dist_list <- trial_great_circle %>%
  map(~spread(., environment2, dist) %>% 
        remove_rownames() %>% 
        column_to_rownames("environment1") %>%
        as.matrix() %>% 
        as.dist() )


## Phenotypic distance metric
# This is described by Bernardo2010

# First calculate the distance based on each environment
ge_mean_D <- S2_MET_BLUEs %>%
  filter(line_name %in% tp, trait %in% traits) %>%
  # Split by trait
  split(.$trait) %>%
  map2(.x = ., .y = tp_vp_env_trait[names(.)], ~filter(.x, environment %in% .y)) %>%
  map(function(df){
    dist1 <- dist_env(x = df, gen.col = "line_name", env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })



# ## Now calculate the line means in each location
# control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
# 
# # Fit a mixed model with fixed G and L effects and random GY and GYL
# temp <- S2_MET_BLUEs %>%
#   filter(line_name %in% tp, trait %in% traits) %>%
#   # Split by trait
#   split(.$trait) %>%
#   map2(.x = ., .y = tp_vp_env_trait[names(.)], ~filter(.x, environment %in% .y)) %>%
#   map(~mutate_at(., vars(line_name, location), as.factor))
# 
# 
# fit <- lmer(formula = value ~ line_name + location + line_name:location + (1|year) + (1|line_name:year) + (1|line_name:year:location), data = temp$GrainYield,
#             control = control)

location_BLUEs <- S2_MET_BLUEs %>%
  filter(trait %in% traits) %>% 
  group_by(trait, line_name, location) %>%
  mutate(value = mean(value)) %>%
  ungroup()

## Correlations among locations
loc_cor_complete <- location_BLUEs %>%
  split(.$trait) %>% 
  map2(.x = ., .y = tp_vp_env_trait[names(.)], ~filter(.x, environment %in% .y)) %>%
  map(~distinct(., line_name, location, value) %>% spread(location, value) %>% as.data.frame() %>% 
        column_to_rownames("line_name") %>% as.matrix() %>% cor(., use = "pairwise.complete.obs"))

# Calculate the distance between locations using all data
gl_mean_D_complete <- location_BLUEs %>%
  filter(line_name %in% tp, trait %in% traits) %>%
  # Split by trait
  split(.$trait) %>%
  map2(.x = ., .y = tp_vp_env_trait[names(.)], ~filter(.x, environment %in% .y)) %>%
  map(~distinct(., line_name, environment, value)) %>% 
  map(function(df){
    dist1 <- dist_env(x = df, gen.col = "line_name", env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })





## Now calculate the line means in each location, excluding 2017 data
location_BLUEs_realistic <- S2_MET_BLUEs %>%
  filter(year != 2017, trait %in% traits) %>%
  group_by(trait, line_name, location) %>%
  summarize(value = mean(value)) %>%
  ungroup() %>%
  full_join(distinct(S2_MET_BLUEs, trait, environment, location, year), .) %>%
  filter(!is.na(line_name))


## Correlations among locations
loc_cor_realistic <- location_BLUEs_realistic %>%
  split(.$trait) %>%
  map2(.x = ., .y = tp_vp_env_trait[names(.)], ~filter(.x, environment %in% .y)) %>%
  map(~distinct(., line_name, location, value)) %>% 
  map(~select(., line_name, location, value) %>% spread(location, value) %>% as.data.frame() %>% 
        column_to_rownames("line_name") %>% as.matrix() %>% cor(., use = "pairwise.complete.obs"))

# Calculate the distance between locations using all data
gl_mean_D_realistic <- location_BLUEs_realistic %>%
  filter(line_name %in% tp) %>%
  # Split by trait
  split(.$trait) %>%
  map2(.x = ., .y = tp_vp_env_trait[names(.)], ~filter(.x, environment %in% .y)) %>%
  map(~distinct(., line_name, environment, value)) %>% 
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
  filter(set == "complete", trait %in% traits) %>% 
  add_row(trait = traits, model = "great_circle_dist", dist = great_circle_dist_list, cov = map(great_circle_dist_list, ~1 - as.matrix(.))) %>% 
  add_row(trait = traits, model = "pheno_dist", dist = ge_mean_D, cov = env_cor_all[-3]) %>% 
  add_row(trait = traits, model = "pheno_loc_dist", dist = gl_mean_D_complete, cov = loc_cor_complete) %>%
  mutate(set = "complete")
  

dist_method_df_realistic <- ec_sim_mat_df1 %>% 
  filter(set == "realistic", trait %in% traits) %>% 
  add_row(trait = traits, model = "great_circle_dist", dist = great_circle_dist_list, cov = map(great_circle_dist_list, ~1 - as.matrix(.))) %>% 
  add_row(trait = traits, model = "pheno_loc_dist", dist = gl_mean_D_realistic, cov = loc_cor_realistic) %>%
  mutate(set = "realistic")


## Combine
## Remove the distance matrices
dist_method_df <- bind_rows(dist_method_df_complete, dist_method_df_realistic)






#### Cluster environments using model-based clustering

# Mininum number of environments in a cluster
min_env <- 2

## First create matrices to use as data sources
## 
## Two procedures:
## 1. For the complete set, cluster all environments / locations simultaneously
## 2. For the realistic set, cluster only the realistic training data

# Complete set
# GCD
gcd_mat_complete <- data_frame(set = "complete", model = "great_circle_dist", trait = names(gcd_mat_list), data = gcd_mat_list) %>%
  mutate(cluster = map(data, ~env_mclust(data = ., min_env = min_env)))

# Pheno dist
pd_mat_complete <- data_frame(set = "complete", model = "pheno_dist", trait = names(ge_mean_D), data = ge_mean_D) %>%
  mutate(data = map(data, cmdscale),
         cluster = map(data, ~env_mclust(data = ., min_env = min_env)))

# Location dist
pld_mat_complete <- data_frame(set = "complete", model = "pheno_location_dist", trait = names(gl_mean_D_complete), data = gl_mean_D_complete) %>%
  mutate(data = map(data, cmdscale),
         cluster = map(data, ~env_mclust(data = ., min_env = min_env)))

## Environmental covariates
ec_mat_complete <- ec_mats %>% 
  filter(set == "complete", trait %in% traits) %>% 
  rename(data = mat) %>%
  left_join(., data_frame(trait = names(tp_vp_env_trait), env = tp_vp_env_trait)) %>% 
  mutate(data = map2(data, env, ~.x[row.names(.x) %in% .y,])) %>%
  mutate(cluster = map(data, ~env_mclust(data = ., min_env = min_env)),
         model = str_replace_all(ec_group, "_", " ") %>% str_to_title() %>% abbreviate(2),
         model = str_c(model, group)) %>%
  select(set, model, trait, data, cluster)


# Combine
cluster_df_complete <- bind_rows(gcd_mat_complete, pd_mat_complete, pld_mat_complete, ec_mat_complete)

  
## Realistic
# GCD
gcd_mat_realistic <- data_frame(set = "realistic", model = "great_circle_dist", trait = names(gcd_mat_list), data = gcd_mat_list,
                               test_env = map(complete_train_env[-3], ~str_subset(., "17"))) %>%
  mutate(cluster = map2(.x = data, .y = test_env, ~env_mclust(data = .x, min_env = min_env, test.env = .y)))

# Location dist
pld_mat_realistic  <- data_frame(set = "realistic", model = "pheno_location_dist", trait = names(gl_mean_D_realistic), data = gl_mean_D_realistic,
                               test_env = map(complete_train_env[-3], ~str_subset(., "17"))) %>%
  mutate(data = map(data, cmdscale),
         test_env = map2(.x = data, .y = test_env, ~intersect(row.names(.x), .y)),
         cluster = map2(.x = data, .y = test_env, ~env_mclust(data = .x, min_env = min_env, test.env = .y)))

## Environmental covariates
ec_mat_realistic  <- ec_mats %>% 
  filter(set == "realistic", trait %in% traits) %>% 
  rename(data = mat) %>%
  left_join(., data_frame(trait = names(tp_vp_env_trait), env = tp_vp_env_trait)) %>% # filter out undesired environments
  mutate(data = map2(data, env, ~.x[row.names(.x) %in% .y,, drop = FALSE])) %>%
  left_join(., map(complete_train_env[-3], ~str_subset(., "17")) %>% data_frame(trait = names(.), test_env = .)) %>%
  mutate(cluster = map2(.x = data, .y = test_env, ~env_mclust(data = .x, min_env = min_env, test.env = .y)),
         model = str_replace_all(ec_group, "_", " ") %>% str_to_title() %>% abbreviate(2),
         model = str_c(model, group)) %>%
  select(set, model, trait, data, cluster)
  

# Combine
cluster_df_realistic <- bind_rows(gcd_mat_realistic, pld_mat_realistic, ec_mat_realistic)


## Combine all
cluster_df <- bind_rows(cluster_df_complete, cluster_df_realistic) %>%
  select(-test_env)



### Test variance components given clustering
cluster_df_tomodel <- cluster_df %>% 
  filter(set == "complete") %>% 
  unnest(cluster) %>% 
  left_join(filter(S2_MET_BLUEs, line_name %in% tp))

## Control for lmer
control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")

cluster_varcomp <- cluster_df_tomodel %>% 
  group_by(model, trait) %>%
  filter(n_distinct(cluster) > 1) %>% # Pass over trait/models that forced a single cluster
  do({
    df <- .
    wts <- df$std_error^2
    
    # Print message
    print(paste(unique(df$trait), unique(df$model)))
    
    # Fit a model
    # Random effects of genotype, cluster, environment in cluster, G x cluster, and g x environment in cluster
    fit <- lmer(value ~ (1|line_name) + (1|cluster) + (1|cluster:environment) + (1|line_name:cluster) + (1|line_name:cluster:environment), 
                data = df, control = control)
    
    # Random anova
    fit_ranova <- ranova(model = fit) %>%
      broom::tidy() %>% 
      filter(!str_detect(term, "none")) %>% 
      mutate(term = str_remove_all(term, "\\(1 \\| |\\)"))
    
    # Combine and tidy, then return
    as.data.frame(VarCorr(fit)) %>% 
      select(term = grp, variance = vcov) %>% 
      left_join(., fit_ranova, by = "term") %>%
      select(term, variance, LRT, df, p.value)
    
  })
    
single_cluster_cases <- cluster_df_tomodel %>% 
  group_by(model, trait) %>%
  filter(n_distinct(cluster) == 1) %>%
  distinct(model, trait)
  

## Plot
g_cluster_varcomp <- cluster_varcomp %>%
  mutate(varprop = variance / sum(variance)) %>%
  ggplot(aes(x = model, y = varprop, fill = term)) +
  geom_col() +
  facet_grid(trait ~ .) +
  theme_presentation2() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1))
  
ggsave(filename = "cluster_varcomp.jpg", plot = g_cluster_varcomp, path = fig_dir, width = 8, height = 9)


# Save a table
cluster_varcomp_table <- cluster_varcomp %>%
  ungroup() %>%
  mutate(varprop = variance / sum(variance), 
         annotation = paste0(formatC(x = varprop, digits = 2), ifelse(p.value < 0.05, "*", "")), 
         annotation = str_remove_all(annotation, "NA")) %>%
  select(trait, model, term, annotation) %>%
  spread(term, annotation)

write_csv(x = cluster_varcomp_table, path = file.path(fig_dir, "cluster_varcomp_table.csv"))


## Of all trait-model combinations, which result in significant G X cluster?
# model       trait              term     variance        LRT df      p.value
# 1    great_circle_dist  GrainYield line_name:cluster 3.242287e+04  61.206548  1 5.139032e-15
# 2    great_circle_dist HeadingDate line_name:cluster 2.764723e-01  10.605145  1 1.127734e-03
# 3             MYEC_All  GrainYield line_name:cluster 3.479472e+04  56.852976  1 4.696428e-14
# 4             MYEC_All HeadingDate line_name:cluster 8.207656e-01  76.461140  1 2.245892e-18
# 5             MYEC_All PlantHeight line_name:cluster 7.876960e-01   7.698212  1 5.527556e-03
# 6            MYEC_IPCA  GrainYield line_name:cluster 2.189498e+04  43.586845  1 4.055561e-11
# 7            MYEC_IPCA HeadingDate line_name:cluster 6.151390e-01  43.976911  1 3.322725e-11
# 8            MYEC_Mean  GrainYield line_name:cluster 2.477421e+04  39.265137  1 3.699888e-10
# 9            MYEC_Mean HeadingDate line_name:cluster 8.828733e-01  88.245794  1 5.780602e-21
# 10           MYEC_Mean PlantHeight line_name:cluster 7.727645e-01  10.012722  1 1.554626e-03
# 11            OYEC_All  GrainYield line_name:cluster 1.649557e+04  19.217038  1 1.166675e-05
# 12            OYEC_All HeadingDate line_name:cluster 4.445435e-01  19.614002  1 9.477208e-06
# 13            OYEC_All PlantHeight line_name:cluster 1.575744e+00  29.782854  1 4.832471e-08
# 14           OYEC_IPCA  GrainYield line_name:cluster 4.357908e+04  42.842169  1 5.933914e-11
# 15           OYEC_IPCA HeadingDate line_name:cluster 4.605247e-01  33.263482  1 8.047986e-09
# 16           OYEC_Mean  GrainYield line_name:cluster 2.195249e+04  38.632631  1 5.115646e-10
# 17           OYEC_Mean HeadingDate line_name:cluster 7.379050e-01  63.731684  1 1.425730e-15
# 18           OYEC_Mean PlantHeight line_name:cluster 5.870112e-01   9.630360  1 1.913871e-03
# 19          pheno_dist  GrainYield line_name:cluster 7.501530e+04 163.011091  1 2.487701e-37
# 20          pheno_dist HeadingDate line_name:cluster 1.165701e+00  87.328155  1 9.192991e-21
# 21          pheno_dist PlantHeight line_name:cluster 3.291683e+00 199.800223  1 2.309020e-45
# 22 pheno_location_dist  GrainYield line_name:cluster 2.789433e+04  40.461515  1 2.005272e-10
# 23 pheno_location_dist HeadingDate line_name:cluster 1.871909e+00 169.755848  1 8.365500e-39
# 24 pheno_location_dist PlantHeight line_name:cluster 9.410029e-01  17.598128  1 2.728568e-05

## Line_name:cluster:



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
env_rank_df <- dist_method_df %>% 
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
save("env_rank_df", "cluster_df", "pred_env_dist_rank","pred_env_rank_random", "pred_env_random", file = save_file)













