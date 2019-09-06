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
library(cluster)
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
# Load the two-way geno/env tables and AMMI results
load(file.path(result_dir, "genotype_environment_phenotypic_analysis.RData"))



# # Create a new data.frame to hold the different datasets
# S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
#   group_by(trait) %>%
#   nest() %>%
#   mutate(data = map(data, droplevels))




## Use fitted values of AMMI model to cluster environments based on the best lines in each
ammi_fitted <- ammi_out %>%
  group_by(set, trait) %>% 
  do({
    df <- .
    
    ## Filter for PC1
    scores <- df$ammi[[1]][c("escores", "gscores")] %>%
      map(filter, PC == "PC1")
    
    # Create matrices of genotype and environment effects
    g_effect <- scores$gscores %>% 
      select(line_name, effect) %>% 
      as.data.frame() %>% 
      column_to_rownames("line_name") %>% 
      as.matrix()
    
    e_effect <- scores$escores %>% 
      select(environment, effect) %>% 
      as.data.frame() %>% 
      column_to_rownames("environment") %>% 
      as.matrix()
    
    
    g_effect1 <- replicate(nrow(e_effect), g_effect[,1])
    e_effect1 <- replicate(nrow(g_effect), e_effect[,1])
    
    main_effect <- g_effect1 + t(e_effect1)
    colnames(main_effect) <- row.names(e_effect1)
    
    # Create interaction scores from the PC
    g_scores <- scores$gscores %>% 
      select(line_name, score) %>%
      as.data.frame() %>% 
      column_to_rownames("line_name") %>% 
      as.matrix()
    
    e_scores <- scores$escores %>% 
      select(environment, score) %>% 
      as.data.frame() %>% 
      column_to_rownames("environment") %>% 
      as.matrix()
    
    int_effect <- tcrossprod(g_scores, e_scores)

    ## Calculate fitted values
    fitted_effect <- main_effect + int_effect
    
    ## Determine the best (and worst) in each environment
    top <- apply(X = fitted_effect, MARGIN = 2, FUN = function(x) row.names(fitted_effect)[order(x, decreasing = TRUE)[1:5]])
    bottom <- apply(X = fitted_effect, MARGIN = 2, FUN = function(x) row.names(fitted_effect)[order(x, decreasing = FALSE)[1:5]])
    
    ## Return the fitted values, the top, and bottom
    fitted_effect1 <- as.data.frame(fitted_effect) %>% rownames_to_column("line_name") %>% gather(environment, fitted, -line_name)
    data_frame(fitted = list(fitted_effect1), top = list(top), bottom = list(bottom))
    
  })
    


## How many mega-environments per trait and set
ammi_fitted %>% 
  mutate(favorable = ifelse(trait == "GrainYield", top, bottom)) %>% 
  summarize(n_cluster = map_dbl(favorable, ~n_distinct(.[1,])))
    
# set           trait       n_cluster
# 1 complete      GrainYield          2
# 2 complete      HeadingDate         2
# 3 complete      PlantHeight         2
# 4 realistic2015 GrainYield          2
# 5 realistic2015 HeadingDate         2
# 6 realistic2015 PlantHeight         1
# 7 realistic2016 GrainYield          3
# 8 realistic2016 HeadingDate         2
# 9 realistic2016 PlantHeight         2
# 10 realistic2017 GrainYield          2
# 11 realistic2017 HeadingDate         1
# 12 realistic2017 PlantHeight         1

## Create clusters
ammi_clusters <- ammi_fitted %>%
  filter(trait %in% traits) %>%
  ungroup() %>%
  mutate(favorable = ifelse(trait == "GrainYield", top, bottom)) %>% 
  mutate(cluster = map(favorable, ~as.data.frame(.[1,]) %>% rownames_to_column("environment") %>% 
                         mutate_at(vars(-environment), funs(cluster = as.numeric(as.factor(.)))) %>% select(environment, cluster)) ) %>%
  ungroup() %>%
  select(set, trait, cluster) %>%
  mutate(model = "AMMI")



## Calculate distance 
ammi_distance <- ammi_out %>%
  group_by(set, trait) %>%
  do(dist = {
    
    row <- .

    e_data <- row$ammi[[1]]$escores %>% 
      filter(PC == "PC1")
    
    ## Calculate distance between scaled environment mean and IPCA score
    e_mean_ipca_dist <- e_data %>%
      select(environment, effect, score) %>% 
      as.data.frame() %>% 
      column_to_rownames("environment") %>% 
      as.matrix() %>% 
      scale() %>% 
      dist()
    
    ## Return
    e_mean_ipca_dist
    
  }) %>% ungroup() %>%
  mutate(model = "AMMI")



# 
# ## What about model-based clustering of environments based on env mean + IPCA score
# ammi_fitted_alt_complete <- ammi_out %>%
#   filter(set == "complete") %>%
#   group_by(set, trait) %>% 
#   do({
#     df <- .
#     
#     ## Filter for PC1
#     scores <- df$ammi[[1]][c("escores", "gscores")] %>%
#       map(filter, PC == "PC1")
# 
#     # Create DF of IPCA and effect
#     e_mat <- scores$escores %>%
#       select(environment, effect, score) %>%
#       as.data.frame() %>%
#       column_to_rownames("environment")
#     
#     env_mclust(data = e_mat, min_env = 2)
#     
#   }) %>% ungroup()
# 
# 
# ## How many mega-environments per trait and set
# ammi_fitted_alt_complete %>% 
#   group_by(trait) %>%
#   summarize(n_cluster = n_distinct(cluster))













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
  filter(line_name %in% tp) %>%
  filter(trait %in% traits) %>%
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
  filter(line_name %in% tp) %>%
  filter(trait %in% traits) %>%
  # Split by trait
  split(.$trait) %>%
  map2(.x = ., .y = tp_vp_env_trait[names(.)], ~filter(.x, environment %in% .y)) %>%
  map(~distinct(., line_name, environment, value)) %>% 
  map(function(df){
    dist1 <- dist_env(x = df, gen.col = "line_name", env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })





## Now calculate the line means in each location, excluding data from each dropped year
location_BLUEs_realistic_df <- distinct(S2_MET_BLUEs, year, trait) %>%
  filter(trait %in% traits) %>%
  mutate(set = paste0("realistic", year), data = list(NULL), cor = list(NULL), dist = list(NULL))

for (i in seq(nrow(location_BLUEs_realistic_df))) {

  location_BLUEs_realistic_df$data[[i]] <- S2_MET_BLUEs %>%
    filter(year != location_BLUEs_realistic_df$year[i], trait %in% location_BLUEs_realistic_df$trait[i]) %>%
    group_by(trait, line_name, location) %>%
    summarize(value = mean(value)) %>%
    ungroup() %>%
    full_join(distinct(S2_MET_BLUEs, trait, environment, location, year), ., by = c("location", "trait")) %>%
    filter(!is.na(line_name))
  
  location_BLUEs_realistic_df$cor[[i]] <- location_BLUEs_realistic_df$data[[i]] %>%
    filter(environment %in% tp_vp_env_trait[[location_BLUEs_realistic_df$trait[i]]]) %>%
    distinct(., line_name, location, value) %>%
    select(., line_name, location, value) %>% 
    spread(location, value) %>% 
    as.data.frame() %>% 
    column_to_rownames("line_name") %>% 
    as.matrix() %>% 
    cor(., use = "pairwise.complete.obs")
  
  location_BLUEs_realistic_df$dist[[i]] <- location_BLUEs_realistic_df$data[[i]] %>%
    filter(line_name %in% tp) %>%
    filter(environment %in% tp_vp_env_trait[[location_BLUEs_realistic_df$trait[i]]]) %>%
    distinct(., line_name, environment, value) %>%
    {
      dist1 <- dist_env(x = ., gen.col = "line_name", env.col = "environment", pheno.col = "value")
      # Replace missing with the mean
      dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
      dist1
    }

}





## Combine the one-year and multi-year ECs
# Create 3 different sets of matrices based on the distance calculations
sim_mat_group <- bind_rows(
  mutate(sim_mat_df, mat_set = "Malosetti"),
  mutate(sim_mat_df1, mat_set = "Jarquin"),
  mutate(sim_mat_df2, mat_set = "MalosettiStand")
) %>% select(-mat)


## Create the model using an abbreviation for the ec_group (multi-year) and the group (All, Mean. IPCA)
ec_sim_mat_df1 <- sim_mat_group %>%
  mutate(model = str_replace_all(ec_group, "multiyear", "multi year") %>% str_to_title() %>% abbreviate(2),
         model = str_c(model, group),
         dist = map(sim_mat, dist)) %>%
  select(mat_set, set, model, trait, cov = sim_mat, dist)


## Combine - scale such that greater values mean "closer" or more similar
dist_method_df_complete <- ec_sim_mat_df1 %>% 
  filter(set == "complete") %>% 
  # filter(trait %in% traits) %>% 
  add_row(trait = unique(.$trait), model = "great_circle_dist", dist = great_circle_dist_list, cov = map(great_circle_dist_list, ~1 - as.matrix(.))) %>% 
  add_row(trait = unique(.$trait), model = "pheno_dist", dist = ge_mean_D, cov = env_cor_all[traits]) %>% 
  add_row(trait = unique(.$trait), model = "pheno_loc_dist", dist = gl_mean_D_complete, cov = loc_cor_complete) %>%
  mutate(set = "complete")
  

dist_method_df_realistic <- ec_sim_mat_df1 %>% 
  # filter(str_detect(set, "realistic")) %>%
  # filter(trait %in% traits) %>% 
  bind_rows(., data_frame(trait = unique(.$trait), model = "great_circle_dist", dist = great_circle_dist_list, cov = map(great_circle_dist_list, ~1 - as.matrix(.))) %>% crossing(set = unique(location_BLUEs_realistic_df$set), .)) %>%
  bind_rows(., location_BLUEs_realistic_df %>% select(set, trait, cov = cor, dist) %>% mutate(model = "pheno_loc_dist")) %>%
  filter(set == "realistic2017")
  

## Combine
## Remove the distance matrices
dist_method_df <- bind_rows(ammi_distance, dist_method_df_complete, dist_method_df_realistic) %>%
  filter(set == "realistic2017")






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

# PAM
gcd_mat_complete_pam <- data_frame(set = "complete", model = "great_circle_dist", trait = names(gcd_mat_list), data = gcd_mat_list) %>%
  mutate(cluster = map(data, ~env_mclust(data = ., min_env = min_env, method = "pam")))


# Pheno dist
pd_mat_complete <- data_frame(set = "complete", model = "pheno_dist", trait = names(ge_mean_D), data = ge_mean_D) %>%
  mutate(data = map(data, cmdscale),
         cluster = map(data, ~env_mclust(data = ., min_env = min_env)))

pd_mat_complete_pam <- data_frame(set = "complete", model = "pheno_dist", trait = names(ge_mean_D), data = ge_mean_D) %>%
  mutate(data = map(data, cmdscale),
         cluster = map(data, ~env_mclust(data = ., min_env = min_env, method = "pam", diss = FALSE)))

# Location dist
pld_mat_complete <- data_frame(set = "complete", model = "pheno_location_dist", trait = names(gl_mean_D_complete), data = gl_mean_D_complete) %>%
  mutate(data = map(data, cmdscale),
         cluster = map(data, ~env_mclust(data = ., min_env = min_env)))

pld_mat_complete_pam <- data_frame(set = "complete", model = "pheno_location_dist", trait = names(gl_mean_D_complete), data = gl_mean_D_complete) %>%
  mutate(data = map(data, cmdscale),
         cluster = map(data, ~env_mclust(data = ., min_env = min_env, method = "pam")))


## Environmental covariates
ec_mats_data <- ec_mats %>% 
  mutate(ec_group = str_replace(ec_group, "multiyear", "multi_year")) %>%
  filter(set == "complete", trait %in% traits) %>% 
  rename(data = mat) %>%
  left_join(., data_frame(trait = names(tp_vp_env_trait), env = tp_vp_env_trait)) %>% 
  mutate(data = map2(data, env, ~.x[row.names(.x) %in% .y,,drop = F]))

ec_mat_complete <- ec_mats_data %>%
  mutate(cluster = map(data, ~env_mclust(data = ., min_env = min_env)),
         model = str_replace_all(ec_group, "_", " ") %>% str_to_title() %>% abbreviate(2),
         model = str_c(model, group)) %>%
  select(set, model, trait, data, cluster)

## Environmental covariates
ec_mat_complete_pam <- ec_mats_data %>%
  mutate(cluster = map(data, ~env_mclust(data = ., min_env = min_env, method = "pam")),
         model = str_replace_all(ec_group, "_", " ") %>% str_to_title() %>% abbreviate(2),
         model = str_c(model, group)) %>%
  select(set, model, trait, data, cluster)


# Combine
cluster_df_complete <- bind_rows(gcd_mat_complete, pd_mat_complete, pld_mat_complete, ec_mat_complete, filter(ammi_clusters, set == "complete"))
cluster_df_complete_pam <- bind_rows(gcd_mat_complete_pam, pd_mat_complete_pam, pld_mat_complete_pam, ec_mat_complete_pam, 
                                     filter(ammi_clusters, set == "complete"))

  





## Realistic
# GCD
gcd_mat_realistic <- distinct(dist_method_df_realistic, set) %>% 
  crossing(., data_frame(model = "great_circle_dist", trait = names(gcd_mat_list), data = gcd_mat_list, test_env = realistic_test_env[traits])) %>% 
  mutate(test_env = map2(set, test_env, ~str_subset(.y, str_extract(.x, "1[0-9]"))),
         cluster = map2(.x = data, .y = test_env, ~env_mclust(data = .x, min_env = min_env, test.env = .y)))

gcd_mat_realistic_pam <- distinct(dist_method_df_realistic, set) %>% 
  crossing(., data_frame(model = "great_circle_dist", trait = names(gcd_mat_list), data = gcd_mat_list, test_env = realistic_test_env[traits])) %>% 
  mutate(test_env = map2(set, test_env, ~str_subset(.y, str_extract(.x, "1[0-9]"))),
         cluster = map2(.x = data, .y = test_env, ~env_mclust(data = .x, min_env = min_env, test.env = .y, method = "pam")))

# Location dist
pld_mat_realistic  <- gcd_mat_realistic %>% 
  distinct(set, trait, test_env) %>% 
  left_join(., location_BLUEs_realistic_df %>% mutate(model = "pheno_location_dist") %>% select(model, set, trait, data = dist)) %>%
  mutate(data = map(data, cmdscale),
         test_env = map2(.x = data, .y = test_env, ~intersect(row.names(.x), .y)),
         cluster = map2(.x = data, .y = test_env, ~env_mclust(data = .x, min_env = min_env, test.env = .y)))

pld_mat_realistic_pam  <- gcd_mat_realistic %>% 
  distinct(set, trait, test_env) %>% 
  left_join(., location_BLUEs_realistic_df %>% mutate(model = "pheno_location_dist") %>% select(model, set, trait, data = dist)) %>%
  mutate(data = map(data, cmdscale),
         test_env = map2(.x = data, .y = test_env, ~intersect(row.names(.x), .y)),
         cluster = map2(.x = data, .y = test_env, ~env_mclust(data = .x, min_env = min_env, test.env = .y, method = "pam")))

## Environmental covariates
ec_mat_realistic  <- ec_mats %>% 
  mutate(ec_group = str_replace(ec_group, "multiyear", "multi_year")) %>%
  filter(str_detect(set, "realistic"), trait %in% traits) %>% 
  rename(data = mat) %>%
  left_join(., data_frame(trait = names(tp_vp_env_trait), env = tp_vp_env_trait)) %>% # filter out undesired environments
  mutate(data = map2(data, env, ~.x[row.names(.x) %in% .y,, drop = FALSE])) %>%
  left_join(., distinct(pld_mat_realistic, set, trait, test_env)) %>%
  mutate(cluster = map2(.x = data, .y = test_env, ~env_mclust(data = .x, min_env = min_env, test.env = .y)),
         model = str_replace_all(ec_group, "_", " ") %>% str_to_title() %>% abbreviate(2),
         model = str_c(model, group)) %>%
  select(set, model, trait, data, cluster, test_env)

## Environmental covariates
ec_mat_realistic_pam  <- ec_mats %>% 
  mutate(ec_group = str_replace(ec_group, "multiyear", "multi_year")) %>%
  filter(str_detect(set, "realistic"), trait %in% traits) %>% 
  rename(data = mat) %>%
  left_join(., data_frame(trait = names(tp_vp_env_trait), env = tp_vp_env_trait)) %>% # filter out undesired environments
  mutate(data = map2(data, env, ~.x[row.names(.x) %in% .y,, drop = FALSE])) %>%
  left_join(., distinct(pld_mat_realistic, set, trait, test_env)) %>%
  mutate(cluster = map2(.x = data, .y = test_env, ~env_mclust(data = .x, min_env = min_env, test.env = .y, method = "pam")),
         model = str_replace_all(ec_group, "_", " ") %>% str_to_title() %>% abbreviate(2),
         model = str_c(model, group)) %>%
  select(set, model, trait, data, cluster, test_env)

  

# Combine
cluster_df_realistic <- bind_rows(gcd_mat_realistic, pld_mat_realistic, ec_mat_realistic) %>%
  # Make sure the test environments are all consistent
  split(.$trait) %>%
  map_df(~mutate(., unique_env = map(data, row.names) %>% subset(., !map_lgl(., is.null)) %>% reduce(intersect) %>% list(),
              unique_test_env = map(test_env, ~.) %>% subset(., !map_lgl(., is.null)) %>% reduce(intersect) %>% list())) %>%
  mutate(test_env = map2(test_env, unique_test_env, intersect),
         cluster = map2(cluster, unique_env, ~filter(.x, environment %in% .y))) %>%
  select(-contains("unique"))
  
cluster_df_realistic_pam <- bind_rows(gcd_mat_realistic_pam, pld_mat_realistic_pam, ec_mat_realistic_pam) %>%
  # Make sure the test environments are all consistent
  split(.$trait) %>%
  map_df(~mutate(., unique_env = map(data, row.names) %>% subset(., !map_lgl(., is.null)) %>% reduce(intersect) %>% list(),
                 unique_test_env = map(test_env, ~.) %>% subset(., !map_lgl(., is.null)) %>% reduce(intersect) %>% list())) %>%
  mutate(test_env = map2(test_env, unique_test_env, intersect),
         cluster = map2(cluster, unique_env, ~filter(.x, environment %in% .y))) %>%
  select(-contains("unique"))


## Combine all
cluster_df <- bind_rows(cluster_df_complete, cluster_df_realistic) %>%
  select(-test_env)

## Combine all
cluster_df_pam <- bind_rows(cluster_df_complete_pam, cluster_df_realistic_pam) %>%
  select(-test_env)






## How many clusters per trait, set, and model?
cluster_summary <- cluster_df %>% 
  unnest(cluster) %>% 
  mutate(model = ifelse(model == "pheno_location_dist", "pheno_loc_dist", model),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr_use)) %>%
  filter(model %in% dist_method_abbr_use) %>%
  group_by(set, trait, model) %>% 
  summarize(nCluster = n_distinct(cluster)) %>% 
  spread(model, nCluster)

## How many clusters per trait, set, and model?
cluster_summary_pam <- cluster_df_pam %>% 
  unnest(cluster) %>% 
  mutate(model = ifelse(model == "pheno_location_dist", "pheno_loc_dist", model),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr_use)) %>%
  filter(model %in% dist_method_abbr_use) %>%
  group_by(set, trait, model) %>% 
  summarize(nCluster = n_distinct(cluster)) %>% 
  spread(model, nCluster)

cluster_summary %>%
  kableExtra::kable() %>%
  kableExtra::kable_styling()

## Write to csv
cluster_summary %>% 
  write_csv(x = ., path = file.path(fig_dir, "number_of_clusters.csv"))

## Write to csv
cluster_summary_pam %>% 
  write_csv(x = ., path = file.path(fig_dir, "number_of_clusters_pam.csv"))



## A function to compare clusters
## This is VERY ad hoc and should not be generalized
compare_clusters <- function(x) {
  comparison <- select(x, model, cluster) %>% 
    crossing(., .)
  
  comparison1 <- map2_dbl(.x = comparison$cluster, .y = comparison$cluster1, ~{
    
    # This calculate the proportion of pairs of environments assigned to the same cluster
    comp1 <- crossing(environment = .x[[1]], environment1 = .y[[1]]) %>%
      filter(environment != environment1) %>%
      mutate(same_cluster_x = NA, same_cluster_y = NA)
    
    for (i in seq(nrow(comp1))) {
      e1 <- comp1$environment[[i]]
      e2 <- comp1$environment1[[i]]
      
      comp1$same_cluster_x[[i]] <- n_distinct(subset(.x, environment %in% c(e1, e2), cluster, drop = T)) == 1
      comp1$same_cluster_y[[i]] <- n_distinct(subset(.y, environment %in% c(e1, e2), cluster, drop = T)) == 1
      
    }
    
    comp2 <- comp1 %>%
      mutate(same_cluster = same_cluster_x == same_cluster_y)
    
    mean(comp2$same_cluster)
    
  })
  
  
  
  ## Find the number of clusters and the mean (and range) of the number 
  ## of environments in each cluster
  model_nCluster <- group_by(x, model) %>%
    do({
      tb <- table(.$cluster[[1]]$cluster)
      data.frame(nCluster = length(tb), nEnv_mean = mean(tb), nEnv_min = min(tb), nEnv_max = max(tb)) %>%
        mutate(annotation = paste0(round(nEnv_mean, 1), "\n(", nEnv_min, ",", nEnv_max, ")"),
               annotation1 = paste0(nCluster, "\n(", round(nEnv_mean, 1), ")"))
      
    }) %>% ungroup()
  
  
  ## Add the mean
  comparison %>%
    select(model, model1) %>%
    mutate(prop_common = comparison1) %>%
    left_join(., model_nCluster, by = "model")
  
}


## What is the proportion of environments assigned to the same cluster across distance method
cluster_compare <- cluster_df %>%
  select(-data) %>% 
  group_by(set, trait) %>% 
  do(compare_clusters(x = .)) %>% 
  ungroup()

## What is the proportion of environments assigned to the same cluster across distance method
cluster_compare_pam <- cluster_df_pam %>%
  select(-data) %>% 
  group_by(set, trait) %>% 
  do(compare_clusters(x = .)) %>% 
  ungroup()


## Filter out
cluster_compare1 <- cluster_compare %>%
  mutate_at(vars(model, model1), funs(ifelse(. == "pheno_location_dist", "pheno_loc_dist", .))) %>%
  mutate_at(vars(model, model1), funs(str_replace_all(., dist_method_abbr))) %>%
  filter_at(vars(model1, model), all_vars(. %in% dist_method_abbr_use)) %>%
  mutate_at(vars(model, model1), funs(factor(., levels = dist_method_abbr_use))) %>%
  arrange(set, trait, model, model1) %>%
  # Designate upper or lower triangle
  split(list(.$set, .$trait)) %>%
  map_df(~mutate(., lower_triangle = duplicated(map2_chr(.x = model, .y = model1, ~paste(sort(c(as.character(.x), as.character(.y))), collapse = "_")))))


cluster_compare1 %>% group_by(set, model) %>% summarize(mean = mean(prop_common))

cluster_compare1 %>% 
  filter_at(vars(model, model1), all_vars(. %in% c("AMMI", "PD", "LocPD"))) %>%
  group_by(set, model, model1) %>% 
  summarize(mean = mean(prop_common))

cluster_compare1 %>% 
  filter_at(vars(model, model1), all_vars(. %in% c("All-EC", "Mean-EC", "IPCA-EC"))) %>%
  group_by(set, model) %>% 
  summarize(mean = mean(prop_common))


## Plot heatmaps
heat_colors <- wesanderson::wes_palette("Zissou1")

cluster_compare_hm_list <- cluster_compare1 %>% 
  filter(!(str_detect(set,"realistic") & model == "AMMI")) %>%
  filter(set %in% c("complete", "realistic2017")) %>%
  filter(trait %in% traits) %>% 
  mutate(set = str_replace_all(set, set_replace),
         # set = str_replace_all(set, "Time-forward", "Leave-one-year-out "),
         set = str_replace_all(set, "2017", ""),
         prop_common = ifelse(model == model1, NA, prop_common),
         annotation_use = ifelse(model == model1, annotation1, round(prop_common, 2))) %>%
  filter(!lower_triangle) %>%
  split(.$set) %>%
  map(~ggplot(data = ., aes(x = model, y = model1, fill = prop_common)) + 
        geom_tile() +  
        geom_text(aes(label = annotation_use), size = 2) + 
        scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[3], high = heat_colors[5], midpoint = 0.5, limits = c(0, 1), na.value = "grey95",
                             name = "Proportion of\noverlap") +
        facet_wrap(~ trait, labeller = labeller(trait = str_add_space)) +
        labs(subtitle = unique(.$set)) +
        theme_presentation2(base_size = 10) +
        theme(axis.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom") )

g_cluster_compare <- plot_grid(plotlist = map(cluster_compare_hm_list, ~. + theme(legend.position = "none")), ncol = 1, 
                               labels = LETTERS[seq_along(cluster_compare_hm_list)], align = "hv")
# g_cluster_compare1 <- plot_grid(g_cluster_compare, get_legend(cluster_compare_hm_list[[1]]), ncol = 1, rel_heights = c(1, 0.1))
g_cluster_compare1 <- plot_grid(g_cluster_compare, get_legend(cluster_compare_hm_list[[1]] + theme(legend.position = "right")), 
                                nrow = 1, rel_widths = c(1, 0.15))


ggsave(filename = "cluster_compare_heatmap.jpg", plot = g_cluster_compare1, path = fig_dir, width = 6.5, height = 5, dpi = 1000)





#### Greyscale version
cluster_compare_hm_list <- cluster_compare1 %>% 
  filter(!(str_detect(set,"realistic") & model == "AMMI")) %>%
  filter(set %in% c("complete", "realistic2017")) %>%
  filter(trait %in% traits) %>% 
  mutate(set = str_replace_all(set, set_replace),
         # set = str_replace_all(set, "Time-forward", "Leave-one-year-out "),
         set = str_replace_all(set, "2017", ""),
         prop_common = ifelse(model == model1, NA, prop_common),
         annotation_use = ifelse(model == model1, annotation1, round(prop_common, 2))) %>%
  filter(!lower_triangle) %>%
  split(.$set) %>%
  map(~ggplot(data = ., aes(x = model, y = model1, fill = prop_common)) + 
        geom_tile() +  
        geom_text(aes(label = annotation_use), size = 2) + 
        scale_fill_gradient2(low = grey.colors(11)[11], mid = grey.colors(11)[5], high = grey.colors(10)[1], 
                             midpoint = 0.5, limits = c(0, 1), na.value = "grey95",
                             name = "Proportion of\noverlap") +
        facet_wrap(~ trait, labeller = labeller(trait = str_add_space)) +
        labs(subtitle = unique(.$set)) +
        theme_presentation2(base_size = 10) +
        theme(axis.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom") )

g_cluster_compare <- plot_grid(plotlist = map(cluster_compare_hm_list, ~. + theme(legend.position = "none")), ncol = 1, 
                               labels = LETTERS[seq_along(cluster_compare_hm_list)], align = "hv")
# g_cluster_compare1 <- plot_grid(g_cluster_compare, get_legend(cluster_compare_hm_list[[1]]), ncol = 1, rel_heights = c(1, 0.1))
g_cluster_compare1 <- plot_grid(g_cluster_compare, get_legend(cluster_compare_hm_list[[1]] + theme(legend.position = "right")), 
                                nrow = 1, rel_widths = c(1, 0.15))


ggsave(filename = "cluster_compare_heatmap_greyscale.jpg", plot = g_cluster_compare1, path = fig_dir, width = 6.5, height = 5, dpi = 1000)










### PAM
## Filter out
cluster_compare1 <- cluster_compare_pam %>%
  mutate_at(vars(model, model1), funs(ifelse(. == "pheno_location_dist", "pheno_loc_dist", .))) %>%
  mutate_at(vars(model, model1), funs(str_replace_all(., dist_method_abbr))) %>%
  filter_at(vars(model1, model), all_vars(. %in% dist_method_abbr_use)) %>%
  mutate_at(vars(model, model1), funs(factor(., levels = dist_method_abbr_use))) %>%
  arrange(set, trait, model, model1) %>%
  # Designate upper or lower triangle
  split(list(.$set, .$trait)) %>%
  map_df(~mutate(., lower_triangle = duplicated(map2_chr(.x = model, .y = model1, ~paste(sort(c(as.character(.x), as.character(.y))), collapse = "_")))))


cluster_compare1 %>% group_by(set, model) %>% summarize(mean = mean(prop_common))

cluster_compare1 %>% 
  filter_at(vars(model, model1), all_vars(. %in% c("AMMI", "PD", "LocPD"))) %>%
  group_by(set, model, model1) %>% 
  summarize(mean = mean(prop_common))

cluster_compare1 %>% 
  filter_at(vars(model, model1), all_vars(. %in% c("All-EC", "Mean-EC", "IPCA-EC"))) %>%
  group_by(set, model) %>% 
  summarize(mean = mean(prop_common))


## Plot heatmaps
heat_colors <- wesanderson::wes_palette("Zissou1")

cluster_compare_hm_list <- cluster_compare1 %>% 
  filter(!(str_detect(set,"realistic") & model == "AMMI")) %>%
  filter(set %in% c("complete", "realistic2017")) %>%
  filter(trait %in% traits) %>% 
  mutate(set = str_replace_all(set, set_replace),
         # set = str_replace_all(set, "Time-forward", "Leave-one-year-out "),
         set = str_replace_all(set, "2017", ""),
         prop_common = ifelse(model == model1, NA, prop_common),
         annotation_use = ifelse(model == model1, annotation1, round(prop_common, 2))) %>%
  filter(!lower_triangle) %>%
  split(.$set) %>%
  map(~ggplot(data = ., aes(x = model, y = model1, fill = prop_common)) + 
        geom_tile() +  
        geom_text(aes(label = annotation_use), size = 2) + 
        scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[3], high = heat_colors[5], midpoint = 0.5, limits = c(0, 1), na.value = "grey95",
                             name = "Proportion of\noverlap") +
        facet_wrap(~ trait, labeller = labeller(trait = str_add_space)) +
        labs(subtitle = unique(.$set)) +
        theme_presentation2(base_size = 10) +
        theme(axis.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom") )

g_cluster_compare <- plot_grid(plotlist = map(cluster_compare_hm_list, ~. + theme(legend.position = "none")), ncol = 1, 
                               labels = LETTERS[seq_along(cluster_compare_hm_list)], align = "hv")
# g_cluster_compare1 <- plot_grid(g_cluster_compare, get_legend(cluster_compare_hm_list[[1]]), ncol = 1, rel_heights = c(1, 0.1))
g_cluster_compare1 <- plot_grid(g_cluster_compare, get_legend(cluster_compare_hm_list[[1]] + theme(legend.position = "right")), 
                                nrow = 1, rel_widths = c(1, 0.15))


ggsave(filename = "cluster_compare_heatmap_pam.jpg", plot = g_cluster_compare1, path = fig_dir, width = 6.5, height = 5, dpi = 1000)












### Are clusters defined more by location or environment?
## Compare the proportion of each year/location in a cluster with the expectation of random
## environment assignments and the cluster assignments 
## What is the proportion of environments assigned to the same cluster across distance method
cluster_grouping <- cluster_df %>% 
  select(-data) %>% 
  group_by(set, trait, model) %>% 
  do({
    row <- .
    
    clus <- row$cluster[[1]] %>%
      mutate(location = str_extract(environment, "[A-Z]{3}"),
             year = str_extract(environment, "[0-9]{2}"))
    
    # List of all locations/years
    locs <- clus$location
    yrs <- clus$year
    
    ## Calculate proportions of locations and year
    clus_props <- clus %>% 
      group_by(cluster) %>% 
      summarize_at(vars(location, year), funs(list(as.data.frame(table(.))))) %>%
      mutate(nEnv = map(location, "Freq") %>% map_dbl(sum),
             nLoc = map_dbl(location, ~n_distinct(.$location)),
             nYr = map_dbl(year, ~n_distinct(.$year))) %>%
      # Calculate expectations
      mutate(expLoc = map_dbl(nEnv, ~mean(replicate(1000, n_distinct(sample(locs, .))))),
             expYr = map_dbl(nEnv, ~mean(replicate(1000, n_distinct(sample(yrs, .)))))) %>%
      # Calculate p value from exact binomial test
      mutate(LocTest = pmap_dbl(list(nLoc, nEnv, expLoc), ~binom.test(x = ..1, n = ..2, p = ..3 / ..2)$p.value),
             YrTest = pmap_dbl(list(nYr, nEnv, expYr), ~binom.test(x = ..1, n = ..2, p = ..3 / ..2)$p.value))
    

    ## Calculate some summary
    summarize_at(clus_props, vars(nEnv:expYr), mean) %>%
      mutate(LocTest = prop.test(x = nLoc, n = nEnv, p = expLoc / nEnv )$p.value,
             YrTest = prop.test(x = nYr, n = nEnv, p = expYr / nEnv)$p.value)
    
    
  }) %>% ungroup()

## Analyze tests
cluster_grouping %>%
  select(set:model, LocTest, YrTest) %>%
  gather(test, p_value, LocTest, YrTest) %>%
  # filter(p_value < alpha)
  arrange(p_value)


## Table
cluster_grouping %>%
  mutate(LocDev = nLoc - expLoc, YrDev = nYr - expYr) %>% 
  select(set:model, LocDev, YrDev) %>% 
  gather(term, deviation, LocDev, YrDev) %>% 
  split(list(.$set, .$trait, .$term)) %>%
  map_df(~arrange(., desc(abs(deviation)))) %>%
  write_csv(x = ., path = file.path(fig_dir, "cluster_location_year_bias.csv"))
  







### Test variance components given clustering
cluster_df_tomodel <- bind_rows(mutate(cluster_df, method = "mclust"), mutate(cluster_df_pam, method = "pam")) %>% 
  filter(set %in% c("complete", "realistic2017")) %>%
  unnest(cluster) %>% 
  left_join(filter(S2_MET_BLUEs, line_name %in% tp))

## Control for lmer
control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")

cluster_varcomp <- cluster_df_tomodel %>% 
  group_by(set, model, trait, method) %>%
  filter(n_distinct(cluster) > 1) %>% # Pass over trait/models that forced a single cluster
  do({
    df <- .
    wts <- df$std_error^2
    
    # Print message
    print(paste(unique(df$trait), unique(df$model)))
    
    # # Fit a model
    # # Random effects of genotype, cluster, environment in cluster, G x cluster, and g x environment in cluster
    # fit <- lmer(value ~ (1|line_name) + (1|cluster) + (1|cluster:environment) + (1|line_name:cluster) + (1|line_name:cluster:environment), 
    #             data = df, control = control)
    
    # ## Alternative model
    # fit <- lmer(value ~ (1|line_name) + (1|cluster) + (1|cluster:location) + (1|cluster:year) + (1|cluster:location:year) + 
    #               (1|line_name:cluster) + (1|line_name:cluster:location) + (1|line_name:cluster:year) + (1|line_name:year) +
    #               (1|line_name:cluster:location:year), 
    #             data = df, control = control)
    
    ## Model as in Atlin2000
    formula <- value ~ 1 + (1|year) + (1|cluster) + (1|cluster:location) + (1|cluster:year) +  (1|cluster:location:year) +
      (1|line_name) + (1|line_name:year) + (1|line_name:cluster) + (1|line_name:cluster:location) + (1|line_name:cluster:year) +
      (1|line_name:cluster:location:year)
    
    fit <- lmer(formula, df, control = control)
    
    ## Calculate heritability before and after clustering
    plot_table <- xtabs(formula = ~ line_name + location + year + cluster, data = df)
    
    ## Harmonic means
    # Locations
    harm_loc <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Year
    harm_year <- apply(X = plot_table, MARGIN = c(2,3), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Reps
    harm_rep <- apply(X = plot_table, MARGIN = c(1,2,3), sum) %>% 
      harm_mean()
    
    ## Clusters
    harm_clust <- n_distinct(df$cluster)
    
    
    ## Calculate heritability
    exp <- "line_name / (line_name + (line_name:cluster / n_c) + (line_name:cluster:location / (n_l * n_c)) + 
(line_name:year / (n_y)) + (line_name:cluster:year / (n_c * n_y)) + (line_name:cluster:location:year / (n_l * n_y * n_c)) +
(Residual / (n_r * n_l * n_y * n_c)))"
    H <- herit(object = fit, exp = exp, n_c = harm_clust, n_y = harm_year, n_r = harm_rep, n_l = harm_loc)
    
    ## Heritability of subdivided data (Yan2016)
    exp <- "(line_name + line_name:cluster) / (line_name + line_name:cluster + n_c * ( (line_name:cluster:location / (n_l)) + 
(line_name:cluster:location:year / (n_l * n_y)) + (Residual / (n_r * n_l * n_y))) + (line_name:year / (n_y)) + 
(line_name:cluster:year / (n_c * n_y)))"
    H1 <- herit(object = fit, exp = exp, n_c = harm_clust, n_y = harm_year, n_r = harm_rep, n_l = harm_loc)
    

    
    # # Random anova
    # fit_ranova <- ranova(model = fit) %>%
    #   broom::tidy() %>% 
    #   filter(!str_detect(term, "none")) %>% 
    #   mutate(term = str_remove_all(term, "\\(1 \\| |\\)"))

    
    ## Return a DF
    data_frame(undivided = list(H), subdivided = list(H1))
    
  }) %>% ungroup()
    
single_cluster_cases <- cluster_df_tomodel %>% 
  group_by(set, model, trait, method) %>%
  filter(n_distinct(cluster) == 1) %>%
  distinct(set, model, trait, method)


## Extract heritability
cluster_herit <- cluster_varcomp %>%
  mutate_at(vars(undivided, subdivided), ~map_dbl(., "heritability")) %>%
  rename(clustered = subdivided, unclustered = undivided) %>%
  mutate(model = ifelse(model  == "pheno_location_dist", "pheno_loc_dist", model),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr_use)) %>%
  filter(model %in% dist_method_abbr_use)

## Calculate ratio
cluster_herit %>%
  mutate(effectiveness = clustered / unclustered) %>%
  select(-unclustered, -clustered) %>%
  spread(model, effectiveness)

# set           trait           method  AMMI     PD LocPD   GCD `All-EC` `Mean-EC` `IPCA-EC`
# 1 complete      GrainYield      mclust  4.93  1.11  0.595 0.534    0.482     0.558     0.764
# 2 complete      GrainYield      pam     4.93  2.00  0.459 0.445    0.408     0.558     0.734
# 3 complete      HeadingDate     mclust  1.14  1.01  0.952 0.874    0.861     0.855     0.930
# 4 complete      HeadingDate     pam     1.14  0.953 1.01  0.879    0.839     0.878     0.842
# 5 complete      HeadingDateAGDD mclust NA     1.01  0.915 0.881   NA        NA        NA    
# 6 complete      HeadingDateAGDD pam    NA     0.951 0.927 0.887   NA        NA        NA    
# 7 complete      PlantHeight     mclust  1.05  1.11  0.716 0.613    0.564     0.657    NA    
# 8 complete      PlantHeight     pam     1.05  1.11  0.504 0.561    0.522     0.749     0.632
# 9 realistic2017 GrainYield      mclust  2.15 NA     0.742 0.654    0.636     0.634     0.607
# 10 realistic2017 GrainYield      pam     2.15 NA     1.13  0.797    0.855     0.870     0.776
# 11 realistic2017 HeadingDate     mclust NA    NA     0.983 0.928    0.927     0.933     0.943
# 12 realistic2017 HeadingDate     pam    NA    NA     0.994 0.919    0.925     0.933     0.928
# 13 realistic2017 HeadingDateAGDD mclust NA    NA     0.978 0.924   NA        NA        NA    
# 14 realistic2017 HeadingDateAGDD pam    NA    NA     0.982 0.925   NA        NA        NA    
# 15 realistic2017 PlantHeight     mclust NA    NA     0.778 0.744    0.728     0.801     0.755
# 16 realistic2017 PlantHeight     pam    NA    NA     0.917 0.745    0.749     0.808     0.853
  

## Plot heritability
g_cluster_herit <- cluster_herit %>%
  mutate(set = str_replace_all(set, set_replace),
         set = str_replace_all(set, "[0-9]{4}", "")) %>%
  filter(!(model == "AMMI" & set == "Time-forward")) %>%
  gather(group, heritability, -set, -model, -trait, -method) %>%
  mutate(group = ifelse(group == "unclustered", paste0("H1: ", group), paste0("H2: ", group))) %>%
  ggplot(aes(x = model, y = heritability, color = method, shape = group)) +
  geom_point(position = position_dodge2(0.5), size = 1.5) +
  # scale_color_manual(values = dist_colors_use, name = "Distance\nmethod", guide = FALSE) +
  scale_shape_discrete(name = NULL, labels = str_to_title, guide = guide_legend(nrow = 2)) +
  facet_grid(trait ~ set, scales = "free_x", space = "free_x", labeller = labeller(trait = str_add_space), switch = "y") +
  xlab("Distance measure") +
  scale_y_continuous(breaks = pretty, name = "Heritability") +
  theme_presentation2(base_size = 10) +
  theme(legend.position = c(0.22, 0.07), legend.margin = margin(), legend.key.height = unit(0.5, "line"), 
        axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside")

ggsave(filename = "cluster_heritability_complete_withAGDD.jpg", plot = g_cluster_herit, path = fig_dir, width = 3.5, height = 6, dpi = 1000)


## Remove AGDD
g_cluster_herit <- cluster_herit %>%
  mutate(set = str_replace_all(set, set_replace),
         set = str_replace_all(set, "[0-9]{4}", "")) %>%
  filter(!(model == "AMMI" & set == "Time-forward")) %>%
  filter(trait %in% traits) %>%
  filter(method == "mclust") %>%
  gather(group, heritability, -set, -model, -trait, -method) %>%
  mutate(group = ifelse(group == "unclustered", paste0("H1: ", group), paste0("H2: ", group))) %>%
  ggplot(aes(x = model, y = heritability, color = model, shape = group)) +
  # ggplot(aes(x = model, y = heritability, color = method, shape = group)) +
  geom_point(position = position_dodge2(0.5), size = 1.5) +
  scale_color_manual(values = dist_colors_use, name = "Distance\nmethod", guide = FALSE) +
  scale_shape_discrete(name = NULL, labels = str_to_title, guide = guide_legend(nrow = 2)) +
  facet_grid(trait ~ set, scales = "free_x", space = "free_x", labeller = labeller(trait = str_add_space), switch = "y") +
  xlab("Distance measure") +
  scale_y_continuous(breaks = pretty, name = "Heritability") +
  theme_presentation2(base_size = 10) +
  theme(legend.position = c(0.22, 0.07), legend.margin = margin(), legend.key.height = unit(0.5, "line"), 
        axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside")

ggsave(filename = "cluster_heritability_complete.jpg", plot = g_cluster_herit, path = fig_dir, width = 3.5, height = 5, dpi = 1000)




cluster_varcomp1 <- cluster_varcomp %>%
  mutate_at(vars(undivided, subdivided), ~map(., "var_comp")) %>%
  mutate(model = ifelse(model  == "pheno_location_dist", "pheno_loc_dist", model),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr_use),
         set = str_to_title(set)) %>%
  filter(model %in% dist_method_abbr_use) %>%
  select(-subdivided) %>%
  unnest() %>%
  ## Calculate variance proportions
  group_by(set, model, trait) %>%
  mutate(varprop = variance / sum(variance)) 



  



## Plot
g_cluster_varcomp <- cluster_varcomp1 %>%
  ggplot(aes(x = model, y = varprop, fill = source)) +
  geom_col() +
  facet_grid(trait ~ set, scales = "free_x", space = "free_x") +
  theme_presentation2() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1))
  


# Save a table
cluster_varcomp_table <- cluster_varcomp1 %>%
  filter(method == "mclust") %>%
  mutate(annotation = list(NULL)) %>%
  mutate(annotation = str_trim(paste0(formatC(x = varprop, digits = 2))),
         annotation = str_remove_all(annotation, "NA")) %>%
  select(set, trait, model, source, annotation) %>%
  spread(source, annotation)

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
  filter(trait %in% traits) %>%
  mutate_at(vars(environment:line_name), as.factor)


## For each prediction environment (the tp+vp envs and just the vp envs), rank the 
## training environments by different distance metrics
pred_envs <- c(tp_vp_env, vp_only_env)
train_envs <- c(tp_vp_env, tp_only_env)


# Summarize the traits available in those environments
val_envs_traits <- S2_MET_BLUEs_use %>%
  filter(environment %in% pred_envs) %>% 
  filter(trait %in% traits) %>%
  group_by(environment) %>% 
  distinct(trait) %>%
  ungroup() %>%
  nest(environment, .key = "val_environments") %>%
  mutate(set = "complete") %>%
  bind_rows(., cluster_df_realistic %>% filter(!map_lgl(test_env, is.null)) %>% subset(., model == "great_circle_dist", c(set, trait, test_env)) %>%
              mutate(val_environments = map(test_env, ~data.frame(environment = .))) %>% select(-test_env))


train_envs_traits <- S2_MET_BLUEs_use %>%
  filter(environment %in% train_envs) %>% 
  group_by(environment) %>% 
  distinct(trait) %>%
  ungroup() %>%
  nest(environment, .key = "train_environments") %>%
  full_join(., val_envs_traits) %>%
  mutate(train_environments = map2(train_environments, val_environments, ~dplyr::setdiff(.x, .y)),
         train_environments = ifelse(set == "complete", val_environments, train_environments)) %>%
  select(-val_environments)
  

## Rank the environments relative to each other
# Do this for both population groups
env_rank_df <- dist_method_df %>% 
  filter(!(str_detect(set, "realistic") & model == "AMMI")) %>%
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



## For each validation environment and model, what proportion of the top 5 environments 
## are from the same location or year?
nE <- 5

rank_prop <- pred_env_dist_rank %>%
  filter(model %in% names(dist_method_abbr), ! mat_set %in% c("Jarquin", "MalosettiStand") ) %>%
  group_by(trait, set, validation_environment, model) %>%
  do({
    eRank <- .$rank[[1]]
    eRank5 <- head(eRank, nE)
    ve <- .$validation_environment[1]
    
    ## Proportion from same year
    nYear <- sum(str_extract(eRank5, "[0-9]{2}") == str_extract(ve, "[0-9]{2}"))
    ## Proportion from same year expected by chance
    pYearChance <- prop.table(table(str_extract(eRank, "[0-9]{2}")))[str_extract(ve, "[0-9]{2}")]
    pYearChance <- ifelse(is.na(pYearChance), 0, pYearChance)
    
    
    ## Proportion from same location 
    nLoc <- sum(str_extract(eRank5, "[A-Z]{3}") == str_extract(ve, "[A-Z]{3}"))
    ## Proportion from same location expected by chance
    pLocChance <- prop.table(table(str_extract(eRank, "[A-Z]{3}")))[str_extract(ve, "[A-Z]{3}")]
    pLocChance <- ifelse(is.na(pLocChance), 0, pLocChance)
    
    ## Return
    data.frame(test = c("year", "location"), obs = c(nYear, nLoc), exp = c(pYearChance, pLocChance), row.names = NULL, stringsAsFactors = FALSE) %>%
      mutate(p_value = map2_dbl(obs, exp, ~binom.test(x = .x, n = nE, p = .y / 5)$p.value), pObs = obs / nE)
    
  })


rank_prop %>%
  group_by(trait, set, model, test) %>% 
  summarize_at(vars(pObs, exp, p_value), mean, na.rm = T) %>% 
  as.data.frame() %>%
  arrange(trait, set, test, desc(pObs)) %>% filter(test == "location")


## Boxplot of heritability from different years






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






# Save this
save_file <- file.path(result_dir, "distance_method_results.RData")
save("env_rank_df", "cluster_df", "cluster_df_pam", "cluster_varcomp", "pred_env_dist_rank","pred_env_rank_random", file = save_file)













