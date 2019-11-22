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


## How many of the top genotypes should be used for clustering?
top_geno <- 1


## Use fitted values of AMMI model to cluster environments based on the best lines in each
ammi_clusters_df <- ammi_out %>%
  filter(trait %in% traits) %>%
  # Add a indicator of whether higher or lower trait values are favorable
  mutate(high_favorable = trait == "GrainYield") %>% 
  group_by(set, trait) %>% 
  do({
    row <- .
    
    ## Filter for PC1
    scores <- row$ammi[[1]][c("escores", "gscores")] %>%
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
    
    ## Determine the most favorable genotypes in each environment
    decreasing <- row$high_favorable
    lines_ordered <- apply(X = fitted_effect, MARGIN = 2, 
                           FUN = function(x) row.names(fitted_effect)[order(x, decreasing = decreasing)[seq(top_geno)]] )
    # Tranpose if dim == null
    lines_ordered <- if (is.null(dim(lines_ordered))) t(lines_ordered) else lines_ordered
    
    ## Assign clusters
    clusters <- apply(X = lines_ordered, MARGIN = 2, FUN = paste0, collapse = "_") %>% 
      as.factor() %>% 
      tibble(environment = names(.), cluster = .) %>% 
      mutate(cluster = as.numeric(cluster))
    
    ## If there are singleton clusters, add the environment to the most similar cluster
    clusters_new <- clusters
    cluster_table <- table(clusters_new$cluster)
    any_singletons <- any(cluster_table == 1)
    while(any_singletons) {
      
      # Find the singleton clusters
      singleton_clusters <- names(which(cluster_table == 1))
      
      # Find the first singleton environment
      singleton_environments <- subset(clusters_new, cluster %in% singleton_clusters, environment, drop = T)[1]
      
      # Add that environment to the closest cluster based on IPCA1 score
      clusters_new_ipca <- subset(clusters_new, environment != singleton_environments) %>%
        left_join(.,  rownames_to_column(as.data.frame(e_scores), "environment"), by = "environment") %>%
        group_by(cluster) %>%
        summarize(score = mean(score))
      
      which_cluster <- clusters_new_ipca %>% 
        mutate(diff = e_scores[singleton_environments,] - score) %>% 
        top_n(x = ., n = 1, wt = -diff) %>%
        pull(cluster)
      
      ## Add the environment to the cluster
      clusters_new$cluster[clusters_new$environment == singleton_environments] <- which_cluster
      # Find singletons
      cluster_table <- table(clusters_new$cluster)
      any_singletons <- any(cluster_table == 1)
      
    }
    
    
    
    ## Return the fitted values, the top, and bottom
    fitted_effect1 <- as.data.frame(fitted_effect) %>% 
      rownames_to_column("line_name") %>% 
      gather(environment, fitted, -line_name)
    tibble(fitted = list(fitted_effect1), clusters = list(clusters_new))
    
  }) %>% ungroup()
    


## How many mega-environments per trait and set
ammi_clusters_df %>%
  unnest(clusters) %>%
  group_by(set, trait, cluster) %>%
  mutate(n_env = n_distinct(environment)) %>%
  group_by(set, trait) %>%
  summarize(n_cluster = n_distinct(cluster),
            min_env_cluster = min(n_env),
            max_env_cluster = max(n_env))
    
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
ammi_clusters <- ammi_clusters_df %>%
  select(set, trait, cluster = clusters) %>%
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




## Create the model using an abbreviation for the ec_group (multi-year) and the group (All, Mean. IPCA)
ec_sim_mat_df1 <- rel_mat_df %>%
  select(-mat) %>%
  mutate(model = str_replace_all(ec_group, "multiyear", "multi year") %>% str_to_title() %>% abbreviate(2),
         model = str_c(model, group),
         dist = map(sim_mat, dist)) %>%
  select(set, model, trait, cov = sim_mat, dist)


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
  bind_rows(., tibble(trait = unique(.$trait), model = "great_circle_dist", dist = great_circle_dist_list, 
                      cov = map(great_circle_dist_list, ~1 - as.matrix(.))) %>% 
              crossing(set = unique(location_BLUEs_realistic_df$set), .)) %>%
  bind_rows(., location_BLUEs_realistic_df %>% select(set, trait, cov = cor, dist) %>% mutate(model = "pheno_loc_dist")) %>%
  filter(set != "complete")
  

## Combine
## Remove the distance matrices
dist_method_df <- bind_rows(ammi_distance, dist_method_df_complete, dist_method_df_realistic) %>%
  # filter(set %in% c("complete", "realistic2017")) %>%
  filter(trait %in% traits)






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





############################
## Use hierarchical clustering instead
############################

# cluster_df_complete <- dist_method_df_complete %>%
#   # filter(model == "pheno_dist", trait == trait[1]) %%>
#   mutate(clust = map(dist, ~agnes(., method = "ward")),
#          K_poss = map(clust, ~seq(2, length(.x$order) - 1)),
#          tree = map2(clust, K_poss, function(c, k) map(k, ~cutree(c, k = .x))),
#          silhouette = map2(tree, dist, function(t, d) map(t, ~silhouette(x = .x, dist = d))),
#          sil_width = map(silhouette, ~map_dbl(., ~mean(.[,"sil_width"]))),
#          # Which trees have clusters with at least 2 environments?
#          good_trees = map(tree, ~map_lgl(.x, ~all(table(.) >= 2))),
#          k = pmap_dbl(list(K_poss, sil_width, good_trees), function(x, y, z) x[which(y == max(y[z]))] ))




# Combine
cluster_df_complete <- bind_rows(gcd_mat_complete, pd_mat_complete, pld_mat_complete, ec_mat_complete, 
                                 filter(ammi_clusters, set == "complete")) %>%
  arrange(trait, model) %>%
  mutate(test_env = map(cluster, "environment"))





## LOYO


# GCD
gcd_mat_realistic <- distinct(dist_method_df_realistic, set) %>% 
  crossing(., tibble(model = "great_circle_dist", trait = names(gcd_mat_list), data = gcd_mat_list)) %>% 
  mutate(test_env = map2(data, set, ~row.names(.x) %>% str_subset(., str_extract(.y, "1[0-9]"))),
         cluster = map2(.x = data, .y = test_env, ~env_mclust(data = .x, min_env = min_env, test.env = .y)))

# Location dist
pld_mat_realistic  <- gcd_mat_realistic %>% 
  distinct(set, trait, test_env) %>% 
  left_join(., location_BLUEs_realistic_df %>% mutate(model = "pheno_location_dist") %>% 
              select(model, set, trait, data = dist)) %>%
  mutate(data = map(data, cmdscale),
         test_env = map2(.x = data, .y = test_env, ~intersect(row.names(.x), .y)),
         cluster = map2(.x = data, .y = test_env, ~env_mclust(data = .x, min_env = min_env, test.env = .y)))


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



# Combine the loyo clustering results
cluster_df_realistic <- bind_rows(gcd_mat_realistic, pld_mat_realistic, ec_mat_realistic) %>%
  # Make sure the test environments are all consistent across sets within traits
  group_by(trait, set) %>%
  do({
    df <- .
    
    ## Make sure test_envs are consistent across sets
    test_env_use <- reduce(df$test_env, intersect)
    ## Make sure clusters are consistent
    cluster_env <- df$cluster %>% map("environment") %>% reduce(intersect)
    # Filter clusters
    cluster_use <- map(df$cluster, ~filter(., environment %in% cluster_env))
    
    ## Make sure all test_envs are clustered
    test_env_use1 <- intersect(test_env_use, cluster_env)
    cluster_use1 <- map(cluster_use, ~filter(., environment %in% c(setdiff(cluster_env, test_env_use1), test_env_use1)))
    
    
    # Repackage and export
    mutate(df, test_env = list(test_env_use1), cluster = cluster_use1) %>%
      select(-trait, -set)
    
  }) %>% ungroup()


## Combine all
cluster_df <- bind_rows(cluster_df_complete, cluster_df_realistic)





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
cluster_df_tomodel <- cluster_df %>%  
  # filter(set %in% c("complete", "realistic2017")) %>%
  unnest(cluster) %>% 
  left_join(filter(S2_MET_BLUEs, line_name %in% tp))

## Control for lmer
control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")

## Option for just ECs
cluster_varcomp_ec <- cluster_df_tomodel %>% 
  filter(str_detect(model, "MYEC")) %>%

# cluster_varcomp <- cluster_df_tomodel %>% 
  group_by(set, model, trait) %>%
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
  group_by(set, model, trait) %>%
  filter(n_distinct(cluster) == 1) %>%
  distinct(set, model, trait)

## Recombine varcomp
cluster_varcomp <- cluster_varcomp %>% 
  anti_join(., select(cluster_varcomp_ec, set:trait)) %>% 
  bind_rows(., cluster_varcomp_ec)


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
 



### Rank environments according to a test environment
### This will be used for prediction
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
# # First for test environments
val_envs_traits <- S2_MET_BLUEs_use %>%
  filter(environment %in% pred_envs) %>% 
  filter(trait %in% traits) %>%
  group_by(environment) %>% 
  distinct(trait) %>%
  ungroup() %>%
  nest(environment, .key = "val_environments") %>%
  mutate(set = "complete") %>%
  bind_rows(., cluster_df_realistic %>% filter(!map_lgl(test_env, is.null)) %>% 
              subset(., model == "great_circle_dist", c(set, trait, test_env)) %>%
              mutate(val_environments = map(test_env, ~data.frame(environment = .))) %>% 
              select(-test_env))


# Now for training environments
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






## Randomly order non-prediction environments
# Number of random samples
n_sample <- 10

set.seed(153)
## For each prediction environment, take the environments in one of the distance
## matrices and randomly sample it
pred_env_rank_random <- pred_env_dist_rank %>%
  group_by(set, trait, validation_environment) %>%
  do({
    df <- .
    smpls <- rerun(.n = n_sample, sample(df$rank[[1]]))
    data_frame(validation_environment = df$validation_environment[1], trait = df$trait[1],
               set = df$set[1], model = str_c("rank_sample", seq_along(smpls)), 
               rank = smpls)
    
  }) %>% ungroup()






# Save this
save_file <- file.path(result_dir, "distance_method_results.RData")
save("env_rank_df", "cluster_df", "cluster_varcomp", "pred_env_dist_rank","pred_env_rank_random", file = save_file)













