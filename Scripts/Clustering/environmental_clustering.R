## Cluster environments in the S2MET project
## 
## Author: Jeff Neyhart
## Last modified: October 04, 2018
## 
## This script will create different distance matrices for use in clustering. The
## clustering algorithm will be constant
## 

# Load packages and the source script
library(tidyverse)
library(readxl)
library(pbr)
library(ggdendro)
library(cowplot)

# The head directory
repo_dir <- getwd()

source(file.path(repo_dir, "source.R"))



## Load the environment covariable distance matrices
load(file.path(result_dir, "environmental_covariable_distance_mat.RData"))

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




# ## Genetic correlation of the TP or all lines
# ## This is deprecated
# 
# # Load the environmental genetic correlation matrices
# load(file.path(result_dir, "environmental_genetic_correlations.RData"))
# 
# # Convert each correlation matrix to a distance matrix
# # Missing data gets a correlation of 0
# # 1 - cor is the distance
# env_cor_all_dist_list <- env_cor_all %>% 
#   map(function(mat) {mat[is.na(mat)] <- 0; return(mat)}) %>%
#   map(~1 - .) %>%
#   map(as.dist)
# 
# env_cor_tp_dist_list <- env_cor_tp %>%
#   map(function(mat) {mat[is.na(mat)] <- 0; return(mat)}) %>%
#   map(~1 - .) %>%
#   map(as.dist)


## Phenotypic distance metric
# This is described by Bernardo2010

# First calculate the distance using all data
ge_mean_D_all <- S2_MET_BLUEs %>%
  # Split by trait
  split(.$trait) %>%
  map(function(df){
    dist1 <- dist_env(x = df, gen.col = "line_name", env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })

# Now calculate it based only on data from the TP
ge_mean_D_tp <- S2_MET_BLUEs %>%
  filter(line_name %in% tp) %>%
  # Split by trait
  split(.$trait) %>%
  map(function(df){
    dist1 <- dist_env(x = df, gen.col = "line_name", env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })


# ## Distance based on SREG
# Deprecated
# 
# ## Iterate over the pairs of environments
# # Fit the SREG model, then calculate the residual sum of squares
# env_pair_list <- combn(x = sort(as.character(unique(pheno$environment))), m = 2, function(env_pair) { 
#   
#   # Subset the phenotype data
#   pheno_sub <- filter(pheno, environment %in% env_pair)
#   
#   # Check the dataset balance
#   pheno_check <- pheno_sub %>%
#     group_by(line_name) %>%
#     summarize(obs = n())
#   
#   ## Check that at least two genotypes are evaluated in both environments
#   if (sum(pheno_check$obs == 2) >= 2) {
# 
#     fit <- sreg(formula = value ~ line_name + environment + line_name:environment, 
#                 data = subset(pheno, environment %in% env_pair), gen.col = "line_name", 
#                 env.col = "environment")
#     rss <- sum(fit$effects$residual^2)
#     mse <- mean(fit$effects$residual^2)
#     
#   } else {
#     rss <- mse <- NA
#   
#   }
#   
#   data.frame(env1 = env_pair[1], env2 = env_pair[2], rss = rss, mse = mse, row.names = NULL, 
#              stringsAsFactors = FALSE)
#   
# }, simplify = FALSE)


# ## PCA of GxE BLUPs
# ## Deprecated
# 
# # Fit a new model with compound symmetry and get the GxE BLUPs
# S2_MET_BLUPs_ge_all <- S2_MET_BLUEs %>%
#   split(.$trait) %>%
#   map(function(df) {
# 
#     # Create an object for the data.frame
#     ## Center and scale the data
#     df <- mutate(df, value = scale(value))
# 
#     # Extract the standard errors
#     wts <- df$std_error^2
# 
#     # lmer control
#     lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore",
#                                 calc.derivs = FALSE)
# 
#     # Fit the model
#     form <- value ~ (1|line_name) + (1|environment) + (1|line_name:environment)
#     fit <- lmer(formula = form, data = df, control = lmer_control, weights = wts)
# 
#     # Extract the model.frame
#     mf <- model.frame(fit)
# 
#     # Extract the BLUPs of each genotype in each environment
#     blups <- ranef(fit)$`line_name:environment` %>%
#       rownames_to_column("term") %>%
#       separate(term, c("line_name", "environment"), sep = ":") %>%
#       rename(value = "(Intercept)")
# 
#     # Extract the variance components
#     var_comp <- as.data.frame(VarCorr(fit))
# 
#     # Return a data.frame
#     data_frame(fit = list(fit), var_comp = list(var_comp), BLUP = list(blups)) })
# 
# 
# # Fit a new model with compound symmetry and get the GxE BLUPs
# # Use only the TP data
# S2_MET_BLUPs_ge_tp <- S2_MET_BLUEs %>%
#   filter(line_name %in% tp_geno) %>%
#   split(.$trait) %>%
#   map(function(df) {
#     
#     # Create an object for the data.frame
#     ## Center and scale the data
#     df <- mutate(df, value = scale(value))
#     
#     # Extract the standard errors
#     wts <- df$std_error^2
#     
#     # lmer control
#     lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore",
#                                 calc.derivs = FALSE)
#     
#     # Fit the model
#     form <- value ~ (1|line_name) + (1|environment) + (1|line_name:environment)
#     fit <- lmer(formula = form, data = df, control = lmer_control, weights = wts)
#     
#     # Extract the model.frame
#     mf <- model.frame(fit)
#     
#     # Extract the BLUPs of each genotype in each environment
#     blups <- ranef(fit)$`line_name:environment` %>%
#       rownames_to_column("term") %>%
#       separate(term, c("line_name", "environment"), sep = ":") %>%
#       rename(value = "(Intercept)")
#     
#     # Extract the variance components
#     var_comp <- as.data.frame(VarCorr(fit))
#     
#     # Return a data.frame
#     data_frame(fit = list(fit), var_comp = list(var_comp), BLUP = list(blups)) })

# # Factor analysis using BLUEs
# ge_mean_FA <- S2_MET_BLUEs_use %>%
#   group_by(trait) %>%
#   mutate(BLUEs = list({
#     data[[1]] %>%
#       select(line_name, environment, value) %>%
#       complete(line_name, environment) %>%
#       group_by(environment) %>%
#       mutate(value = ifelse(is.na(value), mean(value, na.rm = TRUE), value)) %>%
#       spread(environment, value) %>%
#       as.data.frame() %>%
#       remove_rownames() %>%
#       column_to_rownames("line_name") %>%
#       as.matrix() })) %>%
#   do(fa_out = fa(r = .$BLUEs[[1]], nfactors = 2, rotate = "varimax")) %>%
#   # Grab the loadings
#   mutate(delta = list(structure(fa_out$loadings, class = "matrix"))) %>%
#   # Distance matrix
#   mutate(FA = list(dist(delta)))


# # Extract the GxE BLUPs and run PCA
# ge_PCA_dist_all <- S2_MET_BLUPs_ge_all %>%
#   map_df(function(out) {
#     blup_mat <- out$BLUP[[1]] %>%
#       spread(environment, value) %>%
#       as.data.frame() %>%
#       remove_rownames() %>%
#       column_to_rownames("line_name") %>%
#       as.matrix() 
#     
#     # Impute with mean
#     blup_prcomp <- blup_mat %>%
#       apply(MARGIN = 2, FUN = function(env) ifelse(is.na(env), mean(env, na.rm = T), env)) %>%
#       t() %>%
#       prcomp(center = TRUE, scale. = TRUE)
#     
#     # Extract all the PCs
#     blup_pcs <- blup_prcomp$x
#     # Subset the ones to use
#     blup_pcs_use <- blup_pcs[,1:n_PC,drop = FALSE]
#     # Calculate the distance between points using the PCs
#     blup_pc_dist <- dist(blup_pcs_use)
#     
#     # Return a data.frame
#     data_frame(blup_mat = list(blup_mat), blup_prcomp = list(blup_prcomp), 
#                blup_pcs = list(blup_pcs), blup_pc_dist = list(blup_pc_dist))
#     
#   })
#     
# # Copy for the tp
# ge_PCA_dist_tp <- S2_MET_BLUPs_ge_tp %>%
#   map_df(function(out) {
#     blup_mat <- out$BLUP[[1]] %>%
#       spread(environment, value) %>%
#       as.data.frame() %>%
#       remove_rownames() %>%
#       column_to_rownames("line_name") %>%
#       as.matrix() 
#     
#     # Impute with mean
#     blup_prcomp <- blup_mat %>%
#       apply(MARGIN = 2, FUN = function(env) ifelse(is.na(env), mean(env, na.rm = T), env)) %>%
#       t() %>%
#       prcomp(center = TRUE, scale. = TRUE)
#     
#     # Extract all the PCs
#     blup_pcs <- blup_prcomp$x
#     # Subset the ones to use
#     blup_pcs_use <- blup_pcs[,1:n_PC,drop = FALSE]
#     # Calculate the distance between points using the PCs
#     blup_pc_dist <- dist(blup_pcs_use)
#     
#     # Return a data.frame
#     data_frame(blup_mat = list(blup_mat), blup_prcomp = list(blup_prcomp), 
#                blup_pcs = list(blup_pcs), blup_pc_dist = list(blup_pc_dist))
#     
#   })



# ## Environmental covariates
# # Use the PCs of the environmental covariates to form clusters
# EC_one_PCA <- prcomp(x = one_year_env_mat, center = TRUE, scale. = TRUE)
# EC_multi_PCA <- prcomp(x = multi_year_env_mat, center = TRUE, scale = TRUE)
# 
# ## Tidy up
# EC_one_PCA_df <- bind_rows(
#   broom::tidy(x = EC_one_PCA, matrix = "u") %>% rename(label = row),
#   broom::tidy(x = EC_one_PCA, matrix = "v") %>% rename(label = column)
# ) %>% 
#   mutate(PC = str_c("PC", PC),
#          type = ifelse(str_detect(label, "[A-Z]{3}[0-9]{2}"), "env", "cov"))
# 
# EC_multi_PCA_df <- bind_rows(
#   broom::tidy(x = EC_multi_PCA, matrix = "u") %>% rename(label = row),
#   broom::tidy(x = EC_multi_PCA, matrix = "v") %>% rename(label = column)
# ) %>% 
#   mutate(PC = str_c("PC", PC),
#          type = ifelse(str_detect(label, "[A-Z]{3}[0-9]{2}"), "env", "cov"))
#   
# 
# ## Visualize briefly
# EC_one_PCA_df %>% 
#   filter(PC %in% c("PC1", "PC2")) %>% 
#   spread(PC, value) %>% 
#   ggplot(aes(x = PC1, y = PC2, col = type)) +
#   geom_point() + 
#   geom_text(aes(label = label), check_overlap = T) + 
#   theme_classic()
# 
# 
# 
# # Copy the list per trait (the distance is the same for each trait)
# # EC_one_PCA_dist_list <- rerun(length(traits), dist(EC_one_PCA$x[,1:n_PC, drop = FALSE])) %>% 
# #   set_names(traits)
# EC_multi_PCA_dist_list <- rerun(length(traits), dist(EC_multi_PCA$x[,1:n_PC, drop = FALSE])) %>% 
#   set_names(traits)


## Combine the one-year and multi-year ECs
one_year_ec_dist1 <- one_year_ec_dist %>%
  mutate(model = str_c("OYEC_", group))
multi_year_ec_dist1 <- multi_year_ec_dist %>%
  mutate(model = str_c("MYEC_", group))

# Combine and add the other distance matrices
dist_method_df_all <- bind_rows(
  data_frame(trait = traits, model = "great_circle_dist", dist = great_circle_dist_list),
  data_frame(trait = traits, model = "pheno_dist", dist = ge_mean_D_all),
  filter(one_year_ec_dist1, population == "all"),
  filter(multi_year_ec_dist1, population == "all")
)

# Combine and add the other distance matrices
dist_method_df_tp <- bind_rows(
  data_frame(trait = names(ge_mean_D_tp), model = "great_circle_dist", dist = great_circle_dist_list),
  data_frame(trait = names(ge_mean_D_tp), model = "pheno_dist", dist = ge_mean_D_tp),
  filter(one_year_ec_dist1, population == "tp"),
  filter(multi_year_ec_dist1, population == "tp")
)


# 
# ## Combine the lists into a data.frame
# dist_method_df_all <- data_frame(
#   trait = traits, 
#   great_circle_dist = great_circle_dist_list, 
#   # env_cor_dist = env_cor_all_dist_list, 
#   ge_mean_D = ge_mean_D_all, 
#   # ge_PCA_dist = ge_PCA_dist_all$blup_pc_dist,
#   # ec_one_PCA_dist = EC_one_PCA_dist_list,
#   ec_multi_PCA_dist = EC_multi_PCA_dist_list)
# 
# 
# dist_method_df_tp <- data_frame(
#   trait = traits, 
#   great_circle_dist = great_circle_dist_list, 
#   # env_cor_dist = env_cor_tp_dist_list, 
#   ge_mean_D = ge_mean_D_tp, 
#   # ge_PCA_dist = ge_PCA_dist_tp$blup_pc_dist,
#   # ec_one_PCA_dist = EC_one_PCA_dist_list,
#   ec_multi_PCA_dist = EC_multi_PCA_dist_list)


# Tidy the the distance matrix data.frame, then combine the data.frames for
# TP and TP + VP
dist_method_df_all_tidy <- dist_method_df_all %>% 
  mutate(population = "all") %>%
  select(-group)

# For each trait, identify the most common set of environments
dist_method_df_common_env <- dist_method_df_all_tidy %>% 
  group_by(trait) %>% 
  summarize(common_env = list(map(dist, ~row.names(as.matrix(.))) %>% reduce(intersect)))

# Combine with the distance metrics, convert to a matrix, subset the environments,
# then convert back to a dist object
# Then create cluster objects
clust_method_df_all <- full_join(dist_method_df_all_tidy, dist_method_df_common_env, by = "trait") %>% 
  mutate(dist = list(dist, common_env) %>% pmap(~subset_env(dist = .x, envs = .y)),
         cluster = map(dist, ~hclust(., method = "ward.D"))) %>%
  select(-common_env)


# Just the TP
dist_method_df_tp_tidy <- dist_method_df_tp %>% 
  mutate(population = "tp") %>%
  select(-group)

# For each trait, identify the most common set of environments
dist_method_df_common_env <- dist_method_df_tp_tidy %>% 
  group_by(trait) %>% 
  summarize(common_env = list(map(dist, ~row.names(as.matrix(.))) %>% reduce(intersect)))

# Combine with the distance metrics, convert to a matrix, subset the environments,
# then convert back to a dist object
# Then create cluster objects
clust_method_df_tp <- full_join(dist_method_df_tp_tidy, dist_method_df_common_env, by = "trait") %>% 
  mutate(dist = list(dist, common_env) %>% pmap(~subset_env(dist = .x, envs = .y)),
         cluster = map(dist, ~hclust(., method = "ward.D"))) %>%
  select(-common_env)

## Combine the data.frames
## Replicate the clusters and add each cluster k from min_k to max_k
clust_method_df <- bind_rows(clust_method_df_all, clust_method_df_tp)


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
            geom_text(size = 2) + 
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
pred_envs_traits <- S2_MET_BLUEs_use %>%
  filter(environment %in% pred_envs) %>% 
  group_by(environment) %>% 
  distinct(trait) %>%
  ungroup() %>%
  mutate(environment = as.character(environment))

train_envs_traits <- S2_MET_BLUEs_use %>%
  filter(environment %in% train_envs) %>% 
  group_by(environment) %>% 
  distinct(trait) %>%
  ungroup()


## Rank the environments relative to each other
# Do this for both population groups
dist_method_df <- clust_method_df %>%
  split(.$population) %>%
  map(~select(., -cluster) %>%
        mutate(dist = map(dist, as.matrix)))

# Linearize the distance objects and rank
dist_method_df_rank <- dist_method_df %>% 
  map(~mutate(., dist_rank = map(dist, ~as.data.frame(.) %>% 
                           rownames_to_column("env") %>% 
                           filter(env %in% pred_envs) %>%
                           split(.$env) %>% 
                           map(~.[,!names(.) %in% .$env] %>% .[,-1] %>% sort() %>% 
                                 select(., which(names(.) %in% train_envs))) %>% 
                           data_frame(environment = names(.), env_rank = .))) %>%
  select(-dist, -cov) %>%
  unnest() )


# Combine this data with the the trait information for the prediction environments
pred_env_dist_rank <- dist_method_df_rank %>%
  map(~inner_join(pred_envs_traits, .,  by = c("environment", "trait")))




## Randomly sample environments as a control
# Number of random samples
n_sample <- 100

## For each prediction environment, take the environments in one of the distance
## matrices and randomly sample it
pred_env_rank_random <- pred_env_dist_rank %>%
  map(~group_by(., environment, trait) %>%
        do({
          df <- .
          smpls <- rerun(.n = n_sample, sample(df$env_rank[[1]]))
          data_frame(environment = df$environment[1], trait = df$trait[1],
                     population = df$population[1], model = str_c("sample", seq_along(smpls)), 
                     env_rank = smpls)
          
          # # Combine
          # bind_rows(df, sample_df)
          
        }) %>%
        ungroup() )





# Save this
save_file <- file.path(result_dir, "distance_method_results.RData")
save("clust_method_df", "pred_env_dist_rank","pred_env_rank_random", file = save_file)


