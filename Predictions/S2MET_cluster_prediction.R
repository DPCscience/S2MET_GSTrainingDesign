## S2MET Cluster prediction
## 
## This script will compare clustering methods in the effectiveness of genomewide
## prediction
## 

# List of packages
packages <- c("dplyr", "tidyr", "tibble", "purrr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "pbr", "purrrlyr", "rrBLUP", "gws")

# Set the directory of the R packages
package_dir <- NULL
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET/"
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/" 

# Geno, pheno, and enviro data
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/"
pheno_dir <- file.path(proj_dir, "Phenotype_Data/")
env_var_dir <- file.path(proj_dir, "Environmental_Variables")

geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/GBS_Genos"
pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"
env_var_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Environmental_Data"


# Other directories
fig_dir <- file.path(proj_dir, "Figures/")
pred_dir <- file.path(proj_dir, "Predictions")
entry_dir <- file.path(proj_dir, "Plant_Materials")
analysis_dir <- file.path(proj_dir, "Analysis")
result_dir <- file.path(proj_dir, "Results")




# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))
# Load environmental data
load(file.path(env_var_dir, "environmental_data_compiled.RData"))

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))


# Grab the entry names that are not checks
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)

vp <- entry_list %>% 
  filter(Class == "S2C1R") %>% 
  pull(Line)

# Find the tp and vp that are genotypes
tp_geno <- intersect(tp, row.names(s2_imputed_mat))
vp_geno <- intersect(vp, row.names(s2_imputed_mat))

# Define the checks
checks <- entry_list %>% 
  filter(Class == "Check") %>% 
  pull(Line)

entries <- entry_list %>% 
  pull(Line)

# Extract the tp and vp from the G matrix
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno, vp_geno),]

# Filter the BLUEs for the usable genotypes (i.e. those with marker data)
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  # Filter the S2 MET BLUEs for non-irrigated trials
  filter(!grepl(pattern = "BZI|HTM", x = environment)) %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate(line_name = as.factor(line_name)) %>%
  select(environment, line_name:std_error) %>%
  droplevels()

# Detect cores
n_cores <- detectCores()


# Calculate the genomic relationship matrix
A <- A.mat(X = s2_imputed_mat_use, min.MAF = 0, max.missing = 1)



## Load the cluster data.frame
load_file <- file.path(result_dir, "S2_MET_cluster_df.RData")
load(load_file)

# Combine the relevant distance matrices into a single DF
dist_df <- clust_method_df %>%
  select(trait, D, PCA, FA, EC_one_dist, EC_multi_dist) %>%
  ungroup()



## Analysis 1 - Adding more distant environments
## 
## For each distance metric and for each prediction environment, successively add
## environments that are more distant and assess the prediction accuracy in that
## environment. We will use an across-environment model.
## 

# Define a function to convert the distance matrix into a data.frame with
# a column for the environment
dist_to_df <- function(d) d %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("environment")


# Convert each distance object into a data.frame, then gather and unnest()
dist_env_df <- dist_df %>%
  mutate_at(vars(-trait), funs(map), dist_to_df) %>%
  gather(method, data, -trait) %>%
  unnest() %>%
  # Gather again
  gather(environment2, distance, -trait:-environment) %>%
  # Remove env1 == env2
  filter(environment != environment2,
         !is.na(distance)) %>%
  arrange(trait, method, environment)


# Find the prediction environments
# That is, the environments in which the vp is present
pred_env <- S2_MET_BLUEs_use %>%
  group_by(trait, environment) %>%
  filter(sum(line_name %in% vp_geno) > 1) %>%
  distinct(environment, trait)

# Find the training environments
# That is, the environments in which the tp is present
train_env <- S2_MET_BLUEs_use %>%
  group_by(trait, environment) %>%
  filter(sum(line_name %in% tp_geno) > 1) %>%
  distinct(environment, trait)

# # Combine the possible prediction environments with the possible training environments
# # for a single trait
# dist_env_df1 <- left_join(pred_env, train_env, by = "trait") %>%
#   select(trait, pred_env = environment.x, train_env = environment.y) %>%
#   filter(pred_env != train_env) %>%
#   ungroup() %>%
#   mutate_at(vars(contains("env")), as.character) %>%
#   left_join(., dist_env_df, by = c("trait", "pred_env" = "environment", "train_env" = "environment2")) %>%
#   group_by(pred_env) %>%
#   arrange(trait, method, pred_env, distance) %>%
#   group_by(trait, method, pred_env) %>%
#   nest(.key = "env_dist")
# 
# # Group by rows and run the predictions
# add_one_env_pred <- dist_env_df1 %>%
#   mutate(core = sort(rep(seq(n_cores), length.out = nrow(.)))) %>%
#   split(.$core) %>%
#   # Apply by core
#   mclapply(function(core_df) {
#     core_df %>%
#       group_by(trait, method, pred_env) %>%
#       do(pred_results = {
# 
#         # Convert the row to something useable
#         i <- unnest(.)
# 
#         # Add the cumulative mean distance of adding environments
#         i_cummean <- i %>%
#           mutate(cum_mean_dist = cummean(distance))
# 
#         # Add the prediction data
#         i_cummean_pred <- i_cummean %>%
#           left_join(x = ., y = filter(S2_MET_BLUEs_use, line_name %in% vp_geno),
#                     by = c("trait", "pred_env" = "environment")) %>%
#           group_by_at(vars(trait:cum_mean_dist)) %>%
#           nest(line_name:std_error, .key = "pred_data") %>%
#           ungroup()
# 
#         # Add the training data
#         i_cummean_train <- i_cummean %>%
#           left_join(x = ., y = filter(S2_MET_BLUEs_use, line_name %in% tp_geno),
#                     by = c("trait", "train_env" = "environment")) %>%
#           mutate(environment = train_env) %>%
#           group_by_at(vars(trait:cum_mean_dist)) %>%
#           nest(environment, line_name:std_error, .key = "train_data") %>%
#           ungroup() %>%
#           # Accumulate the training data
#           mutate(train_data = accumulate(train_data, bind_rows))
# 
#         # Combine
#         i_tomodel <- full_join(i_cummean_pred, i_cummean_train,
#                                by = c("trait", "method", "pred_env", "train_env", "distance", "cum_mean_dist"))
# 
# 
#         # Vector of prediction accuracies
#         pred_acc_vec <- vector("list", length = nrow(i_tomodel))
# 
#         # Iterate over the list of environments
#         for (r in seq_along(pred_acc_vec)) {
# 
#           # Model frame
#           mf <- model.frame(value ~ line_name + environment + std_error, data = i_tomodel$train_data[[r]])
#           y <- model.response(mf)
#           # Matrix of 1s if r = 1
#           if (n_distinct(mf$environment) == 1) {
#             X <-  model.matrix(~ 1, droplevels(mf))
#           } else {
#             X <- model.matrix(~ 1 + environment, droplevels(mf))
#           }
#           Zlist <- ranef_model_matrix(random = ~ g(line_name), data = mf, vcov = list(line_name = A),
#                                       sparse = TRUE)
# 
#           # R matrix of standard errors
#           R <- solve(Diagonal(x = mf$std_error^2))
# 
#           fit1 <- sommer::mmer(Y = y, X = X, Z = Zlist, R = list(units = R))
# 
#           # Extract PGV
#           pgv <- fit1$u.hat[[1]] %>%
#             as.data.frame() %>%
#             rownames_to_column("line_name") %>%
#             rename(pgv = T1)
# 
#           # Extract the validation data
#           pred_data <- pull(i_tomodel, pred_data)[[r]] %>%
#             mutate(line_name = as.character(line_name))
# 
#           ## Measure accuracy
#           pred_acc_vec[[r]]  <- left_join(pred_data, pgv, by = "line_name") %>%
#             do(boot_cor(x = .$pgv, y = .$value, boot.reps = 1000))
# 
#         }
# 
#         # Return the cummean df with prediction accuracies
#         bind_cols(i_cummean, bind_rows(pred_acc_vec))
# 
#         })
#    }, mc.cores = n_cores)
# 
# # Save
# save_file <- file.path(result_dir, "S2MET_pred_by_env_dist.RData")
# save("add_one_env_pred", file = save_file)
# 



## Random control
## For each environment to predict, randomly sample n environments to use in prediction
##
## First set up the predictions
random_env_df <- left_join(pred_env, train_env, by = "trait") %>%
  select(trait, pred_env = environment.x, train_env = environment.y) %>%
  filter(pred_env != train_env) %>%
  ungroup() %>%
  mutate_at(vars(contains("env")), as.character) %>%
  group_by(trait, pred_env) %>%
  nest(train_env, .key = "train_env")

# For each trait, find the number of possible training environments (n_train)
# Randomly sample n = 1, 2, ..., n_train - 1 environments with replacement 100 times
# Use those environments to predict the validation environment
random_env_df1 <- random_env_df %>%
  group_by(trait) %>%
  mutate(n_train = list(seq(min(map_dbl(train_env, nrow))))) %>%
  unnest(n_train) %>%
  left_join(., random_env_df) %>%
  ungroup()

### Maximum number of training environments
max_train <- max(random_env_df1$n_train)
max_train <- 8

n_iter <- 25

# Add cores and run the predictions
random_env_pred_out <- random_env_df1 %>%
  # Filter n_train to be less than max_train
  filter(n_train <= max_train) %>%
  mutate(core = sort(rep(seq(n_cores), length.out = nrow(.)))) %>%
  split(.$core) %>%
  # Apply by core
  mclapply(function(core_df) {
    core_df %>%
      by_row(function(i) {

        # Number of training environments to sample
        n_train_sample <- i$n_train

        # Grab the validation data
        pred_data <- left_join(i, S2_MET_BLUEs_use, by = c("trait", "pred_env" = "environment")) %>%
          filter(line_name %in% vp_geno) %>%
          select(line_name:std_error) %>%
          mutate(line_name = as.character(line_name))

        # Draw iterations of sampling
        train_data <- rerun(.n = n_iter, {

          unnest(i) %>%
            sample_n(size = n_train_sample, replace = TRUE) %>%
            left_join(., S2_MET_BLUEs_use, by = c("trait", "train_env" = "environment")) %>%
            filter(line_name %in% tp_geno)

        })

        # Run predictions across the iterations
        predictions <- train_data %>%
          map_df(function(td) {

            # Model frame
            mf <- model.frame(value ~ line_name + train_env + std_error, data = td)
            y <- model.response(mf)
            # Matrix of 1s if r = 1
            if (n_distinct(mf$train_env) == 1) {
              X <-  model.matrix(~ 1, droplevels(mf))
            } else {
              X <- model.matrix(~ 1 + train_env, droplevels(mf))
            }

            Zlist <- ranef_model_matrix(random = ~ g(line_name), data = mf, vcov = list(line_name = A),
                                        sparse = TRUE)

            # R matrix of standard errors
            R <- solve(Diagonal(x = mf$std_error^2))

            fit1 <- sommer::mmer(Y = y, X = X, Z = Zlist, silent = TRUE)

            # Extract PGV
            pgv <- fit1$u.hat[[1]] %>%
              as.data.frame() %>%
              rownames_to_column("line_name") %>%
              rename(pgv = T1)

            left_join(pred_data, pgv, by = "line_name") %>%
              do(boot_cor(x = .$pgv, y = .$value, boot.reps = 1000))

          })

        # Return the prediction accuracies
        return(predictions)

      }, .to = "pred_out")

  }, mc.cores = n_cores)

# Save
save_file <- file.path(result_dir, "S2MET_pred_by_random_env.RData")
save("random_env_pred_out", file = save_file)





# ## Prediction
# # For all clusters, drop one environment and use the remaining environments to predict
# # it.
# # Iterate over clusters
# # First do this using clusters created by all entries
# cluster_pred_acc_setup <- clust_df_tomodel %>%
#   filter(between(k, 3, 6)) %>%
#   rename(tr = trait) %>%
#   # Add core
#   mutate(core = sort(rep(seq(1, n_cores), length.out = nrow(.)))) %>%
#   rowwise() %>%
#   
#   # Cut the tree and assign environments to clusters
#   mutate(clust_df = list({
#     cutree(clust, k = k) %>%
#       data_frame(cluster_env = names(.), cluster = .) %>%
#       mutate(cluster = paste("cluster", cluster, sep = "")) })) %>%
#   
#   # Add the clusters to the BLUEs data
#   mutate(phenos_cluster = list({
#     left_join(S2_MET_BLUE_entries, clust_df, c("environment" = "cluster_env")) %>%
#       mutate(cluster = as.factor(cluster)) %>%
#       filter(trait == as.character(tr)) %>%
#       droplevels() })) %>%
#   
#   # Create a df of the the training and testing data for each cluster/environmnent
#   mutate(pred_df = list({
#     clust_df %>% 
#       group_by(cluster, cluster_env) %>% 
#       do(data = filter(phenos_cluster, cluster == .$cluster)) %>%
#       mutate(
#         train_df = list(filter(data, environment != cluster_env, line_name %in% tp)), 
#         test_df = list(filter(data, environment == cluster_env, line_name %in% vp))) %>%
#       select(-data) })) %>%
#   
#   # Remove unnecessary columns
#   select(-clust_df, -phenos_cluster)
#   
# # Run predictions
# cluster_pred_acc <- cluster_pred_acc_setup %>%
#   unnest(pred_df) %>% 
#   group_by(population, tr, method, k, cluster, cluster_env) %>%
#   do(acc = pred_acc_out(train_df = .$train_df[[1]], test_df = .$test_df[[1]])) %>%
#   select(population, trait = tr, method, k, cluster, pred_env = cluster_env, acc)
# 
# 
# # Create a similar data.frame, but for predicting each environment (LOO) using all
# # other environments (non-clustered)
# nocluster_pred_acc_all <- distinct(S2_MET_BLUE_entries, trait, environment) %>%
#   rename(tr = trait, cluster_env = environment) %>% 
#   rowwise() %>%
#   mutate(cluster = "cluster0", method = "no_cluster") %>%
#   
#   # Add the BLUEs data
#   mutate(phenos_cluster = list({
#     S2_MET_BLUE_entries %>%
#       filter(trait == as.character(tr)) %>%
#       droplevels() })) %>%
#   group_by(tr, method, cluster_env, cluster) %>% 
#   do(train_df = {filter(.$phenos_cluster[[1]], environment != .$cluster_env[1], line_name %in% tp)}, 
#      test_df = {filter(.$phenos_cluster[[1]], environment == .$cluster_env[1], line_name %in% vp)})
# 
# nocluster_pred_acc_all_results <- nocluster_pred_acc_all %>%
#   group_by(tr, method, cluster, cluster_env) %>%
#   do(acc = pred_acc_out(train_df = .$train_df[[1]], test_df = .$test_df[[1]])) %>%
#   select(trait = tr, method, cluster, pred_env = cluster_env, acc)
# 
# 
# 
# # Combine
# pred_acc_results <- bind_rows(
#   cluster_pred_acc,
#   nocluster_pred_acc_all_results)
# 
# # Save the data
# save_file <- file.path(pred_dir, "Results/cluster_predictions.RData")
# save("pred_acc_results", file = save_file)
#  


