## 
## Author: Jeff Neyhart
## Last Updated: October 15, 2018
## 



### Run on MSI

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"
source(file.path(repo_dir, "source_MSI.R"))



# ## Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# library(modelr)



## Number of cores
n_core <- 16
n_core <- detectCores()



# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))

## Prepare the BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate(line_name = as.factor(line_name))

# Number of random cluster assignments
n_random <- 25

## Subset the modeling data to exclude the VP
S2_MET_BLUEs_tomodel_train <- S2_MET_BLUEs_tomodel %>%
  filter(line_name %in% tp_geno)



## Data.frame of test environments
test_env_df <- bind_rows(
  data_frame(set = "complete", trait = names(complete_train_env), train_env = complete_train_env, test_env = complete_train_env),
  data_frame(set = "realistic", trait = names(complete_train_env), train_env = realistic_train_env, test_env = realistic_test_env)
)

clusters <- cluster_df %>% 
  left_join(., test_env_df) %>%
  mutate(nClusters = map_dbl(cluster, ~n_distinct(.$cluster)))

cluster_train_test <- clusters %>%
  mutate(data = list(NULL))

## Assign training and testing data to the clusters
for (i in seq(nrow(cluster_train_test))) {
  row <- cluster_train_test[i,]
  
  clus <- row$cluster[[1]]
  tr <- unique(row$trait)
  test_env <-  row$test_env[[1]]
  train_env <- row$train_env[[1]]
  
  # Attach the BLUEs
  data_tomodel <- clus %>% 
    left_join(filter(S2_MET_BLUEs_tomodel, trait == tr), by = "environment")
  
  # Split by cluster
  clus_split <- clus %>% 
    split(.$cluster) %>%
    map_df(~{
      # Create a list of testing environments
      test_env_tomodel <- intersect(test_env, .$environment)
      train_env_tomodel <- intersect(train_env, .$environment)
      
      df_tomodel <- filter(S2_MET_BLUEs_tomodel, trait == tr, environment %in% .$environment)
      
      data_frame(val_environment = test_env_tomodel) %>% 
        mutate(train = map(val_environment, ~filter(df_tomodel, line_name %in% tp_geno, environment %in% train_env_tomodel, environment != .) %>%
                             mutate(environment = as.factor(environment))), 
               test = map(val_environment, ~filter(df_tomodel, line_name %in% vp_geno, environment == .) %>% 
                            mutate(line_name = as.character(line_name)))) %>%
        mutate_at(vars(train, test), ~map(., ~select(., environment, line_name, value, std_error)))
      
    }) %>% mutate(nTrainEnv = map_dbl(train, ~n_distinct(.$environment)))
  
  cluster_train_test$data[[i]] <- left_join(x = clus_split, y = clus, by = c("val_environment" = "environment"))
  
}
  
  


## Split the clusters for parallelization
clusters_split <- cluster_train_test %>%
  unnest(data) %>%
  assign_cores(df = ., n_core = n_core) %>%
  split(.$core)


## Iterate over trait, model, and validation environments
cluster_predictions <- mclapply(X = clusters_split, FUN = function(core_df) {

  # #
  # r = 1
  # core_df <- clusters_split[[r]]
  # #

  results_out <- vector("list", nrow(core_df))

  ## Iterate over rows
  for (i in seq_along(results_out)) {

    # Grab the data
    row <- core_df[i,]

    ##
    ##

   
# cluster_predictions <- cluster_train_test %>%
#   unnest(data) %>%
#   group_by(set, model, trait, val_environment) %>%
#   do(out = {
#     
#     row <- .
#     
#     ##
    
    ## Run the base predictions
    base_pred <- gblup(K = K, train = row$train[[1]], test = row$test[[1]], fit.env = TRUE)$accuracy

    
    # # Get the pool of possible training environments - excluding the current validation environment
    # possible_train_env <- setdiff(filter(cluster_train_test, set == row$set, model == row$model, trait == row$trait)$train_env[[1]], row$val_environment)
    # # Get the data for these possible trainin environments
    # possible_train_data <- subset(S2_MET_BLUEs_tomodel_train, trait == row$trait & environment %in% possible_train_env)
    # 
    # ## If the number of possible training environments is equal to the number of intended training environments,
    # ## skip randomization, since the results will be the same
    # if (length(possible_train_env) == row$nTrainEnv) {
    #   random_pred <- NA
    #   
    # } else {
    # 
    #   ## Randomly sample the same number of environments from this pool
    #   random_train_data <- replicate(n = n_random, sample(possible_train_env, size = row$nTrainEnv), simplify = FALSE) %>% 
    #     map(~filter(possible_train_data, environment %in% .))
    #   # Run predictions
    #   random_pred <- random_train_data %>%
    #     map_dbl(~gblup(K = K, train = ., test = row$test[[1]], fit.env = TRUE)$accuracy)
    #   
    # }

  #   ##
  #       
  #   list(base = base_pred, random = random_pred)
  #   
  # }) %>% ungroup()
  # 
  # ##
    
    #
    random_pred <- NA
  
    
    # Add predictions to the list
    results_out[[i]] <- list(base = base_pred, random = random_pred)

  }

  core_df %>%
    mutate(out = results_out) %>%
    select(-core, -train, -test)

}, mc.cores = n_core) # Close the parallel operation

cluster_predictions <- bind_rows(cluster_predictions)









## Cross-validation predictions
# First generate test/train sets
n_train <- length(vp_geno)
reps <- 50

cv_test_sets <- replicate(n = reps, sort(sample(tp_geno, size = n_train)), simplify = FALSE)
cv_train_sets <- map(cv_test_sets, ~setdiff(tp_geno, .))

## Create a data.frame
cv_sets <- data_frame(cv_rep = seq_along(cv_test_sets), cv_train = cv_train_sets, cv_test = cv_test_sets)


cluster_cv_train_test <- cluster_train_test %>%
  unnest(data) %>%
  left_join(., nest(group_by(S2_MET_BLUEs_tomodel_train, environment, trait)), by = c("val_environment" = "environment", "trait")) %>%
  crossing(., cv_sets) %>%
  mutate(cv_train = map2(.x = train, .y = cv_train, ~filter(.x, line_name %in% .y)),
         cv_test = map2(.x = data, .y = cv_test, ~filter(.x, line_name %in% .y))) %>%
  select(set:val_environment, nTrainEnv, cluster, cv_rep, cv_train, cv_test, vp_test = test)


## Split
clusters_cv_split <- cluster_cv_train_test %>%
  assign_cores(n_core = n_core) %>%
  split(.$core)


## Iterate over trait, model, and validation environments
cluster_cv_predictions <- mclapply(X = clusters_cv_split, FUN = function(core_df) {

  # #
  # r = 1
  # core_df <- clusters_cv_split[[r]]
  # core_df <- clusters_cv_split[[r]] %>% sample_n(size = 200)
  # #
  
  results_out <- vector("list", nrow(core_df))
  
  ## Iterate over rows
  for (i in seq_along(results_out)) {
    
    # Grab the data
    row <- core_df[i,]
    
    ##
    ##


# cluster_cv_predictions <- cluster_cv_train_test %>%
#   unnest(data) %>%
#   filter(model %in% c("great_circle_dist", "pheno_location_dist"), trait == "GrainYield") %>%
#   group_by(set, model, trait, val_environment, cv_rep) %>%
#   do({
# 
#     row <- .

    ##
    
    ## Run the base predictions
    base_pred_cv <- gblup(K = K, train = row$cv_train[[1]], test = row$cv_test[[1]], fit.env = TRUE)$accuracy
    base_pred_vp <- gblup(K = K, train = row$cv_train[[1]], test = row$vp_test[[1]], fit.env = TRUE)$accuracy
    
    ##

#     data.frame(cluster = row$cluster, nTrainEnv = row$nTrainEnv, base_pred_cv, base_vp = base_pred_vp, random = random_pred, 
#                stringsAsFactors = FALSE)
# 
# }) %>% ungroup()

    ##

    # Add predictions to the list
    results_out[[i]] <- list(base_pred_cv, base_vp = base_pred_vp)

    }

  core_df %>%
    mutate(out = results_out) %>%
    select(-core, -cv_train, -cv_test, -vp_test)

}, mc.cores = n_core) # Close the parallel operation



cluster_cv_predictions <- bind_rows(cluster_cv_predictions)





## Save this
save("cluster_predictions", "cluster_cv_predictions", file = file.path(result_dir, "cluster_predictions.RData"))













