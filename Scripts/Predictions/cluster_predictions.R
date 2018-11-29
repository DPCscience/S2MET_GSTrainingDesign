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
  filter(line_name %in% c(tp_geno, vp_geno)) %>% # Only use environments in which both the TP and VP were grown
  mutate(line_name = as.factor(line_name))

# Number of random cluster assignments
n_random <- 25

## Data.frame of test environments
test_env_df <- bind_rows(
  data_frame(set = "complete", trait = names(complete_train_env), test_env = complete_train_env),
  data_frame(set = "realistic", trait = names(complete_train_env), test_env = realistic_test_env)
)

clusters <- cluster_df %>% 
  left_join(., test_env_df) %>%
  mutate(nClusters = map_dbl(cluster, ~n_distinct(.$cluster)))


## Split the clusters for parallelization
clusters_split <- clusters %>%
  assign_cores(df = ., n_core = n_core) %>%
  split(.$core)


## Iterate over trait and model combinations
cluster_predictions <- mclapply(X = clusters_split, FUN = function(core_df) {

  # #
  # r = 1
  # core_df <- clusters_split[[r]]
  # #

  results_out <- vector("list", nrow(core_df))

  ## Iterate over rows
  for (r in seq_along(results_out)) {

    row <- core_df[r,]


# ## Run on a local machine
# 
# cluster_predictions <- clusters %>%
#   group_by(set, trait, model) %>%
#   do({
#     row <- .
    
    clus <- row$cluster[[1]]
    tr <- unique(row$trait)
    test_env <-  row$test_env[[1]]
    
    print(paste(row$trait, row$model))
    
    
    # Attach the BLUEs and split by cluster
    data_tomodel <- clus %>% 
      left_join(filter(S2_MET_BLUEs_tomodel, trait == tr), by = "environment") %>%
      mutate(row = seq(nrow(.)))
    
    data_tomodel_split <- data_tomodel %>%
      split(.$cluster)
    
    # For each cluster, drop one environment, train a model, and predict
    pred_out <- vector("list", length(data_tomodel_split))
    
    for (i in seq_along(pred_out)) {
      
      # Get data for the ith cluster
      df_tomodel <- data_tomodel_split[[i]]
      # Iterate over testing environments and create train/test sets
      test_env_tomodel <- intersect(test_env, df_tomodel$environment)

      train_test <- data_frame(val_environment = test_env_tomodel) %>% 
        mutate(train = map(val_environment, ~filter(df_tomodel, line_name %in% tp_geno, environment != .)), 
               test = map(val_environment, ~filter(df_tomodel, line_name %in% vp_geno, environment == .) %>% mutate(line_name = as.character(line_name))))
      
      
      # Make predictions using the training set
      train_test1 <- train_test %>%
        mutate(pgvs = map(train, ~{
          # mf <- model.frame(value ~ line_name, .)
          mf <- model.frame(value ~ line_name + environment, .)
          y <- model.response(mf)
          
          if (n_distinct(mf$environment) > 1) {
            X <- model.matrix(~ 1 + environment, mf)
          
          } else {
            X <- model.matrix(~ 1, mf)
          }
          
          Z <- model.matrix(~ -1 + line_name, mf)
          
          fit <- mixed.solve(y = y, Z = Z, K = K, X = X)
          fit$u %>% {data.frame(line_name = names(.), pred_value = ., row.names = NULL, stringsAsFactors = FALSE)}
        }))
      
      train_test2 <- train_test1 %>%
        mutate(pgvs = map2(.x = test, .y = pgvs, ~left_join(.x, .y, by = "line_name"))) %>% 
        unnest(pgvs) %>%
        select(trait, val_environment, line_name, value, pred_value)

      pred_acc <- train_test2 %>% 
        group_by(val_environment) %>% 
        summarize(accuracy = cor(value, pred_value)) %>%
        ungroup()
      
      # ### Use a random sample of the same number of environments for prediction
      # n_train_env <- length(pred_envs) - 1
      # 
      # pred_acc_random <- pred_envs %>%
      #   map(~{
      #     val_env <- .
      #     # Randomly sample remaining environments
      #     rand_train_env <- replicate(n_random, sample(x = setdiff(clus$environment, val_env), size = n_train_env), simplify = FALSE)
      #     
      #     # Create train/test sets
      #     train <- map(rand_train_env, ~resample(data = data_tomodel, idx = filter(data_tomodel, environment %in% ., line_name %in% tp_geno)$row))
      #     test <- replicate(length(train), resample(data = data_tomodel, idx = filter(data_tomodel, environment == val_env, line_name %in% vp_geno)$row), simplify = FALSE)
      #       
      #     ## Using the training set, calculate BLUEs over environments
      #     train_test <- map(train, ~as.data.frame(.) %>% geno_means(data = .) %>% 
      #                         mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))))
      #     
      #     # Make predictions using the training set
      #     pgvs <- train_test %>% 
      #       map(~{
      #         mf <- model.frame(value ~ line_name, .)
      #         # mf <- model.frame(value ~ line_name + environment, .)
      #         # mf <- model.frame(value ~ line_name + environment, train_test$train[[1]])
      #         
      #         y <- model.response(mf)
      #         # X <- model.matrix(~ 1 + environment, mf)
      #         X <- model.matrix(~ 1, mf)
      #         
      #         Z <- model.matrix(~ -1 + line_name, mf)
      #         
      #         fit <- mixed.solve(y = y, Z = Z, K = K, X = X)
      #         fit$u %>% {data.frame(line_name = names(.), pred_value = ., row.names = NULL, stringsAsFactors = FALSE)}
      #       })
      #     
      #     # Validate predictions using the testing set
      #     pred_acc <- map2_dbl(.x = test, .y = pgvs, ~left_join(as.data.frame(.x), .y, by = "line_name") %>% {cor(.$value, .$pred_value)})
      #     
      #     # Return predictions
      #     data_frame(environment = val_env, sample = seq(n_random), accuracy = pred_acc)
      #     
      #   })
            
      
      
      # Add predictions to the list
      pred_out[[i]] <- pred_acc # %>% mutate(random_accuracy = pred_acc_random)
      
    }
    
    # Bind the rows of the pred_out list and combine with the cluster df
    # then return
    # left_join(clus, bind_rows(pred_out), by = c("environment" = "val_environment"))
    results_out[[r]] <- left_join(clus, bind_rows(pred_out), by = c("environment" = "val_environment"))
    
  # }) %>% ungroup()
    
  } # Close the row loop

  core_df %>%
    mutate(out = results_out) %>%
    select(-core)

}, mc.cores = n_core) # Close the parallel operation

cluster_predictions <- bind_rows(cluster_predictions)



## Save this
save("cluster_predictions", file = file.path(result_dir, "cluster_predictions.RData"))













