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
  filter(line_name %in% c(tp_geno, vp_geno),
         environment %in% tp_vp_env) %>% # Only use environments in which both the TP and VP were grown
  mutate(line_name = as.factor(line_name))



# Define clusters for each distance method
# Cut the tree to obtain the maximum number of clusters with at least 3 environments in each cluster (arbitrary)
min_env <- c(3, 5, 10)
# Number of random cluster assignments
n_random <- 25

clusters <- env_rank_df %>% 
  crossing(., min_env) %>%
  group_by(set, trait, model, min_env) %>%
  do({
    df <- .
    
    # Convert the distance object to a matrix - subset the relevant environments
    dist_mat <- as.matrix(df$dist[[1]])
    dist_mat1 <- dist_mat %>% 
      subset(row.names(.) %in% tp_vp_env, colnames(.) %in% tp_vp_env) %>%
      as.dist()
    
    clus <- hclust(d = dist_mat1, method = "ward.D")
    minE <- df$min_env
    
    # Find K
    k <- map(seq_along(clus$labels), ~cutree(clus, k = .)) %>% 
      map(table) %>% 
      map_lgl(~all(. >= minE)) %>% 
      which() %>% 
      max()
    
    # Cut the tree
    cluster_df <- cutree(clus, k = k) %>%
      {data.frame(environment = names(.), cluster = ., stringsAsFactors = FALSE, row.names = NULL)}
    
    data_frame(nClust = k, env_cluster = list(cluster_df))
    
  }) %>% ungroup()




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
#   group_by(set, trait, model, min_env) %>%
#   do({
#     row <- .
    

    clus <- row$env_cluster[[1]]
    tr <- unique(row$trait)
    
    
    # Attach the BLUEs and split by cluster
    data_tomodel <- clus %>% 
      left_join(filter(S2_MET_BLUEs_tomodel, trait == tr), by = "environment") %>%
      mutate(row = seq(nrow(.)))
    
    data_tomodel_split <- data_tomodel %>%
      split(.$cluster)
    
    # For each cluster, drop one environment, train a model, and predict
    pred_out <- vector("list", length(data_tomodel_split))
    
    for (i in seq_along(pred_out)) {
      
      df_tomodel <- data_tomodel_split[[i]] %>%
        mutate(row = seq(nrow(.)))
      pred_envs <- unique(df_tomodel$environment)
      train_test <- pred_envs %>%
        map_df(~data_frame(train = list(resample(data = df_tomodel, idx = filter(df_tomodel, environment != ., line_name %in% tp_geno)$row)), 
                           test = list(resample(data = df_tomodel, idx = filter(df_tomodel, environment == ., line_name %in% vp_geno)$row))))
      
      ## Using the training set, calculate BLUEs over environments
      train_test <- train_test %>%
        mutate(train = map(train, ~as.data.frame(.) %>% geno_means(data = .) %>% mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))))
      
      
      # Make predictions using the training set
      pgvs <- train_test$train %>% 
        map(~{
          mf <- model.frame(value ~ line_name, .)
          # mf <- model.frame(value ~ line_name + environment, .)
          # mf <- model.frame(value ~ line_name + environment, train_test$train[[1]])
          
          y <- model.response(mf)
          # X <- model.matrix(~ 1 + environment, mf)
          X <- model.matrix(~ 1, mf)
          
          Z <- model.matrix(~ -1 + line_name, mf)
          
          fit <- mixed.solve(y = y, Z = Z, K = K, X = X)
          fit$u %>% {data.frame(line_name = names(.), pred_value = ., row.names = NULL, stringsAsFactors = FALSE)}
        })
      
      # Validate predictions using the testing set
      pred_acc <- map2_dbl(.x = train_test$test, .y = pgvs, ~left_join(as.data.frame(.x), .y, by = "line_name") %>% {cor(.$value, .$pred_value)})
      
      
      ### Use a random sample of the same number of environments for prediction
      n_train_env <- length(pred_envs) - 1
      
      pred_acc_random <- pred_envs %>%
        map(~{
          val_env <- .
          # Randomly sample remaining environments
          rand_train_env <- replicate(n_random, sample(x = setdiff(clus$environment, val_env), size = n_train_env), simplify = FALSE)
          
          # Create train/test sets
          train <- map(rand_train_env, ~resample(data = data_tomodel, idx = filter(data_tomodel, environment %in% ., line_name %in% tp_geno)$row))
          test <- replicate(length(train), resample(data = data_tomodel, idx = filter(data_tomodel, environment == val_env, line_name %in% vp_geno)$row), simplify = FALSE)
            
          ## Using the training set, calculate BLUEs over environments
          train_test <- map(train, ~as.data.frame(.) %>% geno_means(data = .) %>% 
                              mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno))))
          
          # Make predictions using the training set
          pgvs <- train_test %>% 
            map(~{
              mf <- model.frame(value ~ line_name, .)
              # mf <- model.frame(value ~ line_name + environment, .)
              # mf <- model.frame(value ~ line_name + environment, train_test$train[[1]])
              
              y <- model.response(mf)
              # X <- model.matrix(~ 1 + environment, mf)
              X <- model.matrix(~ 1, mf)
              
              Z <- model.matrix(~ -1 + line_name, mf)
              
              fit <- mixed.solve(y = y, Z = Z, K = K, X = X)
              fit$u %>% {data.frame(line_name = names(.), pred_value = ., row.names = NULL, stringsAsFactors = FALSE)}
            })
          
          # Validate predictions using the testing set
          pred_acc <- map2_dbl(.x = test, .y = pgvs, ~left_join(as.data.frame(.x), .y, by = "line_name") %>% {cor(.$value, .$pred_value)})
          
          # Return predictions
          data_frame(environment = val_env, sample = seq(n_random), accuracy = pred_acc)
          
        })
            
      
      # Add predictions to the list
      pred_out[[i]] <- data.frame(environment = pred_envs, accuracy = pred_acc, row.names = NULL, stringsAsFactors = FALSE) %>%
        as_data_frame() %>%
        mutate(random_accuracy = pred_acc_random)
      
    }
    
    # Bind the rows of the pred_out list and combine with the cluster df
    # then return
    # left_join(clus, bind_rows(pred_out), by = "environment")
    results_out[[r]] <- left_join(clus, bind_rows(pred_out), by = "environment")
    
  # }) %>% ungroup()
    
  } # Close the row loop

  core_df %>%
    mutate(out = results_out) %>%
    select(-core)

}, mc.cores = n_core) # Close the parallel operation

cluster_predictions <- bind_rows(cluster_predictions)



## Save this
save("cluster_predictions", file = file.path(result_dir, "cluster_predictions.RData"))













