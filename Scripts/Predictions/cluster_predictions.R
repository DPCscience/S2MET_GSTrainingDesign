## Predictions based on environmental clustering
## 
## Author: Jeff Neyhart
## Last Updated: October 15, 2018
## 



### Run for MSI

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

clusters <- clust_method_df %>% 
  filter(population == "tp") %>%
  crossing(., min_env) %>%
  group_by(trait, model, min_env) %>%
  do({
    clus <- .$cluster[[1]]
    minE <- .$min_env
    
    # Find K
    k <- map(seq_along(clus$labels), ~cutree(clus, k = .)) %>% 
      map(table) %>% 
      map_lgl(~all(. >= minE)) %>% 
      which() %>% 
      max()
    
    # Cut the tree
    cluster_df <- cutree(clus, k = k) %>%
      {data.frame(environment = names(.), cluster = ., stringsAsFactors = FALSE, row.names = NULL)} %>%
      # Remove unwanted environments
      filter(environment %in% tp_vp_env)
    
    data_frame(env_cluster = list(cluster_df), k = k)
    
  }) %>% ungroup()




## Split the clusters for parallelization
clusters_split <- clusters %>%
  assign_cores(df = ., n_core = n_core) %>%
  split(.$core)


## Iterate over trait and model combinations
cluster_predictions <- mclapply(X = clusters_split, FUN = function(core_df) {
  
  results_out <- vector("list", nrow(core_df))
  
  ## Iterate over rows
  for (r in seq_along(results_out)) {
    
    row <- core_df[r,]
    clus <- row$env_cluster[[1]]
    tr <- unique(row$trait)
    
    # Attach the BLUEs and split by cluster
    data_tomodel <- clus %>% 
      left_join(filter(S2_MET_BLUEs_tomodel, trait == tr), by = "environment") %>%
      split(.$cluster)
    
    # For each cluster, drop one environment, train a model, and predict
    pred_out <- vector("list", length(data_tomodel))
    
    for (i in seq_along(pred_out)) {
      df_tomodel <- data_tomodel[[i]] %>%
        mutate(row = seq(nrow(.)))
      pred_envs <- unique(df_tomodel$environment)
      train_test <- pred_envs %>%
        map_df(~data_frame(train = list(resample(data = df_tomodel, idx = filter(df_tomodel, environment != ., line_name %in% tp_geno)$row)), 
                           test = list(resample(data = df_tomodel, idx = filter(df_tomodel, environment == ., line_name %in% vp_geno)$row))))
      
      # Make predictions using the training set
      pgvs <- train_test$train %>% 
        map(~{
          mf <- model.frame(value ~ line_name + environment, .)
          # mf <- model.frame(value ~ line_name + environment, train_test$train[[1]])
          
          y <- model.response(mf)
          X <- model.matrix(value ~ 1 + environment, droplevels(mf))
          Z <- model.matrix(value ~ -1 + line_name, mf)
          
          fit <- mixed.solve(y = y, Z = Z, K = K, X = X)
          fit$u %>% {data.frame(line_name = names(.), pred_value = ., row.names = NULL, stringsAsFactors = FALSE)}
        })
      
      # Validate predictions using the testing set
      pred_acc <- map2_dbl(.x = train_test$test, .y = pgvs, ~left_join(as.data.frame(.x), .y, by = "line_name") %>% {cor(.$value, .$pred_value)})
      
      # Add predictions to the list
      pred_out[[i]] <- data.frame(environment = pred_envs, accuracy = pred_acc, row.names = NULL, stringsAsFactors = FALSE)
      
    }
    
    # Bind the rows of the pred_out list and combine with the cluster df
    # then return
    results_out[[r]] <- left_join(clus, bind_rows(pred_out), by = "environment")
    
  } # Close the row loop
  
  core_df %>%
    mutate(out = results_out) %>%
    select(-core)
  
}, mc.cores = n_core) # Close the parallel operation

cluster_predictions <- bind_rows(cluster_predictions)



## Save this
save("cluster_predictions", file = file.path(result_dir, "cluster_predictions.RData"))










