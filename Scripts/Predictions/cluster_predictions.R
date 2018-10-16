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


## Number of cores
n_core <- 16
n_core <- detectCores()



# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))

## Prepare the BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate(line_name = as.factor(line_name))

# Define clusters for each distance method
# Cut the tree to obtain the maximum number of clusters with at least 3 environments in each cluster (arbitrary)
min_env <- 3

clusters <- clust_method_df %>% 
  filter(population == "tp") %>%
  group_by(trait, model) %>%
  do({
    clus <- .$cluster[[1]]
    
    # Find K
    k <- map(seq_along(clus$labels), ~cutree(clus, k = .)) %>% 
      map(table) %>% 
      map_lgl(~all(. >= min_env)) %>% 
      which() %>% 
      max()
    
    # Cut the tree
    cluster_df <- cutree(clus, k = k) %>%
      {data.frame(environment = names(.), cluster = ., stringsAsFactors = FALSE, row.names = NULL)}
    
    data_frame(env_cluster = list(cluster_df), k = k)
    
  }) %>% ungroup()


## Merge with the BLUEs and assign clusters


## Iterate over trait and model combinations
cluster_predictions <- clusters %>%
  group_by(trait, model) %>%
  do({
    clus <- .$env_cluster[[1]]
    tr <- unique(.$trait)
    
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
    left_join(clus, bind_rows(pred_out), by = "environment")
    
  })

## Save this
save("cluster_predictions", file = file.path(result_dir, "cluster_predictions.RData"))


# Get the unique distance models, then assign new names
models <- sort(unique(cluster_predictions$model))
models_rename <- setNames(c("GCD", "MYEC_All", "MYEC_Cor1", "MYEC_F1", "MYEC_Cor5", "MYEC_F5", "OYEC_All", "OYEC_Cor1", "OYEC_F1", "OYEC_Cor5", "OYEC_F5", "PD"), models) %>%
  # Order
  .[c(length(.), seq(length(.) - 1))] %>%
  c(., "random" = "Random")

## Calculate the average accuracy
cluster_predictions_summ <- cluster_predictions %>% 
  ungroup() %>% 
  mutate(model = str_replace_all(model, models_rename)) %>%
  group_by(trait, model) %>%
  # group_by(trait, model, cluster) %>%
  # summarize(accuracy = mean(accuracy, na.rm = T)) %>% 
  summarize(accuracy = mean(accuracy, na.rm = T)) 


# Plot
cluster_predictions_summ %>% 
  ggplot(aes(x = trait, y = accuracy, fill = model)) +
  geom_col(position = "dodge") +
  facet_wrap(~trait, ncol = 2, scales = "free_x")


# Model
fits <- cluster_predictions %>% 
  ungroup() %>% 
  mutate(model = str_replace_all(model, models_rename)) %>%
  mutate_at(vars(model, environment, cluster), as.factor) %>%
  split(.$trait) %>%
  map(~lm(accuracy ~ environment + model, data = .))

effects <- fits %>%
  map(~effects::allEffects(.))









