## Random and cluster predictions in the S2MET project
## 
## Author: Jeff Neyhart
## Last Updated: March 12, 2019
## 



## Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))
library(modelr)



# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))


## Prepare the BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate(line_name = as.factor(line_name))

# Number of random cluster assignments
n_random <- 10
# Number of random environments when using all data
n_rand_env <- c(2, 5, 8)
# Number of CV iterations
nCV <- 25
# CV folds
k <- 5




## Data.frame of test environments
test_env_df <- bind_rows(
  data_frame(set = "complete", trait = names(complete_train_env), train_env = complete_train_env, test_env = complete_train_env),
  data_frame(set = "realistic", trait = names(complete_train_env), train_env = realistic_train_env, test_env = realistic_test_env)
)


## Data.frame of the training lines for CV
cv_tp_df <- data.frame(line_name = tp_geno, stringsAsFactors = FALSE)







###### Leave-one-environment-out predictions


## CV00 - predict new genotypes in new environments
### Leave-one-environment-out cross-validation (CV00)
# Generate training and test sets
cv00_train_test <- test_env_df %>%
  group_by(set, trait) %>%
  do({
    row <- .
    
    ## Subset the phenotype data for the particular trait
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
    
    row2 <- row %>% 
      unnest(test_env) %>%
      crossing(., select(row, train_env)) %>%
      mutate(val_environment = test_env,
             test_env = map(test_env, ~.),
             train_env = map2(train_env, test_env, setdiff))
      
    # Create testing and training sets based on that environment and the cv samples
    row2 %>%
      mutate(data = map2(.x = train_env, test_env, ~{
        tr <- .x
        te <- .y
        replicate(n = nCV, crossv_kfold(data = cv_tp_df, k = k), simplify = FALSE) %>%
          map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y)) %>%
          mutate_at(vars(train, test), funs(map(., as.data.frame))) %>%
          mutate(train = map(train, ~filter(pheno_use, environment %in% tr, line_name %in% .[[1]])),
                 test = map(test, ~filter(pheno_use, environment %in% te, line_name %in% c(.[[1]], vp_geno))),
                 nTrainEnv = length(.x)) })) %>%
      unnest(data)
    
    
    
  }) %>% ungroup()





## Predict
cv00_predictions <- cv00_train_test %>%
  group_by(set, trait, val_environment, nTrainEnv, rep) %>%
  do({
    
    df <- .
    
    ## Calculate genotype means for the training sets
    data_adj <- map(df$train, ~geno_means(data = .))
    
    ## Make predictions
    test_val_pred <- map2(.x = data_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]])
    
    ## Add results to the row
    df1 <- df %>%
      mutate(test_out = map(test_val_pred, ~filter(., !line_name %in% vp)),
             val_out = map(test_val_pred, ~filter(., line_name %in% vp))) %>%
      select(.id, test_out, val_out)
    
    ## Calculate test accuracy
    test_acc <- df1 %>%
      unnest(test_out) %>% 
      summarize(cv00 = cor(value, pred_value))
    
    ## Calculate average val accuracy
    val_acc <- df1 %>%
      unnest(val_out) %>%
      group_by(.id) %>%
      summarize(accuracy = cor(value, pred_value)) %>%
      summarize(pocv00 = mean(accuracy))
    
    
    cbind(test_acc, val_acc)
    
  }) %>% ungroup() %>%
  gather(scheme, accuracy, cv00, pocv00)





## Create training and test sets
## POV00 = untested gen1 in untested env
pov00_train_test <- test_env_df %>%
  group_by(set, trait) %>%
  do({
    row <- .
    
    ## Subset the phenotype data for the particular trait
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
    
    # Get a vector of testing environments
    row1 <- row %>% 
      unnest(test_env) %>%
      crossing(., select(row, train_env)) %>%
      mutate(val_environment = test_env,
             test_env = map(test_env, ~.),
             train_env = map2(train_env, test_env, setdiff),
             nTrainEnv = map_dbl(train_env, length)) %>%
      mutate(train = map(train_env, ~filter(pheno_use, environment %in% ., line_name %in% tp_geno)),
             test = map(test_env, ~filter(pheno_use, line_name %in% vp_geno, environment %in% .) %>% 
                          mutate(line_name = as.character(line_name))))
    
    
  }) %>% ungroup()




## Predict
pov00_predictions <- pov00_train_test %>%
  group_by(set, trait, val_environment, nTrainEnv) %>%
  do(pov00 = {
    
    row <- .
    
    ## Calculate genotype means
    data_adj <- geno_means(data = row$train[[1]])
    
    ## Run the base predictions
    gblup(K = K, train = data_adj, test = row$test[[1]], fit.env = FALSE)[[1]]
    # base_pred <- gblup(K = K, train = row$train[[1]], test = row$test[[1]], fit.env = TRUE)$accuracy

    
  }) %>% ungroup() %>%
  unnest(pov00) %>%
  gather(scheme, accuracy, pov00)



## CV0 - predict observed genotypes in unobserved environments
# Generate training and test sets
cv0_train_test <- test_env_df %>%
  group_by(set, trait) %>%
  do({
    row <- .
    
    ## Subset the phenotype data for the particular trait
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
    
    row2 <- row %>% 
      unnest(test_env) %>%
      crossing(., select(row, train_env)) %>%
      mutate(val_environment = test_env,
             test_env = map(test_env, ~.),
             train_env = map2(train_env, test_env, setdiff))
    
    
    # Create testing and training sets based on that environment and the cv samples
    row2 %>%
      mutate(train = map(train_env, ~filter(pheno_use, environment %in% .)),
             test = map(test_env, ~filter(pheno_use, environment %in% .)),
             nTrainEnv = map_dbl(train_env, length))
    
  }) %>% ungroup()



## Predict
cv0_predictions <- cv0_train_test %>%
  group_by(set, trait, val_environment, nTrainEnv) %>%
  do({
    
    row <- .
    
    ## Calculate genotype means for the training sets
    data_adj <- geno_means(data = row$train[[1]])
    
    ## Make predictions
    test_val_pred <- gblup(K = K, train = data_adj, test = row$test[[1]], fit.env = FALSE)[[2]]
    
    ## Calculate accuracy
    test_val_pred %>%
      mutate(scheme = ifelse(line_name %in% tp_geno, "cv0", "pocv0")) %>%
      group_by(scheme) %>%
      summarize(accuracy = cor(value, pred_value))
    
  }) %>% ungroup()


# 
# 
# 
# 
# ## Create training and test sets
# ## POV1 = untested gen1 in tested env
# pov1_train_test <- test_env_df %>%
#   filter(set == "complete") %>%
#   group_by(set, trait) %>%
#   do({
#     row <- .
#     
#     # Get a vector of testing environments
#     row1 <- row %>% 
#       unnest(test_env) %>%
#       rename(val_environment = test_env)
# 
#     # Create testing and training sets based on that environment
#     row1 %>%
#       mutate(train = map2(.x = trait, .y = val_environment, ~filter(S2_MET_BLUEs_tomodel, trait == .x) %>% droplevels() %>%
#                             filter(line_name %in% tp_geno)),
#              test = map2(.x = trait, .y = val_environment, ~filter(S2_MET_BLUEs_tomodel, line_name %in% vp_geno, trait == .x, environment == .y) %>%
#                            mutate(line_name = as.character(line_name))))
# 
#     
#   }) %>% ungroup()
# 
# 
# 
# 
# ## Predict
# pov1_predictions <- pov1_train_test %>%
#   group_by(set, trait, val_environment) %>%
#   do({
#     
#     row <- .
#     
#     ## Calculate genotype means
#     data_adj <- geno_means(data = row$train[[1]])
#     
#     ## Run the base predictions
#     base_pred <- gblup(K = K, train = data_adj, test = row$test[[1]], fit.env = FALSE)
#     # base_pred <- gblup(K = K, train = row$train[[1]], test = row$test[[1]], fit.env = TRUE)$accuracy
#     
#     ## Calculate accuracy and return
#     base_pred$pgv %>%
#       group_by(environment) %>% 
#       summarize(accuracy = cor(value, pred_value))
#     
#     
#   }) %>% ungroup() %>%
#   select(set, trait, val_environment = environment, accuracy) %>%
#   mutate(scheme = "pov1")
# 
# 
# ## Create training and test sets
# ## CV1 = untested gen1 in tested env
# cv1_train_test <- test_env_df %>%
#   filter(set == "complete") %>%
#   group_by(set, trait) %>%
#   do({
#     row <- .
#     
#     # Get a vector of testing environments
#     row1 <- row %>% 
#       unnest(test_env) %>%
#       rename(val_environment = test_env)
#     
#     # Create testing and training sets based on that environment
#     crossing(row1, cv1_cv00_train_test) %>%
#       mutate(train = map2(.x = trait, .y = train, ~filter(S2_MET_BLUEs_tomodel, trait == .x) %>% droplevels() %>%
#                             filter(line_name %in% .y)),
#              test = pmap(list(trait, val_environment, test, val), ~filter(S2_MET_BLUEs_tomodel, line_name %in% c(..3, ..4), trait == ..1, environment %in% ..2) %>% 
#                            mutate(line_name = as.character(line_name))))
#     
#     
#   }) %>% ungroup()
# 
# 
# 
# 
# ## Predict
# cv1_predictions <- cv1_train_test %>%
#   separate(.id, c("rep", ".id"), sep = "_") %>%
#   group_by(set, trait, val_environment, rep) %>%
#   do({
#     
#     df <- .
#     
#     ## Calculate genotype means for the training sets
#     data_adj <- map(df$train, ~geno_means(data = .))
#     
#     ## Make predictions
#     test_val_pred <- map2(.x = data_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)$pgv)
#     
#     ## Add results to the row
#     df1 <- df %>%
#       mutate(test_out = map(test_val_pred, ~filter(., !line_name %in% vp)),
#              val_out = map(test_val_pred, ~filter(., line_name %in% vp))) %>%
#       select(.id, test_out, val_out)
#     
#     ## Calculate test accuracy
#     test_acc <- df1 %>%
#       unnest(test_out) %>% 
#       group_by(environment) %>%
#       summarize(cv1 = cor(value, pred_value))
#     
#     ## Calculate average val accuracy
#     val_acc <- df1 %>%
#       unnest(val_out) %>%
#       group_by(environment, .id) %>%
#       summarize(accuracy = cor(value, pred_value)) %>%
#       summarize(pocv1 = mean(accuracy))
#     
#     full_join(test_acc, val_acc, by = "environment")
#     
#     
#   }) %>% ungroup() %>%
#   gather(scheme, accuracy, cv1, pocv1) %>%
#   select(set, trait, val_environment = environment, rep, scheme, accuracy)








# ## Cross-validation 2
# 
# # Just the TP - predict!
# cv2_predictions <- cv2_train_test %>%
#   group_by(trait, rep) %>%
#   do({
#     
#     df <- .
#     
#     ## Calculate genotype means for the training sets
#     data_adj <- map(df$train, ~geno_means(data = .))
#     
#     ## Make predictions
#     test_val_pred <- map2(.x = data_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)$pgv)
#     
#     ## Add results to the row
#     df1 <- df %>%
#       mutate(test_out = test_val_pred) %>%
#       select(.id, test_out)
#     
#     ## Calculate test accuracy and return
#     df1 %>%
#       unnest(test_out) %>% 
#       group_by(environment) %>%
#       summarize(accuracy = cor(value, pred_value)) %>%
#       mutate(scheme = "cv2")
#     
#   }) %>% ungroup()
# 
# 
# ## Now all information - TP and VP
# pocv2_predictions <- pocv2_train_test %>%
#   group_by(trait, rep) %>%
#   do({
#     
#     df <- .
#     
#     ## Calculate genotype means for the training sets
#     data_adj <- map(df$train, ~geno_means(data = .))
#     
#     ## Make predictions
#     test_val_pred <- map2(.x = data_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)$pgv)
#     
#     ## Add results to the row
#     df1 <- df %>%
#       mutate(test_out = test_val_pred) %>%
#       select(.id, test_out)
#     
#     ## Calculate test accuracy and return
#     df1 %>%
#       unnest(test_out) %>% 
#       group_by(environment) %>%
#       summarize(accuracy = cor(value, pred_value)) %>%
#       mutate(scheme = "pocv2")
#     
#   }) %>% ungroup()












##### Cluster-based Predictions ######


clusters <- cluster_df %>% 
  left_join(., test_env_df) %>%
  mutate(nClusters = map_dbl(cluster, ~n_distinct(.$cluster))) %>%
  ## You can't have AMMI and "realistic"
  filter(!(set == "realistic" & model == "AMMI"))

## Unnest for future use
clusters1 <- unnest(clusters, cluster)

cluster_train_test <- clusters %>%
  mutate(data = list(NULL))








pov00_cluster_train_test <- clusters %>%
  group_by(set, trait, model, nClusters) %>%
  do({
    
    row <- .
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
    
    clus <- unnest(row, cluster)
  
    row1 <- row %>% 
      unnest(test_env) %>% 
      mutate(train_env = row$train_env, 
             train_env = map2(train_env, test_env, setdiff))
    
    ## Get the data for the cluster of the test environment
    row1 %>%
      mutate(data = map2(test_env, train_env, ~{
        cluster_val <- subset(clus, environment == .x, cluster, drop = T)
        cluster_train <- intersect(subset(clus, cluster == cluster_val, environment, drop = T), .y)

        ## Combine into a single DF
        tibble(
          train = list(filter(pheno_use, environment %in% cluster_train, line_name %in% tp_geno)),
          test = list(filter(pheno_use, environment %in% .x, line_name %in% vp_geno)),
          cluster = cluster_val,
          nTrainEnv = length(cluster_train)
        )
        
      })) %>% unnest(data)
        
    
  }) %>% ungroup() %>%
  rename(val_environment = test_env)
  
  


## This is the normal, POV00 predictions (untested genos in untested envs)
pov00_cluster_predictions <- pov00_cluster_train_test %>%
  group_by(set, model, trait, val_environment, nClusters, nTrainEnv, cluster) %>%
  do(pov00 = {

    row <- .

    ## Calculate genotype means
    data_adj <- geno_means(data = row$train[[1]])
    
    ## Run the base predictions
    gblup(K = K, train = data_adj, test = row$test[[1]], fit.env = FALSE)[[1]]
    # base_pred <- gblup(K = K, train = row$train[[1]], test = row$test[[1]], fit.env = TRUE)$accuracy

  }) %>% ungroup() %>%
  unnest(pov00) %>%
  gather(scheme, accuracy, pov00)



## CV00 predictions (with POCV00)
cv00_cluster_train_test <- clusters %>%
  group_by(set, trait, model, nClusters) %>%
  do({
    
    row <- .
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
    
    clus <- unnest(row, cluster)
    
    row1 <- row %>% 
      unnest(test_env) %>% 
      mutate(train_env = row$train_env, 
             train_env = map2(train_env, test_env, setdiff)) %>%
      ## Remove environments not in pheno_use
      filter(test_env %in% clus$environment)
    
    ## Get the data for the cluster of the test environment
    row1 %>%
      mutate(data = map2(test_env, train_env, ~{
        testE <- .x
        cluster_val <- subset(clus, environment == testE, cluster, drop = T)
        cluster_train <- intersect(subset(clus, cluster == cluster_val, environment, drop = T), .y)
        
        ## Generate CV sets
        replicate(n = nCV, expr = crossv_kfold(data = cv_tp_df), simplify = FALSE) %>%
          map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y)) %>%
          mutate_at(vars(train, test), funs(map(., as.data.frame))) %>%
          mutate(train = map(train, ~filter(pheno_use, environment %in% cluster_train, line_name %in% .[[1]])),
                 test = map(test, ~filter(pheno_use, environment %in% testE, line_name %in% c(.[[1]], vp_geno))),
                 cluster = cluster_val,
                 nTrainEnv = length(cluster_train))
        
      })) %>% unnest(data)
    
    
  }) %>% ungroup() %>%
  rename(val_environment = test_env)
    


## Predict!
cv00_cluster_predictions <- cv00_cluster_train_test %>%
  group_by(set, model, trait, val_environment, nClusters, cluster, nTrainEnv, rep) %>%
  do({
    
    df <- .
    train_data <- df$train
    
    ## Calculate genotype means for the training sets
    data_adj <- map(train_data, ~geno_means(data = .))
    
    ## Make predictions
    test_val_pred <- map2(.x = data_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]])
    
    ## Add results to the row
    df1 <- df %>% 
      mutate(pred = test_val_pred) %>% 
      unnest(pred) %>% 
      mutate(scheme = ifelse(line_name %in% tp, "cv00", "pocv00")) %>%
      select(.id, scheme, line_name, value, pred_value)
      
    
    ## Calculate test accuracy
    test_acc <- df1 %>%
      filter(scheme == "cv00") %>% 
      summarize(cv00 = cor(value, pred_value))
    
    ## Calculate average val accuracy
    val_acc <- df1 %>% 
      filter(scheme == "pocv00") %>% 
      group_by(.id) %>% 
      summarize(accuracy = cor(value, pred_value)) %>%
      summarize(pocv00 = mean(accuracy))
    
    
    cbind(test_acc, val_acc)
    
  }) %>% ungroup()



## CV0 - predict observed genotypes in unobserved environments
cv0_cluster_train_test <- clusters %>%
  group_by(set, trait, model, nClusters) %>%
  do({
    
    row <- .
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
    
    clus <- unnest(row, cluster)
    
    row1 <- row %>% 
      unnest(test_env) %>% 
      mutate(train_env = row$train_env, 
             train_env = map2(train_env, test_env, setdiff)) %>%
      ## Remove environments not in pheno_use
      filter(test_env %in% clus$environment)
    
    ## Get the data for the cluster of the test environment
    row1 %>%
      mutate(data = map2(test_env, train_env, ~{
        testE <- .x
        cluster_val <- subset(clus, environment == testE, cluster, drop = T)
        cluster_train <- intersect(subset(clus, cluster == cluster_val, environment, drop = T), .y)
        
        tibble(
          train = list(filter(pheno_use, environment %in% cluster_train)),
          test = list(filter(pheno_use, environment %in% testE)),
          cluster = cluster_val,
          nTrainEnv = length(cluster_train)
        )

      })) %>% unnest(data)
    
    
  }) %>% ungroup() %>%
  rename(val_environment = test_env)


## Predict!
cv0_cluster_predictions <- cv0_cluster_train_test %>%
  group_by(set, model, trait, val_environment, nClusters, cluster, nTrainEnv) %>%
  do({
    
    row <- .
    train_data <- row$train[[1]]
    
    ## Calculate genotype means for the training sets
    data_adj <- geno_means(data = train_data)
    
    ## Make predictions
    test_val_pred <- gblup(K = K, train = data_adj, test = row$test[[1]], fit.env = FALSE)[[2]]
    
    ## Add results to the row
    test_val_pred %>% 
      mutate(scheme = ifelse(line_name %in% tp, "cv0", "pocv0")) %>%
      group_by(scheme) %>%
      summarize(accuracy = cor(value, pred_value))
    
  }) %>% ungroup()



# ## CV1 predictions
# 
# ## Create training and test sets for CV1 and POV1
# ## POV1 = untested gen1 in tested env
# cv1_pov1_cluster_train_test <- bind_rows(pov1_train_test, cv1_train_test) %>%
#   left_join(., clusters) %>%
#   group_by(set, trait, val_environment, nClusters, model, .id) %>%
#   do(data = {
#     
#     row <- .
#     train_data <- row$train[[1]]
# 
#     ## Extract the cluster from each method wherein the validation environment resides
#     cluster_val <- subset(row$cluster[[1]], environment == row$val_environment, cluster, drop = T)
#     ## Other environments in the cluster?
#     cluster_train <- subset(row$cluster[[1]], cluster == cluster_val, environment, drop = T)
#     
#     ## Subset the training data
#     train_data1 <- subset(train_data, environment %in% cluster_train)
#     test1 <- row$test[[1]]
#     
#     cluster_name <- cluster_val
#       
#     list(train = train_data1, test = test1, nTrainEnv = length(cluster_train), cluster = cluster_name)
#     
#   }) %>% ungroup()
# 
# 
# 
# 
# ## Predict
# cv1_pov1_cluster_predictions <- cv1_pov1_cluster_train_test %>%
#   separate(.id, c("rep", ".id"), sep = "_") %>%
#   group_by(set, trait, val_environment, nClusters, model, rep) %>%
#   do({
#     
#     df <- .
#     
#     ## Flow for pov or cv
#     if (all(is.na(df$.id))) {
#     
#       train_data <- df$data[[1]]$train %>% mutate(cluster = df$data[[1]]$cluster)
#       
#       ## Calculate genotype means for the training sets
#       data_adj <- geno_means(data = train_data)
#       
#       base_pred <- gblup(K = K, train = data_adj, test = df$data[[1]]$test, fit.env = FALSE)
#       # base_pred <- gblup(K = K, train = row$train[[1]], test = row$test[[1]], fit.env = TRUE)$accuracy
#       
#       ## Calculate accuracy and return
#       base_pred$pgv %>%
#         group_by(environment) %>% 
#         summarize(accuracy = cor(value, pred_value)) %>%
#         mutate(cluster = df$data[[1]]$cluster, nTrainEnv = df$data[[1]]$nTrainEnv, scheme = "pov1")
#       
#     } else {
#       
#       train_data <- map(df$data, "train") %>% map2(.x = ., .y = map(df$data, "cluster"), ~mutate(.x, cluster = .y))
#       test <- map(df$data, "test")
#       
#       ## Calculate genotype means for the training sets
#       data_adj <- map(train_data, ~geno_means(data = .))
#       
#       ## Make predictions
#       test_val_pred <- map2(.x = data_adj, .y = test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)$pgv)
#       
#       ## Add results to the row
#       df1 <- df %>% 
#         mutate(pred = test_val_pred, cluster = map_dbl(data, "cluster")) %>% 
#         unnest(pred) %>% 
#         mutate(scheme = ifelse(line_name %in% tp, "cv1", "pocv1")) %>%
#         select(environment, .id, scheme, cluster, line_name, value, pred_value)
#       
#       ## Calculate test accuracy
#       test_acc <- df1 %>%
#         filter(scheme == "cv1") %>% 
#         group_by(cluster, environment, scheme) %>% 
#         summarize(accuracy = cor(value, pred_value)) %>%
#         ungroup()
#       
#       ## Calculate average val accuracy
#       val_acc <- df1 %>% 
#         filter(scheme == "pocv1") %>% 
#         group_by(cluster, environment, scheme, .id) %>% 
#         summarize(accuracy = cor(value, pred_value)) %>%
#         summarize(accuracy = mean(accuracy)) %>%
#         ungroup()
#       
#       rbind(test_acc, val_acc) %>%
#         mutate(nTrainEnv = df$data[[1]]$nTrainEnv)
#       
#     }
#     
#     
#   }) %>% ungroup() %>%
#   select(set, model, trait, val_environment = environment, rep, cluster, nClusters, nTrainEnv, scheme, accuracy)




# ## Cross-validation 2 - Fill-in-the-gaps
# ## 
# ## 
# cv2_cluster_train_test <- cv2_train_test %>%
#   left_join(., filter(clusters, set == "complete")) %>%
#   mutate(i = seq(nrow(.))) %>%
#   group_by(set, trait, model, nClusters, rep, .id) %>%
#   do(data = {
#     row <- .
# 
#     train_data <- row$train[[1]]
#     
#     ## For each cluster, get all validation environments
#     clus <- row$cluster[[1]]
#     
#     train_data1 <- left_join(train_data, clus, by = "environment") %>%
#       split(.$cluster)
#     test1 <- left_join(row$test[[1]], clus, by = "environment") %>%
#       split(.$cluster)
#     
#     cluster_name <- names(train_data1)
# 
#     list(train = train_data1, test = test1, nTrainEnv = map_dbl(train_data1, ~n_distinct(.$environment)), cluster = cluster_name)
#     
#   }) %>% ungroup()
# 
# 
# cv2_cluster_predictions <- cv2_cluster_train_test %>%
#   group_by(set, trait, model, nClusters, rep) %>%
#   do({
#     
#     df <- .
#     
#     test <- map(df$data, "test")
#     train_data <- map(df$data, "train") %>% map(~.[names(test[[1]])])
#     
#     ## Calculate genotype means for the training sets
#     data_adj <- map(train_data, ~map(., ~geno_means(data = .)))
#     
#     ## Make predictions
#     test_val_pred <- map2(.x = data_adj, .y = test, ~map2(.x, .y, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)$pgv))
#     
#     ## Add results to the row
#     df1 <- df %>% 
#       mutate(pred = test_val_pred) %>% 
#       unnest(pred) %>% 
#       mutate(nTrainEnv = unlist(map(df$data, "nTrainEnv"))) %>%
#       unnest(pred) %>%
#       mutate(scheme = ifelse(line_name %in% tp, "cv2", "pocv2")) %>%
#       select(environment, scheme, cluster, nTrainEnv, line_name, value, pred_value)
#     
#     ## Calculate accuracy
#     df1 %>%
#       group_by(cluster, environment, scheme) %>%
#       summarize(accuracy = cor(value, pred_value), nTrainEnv = unique(nTrainEnv)) %>%
#       ungroup()
#     
#     
#   }) %>% ungroup() %>%
#   rename(val_environment = environment)
# 
# 
# ## POCV2
# pocv2_cluster_train_test <- pocv2_train_test %>%
#   left_join(., filter(clusters, set == "complete")) %>%
#   mutate(i = seq(nrow(.))) %>%
#   select(-train_env, -test_env, -data) %>%
#   group_by(set, trait, model, nClusters, rep, .id) %>%
#   do(data = {
#     row <- .
#     
#     train_data <- row$train[[1]]
#     
#     ## For each cluster, get all validation environments
#     clus <- row$cluster[[1]]
#     
#     train_data1 <- left_join(train_data, clus, by = "environment") %>%
#       split(.$cluster)
#     test1 <- left_join(row$test[[1]], clus, by = "environment") %>%
#       split(.$cluster)
#     
#     cluster_name <- names(train_data1)
#     
#     list(train = train_data1, test = test1, nTrainEnv = map_dbl(train_data1, ~n_distinct(.$environment)), cluster = cluster_name)
#     
#   }) %>% ungroup()
# 
# 
# pocv2_cluster_predictions <- pocv2_cluster_train_test %>%
#   group_by(set, trait, model, nClusters, rep) %>%
#   do({
#     
#     df <- .
#     
#     test <- map(df$data, "test")
#     train_data <- map(df$data, "train") %>% map(~.[names(test[[1]])])
#     
#     ## Calculate genotype means for the training sets
#     data_adj <- map(train_data, ~map(., ~geno_means(data = .)))
#     
#     ## Make predictions
#     test_val_pred <- map2(.x = data_adj, .y = test, ~map2(.x, .y, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)$pgv))
#     
#     ## Add results to the row
#     df1 <- df %>% 
#       mutate(pred = test_val_pred) %>% 
#       unnest(pred) %>% 
#       mutate(nTrainEnv = unlist(map(df$data, "nTrainEnv"))) %>%
#       unnest(pred) %>%
#       mutate(scheme = ifelse(line_name %in% tp, "cv2", "pocv2")) %>%
#       select(environment, scheme, cluster, nTrainEnv, line_name, value, pred_value)
#     
#     ## Calculate accuracy
#     df1 %>%
#       group_by(cluster, environment, scheme) %>%
#       summarize(accuracy = cor(value, pred_value), nTrainEnv = unique(nTrainEnv)) %>%
#       ungroup()
#     
#     
#   }) %>% ungroup() %>%
#   rename(val_environment = environment)










#### Random environment cross-validation / pov



## Generate pov00 train/test sets
## POV00 = untested geno in untested env
pov00_train_test_random <- test_env_df %>%
  group_by(set, trait) %>%
  do({
    row <- .
    
    ## Subset phenos for the trait
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
    
    # Get a vector of testing environments
    row1 <- row %>% 
      unnest(test_env) %>%
      rename(val_environment = test_env) %>%
      group_by(set, trait, val_environment) %>%
      ## Create lists of random test environment
      do({
        r <- .
        possible_env <- setdiff(unique(pheno_use$environment), r$val_environment)
        
        test_env <- map(n_rand_env, ~replicate(n = nCV, sample(x = possible_env, size = .), simplify = FALSE))
        
        tibble(nTrainEnv = n_rand_env, test_env) %>% 
          unnest(test_env) %>% 
          group_by(nTrainEnv) %>%
          mutate(rep = seq(n())) %>% ungroup()
        
      }) %>% ungroup()
    
    
    # Create testing and training sets based on that environment
    row1 %>% 
      mutate(train = map(test_env, ~filter(pheno_use, environment %in% ., line_name %in% tp_geno) %>% droplevels() %>% 
                           mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))),
             test = map(val_environment, ~filter(pheno_use, line_name %in% vp_geno, environment == .) %>% mutate(line_name = as.character(line_name))))
      
  }) %>% ungroup()



## Predict
pov00_predictions_random <- pov00_train_test_random %>%
  mutate(n = seq(nrow(.))) %>%
  group_by(set, trait, val_environment, nTrainEnv, rep) %>%
  do(accuracy = {
    
    row <- .
  
    ## Calculate genotype means
    data_adj <- geno_means(data = row$train[[1]])
    
    ## Run the base predictions
    gblup(K = K, train = data_adj, test = row$test[[1]], fit.env = FALSE)[[1]]
    # gblup(K = K, train = row$train[[1]], test = row$test[[1]], fit.env = FALSE)$accuracy
    
  }) %>% ungroup() %>%
  unnest(accuracy) %>%
  mutate(scheme = "pov00")






### Leave-one-environment-out cross-validation (CV00)
# Generate training and test sets
cv00_train_test_random  <- test_env_df %>%
  group_by(set, trait) %>%
  do({
    row <- .
    
    ## Subset phenos for the trait
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
    # Unique environments
    unique_env <- unique(pheno_use$environment)
    # Data.frame of TP lines to randomize
    tp_lines_df <- data_frame(line_name = tp_geno)
    
    
    # Get a vector of testing environments
    row1 <- row %>% 
      unnest(test_env) %>%
      rename(val_environment = test_env) %>%
      group_by(set, trait, val_environment) %>%
      ## Create lists of random test environment
      do({
        r <- .
        possible_env <- setdiff(unique_env, r$val_environment)
        test_env <- map(n_rand_env, ~replicate(n = nCV, sample(x = possible_env, size = .), simplify = FALSE))
        
        ## For each test env set, generate a new cross-validation set
        train_cv <- map(test_env, ~map(., ~mutate(crossv_kfold(data = tp_lines_df, k = k), test_env = list(.))))
        
        
        tibble(nTrainEnv = n_rand_env, train_cv) %>% 
          unnest(train_cv) %>% 
          group_by(nTrainEnv) %>%
          mutate(rep = seq(n())) %>% 
          ungroup() %>%
          unnest(train_cv) %>%
          unite(.id, c("rep", ".id"), sep = "_") %>%
          mutate_at(vars(train, test), funs(lapply(., FUN = as.data.frame))) %>%
          mutate_at(vars(train, test), funs(lapply(., pull)))
        
      }) %>% ungroup()
    
    
    # Create testing and training sets based on that environment
    row1 %>% 
      mutate(train = map2(.x = test_env, .y = train, ~filter(pheno_use, environment %in% .x, line_name %in% .y) %>% droplevels() %>% 
                           mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))),
             test = map2(.x = val_environment, .y = test, ~filter(pheno_use, line_name %in% c(.y, vp_geno), environment == .x) %>% 
                           mutate(line_name = as.character(line_name))))
    
  }) %>% ungroup()
    
    
    
## Predict
cv00_predictions_random <- cv00_train_test_random %>%
  separate(.id, c("rep", ".id"), sep = "_") %>%
  group_by(set, trait, val_environment, nTrainEnv, rep) %>%
  do({
    
    df <- .
    
    ## Calculate genotype means for the training sets
    data_adj <- map(df$train, ~geno_means(data = .))
    
    ## Make predictions
    test_val_pred <- map2(.x = data_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]])
    
    ## Add results to the row
    df1 <- df %>%
      mutate(test_out = map(test_val_pred, ~filter(., !line_name %in% vp)),
             val_out = map(test_val_pred, ~filter(., line_name %in% vp))) %>%
      select(.id, test_out, val_out)
    
    ## Calculate test accuracy
    test_acc <- df1 %>%
      unnest(test_out) %>% 
      summarize(cv00 = cor(value, pred_value))
    
    ## Calculate average val accuracy
    val_acc <- df1 %>%
      unnest(val_out) %>%
      group_by(.id) %>%
      summarize(accuracy = cor(value, pred_value)) %>%
      summarize(pocv00 = mean(accuracy))
    
    
    cbind(test_acc, val_acc)
    
  }) %>% ungroup() %>%
  gather(scheme, accuracy, cv00, pocv00)








### CV0 - tested genotyped in untested environments
# Generate training and test sets
cv0_train_test_random  <- test_env_df %>%
  group_by(set, trait) %>%
  do({
    row <- .
    
    ## Subset phenos for the trait
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
    # Unique environments
    unique_env <- unique(pheno_use$environment)
    
    
    # Get a vector of testing environments
    row1 <- row %>% 
      unnest(test_env) %>%
      rename(val_environment = test_env) %>%
      group_by(set, trait, val_environment) %>%
      ## Create lists of random test environment
      do({
        r <- .
        possible_env <- setdiff(unique_env, r$val_environment)
        test_env <- map(n_rand_env, ~replicate(n = nCV, sample(x = possible_env, size = .), simplify = FALSE))
        
        tibble(nTrainEnv = n_rand_env, test_env) %>% 
          unnest(test_env) %>% 
          group_by(nTrainEnv) %>%
          mutate(rep = seq(n())) %>% 
          ungroup()
        
      }) %>% ungroup()
    
    
    # Create testing and training sets based on that environment
    row1 %>% 
      mutate(train = map(test_env, ~filter(pheno_use, environment %in% .) %>% droplevels() %>% 
                            mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))),
             test = map(val_environment, ~filter(pheno_use, environment == .) %>% mutate(line_name = as.character(line_name))))
    
  }) %>% ungroup()


## Predict
cv0_predictions_random <- cv0_train_test_random %>%
  group_by(set, trait, val_environment, nTrainEnv, rep) %>%
  do({
    
    row <- .
    
    ## Calculate genotype means
    data_adj <- geno_means(data = row$train[[1]])
    
    ## Make predictions
    test_val_pred <- gblup(K = K, train = data_adj, test = row$test[[1]], fit.env = FALSE)[[2]]
      
    test_val_pred %>% 
      mutate(scheme = ifelse(line_name %in% tp_geno, "cv0", "pocv0")) %>%
      group_by(scheme) %>%
      summarize(accuracy = cor(value, pred_value))
    
  }) %>% ungroup() 



# ## CV1 - untested genotypes in tested environments
# # Generate training and test sets
# cv1_train_test_random  <- test_env_df %>%
#   filter(set == "complete") %>%
#   group_by(set, trait) %>%
#   do({
#     row <- .
#     
#     ## Subset phenos for the trait
#     pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
#     # Unique environments
#     unique_env <- unique(pheno_use$environment)
#     # Data.frame of TP lines to randomize
#     tp_lines_df <- data_frame(line_name = tp_geno)
#     
#     # Get a vector of testing environments
#     row1 <- row %>% 
#       unnest(test_env) %>%
#       rename(val_environment = test_env) %>%
#       group_by(set, trait, val_environment) %>%
#       ## Create lists of random test environment
#       do({
#         r <- .
#         possible_env <- setdiff(unique_env, r$val_environment)
#         test_env <- map(n_rand_env, ~replicate(n = nCV, sample(x = possible_env, size = .), simplify = FALSE))
#         
#         ## For each test env set, generate a new cross-validation set
#         train_cv <- map(test_env, ~map(., ~mutate(crossv_kfold(data = tp_lines_df, k = k), test_env = list(.))))
#         
#         
#         tibble(nTrainEnv = n_rand_env, train_cv) %>% 
#           unnest(train_cv) %>% 
#           group_by(nTrainEnv) %>%
#           mutate(rep = seq(n())) %>% 
#           ungroup() %>%
#           unnest(train_cv) %>%
#           unite(.id, c("rep", ".id"), sep = "_") %>%
#           mutate_at(vars(train, test), funs(lapply(., FUN = as.data.frame))) %>%
#           mutate_at(vars(train, test), funs(lapply(., pull)))
#         
#       }) %>% ungroup()
#     
#     
#     # Create testing and training sets based on that environment
#     row1 %>% 
#       mutate(train = pmap(list(test_env, train, val_environment), ~filter(pheno_use, environment %in% c(..1, ..3), line_name %in% ..2) %>% droplevels() %>% 
#                             mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))),
#              test = map2(.x = val_environment, .y = test, ~filter(pheno_use, line_name %in% c(.y, vp_geno), environment == .x) %>% 
#                            mutate(line_name = as.character(line_name))))
#     
#   }) %>% ungroup()
# 
# 
# ## Predict
# cv1_predictions_random <- cv1_train_test_random %>%
#   separate(.id, c("rep", ".id"), sep = "_") %>%
#   group_by(set, trait, val_environment, nTrainEnv, rep) %>%
#   do({
#     
#     df <- .
#     
#     ## Calculate genotype means for the training sets
#     data_adj <- map(df$train, ~geno_means(data = .))
#     
#     ## Make predictions
#     test_val_pred <- map2(.x = data_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]])
#     
#     ## Add results to the row
#     df1 <- df %>%
#       mutate(test_out = map(test_val_pred, ~filter(., !line_name %in% vp)),
#              val_out = map(test_val_pred, ~filter(., line_name %in% vp))) %>%
#       select(.id, test_out, val_out)
#     
#     ## Calculate test accuracy
#     test_acc <- df1 %>%
#       unnest(test_out) %>% 
#       summarize(cv1 = cor(value, pred_value))
#     
#     ## Calculate average val accuracy
#     val_acc <- df1 %>%
#       unnest(val_out) %>%
#       group_by(.id) %>%
#       summarize(accuracy = cor(value, pred_value)) %>%
#       summarize(pocv1 = mean(accuracy))
#     
#     
#     cbind(test_acc, val_acc)
#     
#   }) %>% ungroup() 
# 
# 
# ## POV1 - untested genotypes in tested environments
# # Generate training and test sets
# pov1_train_test_random  <- test_env_df %>%
#   filter(set == "complete") %>%
#   group_by(set, trait) %>%
#   do({
#     row <- .
#     
#     ## Subset phenos for the trait
#     pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
#     # Unique environments
#     unique_env <- unique(pheno_use$environment)
#     
#     # Get a vector of testing environments
#     row1 <- row %>% 
#       unnest(test_env) %>%
#       rename(val_environment = test_env) %>%
#       group_by(set, trait, val_environment) %>%
#       ## Create lists of random test environment
#       do({
#         r <- .
#         possible_env <- setdiff(unique_env, r$val_environment)
#         test_env <- map(n_rand_env, ~replicate(n = nCV, sample(x = possible_env, size = .), simplify = FALSE))
#         
#         tibble(nTrainEnv = n_rand_env, test_env) %>% 
#           unnest(test_env) %>% 
#           group_by(nTrainEnv) %>%
#           mutate(rep = seq(n())) %>% 
#           ungroup()
#         
#       }) %>% ungroup()
#     
#     
#     # Create testing and training sets based on that environment
#     row1 %>% 
#       mutate(train = map2(.x = test_env, .y = val_environment, ~filter(pheno_use, environment %in% c(.x, .y), line_name %in% tp_geno) %>% droplevels() %>% 
#                             mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))),
#              test = map(val_environment, ~filter(pheno_use, line_name %in% vp_geno, environment == .) %>% mutate(line_name = as.character(line_name))))
#     
#   }) %>% ungroup()
# 
# 
# ## Predict
# pov1_predictions_random <- pov1_train_test_random %>%
#   group_by(set, trait, val_environment, nTrainEnv, rep) %>%
#   do(accuracy = {
#     
#     row <- .
#     
#     ## Calculate genotype means
#     data_adj <- geno_means(data = row$train[[1]])
#     
#     ## Make predictions
#     gblup(K = K, train = data_adj, test = row$test[[1]], fit.env = FALSE)[[1]]
#     
#   }) %>% ungroup() %>%
#   unnest(accuracy) %>%
#   mutate(scheme = "pov1")
# 
# 
# ## CV2 - predict tested genos in tested envs
# cv2_train_test_random <- test_env_df %>%
#   filter(set == "complete") %>%
#   group_by(set, trait) %>%
#   do({
#     row <- .
#     
#     ## Subset phenos for the trait
#     pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
#     pheno_use_train <- filter(pheno_use, line_name %in% tp_geno)
#     # Unique environments
#     unique_env <- unique(pheno_use$environment)
#     
#     ## For each number of training environments, create cv sets
#     map(n_rand_env, ~replicate(n = nCV * 10, sample(x = unique_env, size = .), simplify = FALSE) %>%
#           map(~crossv_kfold(data = filter(pheno_use_train, environment %in% .), k = k)) %>%
#           map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y)) ) %>%
#       map2_df(.x = ., .y = n_rand_env, ~mutate(.x, nTrainEnv = .y)) %>%
#       unite(.id, c("rep", ".id"), sep = "_") %>%
#       mutate_at(vars(train, test), funs(map(., as.data.frame)))
# 
#     
#   }) %>% ungroup()
# 
# 
# ## Predict
# cv2_predictions_random <- cv2_train_test_random %>%
#   mutate(n = seq(nrow(.))) %>%
#   separate(.id, c("rep", ".id"), sep = "_") %>%
#   group_by(set, trait, nTrainEnv, rep) %>%
#   do({
#     
#     df <- .
#     
#     ## Calculate genotype means for the training sets
#     data_adj <- map(df$train, ~geno_means(data = .))
#     
#     ## Make predictions
#     test_val_pred <- map2(.x = data_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]])
#     
#     ## Add results to the row
#     bind_rows(test_val_pred) %>% 
#       group_by(environment) %>%
#       summarize(accuracy = cor(value, pred_value))
# 
#   }) %>% ungroup() %>%
#   mutate(scheme = "cv2")
# 
# 
# ## POCV2 - predict tested vp in tested environments using complete tp data
# pocv2_train_test_random <- test_env_df %>%
#   filter(set == "complete") %>%
#   group_by(set, trait) %>%
#   do({
#     row <- .
#     
#     ## Subset phenos for the trait
#     pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
#     pheno_use_train <- filter(pheno_use, line_name %in% tp_geno)
#     pheno_use_vp <- filter(pheno_use, line_name %in% vp_geno)
#     # Unique environments
#     unique_env <- unique(pheno_use$environment)
#     
#     ## For each number of training environments, create pocv sets
#     map(n_rand_env, ~replicate(n = nCV * 10, sample(x = unique_env, size = .), simplify = FALSE) %>%
#           map(~crossv_kfold(data = filter(pheno_use_vp, environment %in% .), k = k)) %>%
#           map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y)) ) %>%
#       map2_df(.x = ., .y = n_rand_env, ~mutate(.x, nTrainEnv = .y)) %>%
#       unite(.id, c("rep", ".id"), sep = "_") %>%
#       mutate_at(vars(train, test), funs(map(., as.data.frame))) %>%
#       ## Add available training data
#       mutate(train = map(train, ~bind_rows(subset(pheno_use_train, environment %in% unique(.$environment)), .)))
#     
#     
#   }) %>% ungroup()
# 
# 
# ## Predict
# pocv2_predictions_random <- pocv2_train_test_random %>%
#   separate(.id, c("rep", ".id"), sep = "_") %>%
#   group_by(set, trait, nTrainEnv, rep) %>%
#   do({
#     
#     df <- .
#     
#     ## Calculate genotype means for the training sets
#     data_adj <- map(df$train, ~geno_means(data = .))
#     
#     ## Make predictions
#     test_val_pred <- map2(.x = data_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]])
#     
#     ## Add results to the row
#     bind_rows(test_val_pred) %>% 
#       group_by(environment) %>%
#       summarize(accuracy = cor(value, pred_value))
#     
#   }) %>% ungroup() %>%
#   mutate(scheme = "pocv2")




# List of random predictions
random_pred <- c("cv0_predictions_random", "cv00_predictions_random", "cv1_predictions_random", "cv2_predictions_random", 
                 "pocv2_predictions_random", "pov00_predictions_random", "pov1_predictions_random")



##













##### Random cluster predictions #####




## Number of random cluster assignments
n_random <- nCV


clusters <- cluster_df %>% 
  left_join(., test_env_df) %>%
  mutate(nClusters = map_dbl(cluster, ~n_distinct(.$cluster))) %>%
  ## You can't have AMMI and "realistic"
  filter(!(set == "realistic" & model == "AMMI"))
  

## Training lines for CV
cv_train_df <- data.frame(line_name = tp_geno, stringsAsFactors = FALSE)

## CV00 - predict new genotypes in new environments ##
# For each environment in method, randomly sample n_random sets of the same number of training environments. For each n_random,
# generate nCV cross-validation samples
cv00_cluster_random_train_test <- clusters %>%
  group_by(set, trait, model) %>%
  do({
    row <- .
    
    # Phenotypes to use
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait %in% row$trait)
    
    ## Possible training environments
    train_env <- row$train_env[[1]]
    ## Test environments
    test_env <- row$test_env[[1]]
    
    ## Calculate nTrainEnv depending on the set
    if (row$set == "complete") {
    
      ## Unnest clusters
      clus1 <- unnest(row, cluster) %>%
        group_by(cluster) %>%
        mutate(nTrainEnv = n() - 1) %>%
        ungroup()
      
    } else {
      
      ## Unnest clusters
      clus1 <- unnest(row, cluster) %>%
        group_by(cluster) %>%
        mutate(nTrainEnv = sum(environment %in% train_env)) %>%
        ungroup()
      
    }
    
    clus2 <- clus1 %>%
      filter(environment %in% test_env) %>%
      ## For each environment, sample nTrainEnv environments from the train_env pool
      mutate(random = map2(.x = environment, .y = nTrainEnv, ~replicate(n = n_random, sample(setdiff(train_env, .x), size = .y), simplify = FALSE)),
             cv_sample = map(environment, ~replicate(n = nCV, expr = crossv_kfold(data = cv_train_df, k = k), simplify = FALSE)),
             cv_sample = map(cv_sample, ~map(., ~mutate_at(., vars(-.id), funs(map(., as.data.frame))))))
    
    ## Use the sampled environments and sampled genotypes to create train/test sets
    clus2 %>%
      mutate(train = pmap(list(cv_sample, random, environment), ~{
        .z <- ..3
        map2(.x = ..1, .y = ..2, ~{
          mutate(.x, train = map(train, function(x) subset(pheno_use, environment %in% .y & line_name %in% x[[1]])),
                 test = map(test, function(x) subset(pheno_use, environment == .z & line_name %in% c(x[[1]], vp_geno))))
          
        }) %>%
          map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y))
        
      })) %>%
      select(-random, -cv_sample) %>% 
      unnest(train)
    
  }) %>% ungroup()


## Predict
cv00_cluster_random_predictions <- cv00_cluster_random_train_test %>%
  group_by(set, model, trait, environment, rep, cluster, nTrainEnv) %>%
  do({
    
    df <- .
    
    ## Calculate genotype means for the training sets
    data_adj <- map(df$train, ~geno_means(data = .))
    
    ## Make predictions
    test_val_pred <- map2(.x = data_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]])
    
    ## Add results to the row
    df1 <- df %>% 
      mutate(pred = test_val_pred) %>% 
      unnest(pred) %>% 
      mutate(scheme = ifelse(line_name %in% tp, "cv00", "pocv00")) %>%
      select(environment, .id, scheme, cluster, line_name, value, pred_value)
    
    ## Calculate test accuracy
    test_acc <- df1 %>%
      filter(scheme == "cv00") %>% 
      summarize(cv00 = cor(value, pred_value))
    
    ## Calculate average val accuracy
    val_acc <- df1 %>%
      filter(scheme == "pocv00") %>% 
      group_by(.id) %>%
      summarize(pocv00 = cor(value, pred_value)) %>%
      summarize(pocv00 = mean(pocv00))
    
    cbind(test_acc, val_acc)
    
  }) %>% ungroup() %>%
  rename(val_environment = environment)


## POV00 - predict new genotypes in new environments ##
# For each environment in method, randomly sample n_random sets of the same number of training environments. For each n_random,
# generate nCV cross-validation samples
pov00_cluster_random_train_test <- clusters %>%
  group_by(set, trait, model) %>%
  do({
    row <- .
    
    # Phenotypes to use
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait %in% row$trait)
    
    ## Possible training environments
    train_env <- row$train_env[[1]]
    ## Test environments
    test_env <- row$test_env[[1]]
    
    ## Calculate nTrainEnv depending on the set
    if (row$set == "complete") {
      
      ## Unnest clusters
      clus1 <- unnest(row, cluster) %>%
        group_by(cluster) %>%
        mutate(nTrainEnv = n() - 1) %>%
        ungroup()
      
    } else {
      
      ## Unnest clusters
      clus1 <- unnest(row, cluster) %>%
        group_by(cluster) %>%
        mutate(nTrainEnv = sum(environment %in% train_env)) %>%
        ungroup()
      
    }
    
    clus2 <- clus1 %>%
      filter(environment %in% test_env) %>%
      ## For each environment, sample nTrainEnv environments from the train_env pool
      mutate(random = map2(.x = environment, .y = nTrainEnv, ~replicate(n = n_random, sample(setdiff(train_env, .x), size = .y), simplify = FALSE)))
    
    ## Use the sampled environments and sampled genotypes to create train/test sets
    clus2 %>%
      mutate(data = map2(.x = random, .y = environment, ~tibble(
        train = map(.x, ~subset(pheno_use, environment %in% . & line_name %in% tp_geno)),
        test = list(subset(pheno_use, environment == .y & line_name %in% vp_geno)),
        rep = seq_along(.x)
      ))) %>%
      select(-random) %>% 
      unnest(data)
    
  }) %>% ungroup()

## Predict
pov00_cluster_random_predictions <- pov00_cluster_random_train_test %>%
  group_by(set, model, trait, environment, rep, cluster, nTrainEnv) %>%
  do(pov00 = {
    
    row <- .
    
    ## Calculate genotype means for the training sets
    data_adj <- geno_means(data = row$train[[1]])
    
    ## Make predictions
    gblup(K = K, train = data_adj, test = row$test[[1]], fit.env = FALSE)[[1]]

  }) %>% ungroup() %>%
  rename(val_environment = environment) %>% 
  unnest(pov00) %>%
  gather(scheme, accuracy, pov00)
 


## CV0 and POCV0 - predict observed genotyped in unobserved environment
cv0_cluster_random_train_test <- clusters %>%
  group_by(set, trait, model) %>%
  do({
    row <- .
    
    # Phenotypes to use
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait %in% row$trait)
    
    ## Possible training environments
    train_env <- row$train_env[[1]]
    ## Test environments
    test_env <- row$test_env[[1]]
    
    ## Calculate nTrainEnv depending on the set
    if (row$set == "complete") {
      
      ## Unnest clusters
      clus1 <- unnest(row, cluster) %>%
        group_by(cluster) %>%
        mutate(nTrainEnv = n() - 1) %>%
        ungroup()
      
    } else {
      
      ## Unnest clusters
      clus1 <- unnest(row, cluster) %>%
        group_by(cluster) %>%
        mutate(nTrainEnv = sum(environment %in% train_env)) %>%
        ungroup()
      
    }
    
    clus2 <- clus1 %>%
      filter(environment %in% test_env) %>%
      ## For each environment, sample nTrainEnv environments from the train_env pool
      mutate(random = map2(.x = environment, .y = nTrainEnv, ~replicate(n = n_random, sample(setdiff(train_env, .x), size = .y), simplify = FALSE)))
    
    ## Use the sampled environments and sampled genotypes to create train/test sets
    clus2 %>%
      mutate(data = map2(.x = random, .y = environment, ~tibble(
        train = map(.x, ~subset(pheno_use, environment %in% .)),
        test = list(subset(pheno_use, environment == .y)),
        rep = seq_along(.x)
      ))) %>%
      select(-random) %>% 
      unnest(data)
    
  }) %>% ungroup()


## Predict
cv0_cluster_random_predictions <- cv0_cluster_random_train_test %>%
  group_by(set, model, trait, environment, rep, cluster, nTrainEnv) %>%
  do({
    
    row <- .
    train_data <- row$train[[1]]
    
    ## Calculate genotype means for the training sets
    data_adj <- geno_means(data = train_data)
    
    ## Make predictions
    test_val_pred <- gblup(K = K, train = data_adj, test = row$test[[1]], fit.env = FALSE)[[2]]
    
    test_val_pred %>%
      mutate(scheme = ifelse(line_name %in% tp_geno, "cv0", "pocv0")) %>%
      group_by(scheme) %>%
      summarize(accuracy = cor(value, pred_value))
    
  }) %>% ungroup() %>%
  rename(val_environment = environment)
    
   


# ## CV1 and POCV1 - predict untested genotypes in tested environments
# cv1_pocv1_cluster_random_train_test <- clusters %>%
#   filter(set == "complete") %>%
#   group_by(set, trait, model) %>%
#   do({
#     row <- .
#     
#     # Phenotypes to use
#     pheno_use <- filter(S2_MET_BLUEs_tomodel, trait %in% row$trait)
#     
#     ## Possible training environments
#     train_env <- row$train_env[[1]]
#     ## Test environments
#     test_env <- row$test_env[[1]]
#     
#     ## Unnest clusters
#     clus1 <- unnest(row, cluster) %>%
#       group_by(cluster) %>%
#       mutate(nTrainEnv = n() - 1) %>%
#       ungroup()
#     
#     
#     clus2 <- clus1 %>%
#       ## For each environment, sample nTrainEnv environments from the train_env pool
#       mutate(random = map2(.x = environment, .y = nTrainEnv, ~replicate(n = n_random, c(.x, sample(setdiff(train_env, .x), size = .y)), simplify = FALSE)),
#              cv_sample = map(environment, ~replicate(n = nCV, expr = crossv_kfold(data = cv_train_df, k = k), simplify = FALSE)),
#              cv_sample = map(cv_sample, ~map(., ~mutate_at(., vars(-.id), funs(map(., as.data.frame))))))
#     
#     ## Use the sampled environments and sampled genotypes to create train/test sets
#     clus2 %>%
#       mutate(train = pmap(list(cv_sample, random, environment), ~{
#         .z <- ..3
#         map2(.x = ..1, .y = ..2, ~{
#           mutate(.x, train = map(train, function(x) subset(pheno_use, environment %in% .y & line_name %in% x[[1]])),
#                  test = map(test, function(x) subset(pheno_use, environment == .z & line_name %in% c(x[[1]], vp_geno))))
#           
#         }) %>%
#           map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y))
#         
#       })) %>%
#       select(-random, -cv_sample) %>% 
#       unnest(train)
#     
#     
#   }) %>% ungroup()
# 
# 
# cv1_pocv1_cluster_random_predictions <- cv1_pocv1_cluster_random_train_test %>%
#   group_by(set, model, trait, environment, rep, cluster, nTrainEnv) %>%
#   do({
#     
#     df <- .
#     
#     ## Calculate genotype means for the training sets
#     data_adj <- map(df$train, ~geno_means(data = .))
#     
#     ## Make predictions
#     test_val_pred <- map2(.x = data_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]])
#     
#     ## Add results to the row
#     df1 <- df %>% 
#       mutate(pred = test_val_pred) %>% 
#       unnest(pred) %>% 
#       mutate(scheme = ifelse(line_name %in% tp, "cv00", "pocv00")) %>%
#       select(environment, .id, scheme, cluster, line_name, value, pred_value)
#     
#     ## Calculate test accuracy
#     test_acc <- df1 %>%
#       filter(scheme == "cv1") %>% 
#       summarize(cv00 = cor(value, pred_value))
#     
#     ## Calculate average val accuracy
#     val_acc <- df1 %>%
#       filter(scheme == "pocv1") %>% 
#       group_by(.id) %>%
#       summarize(pocv00 = cor(value, pred_value)) %>%
#       summarize(pocv00 = mean(pocv00))
#     
#     cbind(test_acc, val_acc)
#     
#   }) %>% ungroup() %>%
#   rename(val_environment = environment)
# 
# 
# 
# 
# ## POV1 - predict untested genotypes in tested environments
# pov1_cluster_random_train_test <- clusters %>%
#   filter(set == "complete") %>%
#   group_by(set, trait, model) %>%
#   do({
#     row <- .
#     
#     # Phenotypes to use
#     pheno_use <- filter(S2_MET_BLUEs_tomodel, trait %in% row$trait)
#     
#     ## Possible training environments
#     train_env <- row$train_env[[1]]
#     ## Test environments
#     test_env <- row$test_env[[1]]
#     
#     ## Unnest clusters
#     clus1 <- unnest(row, cluster) %>%
#       group_by(cluster) %>%
#       mutate(nTrainEnv = n() - 1) %>%
#       ungroup()
#     
#     
#     clus2 <- clus1 %>%
#       filter(environment %in% test_env) %>%
#       ## For each environment, sample nTrainEnv environments from the train_env pool
#       mutate(random = map2(.x = environment, .y = nTrainEnv, ~replicate(n = n_random, sample(setdiff(train_env, .x), size = .y), simplify = FALSE)))
#     
#     ## Use the sampled environments and sampled genotypes to create train/test sets
#     clus2 %>%
#       mutate(data = map2(.x = random, .y = environment, ~tibble(
#         train = map(.x, ~subset(pheno_use, environment %in% . & line_name %in% tp_geno)),
#         test = list(subset(pheno_use, environment == .y & line_name %in% vp_geno)),
#         rep = seq_along(.x)
#       ))) %>%
#       select(-random) %>% 
#       unnest(data)
#     
#     
#   }) %>% ungroup()
# 
# 
# pov1_cluster_random_predictions <- pov1_cluster_random_train_test %>% 
#   group_by(set, model, trait, environment, rep, cluster, nTrainEnv) %>%
#   do(pov1 = {
#     
#     row <- .
#     
#     ## Calculate genotype means for the training sets
#     data_adj <- geno_means(data = row$train[[1]])
#     
#     ## Make predictions
#     gblup(K = K, train = data_adj, test = row$test[[1]], fit.env = FALSE)[[1]]
#     
#   }) %>% ungroup() %>%
#   unnest(pov1) %>%
#   rename(val_environment = environment)
# 
# 
# 
# 
# 
# ## CV2 and POCV2 - predict tested genotypes in tested environments
# cv2_cluster_random_train_test <- clusters %>%
#   filter(set == "complete") %>%
#   group_by(set, trait, model) %>%
#   do({
#     row <- .
#     
#     # Phenotypes to use
#     pheno_use <- filter(S2_MET_BLUEs_tomodel, trait %in% row$trait)
#     
#     ## Possible training environments
#     train_env <- row$train_env[[1]]
#     ## Test environments
#     test_env <- row$test_env[[1]]
#     
#     ## Unnest clusters
#     clus1 <- unnest(row, cluster) %>%
#       group_by(cluster) %>%
#       mutate(nTrainEnv = n() - 1) %>%
#       ungroup()
#     
#     
#     clus2 <- clus1 %>%
#       ## For each environment, sample nTrainEnv environments from the train_env pool
#       mutate(random = map2(.x = environment, .y = nTrainEnv, ~replicate(n = n_random, c(.x, sample(setdiff(train_env, .x), size = .y)), simplify = FALSE)),
#              cv_sample = map(environment, ~replicate(n = nCV, expr = crossv_kfold(data = cv_train_df, k = k), simplify = FALSE)),
#              cv_sample = map(cv_sample, ~map(., ~mutate_at(., vars(-.id), funs(map(., as.data.frame))))))
#     
#     ## Use the sampled environments and sampled genotypes to create train/test sets
#     clus2 %>%
#       mutate(train = pmap(list(cv_sample, random, environment), ~{
#         .z <- ..3
#         map2(.x = ..1, .y = ..2, ~{
#           mutate(.x, train = map(train, function(x) subset(pheno_use, environment %in% .y & line_name %in% x[[1]])),
#                  test = map(test, function(x) subset(pheno_use, environment == .z & line_name %in% c(x[[1]], vp_geno))))
#           
#         }) %>%
#           map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y))
#         
#       })) %>%
#       select(-random, -cv_sample) %>% 
#       unnest(train)
#     
#     
#   }) %>% ungroup()
# 
# 
# cv1_pocv1_cluster_random_predictions <- cv1_pocv1_cluster_random_train_test %>%
#   group_by(set, model, trait, environment, rep, cluster, nTrainEnv) %>%
#   do({
#     
#     df <- .
#     
#     ## Calculate genotype means for the training sets
#     data_adj <- map(df$train, ~geno_means(data = .))
#     
#     ## Make predictions
#     test_val_pred <- map2(.x = data_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]])
#     
#     ## Add results to the row
#     df1 <- df %>% 
#       mutate(pred = test_val_pred) %>% 
#       unnest(pred) %>% 
#       mutate(scheme = ifelse(line_name %in% tp, "cv00", "pocv00")) %>%
#       select(environment, .id, scheme, cluster, line_name, value, pred_value)
#     
#     ## Calculate test accuracy
#     test_acc <- df1 %>%
#       filter(scheme == "cv00") %>% 
#       summarize(cv00 = cor(value, pred_value))
#     
#     ## Calculate average val accuracy
#     val_acc <- df1 %>%
#       filter(scheme == "pocv00") %>% 
#       group_by(.id) %>%
#       summarize(pocv00 = cor(value, pred_value)) %>%
#       summarize(pocv00 = mean(pocv00))
#     
#     cbind(test_acc, val_acc)
#     
#   }) %>% ungroup() %>%
#   rename(val_environment = environment)


















pred_list <- cfind(what = "_predictions", class = "data.frame")

# save_file <- file.path(result_dir, "all_predictions.RData")
save_file <- file.path(result_dir, "all_predictions_00.RData")
save(list = pred_list, file = save_file)













