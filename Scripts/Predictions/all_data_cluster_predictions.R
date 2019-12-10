## Random and cluster predictions in the S2MET project
## 
## Author: Jeff Neyhart
## Last Updated: March 12, 2019
## 


# Run the source script
repo_dir <- "/path/to/project/directory/on/supercomputer"

source(file.path(repo_dir, "source_MSI.R"))
library(data.table)
library(modelr)

## Number of cores
n_core <- detectCores()
n_core <- 8




# ## Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# library(data.table)
# library(modelr)



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




## Data.frame of training and test environments per trait
test_env_df <- cluster_df %>%
  group_by(set, trait) %>%
  do({
    df <- .
    test_env1 <- reduce(df$test_env, intersect)
    
    if (unique(df$set) == "complete") {
      train_env1 <- test_env1
    } else {
      train_env1 <- map(df$cluster, "environment") %>% 
        reduce(intersect) %>%
        setdiff(., test_env1)
    }
    
    tibble(test_env = list(test_env1), train_env = list(train_env1))
  }) %>% ungroup()


## Data.frame of the training lines for CV
cv_tp_df <- data.frame(line_name = tp_geno, stringsAsFactors = FALSE)







###### Leave-one-environment-out predictions

loeo_prediction_df <- test_env_df %>%
  mutate(test_env = map2(train_env, test_env, 
                         ~tibble(val_environment = .y, train_env = lapply(X = .y, function(.y1) setdiff(.x, .y1))))) %>%
  unnest(test_env)


# Split
loeo_prediction_df1 <- loeo_prediction_df %>%
  assign_cores(n_core = n_core) %>%
  split(.$core)
  
  


## CV00 - predict new genotypes in new environments
### Leave-one-environment-out cross-validation (CV00)
cv00_predictions <- loeo_prediction_df1 %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
    
    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
    
      ## Filter phenotypes for that trait
      pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
      ve <- row$val_environment
      train_use <- subset(pheno_use, environment != ve)
      test_use <- subset(pheno_use, environment == ve)
      
      
      ## Create CV randomizations
      cv_rand <- replicate(n = nCV, crossv_kfold(data = cv_tp_df, k = k), simplify = FALSE) %>%
        map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y)) %>%
        mutate_at(., vars(train, test), funs(map(., as.data.frame))) %>%
        mutate(train = map(train, ~subset(train_use, line_name %in% .[[1]])),
               test = map(test, ~subset(test_use, line_name %in% c(.[[1]], vp_geno)))) %>%
        nest(train, test, .id)
      
      pred1 <- vector("list", nrow(cv_rand))
      
      # Loop
      for (r in seq(nrow(cv_rand))) {
        
        df <- cv_rand[[2]][[r]]
        
        dat_adj <- map(df$train, ~geno_means(data = .))
        # Predict and validate
        test_val_pred <- map2(.x = dat_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]] %>%
          mutate(scheme = ifelse(line_name %in% tp_geno, "cv00", "pocv00")) ) %>%
          mutate(df, pred = .) %>%
          unnest(pred)
        
        
        cv_acc <- test_val_pred %>%
          subset(scheme == "cv00") %>%
          summarize(cv00 = cor(value, pred_value))
        
        pocv_acc <- test_val_pred %>%
          subset(scheme == "pocv00") %>%
          group_by(.id) %>%
          summarize(pocv00 = cor(value, pred_value)) %>%
          summarize(pocv00 = mean(pocv00))
        
        pred1[[r]] <- cbind(cv_acc, pocv_acc)
        
      }
      
      out[[i]] <- cv_rand %>%
        mutate(pred = pred1) %>%
        unnest(pred) %>%
        select(-data)
      
    }
    
    core_df %>%
      mutate(out = out) %>%
      select(-train_env, -core) %>%
      unnest(out)
    
  }) %>% bind_rows()




## POV00 - predict new genotypes in new environments
### Leave-one-environment-out cross-validation (POV00)
pov00_predictions <- loeo_prediction_df %>%
  group_by(set, trait, val_environment) %>%
  do(pov00 = {
    
    row <- .
    
    ## Filter phenotypes for that trait
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
    ve <- row$val_environment
    train_use <- geno_means(data = subset(pheno_use, environment != ve & line_name %in% tp_geno))
    test_use <- subset(pheno_use, environment == ve & line_name %in% vp_geno)
    
    gblup(K = K, train = train_use, test = test_use, fit.env = FALSE)[[1]]
    
  }) %>% unnest()
  
  


## CV0 - predict observed genotypes in unobserved environments
cv0_predictions <- loeo_prediction_df %>%
  group_by(set, trait, val_environment) %>%
  do({
    
    row <- .
    
    ## Filter phenotypes for that trait
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
    ve <- row$val_environment
    train_use <- geno_means(data = subset(pheno_use, environment != ve))
    test_use <- subset(pheno_use, environment == ve)
    
    gblup(K = K, train = train_use, test = test_use, fit.env = FALSE)[[2]] %>%
      mutate(scheme = ifelse(line_name %in% tp_geno, "cv0", "pocv0")) %>%
      group_by(scheme) %>%
      summarize(accuracy = cor(value, pred_value))
      
    
  }) %>% ungroup()




## List of all predictions
all_pred_list <- c("cv00_predictions", "pov00_predictions", "cv0_predictions")



##### Cluster-based Predictions ######


clusters <- cluster_df %>% 
  select(-test_env) %>%
  inner_join(., test_env_df) %>%
  mutate(nClusters = map_dbl(cluster, ~n_distinct(.$cluster))) %>%
  ## You can't have AMMI and "realistic"
  filter(!(set == "realistic" & model == "AMMI"))

## Unnest for future use
clusters1 <- unnest(clusters, cluster) %>%
  group_by(set, model, trait, nClusters, cluster) %>%
  nest(environment) %>%
  left_join(., select(clusters, set:trait, train_env, test_env))

clusters1_split <- clusters1 %>%
  assign_cores(n_core = n_core, split = TRUE)



## CV00 - predict unobserved genotypes in unobserved environments
cv00_cluster_predictions <- clusters1_split %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {

    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {

          df <- core_df[i,]


# cv00_cluster_predictions <- clusters1 %>%
#   group_by(set, model, trait, nClusters, cluster) %>%
#   do({
#       
#     df <- .
      
      # List of validation environments in the clusters at hand
      ve_use <- intersect(df$test_env[[1]], df$data[[1]][[1]])
      # List of training environments 
      te_use <- intersect(df$train_env[[1]], df$data[[1]][[1]])
      
      ## Filter phenotypes for the trait and for target environments
      pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == df$trait & environment %in% df$data[[1]][[1]])
      
      ve_pred <- vector("list", length(ve_use))
      
      # Iterate over validation environment
      for (r in seq_along(ve_pred)) {
        
        # Subset data
        ve <- ve_use[r]
        ## Training environments
        te <- setdiff(te_use, ve)
          
        train_use <- subset(pheno_use, environment %in% te)
        test_use <- subset(pheno_use, environment == ve)
        
        ## Create CV randomizations for each validation environment
        cv_rand <- replicate(n = nCV, crossv_kfold(data = cv_tp_df, k = k), simplify = FALSE) %>%
          map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y)) %>%
          mutate_at(., vars(train, test), funs(map(., as.data.frame))) %>%
          mutate(train = map(train, ~subset(train_use, line_name %in% .[[1]])),
                 test = map(test, ~subset(test_use, line_name %in% c(.[[1]], vp_geno)))) %>%
          nest(train, test, .id)
        
        
        ## Run adjustments, predictions, and validation
        pred1 <- lapply(cv_rand$data, function(df) {
          dat_adj <- map(df$train, ~geno_means(data = .))
          
          test_val_pred <- map2(.x = dat_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]] %>%
                                  mutate(scheme = ifelse(line_name %in% tp_geno, "cv00", "pocv00")) ) %>%
            mutate(df, pred = .) %>%
            unnest(pred)
          
          cv_acc <- test_val_pred %>%
            subset(scheme == "cv00") %>%
            summarize(cv00 = cor(value, pred_value))
          
          pocv_acc <- test_val_pred %>%
            subset(scheme == "pocv00") %>%
            group_by(.id) %>%
            summarize(pocv00 = cor(value, pred_value)) %>%
            summarize(pocv00 = mean(pocv00))
          
          cbind(cv_acc, pocv_acc)
          
        })
        
        ## Add the accuracies to the list
        ve_pred[[r]] <- cv_rand %>%
          mutate(pred = pred1, nTrainEnv = n_distinct(train_use$environment)) %>%
          unnest(pred) %>%
          select(-data)
        
      }
      
  #     df$data[[1]] %>% 
  #       mutate(predictions = ve_pred) %>% 
  #       rename(val_environment = environment)
  #     
  # })
  #     
      
      
      
      ## Add the cluster results to the out vector
      out[[i]] <- data_frame(environment = ve_use) %>%
        mutate(predictions = ve_pred) %>%
        rename(val_environment = environment) %>%
        unnest()

    }

    # Add the out vector to the core_df
    core_df %>%
      mutate(out = out) %>%
      select(-data, -train_env, -test_env, -core) %>%
      unnest(out)

  }) %>% bind_rows()
      
      
## POV00 - predict unobserved genotypes in unobserved environments
pov00_cluster_predictions <- clusters %>%
  group_by(set, trait, model, nClusters) %>%
  do({
    
    row <- .
    pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
    
    clus <- unnest(row, cluster)
  
    row1 <- row %>% 
      unnest(test_env) %>% 
      mutate(train_env = row$train_env, 
             train_env = map2(train_env, test_env, setdiff)) %>%
      ## Get the data for the cluster of the test environment
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
    
    ## Calculate genotype means
    data_adj <-  map(row1$train, ~geno_means(data = .))
    
    ## Run the base predictions
    preds <- map2_dbl(.x = data_adj, .y = row1$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[1]])
    
    ## Add results and output
    row1 %>% 
      mutate(pov00 = preds) %>% 
      select(set, model, trait, val_environment = test_env, nClusters, cluster, nTrainEnv, pov00)
    
  }) %>% ungroup()
  


## CV0 and POV0 - predict observed genotypes in unobserved environments
cv0_cluster_predictions <- clusters %>%
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
      filter(test_env %in% clus$environment) %>%
      ## Get the data for the cluster of the test environment
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
    
    
    ## Calculate genotype means
    data_adj <-  map(row1$train, ~geno_means(data = .))
    
    ## Run the base predictions
    test_val_pred <- map2(.x = data_adj, .y = row1$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]]) %>%
      map(~mutate(., scheme = ifelse(line_name %in% tp, "cv0", "pocv0")) %>%
            group_by(scheme) %>%
            summarize(accuracy = cor(value, pred_value)) )
    
    ## Add results and output
    row1 %>% 
      mutate(out = test_val_pred) %>% 
      unnest(out) %>%
      rename(val_environment = test_env)
    
  }) %>% ungroup()


cluster_pred_list <- c("cv00_cluster_predictions", "pov00_cluster_predictions", "cv0_cluster_predictions")






##### Random cluster predictions #####

## Number of random cluster assignments
n_random <- nCV


## Training lines for CV
cv_train_df <- data.frame(line_name = tp_geno, stringsAsFactors = FALSE)

## CV00 - predict new genotypes in new environments ##
# For each environment in method, randomly sample n_random sets of the same number of training environments. For each n_random,
# generate nCV cross-validation samples

## CV00 - predict unobserved genotypes in unobserved environments
cv00_cluster_random_predictions <- clusters1_split %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
    
    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {
      
      df <- core_df[i,]
      
      # List of validation environments in the clusters at hand
      ve_use <- intersect(df$test_env[[1]], df$data[[1]][[1]])
      # List of training environments 
      te_use <- intersect(df$train_env[[1]], df$data[[1]][[1]])

      ## Filter phenotypes for the trait and for target environments
      pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == df$trait)
      
      ve_pred <- vector("list", length(ve_use))
      
      # Iterate over validation environment
      for (r in seq_along(ve_pred)) {
        
        # Subset data
        ve <- ve_use[r]
        
        ## Possible training environments
        te <- setdiff(te_use, ve)
        # Number of training environments
        nTrainEnv <- length(te)
        test_use <- subset(pheno_use, environment == ve)
        
        ## Create CV randomizations for each validation environment
        cv_rand <- replicate(n = nCV, crossv_kfold(data = cv_tp_df, k = k), simplify = FALSE) %>%
          map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y, train_env = list(sample(df$train_env[[1]], nTrainEnv)))) %>%
          mutate_at(., vars(train, test), funs(map(., as.data.frame))) %>%
          mutate(train = map2(train, train_env, ~subset(pheno_use, line_name %in% .x[[1]] & environment %in% .y)),
                 test = map(test, ~subset(test_use, line_name %in% c(.[[1]], vp_geno)))) %>%
          nest(train, test, .id, train_env)

        ## Run adjustments, predictions, and validation
        pred1 <- lapply(cv_rand$data, function(df) {
          dat_adj <- map(df$train, ~geno_means(data = .))
          
          test_val_pred <- map2(.x = dat_adj, .y = df$test, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]] %>%
                                  mutate(scheme = ifelse(line_name %in% tp_geno, "cv00", "pocv00")) ) %>%
            mutate(df, pred = .) %>%
            unnest(pred)
          
          cv_acc <- test_val_pred %>%
            subset(scheme == "cv00") %>%
            summarize(cv00 = cor(value, pred_value))
          
          pocv_acc <- test_val_pred %>%
            subset(scheme == "pocv00") %>%
            group_by(.id) %>%
            summarize(pocv00 = cor(value, pred_value)) %>%
            summarize(pocv00 = mean(pocv00))
          
          cbind(cv_acc, pocv_acc)
          
        })
        
        ## Add the accuracies to the list
        ve_pred[[r]] <- cv_rand %>%
          mutate(pred = pred1, nTrainEnv = map_dbl(data, ~length(.$train_env[[1]]))) %>%
          unnest(pred) %>%
          select(-data)
        
      }
      
      #     df$data[[1]] %>% 
      #       mutate(predictions = ve_pred) %>% 
      #       rename(val_environment = environment)
      #     
      # })
      #     
      
      
      
      ## Add the cluster results to the out vector
      out[[i]] <- data_frame(environment = ve_use) %>%
        mutate(predictions = ve_pred) %>%
        rename(val_environment = environment) %>%
        unnest()
      
      
    }
    
    # Add the out vector to the core_df
    core_df %>%
      mutate(out = out) %>%
      select(-data, -train_env, -test_env, -core) %>%
      unnest(out)
    
  }) %>% bind_rows()


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
    


cluster_random_pred_list <- c("cv00_cluster_random_predictions", "pov00_cluster_random_predictions", "cv0_cluster_random_predictions")



pred_list <- c(all_pred_list, cluster_pred_list, cluster_random_pred_list)

# save_file <- file.path(result_dir, "all_predictions.RData")
save_file <- file.path(result_dir, "all_data_cluster_predictions.RData")
save(list = pred_list, file = save_file)













