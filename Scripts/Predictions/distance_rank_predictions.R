## Add-one-environment predictions in the S2MET prediction project
## This script is meant to run on MSI
## 
## Author: Jeff Neyhart
## Last Updated: March 18, 2019
## 


# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions/"
source(file.path(repo_dir, "source_MSI.R"))
library(data.table)


## Number of cores
n_core <- 32
n_core <- detectCores()
n_core <- 4



# ## Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# library(modelr)
# library(data.table)



# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))


## Prepare the BLUEs for modeling
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate(line_name = as.factor(line_name))


# Number of CV iterations
nCV <- 10
# Number of random environment rankings
n_random <- nCV
# CV folds
k <- 5




## Data.frame of test environments
test_env_df <- bind_rows(
  data_frame(set = "complete", trait = names(complete_train_env), train_env = complete_train_env, test_env = complete_train_env),
  data_frame(set = "realistic", trait = names(complete_train_env), train_env = realistic_train_env, test_env = realistic_test_env)
)


## Data.frame of the training lines for CV
cv_tp_df <- data.frame(line_name = tp_geno, stringsAsFactors = FALSE)



##### Distance Rank Prediction #####

# Use the clusters based on the TP-only data
environment_rank_df <- pred_env_dist_rank %>%
  rename(val_environment = validation_environment) %>%
  filter(!mat_set %in% c("Jarquin", "MalosettiStand")) %>%
  filter(model %in% names(dist_method_abbr_use)) %>%
  select(-mat_set) %>%
  mutate(data = list(NULL))


# 
# # ## Create training and test sets
# # pov00_environment_rank_train_test <- environment_rank_df %>%
# #   group_by(set, trait, model, val_environment) %>%
# #   do({
# # 
# #     row <- .
# 
# pov00_environment_rank_train_test <- environment_rank_df %>%
#   assign_cores(n_core) %>%
#   split(.$core) %>%
#   mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
# 
#     out <- vector("list", nrow(core_df))
#     for (i in seq_along(out)) {
# 
#       row <- core_df[i,]
#       
#     
#       ## Filter phenotypes for that trait
#       pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
#       
#       ## Create the ordered list of training environments
#       # Then create training sets
#       # Calculate genotypic means (this will hopefully save space later)
#       train <- accumulate(row$rank[[1]], c) %>%
#         map(~filter(pheno_use, environment %in% ., line_name %in% tp_geno)) %>%
#         map(~group_by(., line_name) %>% summarize(value = mean(value)) %>% ungroup() %>% mutate(environment = "blue", std_error = 0))
#       
#       # # Test sets
#       # tibble(
#       #   train = train,
#       #   test = list(filter(pheno_use, environment == row$val_environment, line_name %in% vp_geno)),
#       #   nTrainEnv = seq_along(train)
#       # )
#       
#       # Test sets
#       out[[i]] <- tibble(
#         train = train,
#         test = list(filter(pheno_use, environment == row$val_environment, line_name %in% vp_geno)),
#         nTrainEnv = seq_along(train)
#       )
#       
#     }
#     core_df %>%
#       mutate(data = out) %>%
#       select(-core) %>%
#       unnest(data)
#     
#     
#   # }) %>% ungroup()
#     
#   }) %>% bind_rows()
# 
# 
# 
# 
# # ## Run the predictions
# # pov00_environment_rank_predictions <- pov00_environment_rank_train_test %>%
# #   group_by(set, model, trait, val_environment, nTrainEnv) %>%
# #   do(pov00 = {
# # 
# #       row <- .
# 
# ## Run the predictions
# pov00_environment_rank_predictions <- pov00_environment_rank_train_test %>%
#   assign_cores(n_core) %>%
#   split(.$core) %>%
#   mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
# 
#     out <- vector("list", nrow(core_df))
#     for (i in seq_along(out)) {
# 
#       row <- core_df[i,]
#       
#       
#       # ## Run the base predictions and return accuracy
#       # gblup(K = K, train = row$train[[1]], test = row$test[[1]], fit.env = FALSE)[[1]]
#       
#       ## Run the base predictions and return accuracy
#       out[[i]] <- gblup(K = K, train = row$train[[1]], test = row$test[[1]], fit.env = FALSE)[[1]]
# 
#     }
# 
#     core_df %>%
#       mutate(pov00 = out) %>%
#       select(-core) %>%
#       unnest(pov00) %>%
#       gather(scheme, accuracy, pov00)
# 
#     
#   # }) %>% ungroup() %>%
#   # unnest(pov00) %>%
#   # gather(scheme, accuracy, pov00)
#     
#   }) %>% bind_rows()
# 
# 
# 
# # ## Random environment rankings
# # pov00_environment_rank_random_df <- environment_rank_df %>%
# #   filter(model == "great_circle_dist") %>%
# #   group_by(set, trait, val_environment) %>%
# #   do({
# # 
# #     row <- .
# 
# 
# ## Random environment rankings
# pov00_environment_rank_random_df <- environment_rank_df %>%
#   filter(model == "great_circle_dist") %>%
#   assign_cores(n_core) %>%
#   split(.$core) %>%
#   mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
# 
#     out <- vector("list", nrow(core_df))
#     for (i in seq_along(out)) {
# 
#       row <- core_df[i,]
#     
#       ## Filter phenotypes for that trait
#       pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
#       
#       ## Create random orders of environments
#       ## Then create training sets
#       ## Then adjust
#       train <- replicate(n = n_random, sample(row$rank[[1]]), simplify = FALSE) %>%
#         map(~accumulate(., c) %>%
#               map(~filter(pheno_use, environment %in% ., line_name %in% tp_geno)) %>%
#               map(~group_by(., line_name) %>% summarize(value = mean(value)) %>% ungroup() %>% mutate(environment = "blue", std_error = 0)) )
#       
#       # train %>% 
#       #   map(~tibble(train = ., nTrainEnv = seq_along(.))) %>% 
#       #   map2_df(.x = ., .y = seq_along(.), ~mutate(.x, model = paste0("sample", .y))) %>%
#       #   mutate(test = list(filter(pheno_use, environment == row$val_environment, line_name %in% vp_geno)))
#       
#       out[[i]] <- train %>%
#         map(~tibble(train = ., nTrainEnv = seq_along(.))) %>% 
#         map2_df(.x = ., .y = seq_along(.), ~mutate(.x, model = paste0("sample", .y))) %>%
#         mutate(test = list(filter(pheno_use, environment == row$val_environment, line_name %in% vp_geno)))
# 
#     }
# 
#     core_df %>%
#       mutate(data = out) %>%
#       select(-core) %>%
#       unnest(data) %>%
#       select(-model) %>% # Rename model1 to avoid duplicates
#       rename(model = model1)
#     
#   # }) %>% ungroup()
#     
#   }) %>% bind_rows()
# 
# 
# 
# 
# # ## Run the predictions
# # pov00_environment_rank_random_predictions <- pov00_environment_rank_random_df %>%
# #   group_by(set, model, trait, val_environment, nTrainEnv) %>%
# #   do(pov00 = {
# # 
# #     row <- .
#     
# ## Run the predictions
# pov00_environment_rank_random_predictions <- pov00_environment_rank_random_df %>%
#   assign_cores(n_core) %>%
#   split(.$core) %>%
#   mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
# 
#     out <- vector("list", nrow(core_df))
#     for (i in seq_along(out)) {
# 
#       row <- core_df[i,]
#     
#     # gblup(K = K, train = row$train[[1]], test = row$test[[1]], fit.env = FALSE)[[1]]
#     
#       ## Run the base predictions and return accuracy
#       out[[i]] <- gblup(K = K, train = row$train[[1]], test = row$test[[1]], fit.env = FALSE)[[1]]
# 
#     }
# 
#     core_df %>%
#       mutate(pov00 = out) %>%
#       select(-core, -train, -test) %>%
#       unnest(pov00) %>%
#       gather(scheme, accuracy, pov00)
#     
#     
#   # }) %>% ungroup() %>%
#   # unnest(pov00) %>%
#   # gather(scheme, accuracy, pov00)
#     
#   }) %>% bind_rows()
# 
# 
# 



# ## Cross-validation (CV00) with adding environments
# ## Create training and test sets
# cv00_environment_rank_train_test <- environment_rank_df %>%
#   # slice(1:5) %>%
#   group_by(set, trait, model, val_environment) %>%
#   do({
# 
#     row <- .


## Cross-validation (CV00) with adding environments
## Create training and test sets
cv00_environment_rank_train_test <- environment_rank_df %>%
  assign_cores(n_core) %>%
  split(.$core) %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {

    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {

      row <- core_df[i,]
    
      ## Filter phenotypes for that trait
      pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
      ve <- row$val_environment
      
      ## Create CV randomizations
      cv_rand <- replicate(n = nCV, crossv_kfold(data = cv_tp_df, k = k), simplify = FALSE) %>%
        map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y)) %>%
        mutate_at(., vars(train, test), funs(map(., as.data.frame))) %>%
        crossing(., tibble(train_env = accumulate(row$rank[[1]], c)))
      
      # # Then create training sets
      # cv_rand %>%
      #   mutate(train = map2(train, train_env, ~filter(pheno_use, environment %in% .y, line_name %in% .x[[1]])),
      #          # train1 = map(train, ~geno_means(data = .)),
      #          train = map(train, ~group_by(., line_name) %>% summarize(value = mean(value)) %>% ungroup() %>% mutate(environment = "blue", std_error = 0)),
      #          test = map(test, ~filter(pheno_use, environment %in% row$val_environment, line_name %in% c(.[[1]], vp_geno))),
      #          nTrainEnv = map_dbl(train_env, length))
 
      
      # Then create training sets
      cv_rand1 <- mutate_at(cv_rand, vars(train, test), ~list(NULL)) %>%
        mutate(nTrainEnv = 0)
      train1 <- test1 <- vector("list", nrow(cv_rand1))
      nTrainEnv1 <- numeric(length = nrow(cv_rand1))
      
      # Loop
      for (r in seq(nrow(cv_rand))) {
        train_r <- subset(pheno_use, environment %in% cv_rand[[5]][[r]] & line_name %in% cv_rand[[1]][[r]][[1]])
        train_r1 <- as.data.frame(as.data.table(train_r)[, lapply(.SD, mean), by = line_name, .SDcols = "value"])
        train1[[r]] <- mutate(train_r1, environment = "blue", std_error = 0)
        
        test1[[r]] <- subset(pheno_use, environment %in% ve & line_name %in% c(cv_rand[[2]][[r]][[1]], vp_geno))
        nTrainEnv1[r] <- length(cv_rand[[5]][[r]])
        
      }
      
      # mutate(cv_rand1, train = train1, test = test1, nTrainEnv = nTrainEnv1)
      
      
  out[[i]] <- mutate(cv_rand1, train = train1, test = test1, nTrainEnv = nTrainEnv1)
  }
  core_df %>%
    mutate(data = out) %>%
    select(-core) %>%
    unnest(data)

  }) %>% bind_rows()
    
      
  # }) %>% ungroup()
  


# # CV00 predictions
# cv00_environment_rank_predictions <- cv00_environment_rank_train_test %>%
#   # Must nest .ids
#   group_by(trait, set, val_environment, model, rep, nTrainEnv) %>% 
#   do({
#     
#     df <- .


# CV00 predictions
cv00_environment_rank_predictions <- cv00_environment_rank_train_test %>%
  # Must nest .ids
  group_by(trait, set, val_environment, model, rep, nTrainEnv) %>%
  nest() %>%
  assign_cores(n_core) %>%
  split(.$core) %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {

    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {

      datai <- core_df$data[[i]]
      df <- core_df[i,]

      
      ## Run predictions
      test_val_pred <- map2(.x = datai[[1]], .y = datai[[2]], ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]]) %>%
        map2_df(.x = ., .y = datai$.id, ~mutate(.x, .id = .y)) %>%
        mutate(scheme = ifelse(line_name %in% tp_geno, "cv00", "pocv00"))
      
      cv_acc <- test_val_pred %>%
        filter(scheme == "cv00") %>%
        summarize(cv00 = cor(value, pred_value))
      
      pocv_acc <- test_val_pred %>%
        filter(scheme == "pocv00") %>%
        group_by(.id) %>%
        summarize(pocv00 = cor(value, pred_value)) %>%
        summarize(pocv00 = mean(pocv00))
      
      # cbind(cv_acc, pocv_acc)
      
      out[[i]] <- cbind(cv_acc, pocv_acc)
      
    }

    core_df %>%
      mutate(data = out) %>%
      select(-core) %>%
      unnest(data) %>%
      gather(scheme, accuracy, cv00, pocv00)
    
    
    # }) %>% ungroup() %>%
    # gather(scheme, accuracy, cv00, pov00)
    
  })



# ## Cross-validation (CV00) with adding environments
# ## Create training and test sets
# cv00_environment_rank_random_df <- environment_rank_df %>%
#   filter(model == "great_circle_dist") %>%
#   group_by(set, trait, val_environment) %>%
#   do({
# 
#     row <- .


## Random environment rankings
cv00_environment_rank_random_df <- environment_rank_df %>%
  filter(model == "great_circle_dist") %>%
  assign_cores(n_core) %>%
  split(.$core) %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {

    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {

      row <- core_df[i,]
      
      ## Filter phenotypes for that trait
      pheno_use <- subset(S2_MET_BLUEs_tomodel, trait == row$trait)
      
      ## Create random orders of environments
      ## Then create training sets
      ## Then adjust
      train_env <- replicate(n = n_random, sample(row$rank[[1]]), simplify = FALSE) %>%
        map(~accumulate(., c)) %>%
        map2_df(., seq_along(.), ~tibble(train_env = .x, rep = .y))
      
      
      ## Create CV randomizations for each train_env randomization
      cv_rand <- replicate(n = nCV, crossv_kfold(data = cv_tp_df, k = k), simplify = FALSE) %>%
        map2_df(.x = ., .y = seq_along(.), ~mutate(.x, rep = .y)) %>%
        mutate_at(., vars(train, test), funs(map(., as.data.frame))) %>%
        full_join(., train_env, by = "rep") %>%
        mutate(model = paste0("sample", rep))
      
      # Then create training sets
      cv_rand1 <- mutate_at(cv_rand, vars(train, test), ~list(NULL)) %>%
        mutate(nTrainEnv = 0)
      train1 <- test1 <- vector("list", nrow(cv_rand1))
      nTrainEnv1 <- numeric(length = nrow(cv_rand1))
      
      # Loop
      for (r in seq(nrow(cv_rand))) {
        train_r <- subset(pheno_use, environment %in% cv_rand[[5]][[r]] & line_name %in% cv_rand[[1]][[r]][[1]])
        train_r1 <- as.data.frame(as.data.table(train_r)[, lapply(.SD, mean), by = line_name, .SDcols = "value"])
        train1[[r]] <- mutate(train_r1, environment = "blue", std_error = 0)
        
        test1[[r]] <- subset(pheno_use, environment %in% ve & line_name %in% c(cv_rand[[2]][[r]][[1]], vp_geno))
        nTrainEnv1[r] <- length(cv_rand[[5]][[r]])
        
      }
      
      # mutate(cv_rand1, train = train1, test = test1, nTrainEnv = nTrainEnv1)
      

      out[[i]] <- mutate(cv_rand1, train = train1, test = test1, nTrainEnv = nTrainEnv1)
    }

    core_df %>%
      mutate(data = out) %>%
      select(-core) %>%
      unnest(data) %>%
      select(-model) %>% # Rename model1 to avoid duplicates
      rename(model = model1)

   # }) %>% ungroup()
    
  }) %>% bind_rows()


## Predictions
cv00_environment_rank_random_predictions <- cv00_environment_rank_random_df %>%
  # Must nest .ids
  group_by(trait, set, val_environment, model, rep, nTrainEnv) %>% 
  nest() %>%
  assign_cores(n_core) %>% 
  split(.$core) %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
    
    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {
      
      datai <- core_df$data[[i]]
      df <- core_df[i,]
      
      
      ## Run predictions
      test_val_pred <- map2(.x = datai[[1]], .y = datai[[2]], ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]]) %>%
        map2_df(.x = ., .y = datai$.id, ~mutate(.x, .id = .y)) %>%
        mutate(scheme = ifelse(line_name %in% tp_geno, "cv00", "pocv00"))
      
      cv_acc <- test_val_pred %>%
        filter(scheme == "cv00") %>%
        summarize(cv00 = cor(value, pred_value))
      
      pocv_acc <- test_val_pred %>%
        filter(scheme == "pocv00") %>%
        group_by(.id) %>%
        summarize(pocv00 = cor(value, pred_value)) %>%
        summarize(pocv00 = mean(pocv00))
      
      # cbind(cv_acc, pocv_acc)
      
      out[[i]] <- cbind(cv_acc, pocv_acc)
      
    }
    
    core_df %>%
      mutate(data = out) %>%
      select(-core) %>%
      unnest(data) %>%
      gather(scheme, accuracy, cv00, pocv00)
    
    
    # }) %>% ungroup() %>%
    # unnest(pov00) %>%
    # gather(scheme, accuracy, pov00)
    
  })




pred_list <- ls(pattern = "_predictions")

# Save
save_file <- file.path(result_dir, "distance_rank_predictions.RData")
save(list = pred_list, file = save_file)



