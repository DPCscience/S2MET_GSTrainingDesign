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
n_core <- 8



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
nCV <- 25
# Number of random environment rankings
n_random <- nCV
# CV folds
k <- 5




## Data.frame of test environments
test_env_df <- bind_rows(
  data_frame(set = "complete", trait = names(complete_train_env), train_env = complete_train_env, test_env = complete_train_env),
  data_frame(set = "realistic2017", trait = names(complete_train_env), train_env = realistic_train_env, test_env = realistic_test_env)
)


## Data.frame of the training lines for CV
cv_tp_df <- data.frame(line_name = tp_geno, stringsAsFactors = FALSE)



##### Distance Rank Prediction #####

# Use the clusters based on the TP-only data
environment_rank_df <- pred_env_dist_rank %>%
  rename(val_environment = validation_environment) %>%
  filter(!mat_set %in% c("Jarquin", "MalosettiStand")) %>%
  filter(model %in% names(dist_method_abbr_use)) %>%
  filter(set %in% c("complete", "realistic2017")) %>%
  select(-mat_set) %>%
  mutate(data = list(NULL))

## Assign cores and split
environment_rank_df1 <- environment_rank_df %>%
  assign_cores(n_core) %>%
  split(.$core)



pov00_environment_rank_predictions <- environment_rank_df1 %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
  # lapply(X = ., FUN = function(core_df) {
    
    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {

      row <- core_df[i,]

      ## Filter phenotypes for that trait
      pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
      ve <- row$val_environment
      # Rank of environments
      eRank <- row$rank[[1]]
      
      ## Create the ordered list of training environments
      # Then create training sets
      # Calculate genotypic means (this will hopefully save space later)
      train <- accumulate(eRank, c) %>%
        map(~subset(pheno_use, environment %in% . & line_name %in% tp_geno)) %>%
        map(~as.data.frame(as.data.table(.)[, lapply(.SD, mean), by = line_name, .SDcols = "value"])) %>%
        map(~mutate(., environment = "blue", std_error = 0))
      
      # Test set
      test <- subset(pheno_use, environment == ve & line_name %in% vp_geno)
      
      # Empty vector to store prediction accuracy estimates
      pred1 <- vector("numeric", length(train))
      
      # Loop
      for (r in seq_along(pred1)) {

        ## Predict
        pred1[[r]] <- gblup(K = K, train = train[[r]], test = test, fit.env = FALSE)[[1]]

      }
      
      ## output data.frame 
      
      out[[i]] <- data.frame(nTrainEnv = seq_along(eRank), scheme = "pov00", accuracy = pred1)
      
    }
    
    
    core_df %>%
      mutate(data = out) %>%
      select(-core) %>%
      unnest(data)

  })




## Random environment rankings
pov00_environment_rank_random_predictons <- environment_rank_df %>%
  filter(model == "great_circle_dist") %>%
  assign_cores(n_core) %>%
  split(.$core) %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
  # lapply(X = ., FUN = function(core_df) {

    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {

      row <- core_df[i,]

      ## Filter phenotypes for that trait
      pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
      ve <- row$val_environment
      # Rank of environments
      eRank <- row$rank[[1]]
      
      ## Create the ordered list of training environments
      # Then create training sets
      # Calculate genotypic means (this will hopefully save space later)
      train <- replicate(n = n_random, sample(eRank), simplify = FALSE) %>%
        map(~accumulate(., c) %>%
              map(~subset(pheno_use, environment %in% . & line_name %in% tp_geno)) %>%
              map(~as.data.frame(as.data.table(.)[, lapply(.SD, mean), by = line_name, .SDcols = "value"])) %>%
              map(~mutate(., environment = "blue", std_error = 0)) )
      
      # Test set
      test <- subset(pheno_use, environment == ve & line_name %in% vp_geno)
      
      # Empty vector to store prediction accuracy estimates
      pred1 <- vector("list", length(train))
      
      # Loop
      for (r in seq_along(pred1)) {
        train_r <- train[[r]]
        
        ## Predict
        pred_r <- map_dbl(train_r, ~gblup(K = K, train = ., test = test, fit.env = FALSE)[[1]])
        # Data.frame
        pred1[[r]] <- data.frame(nTrainEnv = seq_along(pred_r), rep = r, scheme = "pov00", accuracy = pred_r)
        
      }
      
      ## output data.frame 
      
      out[[i]] <- do.call("rbind", pred1)
      
    }
    
    
    core_df %>%
      mutate(data = out) %>%
      select(-core) %>%
      unnest(data)

  })









## Cross-validation (CV00) with adding environments
## Create training and test sets
## Combine test/train set creation and prediction
cv00_environment_rank_predictions <- environment_rank_df1 %>%
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
      cv_rand <- cv_rand %>% 
        mutate(nTrainEnv = map_dbl(train_env, length)) %>% 
        group_by(rep, nTrainEnv) %>% 
        nest() %>%
        mutate(pred = list(NULL))
      # train1 <- test1 <- vector("list", nrow(cv_rand1))
      pred1 <- vector("list", nrow(cv_rand))
      
      # Loop
      for (r in seq(nrow(cv_rand))) {
        dat <- cv_rand[[3]][[r]]
        
        train_use <- lapply(dat[[1]], "[[", "line_name") %>% 
          map(~subset(pheno_use, environment %in% dat[[4]][[1]] & line_name %in% .)) %>%
          map(~as.data.frame(as.data.table(.)[, lapply(.SD, mean), by = line_name, .SDcols = "value"])) %>%
          map(~mutate(., environment = "blue", std_error = 0))
        
        test_use <- lapply(X = dat[[2]], "[[", "line_name") %>%
          map(~subset(pheno_use, environment %in% ve & line_name %in% c(., vp_geno)))
        
        ## Predict
        test_val_pred <- map2(.x = train_use, .y = test_use, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]]) %>%
          map2_df(.x = ., .y = dat[[3]], ~mutate(.x, .id = .y)) %>%
          mutate(scheme = ifelse(line_name %in% tp_geno, "cv00", "pocv00"))
        
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
        select(-data) %>% 
        unnest(pred)
      
    }
    
    
  core_df %>%
    mutate(data = out) %>%
    select(-core) %>%
    unnest(data)

  })
    
      
  # }) %>% ungroup()
  



## Random environment rankings
cv00_environment_rank_random_predictions <- environment_rank_df %>%
  filter(model == "great_circle_dist") %>%
  assign_cores(n_core) %>%
  split(.$core) %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {

    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {

      row <- core_df[i,]
      
      ## Filter phenotypes for that trait
      pheno_use <- subset(S2_MET_BLUEs_tomodel, trait == row$trait)
      ve <- row$val_environment
      
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
      cv_rand <- cv_rand %>% 
        mutate(nTrainEnv = map_dbl(train_env, length)) %>% 
        group_by(model, rep, nTrainEnv) %>% 
        nest() %>%
        mutate(pred = list(NULL))
      
      # train1 <- test1 <- vector("list", nrow(cv_rand1))
      pred1 <- vector("list", nrow(cv_rand))
      
      # Loop
      for (r in seq(nrow(cv_rand))) {
        dat <- cv_rand[[4]][[r]]
        
        train_use <- lapply(dat[[1]], "[[", "line_name") %>% 
          map(~subset(pheno_use, environment %in% dat[[4]][[1]] & line_name %in% .)) %>%
          map(~as.data.frame(as.data.table(.)[, lapply(.SD, mean), by = line_name, .SDcols = "value"])) %>%
          map(~mutate(., environment = "blue", std_error = 0))
        
        test_use <- lapply(X = dat[[2]], "[[", "line_name") %>%
          map(~subset(pheno_use, environment %in% ve & line_name %in% c(., vp_geno)))
        
        ## Predict
        test_val_pred <- map2(.x = train_use, .y = test_use, ~gblup(K = K, train = .x, test = .y, fit.env = FALSE)[[2]]) %>%
          map2_df(.x = ., .y = dat[[3]], ~mutate(.x, .id = .y)) %>%
          mutate(scheme = ifelse(line_name %in% tp_geno, "cv00", "pocv00"))
        
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
        select(-data) %>% 
        unnest(pred)
      
    }

    core_df %>%
      mutate(data = out) %>%
      select(-core) %>%
      unnest(data) %>%
      select(-model) %>% # Rename model1 to avoid duplicates
      rename(model = model1)

   # }) %>% ungroup()
    
  })









## Cross-validation (CV0) with adding environments
## Use all lines to predict in new environments
## Create training and test sets
## Combine test/train set creation and prediction
cv0_pocv0_environment_rank_predictions <- environment_rank_df1 %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
  # lapply(X = ., FUN = function(core_df) {
      
    
    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      
      ## Filter phenotypes for that trait
      pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
      ve <- row$val_environment
      # Rank of environments
      eRank <- row$rank[[1]]
      
      ## Create the ordered list of training environments
      # Then create training sets
      # Calculate genotypic means (this will hopefully save space later)
      train <- accumulate(eRank, c) %>%
        map(~subset(pheno_use, environment %in% . & line_name %in% c(tp_geno, vp_geno))) %>%
        map(~as.data.frame(as.data.table(.)[, lapply(.SD, mean), by = line_name, .SDcols = "value"])) %>%
        map(~mutate(., environment = "blue", std_error = 0))


      
      # Test set
      test <- subset(pheno_use, environment == ve & line_name %in% c(tp_geno, vp_geno))
      
      # Empty vector to store prediction accuracy estimates
      test_val_pred <- map(train, ~gblup(K = K, train = ., test = test, fit.env = FALSE)[[2]]) %>%
        map2_df(., seq_along(.), ~mutate(.x, nTrainEnv = .y, scheme = ifelse(line_name %in% tp_geno, "cv0", "pocv0")))
      
    
      ## output data.frame 
      out[[i]] <- test_val_pred %>% 
        group_by(nTrainEnv, scheme) %>% 
        summarize(accuracy = cor(value, pred_value)) %>%
        ungroup()
      
    }
    
    
    core_df %>%
      mutate(data = out) %>%
      select(-core) %>%
      unnest(data)
    
  })




## Random rankings
cv0_pocv0_environment_rank_random_predictions <- environment_rank_df %>%
  filter(model == "great_circle_dist") %>%
  assign_cores(n_core) %>%
  split(.$core) %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
  # lapply(X = ., FUN = function(core_df) {
    
    
    out <- vector("list", nrow(core_df))
    for (i in seq_along(out)) {
      
      row <- core_df[i,]
      
      ## Filter phenotypes for that trait
      pheno_use <- filter(S2_MET_BLUEs_tomodel, trait == row$trait)
      ve <- row$val_environment
      # Rank of environments
      eRank <- row$rank[[1]]
      
      ## Create the ordered list of training environments
      # Then create training sets
      # Calculate genotypic means (this will hopefully save space later)
      train <- replicate(n = n_random, sample(eRank), simplify = FALSE) %>%
        map(~accumulate(., c) %>%
              map(~subset(pheno_use, environment %in% . & line_name %in% c(tp_geno, vp_geno))) %>%
              map(~as.data.frame(as.data.table(.)[, lapply(.SD, mean), by = line_name, .SDcols = "value"])) %>%
              map(~mutate(., environment = "blue", std_error = 0)) )
      
      # Test set
      test <- subset(pheno_use, environment == ve & line_name %in% c(tp_geno, vp_geno))
      
      # Empty vector to store prediction accuracy estimates
      pred1 <- vector("list", length(train))
      
      # Loop
      for (r in seq_along(pred1)) {
        train_r <- train[[r]]
        
        # Prediction accuracy estimates
        test_val_pred <- map(train_r, ~gblup(K = K, train = ., test = test, fit.env = FALSE)[[2]]) %>%
          map2_df(., seq_along(.), ~mutate(.x, nTrainEnv = .y, scheme = ifelse(line_name %in% tp_geno, "cv0", "pocv0")))
        
        # Data.frame
        pred1[[r]] <- test_val_pred %>% 
          group_by(nTrainEnv, scheme) %>% 
          summarize(accuracy = cor(value, pred_value)) %>%
          ungroup()
        
      }
      
      ## output data.frame 
      out[[i]] <- map2_df(pred1, seq_along(pred1), ~mutate(.x, rep = .y))
      
    }
    
    
    core_df %>%
      mutate(data = out) %>%
      select(-core) %>%
      unnest(data)
    
  })









pred_list <- ls(pattern = "_predictions")

# Save
save_file <- file.path(result_dir, "distance_rank_predictions.RData")
save(list = pred_list, file = save_file)



