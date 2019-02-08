## S2MET cross-validation
## 
## This script will generate cross-validation samples and then test different models for prediction
## accuracy.
## 
## Author: Jeff Neyhart
## Last modified: January 17, 2019
## 

# ## Parent-offspring cross-validation
# 
# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# 
# # Other packages
# library(modelr)



# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions/"
source(file.path(repo_dir, "source_MSI.R"))

## Number of cores
n_core <- 32
n_core <- detectCores()


load(file.path(result_dir, "distance_method_results.RData"))



# Number of folds
k <- 5

# Number of CV iterations
nCV <- 50






## Relationship matrix for CV and POV
K_cv <- K[tp_geno, tp_geno]
K_pov <- K


# ## List of E correlation matrices
env_cor_mats <- env_rank_df %>%
  filter((str_detect(model, "MYEC") & mat_set == "Jarquin") | model %in% c("pheno_loc_dist", "pheno_dist")) %>%
  select(-dist, -env_rank)




### Parent-offspring validation

## Data to use for CV
pov_data <- S2_MET_BLUEs %>% 
  filter(environment %in% tp_vp_env) %>%
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.)))


## POV1 - predict the VP in tested environments
pov1_rand <- pov_data %>%
  group_by(trait) %>%
  do(data_frame(train = list(filter(., line_name %in% tp_geno)), test = list(filter(., line_name %in% vp_geno)))) %>% 
  ungroup() %>%
  mutate(.id = "01")


## POV0 - predict the VP in untested environments using the TP and VP from tested environments
# LOEO
pov0_rand_loeo <- pov_data %>%
  group_by(trait) %>%
  do({
    df <- .
    
    distinct(df, environment) %>%
      mutate(train = map(environment, ~filter(df, line_name %in% c(tp_geno, vp_geno), environment != .)),
             test = map(environment, ~filter(df, line_name %in% vp_geno, environment == .))) }) %>%
  ungroup() %>%
  mutate(.id = "01")

pov0_rand_future <- pov_data %>%
  group_by(trait) %>%
  do(data_frame(train = list(filter(., line_name %in% c(tp_geno, vp_geno), year != 2017)), 
                test = list(filter(., line_name %in% vp_geno, year == 2017)))) %>% 
  ungroup() %>%
  mutate(.id = "01")
    



## POV00 - predict the VP in untested environments
# LOEO - leave one environment out
pov00_rand_loeo <- pov_data %>%
  group_by(trait) %>%
  do({
    df <- .
    
    distinct(df, environment) %>%
      mutate(train = map(environment, ~filter(df, line_name %in% tp_geno, environment != .)),
             test = map(environment, ~filter(df, line_name %in% vp_geno, environment == .))) }) %>%
  ungroup() %>%
  mutate(.id = "01")

pov00_rand_future <- pov_data %>%
  group_by(trait) %>%
  do(data_frame(train = list(filter(., line_name %in% tp_geno, year != 2017)), 
                test = list(filter(., line_name %in% vp_geno, year == 2017)))) %>% 
  ungroup() %>%
  mutate(.id = "01")





## Predict
pov1_prediction <- pov1_rand %>%
  group_by(trait) %>%
  do({
    row <- .
    tr_env <- tp_vp_env_trait[[row$trait[[1]]]]
    
    
    ## Extract train/test
    train <- as.data.frame(row$train[[1]]) %>%
      mutate(environment = factor(environment, levels = tr_env))
    test <-  as.data.frame(row$test[[1]])
    
    
    ## Model 4
    model4_out <- model4(train = train, test = test, Kg = K_pov)
    ## Model 5
    # model5_out <- model5(train = train, test = test, Kg = K_pov)
    
    ## Model 5 PD uses the phenotypic correlation as the E covariance matrix
    # Calculate the correlation between environments
    Ke <- subset(env_cor_mats, model == "pheno_dist" & trait == row$trait[1], cov, drop = T)[[1]] 
    Ke <- Ke[tr_env, tr_env]
    
    model5_PD_out <- model5(train = train, test = test, Kg = K_pov, Ke = Ke)
    
    data_frame(model = c("M4", "M5_PD"), 
               prediction = list(model4_out, model5_PD_out))
  
  }) %>% ungroup()




## Predict
pov0_loeo_prediction <- pov0_rand_loeo %>%
  group_by(trait, environment) %>%
  do({
    row <- .
    tr_env <- tp_vp_env_trait[[row$trait[[1]]]]
    
    
    ## Extract train/test
    train <- as.data.frame(row$train[[1]]) %>%
      mutate(environment = factor(environment, levels = tr_env))
    test <-  as.data.frame(row$test[[1]])
    
    
    ## Model 4
    model4_out <- model4(train = train, test = test, Kg = K_pov)
    ## Model 5
    # model5_out <- model5(train = train, test = test, Kg = K_pov)
    
    ## Model 5 PD uses the phenotypic correlation as the E covariance matrix
    # Calculate the correlation between environments
    Ke <- subset(env_cor_mats, model == "pheno_dist" & trait == row$trait[1], cov, drop = T)[[1]] 
    Ke <- Ke[tr_env, tr_env]
    
    model5_PD_out <- model5(train = train, test = test, Kg = K_pov, Ke = Ke)
    
    data_frame(model = c("M4", "M5_PD"), 
               prediction = list(model4_out, model5_PD_out))
    
  }) %>% ungroup()


## Predict
pov0_future_prediction <- pov0_rand_future %>%
  group_by(trait) %>%
  do({
    row <- .
    tr_env <- tp_vp_env_trait[[row$trait[[1]]]]
    
    
    ## Extract train/test
    train <- as.data.frame(row$train[[1]]) %>%
      mutate(environment = factor(environment, levels = tr_env))
    test <-  as.data.frame(row$test[[1]])
    
    
    ## Model 4
    model4_out <- model4(train = train, test = test, Kg = K_pov)
    ## Model 5
    # model5_out <- model5(train = train, test = test, Kg = K_pov)
    
    ## Model 5 PD uses the phenotypic correlation as the E covariance matrix
    # Calculate the correlation between environments
    Ke <- subset(env_cor_mats, model == "pheno_dist" & trait == row$trait[1], cov, drop = T)[[1]] 
    Ke <- Ke[tr_env, tr_env]
    
    model5_PD_out <- model5(train = train, test = test, Kg = K_pov, Ke = Ke)
    
    data_frame(model = c("M4", "M5_PD"), 
               prediction = list(model4_out, model5_PD_out))
    
  }) %>% ungroup()



## Predict
pov00_loeo_prediction <- pov00_rand_loeo %>%
  group_by(trait, environment) %>%
  do({
    row <- .
    tr_env <- tp_vp_env_trait[[row$trait[[1]]]]
    
    
    ## Extract train/test
    train <- as.data.frame(row$train[[1]]) %>%
      mutate(environment = factor(environment, levels = tr_env))
    test <-  as.data.frame(row$test[[1]])
    
    
    ## Model 4
    model4_out <- model4(train = train, test = test, Kg = K_pov)
    ## Model 5
    # model5_out <- model5(train = train, test = test, Kg = K_pov)
    
    ## Model 5 PD uses the phenotypic correlation as the E covariance matrix
    # Calculate the correlation between environments
    Ke <- subset(env_cor_mats, model == "pheno_dist" & trait == row$trait[1], cov, drop = T)[[1]] 
    Ke <- Ke[tr_env, tr_env]
    
    model5_PD_out <- model5(train = train, test = test, Kg = K_pov, Ke = Ke)
    
    data_frame(model = c("M4", "M5_PD"), 
               prediction = list(model4_out, model5_PD_out))
    
  }) %>% ungroup()


## Predict
pov00_future_prediction <- pov00_rand_future %>%
  group_by(trait) %>%
  do({
    row <- .
    tr_env <- tp_vp_env_trait[[row$trait[[1]]]]
    
    
    ## Extract train/test
    train <- as.data.frame(row$train[[1]]) %>%
      mutate(environment = factor(environment, levels = tr_env))
    test <-  as.data.frame(row$test[[1]])
    
    
    ## Model 4
    model4_out <- model4(train = train, test = test, Kg = K_pov)
    ## Model 5
    # model5_out <- model5(train = train, test = test, Kg = K_pov)
    
    ## Model 5 PD uses the phenotypic correlation as the E covariance matrix
    # Calculate the correlation between environments
    Ke <- subset(env_cor_mats, model == "pheno_dist" & trait == row$trait[1], cov, drop = T)[[1]] 
    Ke <- Ke[tr_env, tr_env]
    
    model5_PD_out <- model5(train = train, test = test, Kg = K_pov, Ke = Ke)
    
    data_frame(model = c("M4", "M5_PD"), 
               prediction = list(model4_out, model5_PD_out))
    
  }) %>% ungroup()



## Save
save("pov1_prediction", "pov0_loeo_prediction", "pov0_future_prediction", "pov00_loeo_prediction", "pov00_future_prediction", 
     file = file.path(result_dir, "cross_validation_POV.RData"))

