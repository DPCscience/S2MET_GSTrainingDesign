## Test Predictions in the S2MET Prediction Project
## 
## Author: Jeff Neyhart
## Last Updated: March 15, 2018
## 
## This script will conduct some preliminary predictions using the S2MET data
## and potentially some clustering strategy
## 

# Load packages and the source script
library(tidyverse)
library(broom)
library(stringr)
library(readxl)
library(modelr)
library(pbr)
library(rrBLUP)
library(ggridges)

# The head directory
repo_dir <- getwd()

source(file.path(repo_dir, "source.R"))





## Question 1
## What is the relationship between the level of GxE in a "cluster" and the mean
## prediction accuracy in that cluster?
## 

# Test with grain yield
# Remove environments without both the TP and the VP
pheno_tomodel <- S2_MET_BLUEs %>% 
  # filter(trait == "GrainYield") %>%
  filter(trait == "HeadingDate") %>%
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) > 1,
         sum(line_name %in% vp_geno)) %>%
  ungroup() %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  mutate_at(vars(environment:line_name), as.factor)

# Randomly sample 2 environments
envs <- unique(pheno_tomodel$environment)
sample_envs_list <- sample(x = combn(x = envs, 2, simplify = FALSE), size = 250)

# Iterate over the sample environments and calculate GxE
# Then iterate over the environments in each pair and make predictions
sample_envs_results <- sample_envs_list %>%
  map_df(function(sample_envs) {
    
    # print(str_c(sample_envs, collapse = "_"))
    
    df <- pheno_tomodel %>%
      filter(environment %in% sample_envs)
    
    # Calculate variance
    var_comp <- calc_variance(data = df)
    
    # Save the gxe estimate
    sample_gxe <- data.frame(environment1 = sample_envs[1], environment2 = sample_envs[2], 
                             var_gxe = subset(var_comp, var_comp == "line_name:environment", variance, drop = TRUE))
    
    # Iterate over the sample environments, using each as a validation environment
    sample_pred <- map_df(sample_envs, function(val_env) {
      # The remaining environments are used for prediction
      pred_envs <- setdiff(sample_envs, val_env)
      
      # Create the training df
      train_df <- df %>% 
        filter(environment %in% pred_envs,
               line_name %in% tp_geno)
      
      # Create the validation df
      val_df <- df %>% 
        filter(environment %in% val_env,
               line_name %in% vp_geno)
      
      # Subset the heritability of that environment
      val_herit <- stage_one_data %>%
        filter(environment == val_env,
               trait == unique(df$trait)) %>%
        pull(heritability) %>%
        tail(1)
      
      # Create model matrices
      mf <- model.frame(value ~ line_name + environment, train_df)
      y <- model.response(mf)
      Z <- model.matrix(~ line_name, mf)

      # Fit the model
      pred <- mixed.solve(y = y, Z = Z, K = K)
      
      # Combine with the phenotypic data in the validation environment
      pred_vp_pheno <- pred$u %>% 
        data.frame(line_name = names(.), pred_value = ., row.names = NULL, stringsAsFactors = FALSE) %>% 
        left_join(val_df, ., by = "line_name")
      
      # Calculate and return accuracy
      data.frame(train_environment = pred_envs, val_environment = val_env,
                 accuracy = cor(pred_vp_pheno$pred_value, pred_vp_pheno$value),
                 val_environment_herit = val_herit,
                 stringsAsFactors = FALSE, row.names = NULL)
      
    })
    
    # Adjust the accuracies and take the mean
    sample_pred_adj <- sample_pred %>%
      mutate(pred_ability = accuracy / sqrt(val_environment_herit))
    
    # Take the mean
    sample_pred_adj_mean <- mean(sample_pred_adj$accuracy)
    sample_pred_abil_mean <- mean(sample_pred_adj$pred_ability)
    
    # Return a nice data.frame
    as_data_frame(sample_gxe) %>% 
      mutate(pred_acc = sample_pred_adj_mean,
             pred_ability = sample_pred_abil_mean, sample_pred = list(sample_pred_adj))
    
  })


# Plot
sample_envs_results %>% 
  ggplot(aes(x = var_gxe, y = pred_acc)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)
















