## S2MET Simple Genomic Predictions
## 
## Predict each environment only using the data for that environment

library(tidyverse)
library(stringr)
library(readxl)
library(sommer)
library(modelr)
library(gws)
library(pbr)
library(purrrlyr)

source("source.R")

# A matrix
A <- A.mat(X = s2_imputed_mat_use, min.MAF = 0, max.missing = 1)

## Filter out environments with no training data or environments with no
## vp data
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  group_by(trait, environment) %>% 
  filter(sum(line_name %in% tp_geno) > 1, 
         sum(line_name %in% vp_geno) > 1) %>%
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  droplevels()


# For each trait and environment, create training and prediction datasets
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs_use %>% 
  group_by(environment, trait) %>% 
  do(train_dat = filter(., line_name %in% tp_geno),
     test_dat = filter(., line_name %in% vp_geno))


# Group by trait and environment and run predictions
intra_env_pred <- S2_MET_BLUEs_tomodel %>%
  by_row(function(i) {
    
    dat <- i$train_dat[[1]]
    test_dat <- i$test_dat[[1]] %>%
      mutate(line_name = as.character(line_name))
    
    # Create a model.frame
    pred_mf <- model.frame(value ~ line_name + std_error, data = dat, subset = line_name %in% tp_geno)
    y <- model.response(pred_mf)
    X <- model.matrix(~ 1, pred_mf)
    Z <- model.matrix(~ line_name, pred_mf) %>%
      structure(dimnames = list(NULL, levels(pred_mf$line_name)))
    
    Zlist <- list(line_name = list(Z = Z, K = A))
    
    # R matrix
    wts <- pred_mf$std_error^2
    R <- solve(diag(wts))
    
    # Fit with and without the R matrix
    fit <- mmer(Y = y, X = X, Z = Zlist, R = list(res = R), silent = TRUE)

    # Extract the predictions
    pgv <- fit$u.hat$line_name %>%
      as.data.frame() %>% 
      rownames_to_column("line_name") %>% 
      rename(pred = T1)
    
    # Return the PGVs for the vp
    left_join(test_dat, pgv, by = "line_name") %>% 
      dplyr::select(line_name, value, pred)
    
  }, .to = "pred_out")



# Combine with the phenos and measure accuracy
intra_env_pred_acc <- intra_env_pred %>% 
  unnest(pred_out) %>%
  group_by(environment, trait) %>%
  do(boot_cor(x = .$value, y = .$pred, boot.reps = 1000))


# Save
save_file <- file.path(result_dir, "S2MET_intra_env_pred.RData")
save("intra_env_pred_acc", file = save_file)


