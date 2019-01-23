## S2MET cross-validation
## 
## # Regular cross-validation
## # CV1 and CV2
## 
## 
## 

## Source the base script
pred_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions/Scripts/Predictions/CrossValidation/"
pred_dir <- "C:/Users/jln54/GoogleDrive/BarleyLab/Projects/S2MET_Predictions/Scripts/Predictions/CrossValidation/"
source(file.path(pred_dir, "cross_validation_base.R"))


## Run cross-validation


## Use parallelization
cv12_prediction <- CV12 %>%
  assign_cores(n_core = n_core) %>% 
  split(.$core) %>%
  mclapply(X = ., FUN = function(core_df) {
    
    results_out <- vector("list", nrow(core_df))
    # Iterate over core_df
    for (i in seq_along(results_out)) {
      
      row <- core_df[i,]
      tr_env <- tp_vp_env_trait[[row$trait[[1]]]]
      
      
      ## Extract train/test
      train <- as.data.frame(row$train[[1]]) %>%
        mutate(environment = factor(environment, levels = tr_env))
      test <-  as.data.frame(row$test[[1]])
      
      
      ## Model 4
      model4_out <- model4(train = train, test = test, Kg = K_cv)
      ## Model 5
      # model5_out <- model5(train = train, test = test, Kg = K_cv)
      
      ## Model 5 PD uses the phenotypic correlation as the E covariance matrix
      # Calculate the correlation between environments
      Ke <- subset(env_cor_mats, model == "pheno_dist" & trait == row$trait[1], cov, drop = T)[[1]] 
      Ke <- Ke[tr_env, tr_env]
      
      model5_PD_out <- model5(train = train, test = test, Kg = K_cv, Ke = Ke)
      
      results_out[[i]] <- data_frame(model = c("M4", "M5_PD"), 
                                     prediction = list(model4_out, model5_PD_out))
      
      # ## Model 5 EC uses environmental covariates as the E covariance matrix
      # Ke <- subset(env_cor_mats, model == "MYEC_IPCA" & trait == row$trait[1] & set == "complete", cov, drop = T)[[1]] 
      # Ke <- Ke[unique(as.data.frame(row$train[[1]])$environment), unique(as.data.frame(row$train[[1]])$environment)]
      # 
      # model5_EC_IPCA_out <- model5(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv, Ke = Ke)
      # 
      # 
      # results_out[[i]] <- data_frame(model = c("M4", "M5", "M5_PD", "M5_EC"), 
      #                                prediction = list(model4_out, model5_out, model5_PD_out, model5_EC_IPCA_out))
      
      
    }
    
    core_df %>% 
      mutate(results = results_out) %>%
      select(cv, trait, .id, results) %>%
      unnest(results)
      
  })

cv12_prediction <- bind_rows(cv12_prediction)



## Save
save("cv12_prediction", file = file.path(result_dir, "cross_validation_CV12.RData"))






