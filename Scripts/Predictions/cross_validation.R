## S2MET cross-validation
## 
## This script will generate cross-validation samples and then test different models for prediction
## accuracy.
## 
## Author: Jeff Neyhart
## Last modified: January 17, 2019
## 
## 

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Predictions/"
source(file.path(repo_dir, "source_MSI.R"))

## Number of cores
n_core <- 32
n_core <- detectCores()


# ## Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# 
# # Other packages
# library(modelr)

load(file.path(result_dir, "distance_method_results.RData"))


## Some parameters
# The proportion of lines to use for testing (in CV schemes)
# pTest <- 0.2
pTest <- length(vp_geno) / length(tp_geno)
pTrain <- 1 - pTest

# Number of CV iterations
nCV <- 10






## Relationship matrix for CV and POV
K_cv <- K[tp_geno, tp_geno]
K_pov <- K


## Data to use for CV
cv_data <- S2_MET_BLUEs %>%
  # filter(environment %in% sample(tp_vp_env, 10)) %>%
  filter(line_name %in% tp_geno,
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.)))


## Generate CV randomizations
# First generate training/test lines
set.seed(1171039)
train_test <- data_frame(.id = str_pad(as.character(seq(nCV)), width = 2, side = "left", pad = 0), train = replicate(n = nCV, sample_frac(tbl = distinct(cv_data, line_name), size = pTrain)$line_name, simplify = FALSE)) %>%
  mutate(test = map(train, ~setdiff(tp_geno, .))) %>%
  crossing(trait = traits, .)

# CV1 - predict totally untested genotypes using all environments
cv1_rand <- train_test %>%
  mutate_at(vars(train, test), funs(map2(.x = ., .y = trait, ~filter(cv_data, line_name %in% .x, trait == .y))))

# CV2 - predict partially tested genotypes
cv2_rand <- cv_data %>%
  group_by(trait) %>%
  do(crossv_mc(data = ., n = nCV, test = pTest)) %>%
  ungroup()

# CV0 - predict totally tested genotypes in unobserved environments
cv0_rand <- cv_data %>%
  group_by(trait) %>%
  do({
    df <- .
    
    distinct(df, environment) %>% 
      mutate(train = map(environment, ~filter(df, environment != .)), 
             test = map(environment, ~filter(df, environment == .)))
  }) %>% ungroup() %>%
  mutate(.id = "01")

# CV00 - predict totally untested genotypes in totally untested environments
cv00_rand <- full_join(train_test, select(cv0_rand, -.id), by = "trait") %>%
  mutate(train = map2(train.x, train.y, ~filter(.y, line_name %in% .x)),
         test = map2(test.x, test.y, ~filter(.y, line_name %in% .x))) %>%
  select(trait, environment, .id, train, test) %>%
  arrange(trait, environment)


## Combine CV1 and CV2
CV12 <- bind_rows(
  mutate(cv1_rand, cv = "cv1"),
  mutate(cv2_rand, cv = "cv2")
)

CV_zero <- bind_rows(
  mutate(cv0_rand, cv = "cv0"),
  mutate(cv00_rand, cv = "cv00")
)



# ## Predict
# cv12_prediction <- CV12 %>%
#   group_by(cv, trait, .id) %>%
#   do({
#     row <- .
#     
#     ## Model 2
#     model2_out <- model2(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv)
#     ## Model 3
#     model3_out <- model3(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv)
#     
#     data_frame(model = c("M2", "M3"), prediction = list(model2_out, model3_out))
# 
#   })
# 
# 
# 
# cv_zero_prediction <- CV_zero %>%
#   group_by(cv, trait, environment, .id) %>%
#   do({
#     row <- .
#     
#     ## Model 2
#     model2_out <- model2(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv)
#     ## Model 3
#     model3_out <- model3(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv)
#     
#     data_frame(model = c("M2", "M3"), prediction = list(model2_out, model3_out))
#     
#   })




## Use parallelization
cv12_prediction <- CV12 %>%
  assign_cores(n_core = n_core) %>% 
  split(.$core) %>%
  mclapply(X = ., FUN = function(core_df) {
    
    results_out <- vector("list", nrow(core_df))
    # Iterate over core_df
    for (i in seq_along(results_out)) {
      
      row <- core_df[i,]
    
      # ## Model 2
      # model2_out <- model2(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv)
      # ## Model 3
      # model3_out <- model3(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv)
      # 
      # results_out[[i]] <- data_frame(model = c("M2", "M3"), prediction = list(model2_out, model3_out))
      # 
      
      ## Model 4
      model4_out <- model4(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv)
      ## Model 5
      model5_out <- model5(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv)
    
      results_out[[i]] <- data_frame(model = c("M4", "M5"), prediction = list(model4_out, model5_out))
      
      
    }
    
    core_df %>% 
      mutate(results = results_out) %>%
      select(cv, trait, .id, results) %>%
      unnest(results)
      
  })

cv12_prediction <- bind_rows(cv12_prediction)

## CV-zero
cv_zero_prediction <- CV_zero %>%
  assign_cores(n_core = n_core) %>% 
  split(.$core) %>%
  mclapply(X = ., FUN = function(core_df) {
    
    results_out <- vector("list", nrow(core_df))
    # Iterate over core_df
    for (i in seq_along(results_out)) {
      
      row <- core_df[i,]
      
      # ## Model 2
      # model2_out <- model2(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv)
      # ## Model 3
      # model3_out <- model3(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv)
      # 
      # results_out[[i]] <- data_frame(model = c("M2", "M3"), prediction = list(model2_out, model3_out))
      # 
      
      ## Model 4
      model4_out <- model4(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv)
      ## Model 5
      model5_out <- model5(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_cv)
      
      results_out[[i]] <- data_frame(model = c("M4", "M5"), prediction = list(model4_out, model5_out))
      
    }
    
    core_df %>% 
      mutate(results = results_out) %>%
      select(cv, trait, .id, results) %>% 
      unnest(results)
    
  })

cv_zero_prediction <- bind_rows(cv_zero_prediction)



## Combine
cv_predictions <- bind_rows(cv12_prediction, cv_zero_prediction)








### Parent-offspring cross-validation

## Data to use for CV
pocv_data <- S2_MET_BLUEs %>% 
  filter(environment %in% unique(cv_data$environment)) %>%
  filter(line_name %in% c(tp_geno, vp_geno),
         trait %in% traits,
         environment %in% tp_vp_env) %>%
  mutate(id = seq(nrow(.)))


## Generate CV randomizations
# First generate training/test lines
train_test1 <- train_test %>%
  mutate(test = list(vp_geno))


# POCV1 - predict the VP using all environments
pocv1_rand <- train_test1 %>%
  mutate_at(vars(train, test), funs(map2(.x = ., .y = trait, ~filter(pocv_data, line_name %in% .x, trait == .y))))

# POCV2 - predict partially tested VP
pocv2_rand <- pocv_data %>%
  group_by(trait) %>%
  do(crossv_mc(data = ., n = nCV, test = pTest)) %>%
  ungroup()

# POCV0 - predict totally tested VP in unobserved environments
pocv0_rand <- pocv_data %>%
  group_by(trait) %>%
  do({
    df <- .
    
    distinct(df, environment) %>% 
      mutate(train = map(environment, ~filter(df, environment != .)), 
             test = map(environment, ~filter(df, environment == ., line_name %in% vp_geno)))
  }) %>% ungroup() %>%
  mutate(.id = "01")

# POCV00 - predict totally untested VP in totally untested environments
pocv00_rand <- full_join(train_test1, select(pocv0_rand, -.id), by = "trait") %>%
  mutate(train = map2(train.x, train.y, ~filter(.y, line_name %in% .x)), # Take the available training data and subset the TP
         test = map2(test.x, test.y, ~filter(.y, line_name %in% .x))) %>% # Take the available testing data and subset the VP
  select(trait, environment, .id, train, test) %>%
  arrange(trait, environment)

## Combine POCV1 and POCV2
POCV12 <- bind_rows(
  mutate(pocv1_rand, cv = "pocv1"),
  mutate(pocv2_rand, cv = "pocv2")
)

POCV_zero <- bind_rows(
  mutate(pocv0_rand, cv = "pocv0"),
  mutate(pocv00_rand, cv = "pocv00")
)


# ## Predict
# pocv12_prediction <- POCV12 %>%
#   group_by(cv, trait, .id) %>%
#   do({
#     row <- .
# 
#     ## Model 2
#     model2_out <- model2(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_pov)
#     ## Model 3
#     model3_out <- model3(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_pov)
# 
#     data_frame(model = c("M2", "M3"), prediction = list(model2_out, model3_out))
# 
#   })
# 
# 
# 
# pocv_zero_prediction <- POCV_zero %>%
#   group_by(cv, trait, environment, .id) %>%
#   do({
#     row <- .
# 
#     ## Model 2
#     model2_out <- model2(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_pov)
#     ## Model 3
#     model3_out <- model3(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_pov)
# 
#     data_frame(model = c("M2", "M3"), prediction = list(model2_out, model3_out))
# 
#   })






## Use parallelization
pocv12_prediction <- POCV12 %>%
  assign_cores(n_core = n_core) %>% 
  split(.$core) %>%
  mclapply(X = ., FUN = function(core_df) {
    
    results_out <- vector("list", nrow(core_df))
    # Iterate over core_df
    for (i in seq_along(results_out)) {
      
      row <- core_df[i,]
      
      # ## Model 2
      # model2_out <- model2(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_pov)
      # ## Model 3
      # model3_out <- model3(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_pov)
      # 
      # results_out[[i]] <- data_frame(model = c("M2", "M3"), prediction = list(model2_out, model3_out))
      # 
      
      ## Model 4
      model4_out <- model4(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_pov)
      ## Model 5
      model5_out <- model5(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_pov)
      
      results_out[[i]] <- data_frame(model = c("M4", "M5"), prediction = list(model4_out, model5_out))
      
    }
    
    core_df %>% 
      mutate(results = results_out) %>%
      select(cv, trait, .id, results)
    
  })

pocv12_prediction <- bind_rows(pocv12_prediction)

## CV-zero
pocv_zero_prediction <- POCV_zero %>%
  assign_cores(n_core = n_core) %>% 
  split(.$core) %>%
  mclapply(X = ., FUN = function(core_df) {
    
    results_out <- vector("list", nrow(core_df))
    # Iterate over core_df
    for (i in seq_along(results_out)) {
      
      row <- core_df[i,]
      
      # ## Model 2
      # model2_out <- model2(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_pov)
      # ## Model 3
      # model3_out <- model3(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_pov)
      # 
      # results_out[[i]] <- data_frame(model = c("M2", "M3"), prediction = list(model2_out, model3_out))
      # 
      
      ## Model 4
      model4_out <- model4(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_pov)
      ## Model 5
      model5_out <- model5(train = as.data.frame(row$train[[1]]), test = as.data.frame(row$test[[1]]), Kg = K_pov)
      
      results_out[[i]] <- data_frame(model = c("M4", "M5"), prediction = list(model4_out, model5_out))
      
    }
    
    core_df %>% 
      mutate(results = results_out) %>%
      select(cv, trait, .id, results) %>% 
      unnest(results)
    
  })

pocv_zero_prediction <- bind_rows(pocv_zero_prediction)



## Combine
pocv_predictions <- bind_rows(pocv12_prediction, pocv_zero_prediction)



## Save
save("cv_predictions", "pocv_predictions", file = file.path(result_dir, "cross_validation.RData"))






