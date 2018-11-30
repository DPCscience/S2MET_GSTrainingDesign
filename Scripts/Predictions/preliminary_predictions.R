## Test Predictions in the S2MET Prediction Project
## 
## Author: Jeff Neyhart
## Last Updated: April 4, 2018
## 
## This script will conduct some preliminary predictions using the S2MET data. 
## This will mainly include intra-environment predictions and normalization 
## based on heritability.
## 

## Run on a local machine
library(modelr)

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
  filter(line_name %in% c(tp_geno, vp_geno),
         environment %in% tp_vp_env) %>%
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







### Question 2 - What is the intra-environmental prediction accuracy (i.e. what is
### the accuracy of an environment to predict itself?)
### We may also ask what the accuracy is of a random sample of environments to predict
### that environment.

# Subset the dataset for only environments in which both the TP and VP were phenotyped
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         environment %in% tp_vp_env) %>%
  mutate_at(vars(environment:line_name), as.factor)

# Group by traits and environments and create training and test datasets
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs_use %>% 
  mutate(pop = if_else(line_name %in% tp_geno, "train", "test"),
         env = environment) %>% # Copy environment for nesting 
  group_by(trait, environment, pop) %>%
  nest(env, line_name, value, std_error) %>%
  spread(pop, data)

# Iterate over traits and environments
intra_environment_predictions <- S2_MET_BLUEs_tomodel %>%
  group_by(trait, environment) %>%
  do({
    df <- .
    train <- df$train[[1]]
    test <- df$test[[1]]
    
    # cat(.$trait, as.character(.$environment), "\n")
    
    # Run predictions
    pred_out <- gblup(K = K, train = train, test = test, bootrep = 1000)
    # Convert to data.frame and return
    data_frame(accuracy = list(pred_out$boot), pgv = list(pred_out$pgv))
    
  })

# Add the heritability information
intra_environment_predictions1 <- stage_one_data %>% 
  group_by(environment, trait) %>% 
  distinct(environment, trait, heritability) %>% 
  slice(1) %>% 
  right_join(., intra_environment_predictions, by = c("environment", "trait")) %>%
  ungroup()


## Create scatterplots of predicted versus observed values
# First create the base plot
acc_plots_paginate <- intra_environment_predictions %>%
  unnest(pgv) %>% 
  ggplot(aes(x = value, y = pred_value)) + 
  geom_point() + 
  theme_bw(base_size = 8)

acc_plots_paginate1 <- acc_plots_paginate + 
  facet_wrap_paginate(trait ~ environment, nrow = 5, ncol = 5, scales = "free",
                      labeller = labeller(.multi_line = FALSE), page = 1)

acc_plots_paginate2 <- acc_plots_paginate + 
  facet_wrap_paginate(trait ~ environment, nrow = 5, ncol = 5, scales = "free",
                      labeller = labeller(.multi_line = FALSE), page = 2)

acc_plots_paginate3 <- acc_plots_paginate + 
  facet_wrap_paginate(trait ~ environment, nrow = 5, ncol = 5, scales = "free",
                      labeller = labeller(.multi_line = FALSE), page = 3)

acc_plots_paginate4 <- acc_plots_paginate + 
  facet_wrap_paginate(trait ~ environment, nrow = 5, ncol = 5, scales = "free",
                      labeller = labeller(.multi_line = FALSE), page = 4)



# Plot using barplots
intra_environment_predictions_toplot <- intra_environment_predictions1 %>%
  unnest(accuracy) %>%
  group_by(trait, environment) %>%
  mutate(accuracy_adj = cor / sqrt(heritability),
         significant = between(x = 0, left = ci_lower, right = ci_upper), 
         significant_ann = if_else(significant, "*", ""))

g_plot <- intra_environment_predictions_toplot %>% 
  ggplot(aes(x = environment, y = cor)) +
  geom_col(fill = umn_palette(2, 3)[3]) + 
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.5) + 
  ylim(c(-0.5), 1) +
  ylab(expression(r[MG])) +
  xlab("Environment") +
  labs(title = "Intra-environment Predictions") +
  facet_grid(trait ~ .) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Save
save_file <- file.path(fig_dir, "intra_environment_predictions.jpg")
ggsave(filename = save_file, plot = g_plot, width = 10, height = 6, dpi = 1000)


## Sort on grain yield accuracy
environment_order <- intra_environment_predictions_toplot %>% 
  filter(trait == "GrainYield") %>% 
  arrange(cor) %>% 
  pull(environment) %>%
  c(., setdiff(tp_vp_env, .))

intra_environment_predictions_toplot1 <- intra_environment_predictions_toplot %>% 
  ungroup() %>% 
  mutate(environment = factor(environment, levels = environment_order))

g_plot <- intra_environment_predictions_toplot1 %>% 
  ggplot(aes(x = environment, y = cor)) +
  geom_col(fill = umn_palette(2, 3)[3]) + 
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.5) + 
  ylim(c(-0.5), 1) +
  ylab(expression(r[MG])) +
  xlab("Environment") +
  labs(title = "Intra-environment Predictions") +
  facet_grid(trait ~ .) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Save
save_file <- file.path(fig_dir, "intra_environment_predictions_sorted.jpg")
ggsave(filename = save_file, plot = g_plot, width = 10, height = 6, dpi = 1000)






### Question 3 - what is the relationship between the intra-environment prediction
### accuracy and the environment heritability?
# subset the heritability information
env_heritability <- stage_one_data %>% 
  group_by(environment, trait) %>% 
  select(trait, environment, heritability:ci_upper) %>% 
  slice(1) %>% 
  ungroup() %>%
  rename(estimate = heritability) %>%
  rename_at(vars(estimate:ci_upper), ~str_c("heritability_", .))
  
env_accuracy <- intra_environment_predictions_toplot %>% 
  select(environment, trait, estimate = cor, se:ci_upper) %>% 
  ungroup() %>% 
  rename_at(vars(estimate:ci_upper), ~str_c("accuracy_", .))

# Combine
env_measures <- left_join(env_accuracy, env_heritability, by = c("environment", "trait"))

## Plot
g_acc_herit <- env_measures %>% 
  ggplot(aes(x = (heritability_estimate), y = accuracy_estimate)) + 
  geom_segment(aes(x = heritability_ci_lower, xend = heritability_ci_upper, yend = accuracy_estimate),
               color = umn_palette(3, 4)[3]) + 
  geom_segment(aes(xend = heritability_estimate, y = accuracy_ci_lower, yend = accuracy_ci_upper),
               color = umn_palette(3, 4)[4]) + 
  geom_point() + 
  # geom_smooth(method = "lm", se = FALSE) +
  ylab("Prediction Accuracy") +
  xlab("Heritability") +
  facet_wrap(~trait, ncol = 2, scales = "free_x") +
  theme_bw()

# Save
save_file <- file.path(fig_dir, "accuracy_and_heritability.jpg")
ggsave(filename = save_file, plot = g_acc_herit, width = 6, height = 6, dpi = 1000)







## What is the accuracy to predict each environment using data from all other environments (i.e. leave-one-env-out)?
## 
## 


# Modify the BLUEs for predictions
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno),
         environment %in% tp_vp_env,
         trait %in% traits) %>%
  mutate_at(vars(environment:line_name), as.factor) %>%
  # group_by(trait, environment) %>%
  # mutate(value = scale(value)) %>%
  ungroup()


## Create training and test sets
## Complete data
train_test_complete <- S2_MET_BLUEs_use %>%
  group_by(trait) %>%
  do({
    df <- droplevels(.)
    
    # Get a vector of testing environments
    val_environment <- distinct(df, environment)
    
    # Create testing and training sets based on that environment
    val_environment %>% 
      mutate(train = map(environment, ~filter(df, environment != .) %>% droplevels() %>% filter(line_name %in% tp_geno)),
             test = map(environment, ~filter(df, line_name %in% vp_geno, environment == .) %>% mutate(line_name = as.character(line_name))))
    
  }) %>% ungroup()

# Realistic data
train_test_realistic <- S2_MET_BLUEs_use %>%
  group_by(trait) %>%
  do({
    df <- droplevels(.)
    
    val_environment <- distinct(df, environment) %>%
      filter(str_detect(environment, "17"))
    
    # Get a vector of testing environments
    all_val_environment <- as.character(val_environment$environment)
      
    # Create testing and training sets based on that environment
    val_environment %>% 
      mutate(train = map(environment, ~filter(df, environment %in% all_val_environment) %>% droplevels() %>% filter(line_name %in% tp_geno)),
             test = map(environment, ~filter(df, line_name %in% vp_geno, environment == .) %>% mutate(line_name = as.character(line_name))))
    
  }) %>% ungroup()


## Combine
train_test <- bind_rows(
  mutate(train_test_complete, set = "complete"),
  mutate(train_test_realistic, set = "realistic")
)
  
  


### Leave-one-environment-out predictions using the genotype mean over training environments
environment_loeo_predictions_geno_mean <- train_test %>%
  group_by(set, trait, environment) %>%
  do(predictions = gblup(K = K, train = .$train[[1]], test = .$test[[1]], fit.env = TRUE)) %>%
  ungroup()


## Add the number of training environments to the results
environment_loeo_predictions_geno_mean <- environment_loeo_predictions_geno_mean %>% 
  left_join(train_test) %>% 
  mutate(nTrainEnv = map_dbl(train, ~n_distinct(.$environment))) %>% 
  select(set, trait, environment, nTrainEnv, predictions)
  



## Save
save_file <- file.path(result_dir, "all_data_environmental_predictions.RData")
save("environment_loeo_predictions_geno_mean", file = save_file)








