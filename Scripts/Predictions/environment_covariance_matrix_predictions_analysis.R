## Analysis of predictions using different covariance matrices
##
##
## Author: Jeff Neyhart
## Last updated: October 10, 2018
##

# # Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load the results
load(file.path(result_dir, "env_cov_mat_predictions.RData"))




## Summarize predictive ability
predictions_df <- predictions_out$ge %>%
  unnest(predictions) %>%
  select(environment, trait, model, line_name, value, pred_value_ge = pred_value) %>%
  left_join(., unnest(predictions_out$g) %>% select(environment, trait, line_name, pred_value_g = pred_value))

prediction_acc_summ <- predictions_df %>% 
  group_by(trait, model, environment) %>% 
  summarize_at(vars(contains("pred_value")), funs(acc = cor(., value))) %>%
  # Difference from the g model
  mutate(advantage = pred_value_ge_acc - pred_value_g_acc)

# Plot



## Calculate accuracy for all environments (it's LOEO, after all)
prediction_acc_summ_all <- predictions_df %>% 
  group_by(trait, model) %>% 
  summarize_at(vars(contains("pred_value")), funs(acc = cor(., value))) %>%
  # Difference from the g model
  mutate(advantage = pred_value_ge_acc - pred_value_g_acc)

# Plot
prediction_acc_summ_all %>%
  ggplot(aes(x = model, y = advantage)) +
  geom_col()

# Plot the correlation
predictions_df %>% 
  ggplot(aes(x = pred_value_ge, y = value, color = environment)) +
  # geom_abline(intercept = 0, slope = 1) + 
  geom_point() + 
  facet_wrap(~ trait + model, scales = "free_x")














