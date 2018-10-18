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



## MC predictions
mc_pred_tidy <- environment_mc_predictions %>%
  unnest() %>% unnest() %>% 
  select(-trait1, -trait2) %>%
  mutate_at(vars(trait, model, environment), as.factor)

## Summarize the correlations across all environments for each iteration
mc_pred_summ <- mc_pred_tidy %>% 
  group_by(trait, model, environment, iter) %>% 
  summarize(accuracy = cor(value, pgv))

# Now take the mean over iterations
mc_pred_summ1 <- mc_pred_summ %>% 
  summarize(accuracy = mean(accuracy))

# Plot for each model and trait
g_model_acc <- mc_pred_summ1 %>%
  ggplot(aes(x = trait, y = accuracy, fill = model)) +
  geom_boxplot(position = "dodge", alpha = 0.5) +
  xlab("Trait") + 
  ylab("Prediction accuracy") +
  scale_fill_discrete(name = "Model") +
  theme_classic()

ggsave(filename = "environmental_cov_model_accuracy.jpg", plot = g_model_acc, path = fig_dir, width = 6, height = 4, dpi = 1000)
