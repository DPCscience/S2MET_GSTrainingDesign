## S2MET
## 
## Script for analyzing prediction results
## 
## Author: Jeff Neyhart
## Last modified: October 18, 2018
## 


# # Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# # Load some packages
library(lubridate)
library(ggforce)
library(ggridges)
library(gridExtra)

# Load the environmental distance df
load(file.path(result_dir, "distance_method_results.RData"))
# Load the LOEO prediction results
load(file.path(result_dir, "all_data_environmental_predictions.RData"))

## Significant level
alpha <- 0.05


## Create a factor for the distance methods
dist_method_replace <- c("great_circle_dist" = "Great Circle Distance", "pheno_dist" = "Phenotypic Distance", "OYEC_AllEC" = "One Year ECs",
                         "OYEC_Top5EC_Cor" = "One Year EC Top 5", "MYEC_AllEC" = "Ten Year ECs", "MYEC_Top5EC_Cor" = "Ten Year EC Top 5")
dist_method_abbr <- abbreviate(dist_method_replace)

colors <- umn_palette(3)
dist_colors <- c(setNames(colors[c(1:2, 3, 8, 4, 9)], dist_method_abbr), "Random" = "grey75")


#### Leave-one-environment-out predictions ####

# Calculate accuracy
loeo_accuracy <- environment_loeo_predictions_mean %>% 
  unnest() %>%
  group_by(trait, environment) %>%
  do(bootstrap(x = .$value, y = .$pgv, fun = "cor", alpha = alpha)) %>%
  ungroup()

# Calculate the mean per trait
loeo_accuracy %>% 
  group_by(trait) %>% 
  summarize(accuracy = mean(base))

# Plot per trait
g_loeo <- loeo_accuracy %>% 
  ggplot(aes(x = trait, y = base, fill = trait, color = trait)) +
  geom_boxplot(alpha = 0.5) + 
  geom_jitter(width = 0.25) + 
  ylab("Prediction accuracy") +
  xlab("Trait") +
  scale_y_continuous(breaks = pretty) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(filename = "loeo_predictions.jpg", plot = g_loeo, path = fig_dir, width = 5, height = 4, dpi = 1000)



#### Environmental distance predictions ####

# Load data
load(file.path(result_dir, "cluster_predictions_tp.RData"))


## Cumulative predictions


# Rename the distance methods
cumulative_pred_results <- cluster_pred_out %>% 
  unnest() %>%
  rename(dist_method = model) %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace),
         dist_method_abbr = factor(abbreviate(dist_method), levels = dist_method_abbr)) %>%
  left_join(., select(loeo_accuracy, trait, environment, max_accuracy = base), by = c("environment", "trait"))

# Split by trait
cumulative_pred_results_split <- cumulative_pred_results %>% 
  split(.$trait) %>% 
  map(~mutate(., environment = factor(environment, levels = unique(.$environment[order(.$max_accuracy, decreasing = TRUE)]))))

# Plot
g_plotlist <- cumulative_pred_results_split %>%
  map(~ggplot(data = ., aes(x = n_e, y = accuracy, color = dist_method_abbr)) + 
        geom_hline(aes(yintercept = max_accuracy, lty = "Accuracy\nUsing\nAll Data")) +
        # geom_point() + 
        geom_line() + 
        scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
        scale_x_continuous(breaks = pretty) +
        scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
        facet_wrap(~ environment, ncol = 5) + 
        xlab("Number of training environments") +
        ylab("Prediction accuracy") +
        theme_minimal() )

# Save
for (i in seq_along(g_plotlist)) {
  ggsave(filename = paste0("cumulative_environment_prediction_", names(g_plotlist)[i], ".jpg"), plot = g_plotlist[[i]],
         height = 4, width = 8, dpi = 1000, path = fig_dir)
}



## For each trait and environment, calculate the difference between the accuracy
## using the nth training environment and the max accuracy
## Then summarize for each trait
cumulative_pred_diff <- cumulative_pred_results %>% 
  mutate(diff_accuracy = accuracy - max_accuracy) %>%
  group_by(trait, dist_method_abbr, n_e) %>% 
  summarize_at(vars(diff_accuracy), funs(mean, sd, n())) %>%
  ungroup() %>%
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)

# Plot - no random lines
g_diff_accuracy_norandom <- cumulative_pred_diff %>% 
  ggplot(aes(x = n_e, y = mean, color = dist_method_abbr)) + 
  geom_hline(aes(yintercept = 0, lty = "Accuracy\nUsing\nAll Data")) +
  geom_line() + 
  # geom_ribbon(aes(ymin = lower, ymax = upper, fill = dist_method_abbr), alpha = 0.15) +
  scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
  scale_x_continuous(breaks = pretty) +
  scale_y_continuous(breaks = pretty) +
  facet_grid(trait ~ ., switch = "y", scales = "free_y") + 
  xlab("Number of training environments") +
  ylab("Prediction accuracy (relative to using all data)") +
  theme_minimal()

# Save
ggsave(filename = "cumulative_environment_prediction_relative_norandom.jpg", plot = g_diff_accuracy_norandom,
       height = 6, width = 4, dpi = 1000, path = fig_dir)


## Choose an arbitrary number of environments and show the average accuracy advantage
cumulative_pred_diff_example <- cumulative_pred_diff %>% 
  filter(n_e %in% c(3, 5, 10)) %>% 
  select(trait:n_e, mean) %>% 
  mutate(mean = round(mean, 3)) %>% 
  spread(dist_method_abbr, mean) %>% 
  arrange(n_e)

# Save as an image
t_cumulative_pred_diff <- grid.arrange(tableGrob(cumulative_pred_diff_example, rows = NULL, theme = ttheme_minimal()))

ggsave(filename = "cumulative_environment_prediction_example.jpg", plot = t_cumulative_pred_diff, path = fig_dir,
       height = 5, width = 7, dpi = 1000)






## Window predictions


# Rename the distance methods
window_pred_results <- cluster_pred_out_window %>% 
  unnest() %>%
  rename(dist_method = model) %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace),
         dist_method_abbr = factor(abbreviate(dist_method), levels = dist_method_abbr)) %>%
  left_join(., select(loeo_accuracy, trait, environment, max_accuracy = base), by = c("environment", "trait"))

# Split by trait
window_pred_results_split <- window_pred_results %>% 
  split(.$trait) %>% 
  map(~mutate(., environment = factor(environment, levels = unique(.$environment[order(.$max_accuracy, decreasing = TRUE)]))))

# Plot
g_plotlist <- window_pred_results_split %>%
  map(~ggplot(data = ., aes(x = window, y = accuracy, color = dist_method_abbr)) + 
        geom_hline(aes(yintercept = max_accuracy, lty = "Accuracy\nUsing\nAll Data")) +
        # geom_point() + 
        geom_line() + 
        scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
        scale_x_continuous(breaks = pretty) +
        scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
        facet_wrap(~ environment, ncol = 5) + 
        xlab("Sliding window number") +
        ylab("Prediction accuracy") +
        theme_minimal() )

# Save
for (i in seq_along(g_plotlist)) {
  ggsave(filename = paste0("window_environment_prediction_", names(g_plotlist)[i], ".jpg"), plot = g_plotlist[[i]],
         height = 4, width = 8, dpi = 1000, path = fig_dir)
}



# For each trait and environment, calculate the difference between the accuracy
## using the nth training environment and the max accuracy
## Then summarize for each trait
window_pred_diff <- window_pred_results %>% 
  mutate(diff_accuracy = accuracy - max_accuracy) %>%
  group_by(trait, dist_method_abbr, window) %>% 
  summarize_at(vars(diff_accuracy), funs(mean, sd, n())) %>%
  ungroup() %>%
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)

# Plot - no random lines
g_diff_accuracy_norandom <- window_pred_diff %>% 
  ggplot(aes(x = window, y = mean, color = dist_method_abbr)) + 
  geom_hline(aes(yintercept = 0, lty = "Accuracy\nUsing\nAll Data")) +
  geom_line() +
  # geom_point(size = 0.5) + geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  # geom_ribbon(aes(ymin = lower, ymax = upper, fill = dist_method_abbr), alpha = 0.15) +
  scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
  scale_x_continuous(breaks = pretty) +
  scale_y_continuous(breaks = pretty) +
  facet_grid(trait ~ ., switch = "y", scales = "free_y") + 
  xlab("Sliding window number") +
  ylab("Prediction accuracy (relative to using all data)") +
  theme_minimal()

# Save
ggsave(filename = "window_environment_prediction_relative_norandom.jpg", plot = g_diff_accuracy_norandom,
       height = 6, width = 4, dpi = 1000, path = fig_dir)













#### Environment covariance matrix prediction ####

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












#### Environmental cluster predictions ####











