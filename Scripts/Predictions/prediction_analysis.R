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
                         "OYEC_Top5EC_Cor" = "One Year EC Top 5", "MYEC_AllEC" = "Ten Year ECs", "MYEC_Top5EC_Cor" = "Ten Year EC Top 5",
                         "sample" = "Random")
dist_method_abbr <- abbreviate(dist_method_replace)

colors <- umn_palette(3)
dist_colors <- setNames(c(colors[c(1:2, 3, 8, 4, 9)], "grey75"), dist_method_abbr)










#### Leave-one-environment-out predictions ####


# Calculate accuracy
# loeo_accuracy <- environment_loeo_predictions_mean %>% 
loeo_accuracy <- environment_loeo_predictions_geno_means %>%
  filter(trait %in% traits) %>%
  mutate(results = map(predictions, "boot")) %>% 
  unnest(results) %>%
  rename(environment = testEnv)

# Calculate the mean per trait
(loeo_mean <- loeo_accuracy %>% 
  group_by(trait) %>% 
  summarize(accuracy = mean(base)) %>%
  ungroup())

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
  unnest(out) %>%
  rename(dist_method = model) %>%
  mutate(iter = parse_number(dist_method),
         iter = ifelse(str_detect(dist_method, "sample"), iter, 1),
         dist_method = ifelse(str_detect(dist_method, "sample"), "sample", dist_method),
         dist_method = str_replace_all(dist_method, dist_method_replace),
         dist_method_abbr = factor(abbreviate(dist_method), levels = dist_method_abbr)) %>%
  left_join(., select(loeo_accuracy, trait, environment, max_accuracy = base), by = c("environment", "trait")) %>%
  group_by(environment, trait, dist_method, iter) %>%
  mutate(max_accuracy_adding = accuracy[n_e == max(n_e)]) %>%
  ungroup()


## Summarize the accuracies and create CIs for the random method
cumulative_pred_results_summ <- cumulative_pred_results %>%
  group_by(environment, trait, dist_method, dist_method_abbr, n_e, max_accuracy, max_accuracy_adding) %>% 
  summarize_at(vars(accuracy), funs(mean, sd, n())) %>% 
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat) %>%
  ungroup()


# Split by trait
cumulative_pred_results_split <- cumulative_pred_results_summ %>% 
  split(.$trait) %>% 
  map(~mutate(., environment = factor(environment, levels = unique(.$environment[order(.$max_accuracy_adding, decreasing = TRUE)]))))

# Plot
g_plotlist <- cumulative_pred_results_split %>%
  map(~ggplot(data = ., aes(x = n_e, y = mean, color = dist_method_abbr)) + 
        geom_hline(aes(yintercept = max_accuracy_adding, lty = "Accuracy\nUsing\nAll Data")) +
        # geom_point() + 
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
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
         height = 10, width = 8, dpi = 1000, path = fig_dir)
}



## For each trait and environment, calculate the difference between the accuracy
## using the nth training environment and the max accuracy
## Then summarize for each trait
cumulative_pred_diff <- cumulative_pred_results %>% 
  mutate(diff_accuracy = accuracy - max_accuracy_adding) %>%
  group_by(trait, dist_method_abbr, n_e, iter) %>% 
  summarize(diff_accuracy = mean(diff_accuracy)) %>% ## Take the mean over all environments for a method/iteration
  summarize_at(vars(diff_accuracy), funs(mean, sd, n())) %>%
  ungroup() %>%
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)

# Plot - Random
g_diff_accuracy <- cumulative_pred_diff %>% 
  ggplot(aes(x = n_e, y = mean, color = dist_method_abbr)) + 
  geom_hline(aes(yintercept = 0, lty = "Accuracy\nUsing\nAll Data")) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
  scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
  scale_x_continuous(breaks = pretty) +
  scale_y_continuous(breaks = pretty) +
  facet_grid(trait ~ ., switch = "y", scales = "free_y") + 
  xlab("Number of training environments") +
  ylab("Prediction accuracy (relative to using all data)") +
  theme_presentation() +
  theme(strip.placement = "outside")

# Save
ggsave(filename = "cumulative_environment_prediction_relative.jpg", plot = g_diff_accuracy,
       height = 8, width = 7, dpi = 1000, path = fig_dir)


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


## For each distance method, find the average number of environments in which the accuracy is maximized
cumulative_pred_nE <- cumulative_pred_results %>% 
  group_by(trait, dist_method_abbr, iter, environment) %>%
  top_n(x = ., n = 1, wt = accuracy) %>% summarize(n_e = mean(n_e))

cumulative_pred_nE_summ <- cumulative_pred_nE %>%
  summarize(n_e = mean(n_e)) %>%
  summarize_at(vars(n_e), funs(mean, sd, n())) %>%
  ungroup() %>%
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)

# Plot
g_n_train <- cumulative_pred_nE_summ %>% 
  ggplot(aes(x = trait, y = mean, color = dist_method_abbr)) +
  geom_point(position = position_dodge(0.6), size = 3) +
  geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6)) +
  scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_y_continuous(breaks = pretty) +
  ylab("Peak accuracy training environments") +
  theme_presentation() +
  theme(axis.title.x = element_blank())
  
ggsave(filename = "cumulative_max_accuracy_nE.jpg", plot = g_n_train, path = fig_dir,
       height = 5, width = 8, dpi = 1000)


## Using the phenotypic distance measure, calculate the average distance at which
## the maximum prediction accuracy is achieved

# First calculate the mean distance for each training set
average_pheno_distance <- pred_env_dist_rank$tp %>% 
  filter(model == "pheno_dist") %>% 
  mutate(env_rank = map(env_rank, ~unlist(.) %>% data_frame(pred_environment = names(.), distance = .) %>% mutate(n_e = seq(nrow(.))))) %>% 
  unnest() %>%
  group_by(trait, environment) %>% 
  mutate(distance = cummean(distance)) %>%
  ungroup() %>%
  select(-model)


cumulative_pred_nE_summ_avg <- cumulative_pred_nE %>% 
  left_join(., average_pheno_distance) %>%
  summarize(distance = mean(distance)) %>%
  summarize_at(vars(distance), funs(mean, sd, n())) %>%
  ungroup() %>%
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)


# Plot
g_dist<- cumulative_pred_nE_summ_avg %>% 
  ggplot(aes(x = trait, y = mean, color = dist_method_abbr)) +
  geom_point(position = position_dodge(0.6), size = 3) +
  geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6)) +
  scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_y_continuous(breaks = function(x) pretty(x, n = 3)) +
  ylab("Peak accuracy training environments") +
  facet_wrap(~trait, scales = "free", nrow = 1) +
  theme_presentation() +
  theme(axis.title.x = element_blank(), strip.text = element_blank())

ggsave(filename = "cumulative_max_accuracy_pheno_dist.jpg", plot = g_dist, path = fig_dir,
       height = 5, width = 10, dpi = 1000)











## Window predictions


# Rename the distance methods
window_pred_results <- cluster_pred_out_window %>% 
  mutate(out = map(out, ~mutate(., window = seq(nrow(.))))) %>%
  unnest(out) %>%
  rename(dist_method = model) %>%
  mutate(iter = parse_number(dist_method),
         iter = ifelse(str_detect(dist_method, "sample"), iter, 1),
         dist_method = ifelse(str_detect(dist_method, "sample"), "sample", dist_method),
         dist_method = str_replace_all(dist_method, dist_method_replace),
         dist_method_abbr = factor(abbreviate(dist_method), levels = dist_method_abbr)) %>%
  left_join(., select(loeo_accuracy, trait, environment, max_accuracy = base), by = c("environment", "trait")) 


## Summarize the accuracies and create CIs for the random method
window_pred_results_summ <- window_pred_results %>%
  group_by(environment, trait, dist_method, dist_method_abbr, window, max_accuracy) %>% 
  summarize_at(vars(accuracy), funs(mean, sd, n())) %>% 
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat) %>%
  ungroup()


# Split by trait
window_pred_results_split <- window_pred_results_summ %>% 
  split(.$trait) %>% 
  map(~mutate(., environment = factor(environment, levels = unique(.$environment[order(.$max_accuracy, decreasing = TRUE)]))))

# Plot
g_plotlist <- window_pred_results_split %>%
  map(~ggplot(data = ., aes(x = window, y = mean, color = dist_method_abbr)) + 
        geom_hline(aes(yintercept = max_accuracy, lty = "Accuracy\nUsing\nAll Data")) +
        # geom_point() + 
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
        geom_line() + 
        scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
        scale_x_continuous(breaks = pretty) +
        scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
        facet_wrap(~ environment, ncol = 5) + 
        xlab("Window number") +
        ylab("Prediction accuracy") +
        theme_minimal() )

# Save
for (i in seq_along(g_plotlist)) {
  ggsave(filename = paste0("window_environment_prediction_", names(g_plotlist)[i], ".jpg"), plot = g_plotlist[[i]],
         height = 10, width = 8, dpi = 1000, path = fig_dir)
}






# For each trait and environment, calculate the difference between the accuracy
## using the nth training environment and the max accuracy
## Then summarize for each trait
window_pred_diff <- window_pred_results %>% 
  mutate(diff_accuracy = accuracy - max_accuracy) %>%
  group_by(trait, dist_method_abbr, window, iter) %>% 
  summarize(diff_accuracy = mean(diff_accuracy)) %>% ## Take the mean over all environments for a method/iteration
  summarize_at(vars(diff_accuracy), funs(mean, sd, n())) %>%
  ungroup() %>%
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)

# Plot - Random
g_diff_accuracy <- window_pred_diff %>% 
  ggplot(aes(x = window, y = mean, color = dist_method_abbr)) + 
  geom_hline(aes(yintercept = 0, lty = "Accuracy\nUsing\nAll Data")) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
  scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
  scale_x_continuous(breaks = pretty) +
  scale_y_continuous(breaks = pretty) +
  facet_grid(trait ~ ., switch = "y", scales = "free_y") + 
  xlab("Window number") +
  ylab("Prediction accuracy (relative to using all data)") +
  theme_presentation() +
  theme(strip.placement = "outside")

# Save
ggsave(filename = "window_environment_prediction_relative.jpg", plot = g_diff_accuracy,
       height = 8, width = 7, dpi = 1000, path = fig_dir)


## Fit smoothing lines
g_diff_accuracy_smooth <- window_pred_diff %>% 
  ggplot(aes(x = window, y = mean, color = dist_method_abbr)) + 
  geom_hline(aes(yintercept = 0, lty = "Accuracy\nUsing\nAll Data")) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
  scale_x_continuous(breaks = pretty) +
  scale_y_continuous(breaks = pretty) +
  facet_grid(trait ~ ., switch = "y", scales = "free_y") + 
  xlab("Window number") +
  ylab("Prediction accuracy (relative to using all data)") +
  theme_presentation() +
  theme(strip.placement = "outside")

# Save
ggsave(filename = "window_environment_prediction_relative_smooth.jpg", plot = g_diff_accuracy_smooth,
       height = 8, width = 7, dpi = 1000, path = fig_dir)



## For each distance method, find the average window number in which the accuracy is maximized
window_pred_wn <- window_pred_results %>% 
  group_by(trait, dist_method_abbr, iter, environment) %>%
  top_n(x = ., n = 1, wt = accuracy) %>% summarize(window = mean(window))

window_pred_wn_summ <- window_pred_wn %>%
  summarize(window = mean(window)) %>%
  summarize_at(vars(window), funs(mean, sd, n())) %>%
  ungroup() %>%
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)

# Plot
g_window_train <- window_pred_wn_summ %>% 
  ggplot(aes(x = trait, y = mean, color = dist_method_abbr)) +
  geom_point(position = position_dodge(0.6), size = 3) +
  geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6)) +
  scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_y_continuous(breaks = pretty) +
  ylab("Peak accuracy training environments") +
  theme_presentation() +
  theme(axis.title.x = element_blank())

ggsave(filename = "window_max_accuracy_wn.jpg", plot = g_window_train, path = fig_dir,
       height = 5, width = 8, dpi = 1000)

















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

# Load the cluster results
load(file.path(result_dir, "cluster_predictions.RData"))


## Unnest and tidy
cluster_predictions_out <- cluster_predictions %>% 
  # unnest(out) %>%
  filter(model %in% names(dist_method_replace)) %>%
  filter(trait %in% traits) %>%
  mutate(dist_method = str_replace_all(model, dist_method_replace),
         dist_method_abbr = factor(str_replace_all(dist_method, dist_method_abbr), levels = dist_method_abbr)) %>%
  mutate_at(vars(min_env, cluster, environment), as.factor)


# Fit a model per trait
cluster_predictions_model <- cluster_predictions_out %>%
  group_by(trait) %>%
  do(fit = lm(accuracy ~ environment + dist_method_abbr + min_env + min_env:dist_method_abbr + cluster:dist_method_abbr, data = .)) %>%
  ungroup() %>%
  mutate(dist_method_effect = map(fit, ~as.data.frame(effects::Effect(focal.predictors = "dist_method_abbr", .))),
         dist_method_min_env_effect = map(fit, ~as.data.frame(effects::allEffects(mod = .)$`dist_method_abbr:min_env`)))

## Anovas
cluster_predictions_model$fit %>% map(anova)





## Plot
g_cluster_pred <- cluster_predictions_model %>% 
  unnest(dist_method_effect) %>% 
  rename(model = dist_method_abbr) %>%
  mutate(model = factor(model, levels = dist_method_abbr)) %>%
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) +
  geom_hline(data = loeo_mean, aes(yintercept = accuracy, lty = "Accuracy using\nall data"), color = "grey75") +
  geom_pointrange() +
  scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty) + 
  facet_grid(~ trait) +
  ylab("Prediction accuracy") + 
  theme_presentation2(base_size = 16) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.box.margin = margin())
  
ggsave(filename = "cluster_prediction_models.jpg", plot = g_cluster_pred, path = fig_dir, width = 7, height = 5, dpi = 1000)



## Plot distance by minimum number of environments
g_cluster_pred_minenv <- cluster_predictions_model %>% 
  unnest(dist_method_min_env_effect) %>% 
  rename(model = dist_method_abbr) %>%
  mutate(model = factor(model, levels = dist_method_abbr), min_env = factor(min_env, levels = c(3, 5, 10))) %>%
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) +
  geom_hline(data = loeo_mean, aes(yintercept = accuracy, lty = "Accuracy using\nall data"), color = "grey75") +
  geom_pointrange() +
  scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_linetype_manual(values = 2, name = NULL) +
  # scale_y_continuous(breaks = pretty) + 
  facet_grid(trait ~ min_env, scales = "free_y") +
  ylab("Prediction accuracy") + 
  theme_presentation2(base_size = 16) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.box.margin = margin())

ggsave(filename = "cluster_prediction_models_minenv.jpg", plot = g_cluster_pred_minenv, path = fig_dir, width = 7, height = 7, dpi = 1000)









