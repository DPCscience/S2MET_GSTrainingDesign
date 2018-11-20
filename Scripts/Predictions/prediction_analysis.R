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
library(effects) # For ls means / marginal means
library(ggforce)
library(ggridges)
library(gridExtra)

# Load the environmental distance df
load(file.path(result_dir, "distance_method_results.RData"))
# Load the LOEO prediction results
load(file.path(result_dir, "all_data_environmental_predictions.RData"))

# Load data
load(file.path(result_dir, "distance_rank_predictions.RData"))

## Significant level
alpha <- 0.05


## Unique models
unique(cluster_pred_out$model)


## Create a factor for the distance methods
dist_method_replace <- c("pheno_dist" = "Phenotypic Distance", "pheno_loc_dist" = "Location Phenotypic Distance",  "great_circle_dist" = "Great Circle Distance", 
                         "OYEC_All" = "One Year All ECs", "OYEC_Mean" = "One Year Mean Cor EC", "OYEC_IPCA" = "One Year IPCA Cor EC", 
                         "MYEC_All" = "Multi Year All ECs", "MYEC_Mean" = "Multi Year Mean Cor EC", "MYEC_IPCA" = "Multi Year IPCA Cor EC",
                         "sample" = "Random")
dist_method_abbr <- abbreviate(dist_method_replace)

colors <- umn_palette(3)
colors_use <- c(colors[2], "#FFDE7A", colors[c(1, 3, 8, 4, 9, 5, 10)], "grey75")
dist_colors <- setNames(colors_use, dist_method_abbr)










#### Leave-one-environment-out predictions ####


# Calculate accuracy
# loeo_accuracy <- environment_loeo_predictions_mean %>% 
loeo_accuracy <- environment_loeo_predictions_geno_means %>%
  filter(trait %in% traits) %>%
  mutate(results = purrr::map(predictions, "boot")) %>% 
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



## Predictions when adding one environment at a time

# Convert the accuracies to z scores
# Adjust the names of the distance models
cluster_pred_out1 <- cluster_pred_out %>%
  unnest(out) %>% # Comment this if the data was generated on a local machine
  mutate(zscore = ztrans(accuracy),
         nEnv = as.factor(n_e),
         environment = as.factor(validation_environment),
         model = ifelse(str_detect(model, "sample"), "sample", model),
         model = abbreviate(str_replace_all(model, dist_method_replace)),
         model = factor(model, levels = dist_method_abbr)) %>%
  group_by(trait, set, environment, model, nEnv, n_e) %>% 
  summarize(zscore = mean(zscore)) %>%
  ungroup()

## Fit a model per triat
## The accuracy is a function of:
## 1. Environment
## 2. Model
## 3. Number of training environments
## 4. Interactions

cluster_pred_fit <- cluster_pred_out1 %>%
  # filter(model != "Rndm") %>%
  group_by(set, trait) %>% 
  filter(n_e < max(n_e)) %>% droplevels() %>% # Remove the result when using all data, since no difference will exist among models
  do(fit = lm(zscore ~ environment + model + nEnv + environment:model + environment:nEnv + model:nEnv, data = .)) %>%
  ungroup() %>%
  mutate(model_effects = map(fit, ~Effect(focal.predictors = "model", .)),
         environment_effects = map(fit, ~Effect(focal.predictors = "environment", .)),
         nEnv_effects = map(fit, ~Effect(focal.predictors = "nEnv", .)),
         model_nEnv_effects = map(fit, ~Effect(focal.predictors = c("model", "nEnv"), .)),
         model_environment_effects = map(fit, ~Effect(focal.predictors = c("model", "environment"), .)))
     

# Quick anova scan
cluster_pred_fit$fit %>% map(anova)

## Output the variance proportion table
cluster_pred_fit %>% 
  mutate(anova = map(fit, ~anova(.) %>% broom::tidy())) %>% 
  unnest(anova) %>% 
  group_by(set, trait) %>% 
  mutate(varprop = sumsq / sum(sumsq)) %>%
  select(trait, term, varprop) %>% 
  spread(term, varprop)


## Collect the effects
cluster_pred_fit_effects <- cluster_pred_fit %>% 
  filter(! str_detect(trait, "AGDD")) %>%
  mutate_at(vars(contains("effects")), ~map(., as.data.frame))






## Fit a model when using data from all environments
cluster_pred_fit_effects_all <- cluster_pred_out1 %>%
  group_by(set, trait) %>% 
  filter(n_e == max(n_e)) %>% droplevels() %>% # Remove the result when using all data, since no difference will exist among models
  do(fit = lm(zscore ~ environment, data = ., contrasts = list(environment = "contr.sum"))) %>%
  ungroup() %>%
  mutate(environment_effects = map(fit, ~Effect(focal.predictors = "environment", .) %>% as.data.frame),
         mean_effects = zexp(map_dbl(fit, ~coef(.)[1])))




## Plot the effect of model
g_model_effect <- cluster_pred_fit_effects %>% 
  filter(set == "complete") %>%
  unnest(model_effects) %>% 
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(model = factor(model, levels = dist_method_abbr)) %>%
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) + 
  geom_point() + 
  geom_errorbar(width = 0.5) +
  facet_wrap(~trait, scales = "free") +
  scale_color_manual(values = dist_colors) +
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  theme_presentation2() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "bottom")

# Save
ggsave(filename = "cumulative_pred_model_effect.jpg", plot = g_model_effect, path = fig_dir, height = 5, width = 8, dpi = 1000)



## Plot the effect of model at an early number of training environments
## Highlight 1, 5, and 10 environments
nEnv_select <- c("1", "5", "10")
  
g_model_effect_nenv <-  cluster_pred_fit_effects %>% 
  filter(set == "complete") %>%
  left_join(., select(cluster_pred_fit_effects_all, set, trait, mean_effects)) %>% # Add the accuracy when adding using all environments
  unnest(model_nEnv_effects) %>% 
  filter(nEnv %in% nEnv_select) %>%
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(model = factor(model, levels = dist_method_abbr),
         nEnv = factor(nEnv, levels = nEnv_select)) %>%
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) + 
  geom_hline(aes(yintercept = mean_effects, lty = "Accuracy using\nall data" )) +
  geom_point() + 
  geom_errorbar(width = 0.5) +
  facet_grid(trait ~ nEnv, scales = "free") +
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 3)) +
  scale_linetype_manual(values = 2, name = NULL) + 
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  theme_presentation2() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "bottom")

# Save
ggsave(filename = "cumulative_pred_model_effect_envselect.jpg", plot = g_model_effect_nenv, path = fig_dir, height = 8, width = 8, dpi = 1000)





## Plot the effect of number of environment
g_nenv_effect <- cluster_pred_fit_effects %>% 
  filter(set == "complete") %>%
  unnest(nEnv_effects) %>%
  filter(trait == "PlantHeight") %>%
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(nEnv = parse_number(nEnv)) %>%
  ggplot(aes(x = nEnv, y = fit, ymin = lower, ymax = upper)) + 
  geom_point() + 
  geom_line() +
  # geom_errorbar(width = 0.5) +
  geom_ribbon(color = "grey75", alpha = 0.2) +
  facet_wrap(~trait, scales = "free") +
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  xlab("Number of training environments") +
  theme_presentation2()

# Save
ggsave(filename = "cumulative_pred_nEnv_effect.jpg", plot = g_nenv_effect, path = fig_dir, height = 5, width = 8, dpi = 1000)


## Plot the interaction effect of model and number of environment
g_model_nenv_effect <- cluster_pred_fit_effects %>% 
  filter(set == "complete") %>%
  unnest(model_nEnv_effects) %>% 
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(nEnv = parse_number(nEnv),
         model = factor(model, levels = dist_method_abbr)) %>%
  left_join(., select(cluster_pred_fit_effects_all, set, trait, mean_effects)) %>%
  ggplot(aes(x = nEnv, y = fit, ymin = lower, ymax = upper, color = model, fill = model)) + 
  geom_hline(aes(yintercept = mean_effects, lty = "Accuracy using\nall data")) +
  geom_point() + 
  geom_line() +
  # geom_errorbar(width = 0.5) +
  # geom_ribbon(color = "grey75", alpha = 0.2) +
  facet_wrap(~trait, scales = "free") +
  scale_color_manual(values = dist_colors, name = "Model") +
  scale_fill_manual(values = dist_colors, name = "Model") +
  scale_y_continuous(breaks = pretty) +
  scale_linetype_manual(values = 2, name = NULL) +
  ylab("Accuracy") +
  xlab("Number of training environments") +
  theme_presentation2() +
  theme(legend.position = "bottom")

# Save
ggsave(filename = "cumulative_pred_model_nEnv_effect.jpg", plot = g_model_nenv_effect, path = fig_dir, height = 5, width = 8, dpi = 1000)





## Realistic results


## Plot the effect of model
g_model_effect <- cluster_pred_fit_effects %>% 
  filter(set == "realistic") %>%
  unnest(model_effects) %>% 
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(model = factor(model, levels = dist_method_abbr)) %>%
  # left_join(., select(cluster_pred_fit_effects_all, set, trait, mean_effects)) %>% # Add the accuracy when adding using all environments
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) + 
  # geom_hline(aes(yintercept = mean_effects, lty = "Accuracy using\nall data" )) +
  geom_point() + 
  geom_errorbar(width = 0.5) +
  facet_wrap(~trait) +
  scale_color_manual(values = dist_colors, name = "Model") +
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  theme_presentation2() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "bottom")

# Save
ggsave(filename = "cumulative_pred_model_effect_realistic.jpg", plot = g_model_effect, path = fig_dir, height = 5, width = 8, dpi = 1000)



## Plot the effect of model at an early number of training environments
## Highlight 1, 5, and 10 environments
nEnv_select <- c("1", "5", "10")

g_model_effect_nenv <-  cluster_pred_fit_effects %>% 
  filter(set == "realistic") %>%
  left_join(., select(cluster_pred_fit_effects_all, set, trait, mean_effects)) %>% # Add the accuracy when adding using all environments
  unnest(model_nEnv_effects) %>% 
  filter(nEnv %in% nEnv_select) %>%
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(model = factor(model, levels = dist_method_abbr),
         nEnv = factor(nEnv, levels = nEnv_select)) %>%
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) + 
  geom_hline(aes(yintercept = mean_effects, lty = "Accuracy using\nall data" )) +
  geom_point() + 
  geom_errorbar(width = 0.5) +
  facet_grid(trait ~ nEnv) +
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 3)) +
  scale_linetype_manual(values = 2, name = NULL) + 
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  theme_presentation2() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "bottom")

# Save
ggsave(filename = "cumulative_pred_model_effect_envselect_realistic.jpg", plot = g_model_effect_nenv, path = fig_dir, height = 8, width = 8, dpi = 1000)





## Plot the effect of number of environment
g_nenv_effect <- cluster_pred_fit_effects %>%
  filter(set == "realistic") %>%
  unnest(nEnv_effects) %>% 
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(nEnv = parse_number(nEnv)) %>%
  ggplot(aes(x = nEnv, y = fit, ymin = lower, ymax = upper)) + 
  geom_point() + 
  geom_line() +
  # geom_errorbar(width = 0.5) +
  geom_ribbon(color = "grey75", alpha = 0.2) +
  facet_wrap(~trait, scales = "free") +
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  xlab("Number of training environments") +
  theme_presentation2()

# Save
ggsave(filename = "cumulative_pred_nEnv_effect_realistic.jpg", plot = g_nenv_effect, path = fig_dir, height = 5, width = 8, dpi = 1000)


## Plot the interaction effect of model and number of environment
g_model_nenv_effect <- cluster_pred_fit_effects %>% 
  filter(set == "realistic") %>%
  unnest(model_nEnv_effects) %>% 
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(nEnv = parse_number(nEnv)) %>%
  ggplot(aes(x = nEnv, y = fit, ymin = lower, ymax = upper, color = model, fill = model)) + 
  geom_point() + 
  geom_line() +
  # geom_errorbar(width = 0.5) +
  # geom_ribbon(color = "grey75", alpha = 0.2) +
  facet_wrap(~trait) +
  scale_color_manual(values = dist_colors, name = "Model") +
  scale_fill_manual(values = dist_colors, name = "Model") +
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  xlab("Number of training environments") +
  theme_presentation2() +
  theme(legend.position = "bottom")

# Save
ggsave(filename = "cumulative_pred_model_nEnv_effect_realistic.jpg", plot = g_model_nenv_effect, path = fig_dir, height = 5, width = 8, dpi = 1000)


## Plot the model x environment means
g_model_env_effect <- cluster_pred_fit_effects %>% 
  filter(set == "realistic") %>%
  unnest(model_environment_effects) %>% 
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  ggplot(aes(x = environment, y = fit, ymin = lower, ymax = upper, color = model)) + 
  geom_point(position = position_dodge(0.9)) + 
  geom_errorbar(position = position_dodge(0.9), width = 0.5) +
  facet_wrap(~trait, scales = "free") +
  scale_color_manual(values = dist_colors, name = "Model") +
  scale_fill_manual(values = dist_colors, name = "Model") +
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  xlab("Environment") +
  theme_presentation2() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

# Save
ggsave(filename = "cumulative_pred_model_environment_effect_realistic.jpg", plot = g_model_env_effect, path = fig_dir,
       height = 5, width = 10, dpi = 1000)















#### Window predictions

# Convert the accuracies to z scores
cluster_pred_out_window1 <- cluster_pred_out_window %>% 
  unnest(out) %>%
  group_by(trait, set, validation_environment, model) %>% mutate(window = seq(n())) %>% ungroup() %>%
  mutate(model = ifelse(str_detect(model, "sample"), "sample", model),
         zscore = ztrans(accuracy),
         # window = as.factor(window), # As a numeric, estimate a regression slope
         environment = as.factor(validation_environment),
         model = abbreviate(str_replace_all(model, dist_method_replace)),
         model = factor(model, levels = dist_method_abbr)) %>%
  group_by(trait, set, environment, model, window) %>%
  summarize(zscore = mean(zscore)) %>%
  ungroup()

## Fit a model per triat
## The accuracy is a function of:
## 1. Environment
## 2. Model
## 3. Number of training environments
## 4. Interactions

cluster_pred_window_fit <- cluster_pred_out_window1 %>%
  group_by(set, trait) %>% 
  do(fit = lm(zscore ~ environment + model + window + environment:window + model:window, data = .)) %>%
  ungroup() %>%
  mutate(model_effects = map(fit, ~Effect(focal.predictors = "model", .)),
         environment_effects = map(fit, ~Effect(focal.predictors = "environment", .)),
         window_effects = map(fit, ~Effect(focal.predictors = "window", .)),
         model_window_effects = map(fit, ~Effect(focal.predictors = c("model", "window"), .)),
         model_environment_window_effects = map(fit, ~Effect(focal.predictors = c("environment", "model", "window"), .)))

# Quick anova scan
cluster_pred_window_fit$fit %>% map(anova)

## Output the variance proportion table
cluster_pred_window_fit %>% 
  mutate(anova = map(fit, ~anova(.) %>% broom::tidy())) %>% 
  unnest(anova) %>% 
  group_by(set, trait) %>% 
  mutate(varprop = sumsq / sum(sumsq)) %>%
  select(trait, term, varprop) %>% 
  spread(term, varprop)

# set       trait           environment `environment:window`    model `model:window` Residuals   window
# 1 complete  GrainYield            0.821              0.00904 0.000179        0.00302    0.160  0.00678 
# 2 complete  HeadingDate           0.903              0.00716 0.000438        0.00280    0.0864 0.000525
# 3 complete  HeadingDateAGDD       0.903              0.00640 0.000618        0.00183    0.0874 0.000723
# 4 complete  PlantHeight           0.763              0.0231  0.000495        0.00850    0.201  0.00454 
# 5 realistic GrainYield            0.662              0.0727  0.000818        0.00783    0.253  0.00359 
# 6 realistic HeadingDate           0.890              0.0136  0.000900        0.00790    0.0874 0.000130
# 7 realistic HeadingDateAGDD       0.898              0.0144  0.000550        0.00599    0.0793 0.00162 
# 8 realistic PlantHeight           0.765              0.0438  0.00156         0.0129     0.166  0.0107




## Collect the effects
cluster_pred_window_fit_effects <- cluster_pred_window_fit %>% 
  filter(! str_detect(trait, "AGDD")) %>%
  mutate_at(vars(contains("effects")), ~map(., as.data.frame))







## Plot the effect of model at an early number of training environments
g_model_effect1 <- cluster_pred_window_fit_effects %>% 
  filter(set == "complete") %>%
  left_join(., select(cluster_pred_fit_effects_all, set, trait, mean_effects)) %>% # Add results when using all data
  unnest(model_window_effects) %>% 
  filter(window == 1) %>%
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(model = factor(model, levels = dist_method_abbr)) %>%
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) + 
  geom_hline(aes(yintercept = mean_effects, lty = "Accuracy using\nall data")) +
  geom_point() + 
  geom_errorbar(width = 0.5) +
  facet_wrap(~trait) +
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 3)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  theme_presentation2() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "bottom")

# Save
ggsave(filename = "window_pred_model_effect_window1.jpg", plot = g_model_effect1, path = fig_dir, height = 5, width = 8, dpi = 1000)




## Plot the effect of window
g_window_effect <- cluster_pred_window_fit_effects %>% 
  filter(set == "complete") %>%
  unnest(window_effects) %>% 
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  ggplot(aes(x = window, y = fit, ymin = lower, ymax = upper)) + 
  # geom_point() + 
  geom_line() +
  # geom_errorbar(width = 0.5) +
  geom_ribbon(color = "grey75", alpha = 0.2) +
  facet_wrap(~trait) +
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  xlab("Sliding window interval") +
  theme_presentation2()

# Save
ggsave(filename = "window_pred_window_effect.jpg", plot = g_window_effect, path = fig_dir, height = 5, width = 8, dpi = 1000)

## Overall, what is the percent change in prediction accuracy over the sliding windows?
cluster_pred_window_fit_effects %>% 
  filter(set == "complete") %>%
  unnest(window_effects) %>% 
  filter(window %in% c(min(window), max(window))) %>% 
  group_by(trait) %>%
  summarize(acc_per_change = (fit[1] - fit[2]) / fit[1])


# trait       acc_per_change
# 1 GrainYield          0.156 
# 2 HeadingDate         0.0314
# 3 PlantHeight         0.0820


## Plot the interaction effect of model and number of environment
g_model_window_effect <- cluster_pred_window_fit_effects %>% 
  filter(set == "complete") %>%
  unnest(model_window_effects) %>% 
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(model = factor(model, levels = dist_method_abbr)) %>%
  ggplot(aes(x = window, y = fit, ymin = lower, ymax = upper, fill = model)) + 
  # geom_point() + 
  geom_line() +
  # geom_errorbar(width = 0.5) +
  geom_ribbon(color = "grey75", alpha = 0.2) +
  facet_wrap(~trait, scales = "free") +
  scale_color_manual(values = dist_colors, name = "Model") +
  scale_fill_manual(values = dist_colors, name = "Model") +
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  xlab("Sliding window interval") +
  theme_presentation2() +
  theme(legend.position = "bottom")

# Save
ggsave(filename = "window_pred_model_window_effect.jpg", plot = g_model_window_effect, path = fig_dir, height = 5, width = 8, dpi = 1000)


## Overall, what is the percent change in prediction accuracy over the sliding windows, for each model?
cluster_pred_window_fit_effects %>% 
  filter(set == "complete") %>%
  unnest(model_window_effects) %>% 
  filter(window %in% c(min(window), max(window))) %>% 
  group_by(trait, model) %>%
  summarize(acc_per_change = (fit[1] - fit[2]) / fit[1]) %>%
  spread(model, acc_per_change)

# trait          GrCD  LcPD  MYAE   MYICE   MYMCE    OYAE   OYICE   OYMCE   PhnD
# 1 GrainYield   0.205  0.276 0.106  0.174   0.114  0.198   -0.0700  0.158  0.340 
# 2 HeadingDate -0.0577 0.115 0.123 -0.0619 -0.0813 0.00948  0.0330  0.116  0.0822
# 3 PlantHeight  0.0877 0.154 0.106  0.107   0.0461 0.0641  -0.0662 -0.0738 0.333








#### Realistic results


## Plot the effect of model at an early number of training environments
g_model_effect1 <- cluster_pred_window_fit_effects %>% 
  filter(set == "realistic") %>%
  left_join(., select(cluster_pred_fit_effects_all, set, trait, mean_effects)) %>% # Add results when using all data
  unnest(model_window_effects) %>% 
  filter(window == 1) %>%
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(model = factor(model, levels = dist_method_abbr)) %>%
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) + 
  geom_hline(aes(yintercept = mean_effects, lty = "Accuracy using\nall data")) +
  geom_point() + 
  geom_errorbar(width = 0.5) +
  facet_wrap(~trait) +
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 2)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty, limits = c(0.25, 0.45)) +
  ylab("Accuracy") +
  theme_presentation2() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "bottom")

# Save
ggsave(filename = "window_pred_model_effect_window1_realistic.jpg", plot = g_model_effect1, path = fig_dir, height = 5, width = 8, dpi = 1000)




## Plot the effect of window
g_window_effect <- cluster_pred_window_fit_effects %>% 
  filter(set == "realistic") %>%
  unnest(window_effects) %>% 
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  ggplot(aes(x = window, y = fit, ymin = lower, ymax = upper)) + 
  # geom_point() + 
  geom_line() +
  # geom_errorbar(width = 0.5) +
  geom_ribbon(color = "grey75", alpha = 0.2) +
  facet_wrap(~trait) +
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  xlab("Sliding window interval") +
  theme_presentation2()

# Save
ggsave(filename = "window_pred_window_effect_realistic.jpg", plot = g_window_effect, path = fig_dir, height = 5, width = 8, dpi = 1000)

## Overall, what is the percent change in prediction accuracy over the sliding windows?
cluster_pred_window_fit_effects %>% 
  filter(set == "realistic") %>%
  unnest(window_effects) %>% 
  filter(window %in% c(min(window), max(window))) %>% 
  group_by(trait) %>%
  summarize(acc_per_change = (fit[1] - fit[2]) / fit[1])


# trait       acc_per_change
# 1 GrainYield          0.136 
# 2 HeadingDate         0.0295
# 3 PlantHeight         0.259 


## Plot the interaction effect of model and number of environment
g_model_window_effect <- cluster_pred_window_fit_effects %>% 
  filter(set == "realistic") %>%
  unnest(model_window_effects) %>% 
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(model = factor(model, levels = dist_method_abbr)) %>%
  ggplot(aes(x = window, y = fit, ymin = lower, ymax = upper, fill = model)) + 
  # geom_point() + 
  geom_line() +
  # geom_errorbar(width = 0.5) +
  geom_ribbon(color = "grey75", alpha = 0.2) +
  facet_wrap(~trait, scales = "free") +
  scale_color_manual(values = dist_colors, name = "Model") +
  scale_fill_manual(values = dist_colors, name = "Model") +
  scale_y_continuous(breaks = pretty) +
  ylab("Accuracy") +
  xlab("Sliding window interval") +
  theme_presentation2() +
  theme(legend.position = "bottom")

# Save
ggsave(filename = "window_pred_model_window_effect_realistic.jpg", plot = g_model_window_effect, path = fig_dir, 
       height = 5, width = 8, dpi = 1000)


## Overall, what is the percent change in prediction accuracy over the sliding windows, for each model?
cluster_pred_window_fit_effects %>% 
  filter(set == "realistic") %>%
  unnest(model_window_effects) %>% 
  filter(window %in% c(min(window), max(window))) %>% 
  group_by(trait, model) %>%
  summarize(acc_per_change = (fit[1] - fit[2]) / fit[1]) %>%
  spread(model, acc_per_change)

# trait         GrCD   LcPD  MYAE  MYICE  MYMCE
# 1 GrainYield   0.127  0.494 0.184 -0.292  0.267
# 2 HeadingDate -0.313  0.222 0.205  0.217 -0.269
# 3 PlantHeight  0.290 -0.168 0.742  0.220  0.270








#### Environmental cluster predictions ####

# Load the cluster results
load(file.path(result_dir, "cluster_predictions.RData"))


## Unnest and tidy
cluster_predictions_out <- cluster_predictions %>% 
  # unnest(out) %>%
  filter(model %in% names(dist_method_replace)) %>%
  filter(trait %in% traits) %>%
  mutate(dist_method = str_replace_all(model, dist_method_replace),
         model = factor(abbreviate(dist_method), levels = dist_method_abbr),
         zscore = ztrans(accuracy)) %>% 
  mutate_at(vars(min_env, cluster, environment), as.factor)



## Model formula
## Fit the fixed effects of:
## 1. model
## 2. minimum number of training environments in a cluster
## 3. the interaction
## Fit the random effects of:
## 1. environment
## 2. environment nested in cluster
## 3. interaction of model and cluster
## 

# Fit a model per trait
cluster_predictions_model <- distinct(cluster_predictions_out, set, trait) %>%
  mutate(out = list(NULL))
  
for (i in seq(nrow(cluster_predictions_model))) {
  tr <- cluster_predictions_model$trait[i]
  st <- cluster_predictions_model$set[i]
  
  fit <- lmer(formula = zscore ~ model + min_env + model:min_env + (1|environment) + (1|environment:cluster) + (1|model:cluster), 
              data = cluster_predictions_out, subset = trait == tr & set == st)
  
  out <- data_frame(fitted = list(fit),
                    model_effects = list(Effect("model", fit)),
                    min_env_effects = list(Effect("min_env", fit)),
                    model_min_env_effects = list(Effect(c("model", "min_env"), fit)),
                    fixef_test = list(broom::tidy(anova(fit))),
                    ranef_test = list(broom::tidy(ranova(fit)))) %>%
    mutate_at(vars(contains("effects")), ~map(., as.data.frame))
  
  cluster_predictions_model$out[[i]] <- out
  
}

cluster_predictions_model1 <- unnest(cluster_predictions_model)
  

## ANOVAs
unnest(cluster_predictions_model1, fixef_test)


## RANOVAs
unnest(cluster_predictions_model1, ranef_test) %>% filter(!str_detect(term, "none"))






## Plot the effect of model
g_cluster_model_effect <- cluster_predictions_model1 %>% 
  filter(set == "complete") %>%
  left_join(., select(cluster_pred_fit_effects_all, set, trait, mean_effects)) %>%
  unnest(model_effects) %>%
  mutate(model = factor(model, levels = dist_method_abbr)) %>%
  mutate_at(vars(fit, lower, upper), zexp) %>%
  ## Plot
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) +
  geom_hline(aes(yintercept = mean_effects, lty = "Accuracy using\nall data"), color = "grey75") +
  geom_point() +
  geom_errorbar(width = 0.5) + 
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 3)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty) + 
  facet_grid(~ trait) +
  ylab("Prediction accuracy") + 
  theme_presentation2(base_size = 16) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.box.margin = margin())

## Save
ggsave(filename = "cluster_model_effect.jpg", plot = g_cluster_model_effect, path = fig_dir, width = 8, height = 5, dpi = 1000)


## Plot the effect of model conditional on min_env
g_cluster_model_effect1 <- cluster_predictions_model1 %>% 
  filter(set == "complete") %>%
  left_join(., select(cluster_pred_fit_effects_all, set, trait, mean_effects)) %>%
  unnest(model_min_env_effects) %>%
  mutate(model = factor(model, levels = dist_method_abbr),
         min_env = factor(min_env, levels = levels(cluster_predictions_out$min_env))) %>%
  mutate_at(vars(fit, lower, upper), zexp) %>%
  ## Plot
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) +
  geom_hline(aes(yintercept = mean_effects, lty = "Accuracy using\nall data")) +
  geom_point() +
  geom_errorbar(width = 0.5) + 
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 3)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty) + 
  facet_grid(trait ~ min_env, scales = "free_y") +
  ylab("Prediction accuracy") + 
  theme_presentation2(base_size = 16) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.box.margin = margin())


## Save
ggsave(filename = "cluster_model_min_env_effect.jpg", plot = g_cluster_model_effect1, path = fig_dir, width = 8, height = 8, dpi = 1000)








### Realistic scenario
## Plot the effect of model
g_cluster_model_effect <- cluster_predictions_model1 %>% 
  filter(set == "realistic") %>%
  left_join(., select(cluster_pred_fit_effects_all, set, trait, mean_effects)) %>%
  unnest(model_effects) %>%
  mutate(model = factor(model, levels = dist_method_abbr)) %>%
  mutate_at(vars(fit, lower, upper), zexp) %>%
  ## Plot
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) +
  geom_hline(aes(yintercept = mean_effects, lty = "Accuracy using\nall data"), color = "grey75") +
  geom_point() +
  geom_errorbar(width = 0.5) + 
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 2)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty) + 
  facet_grid(~ trait) +
  ylab("Prediction accuracy") + 
  theme_presentation2(base_size = 16) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.box.margin = margin())

## Save
ggsave(filename = "cluster_model_effect_realistic.jpg", plot = g_cluster_model_effect, path = fig_dir, width = 8, height = 5, dpi = 1000)


## Plot the effect of model conditional on min_env
g_cluster_model_effect1 <- cluster_predictions_model1 %>% 
  filter(set == "realistic") %>%
  left_join(., select(cluster_pred_fit_effects_all, set, trait, mean_effects)) %>%
  unnest(model_min_env_effects) %>%
  mutate(model = factor(model, levels = dist_method_abbr),
         min_env = factor(min_env, levels = levels(cluster_predictions_out$min_env))) %>%
  mutate_at(vars(fit, lower, upper), zexp) %>%
  ## Plot
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) +
  geom_hline(aes(yintercept = mean_effects, lty = "Accuracy using\nall data")) +
  geom_point() +
  geom_errorbar(width = 0.5) + 
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 2)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty) + 
  facet_grid(trait ~ min_env, scales = "free_y") +
  ylab("Prediction accuracy") + 
  theme_presentation2(base_size = 16) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.box.margin = margin())


## Save
ggsave(filename = "cluster_model_min_env_effect_realistic.jpg", plot = g_cluster_model_effect1, path = fig_dir, width = 8, height = 8, dpi = 1000)














################
################
## Appendix ####
################
################



# # Rename the distance methods
# cumulative_pred_results <- cluster_pred_out %>% 
#   unnest(out) %>%
#   rename(dist_method = model) %>%
#   mutate(iter = parse_number(dist_method),
#          iter = ifelse(str_detect(dist_method, "sample"), iter, 1),
#          dist_method = ifelse(str_detect(dist_method, "sample"), "sample", dist_method),
#          dist_method = str_replace_all(dist_method, dist_method_replace),
#          dist_method_abbr = factor(abbreviate(dist_method), levels = dist_method_abbr)) %>%
#   left_join(., select(loeo_accuracy, trait, environment, max_accuracy = base), by = c("environment", "trait")) %>%
#   group_by(environment, trait, dist_method, iter) %>%
#   mutate(max_accuracy_adding = accuracy[n_e == max(n_e)]) %>%
#   ungroup()
# 
# 
# ## Summarize the accuracies and create CIs for the random method
# cumulative_pred_results_summ <- cumulative_pred_results %>%
#   group_by(environment, trait, dist_method, dist_method_abbr, n_e, max_accuracy, max_accuracy_adding) %>% 
#   summarize_at(vars(accuracy), funs(mean, sd, n())) %>% 
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat) %>%
#   ungroup()
# 
# 
# # Split by trait
# cumulative_pred_results_split <- cumulative_pred_results_summ %>% 
#   split(.$trait) %>% 
#   purrr::map(~mutate(., environment = factor(environment, levels = unique(.$environment[order(.$max_accuracy_adding, decreasing = TRUE)]))))
# 
# # Plot
# g_plotlist <- cumulative_pred_results_split %>%
#   purrr::map(~ggplot(data = ., aes(x = n_e, y = mean, color = dist_method_abbr)) + 
#                geom_hline(aes(yintercept = max_accuracy_adding, lty = "Accuracy\nUsing\nAll Data")) +
#                # geom_point() + 
#                geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
#                geom_line() + 
#                scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#                scale_x_continuous(breaks = pretty) +
#                scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
#                facet_wrap(~ environment, ncol = 5) + 
#                xlab("Number of training environments") +
#                ylab("Prediction accuracy") +
#                theme_minimal() )
# 
# # Save
# for (i in seq_along(g_plotlist)) {
#   ggsave(filename = paste0("cumulative_environment_prediction_", names(g_plotlist)[i], ".jpg"), plot = g_plotlist[[i]],
#          height = 10, width = 8, dpi = 1000, path = fig_dir)
# }
# 
# 
# ## Highlight STP16
# stp_subset <- cumulative_pred_results_summ %>%
#   filter(environment == "STP16", trait != "HeadingDateAGDD") %>%
#   filter(dist_method_abbr %in% c("GrCD", "PhnD", "TYEC", "Rndm"))
# 
# g_cumulative_pred_example <- stp_subset %>%
#   # filter(dist_method_abbr %in% c("GrCD", "Rndm")) %>%
#   ggplot(data = ., aes(x = n_e, y = mean, color = dist_method_abbr)) +
#   geom_hline(aes(yintercept = max_accuracy_adding, lty = "Accuracy using\nall environments")) +
#   geom_point() +
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
#   geom_line() +
#   geom_blank(data = stp_subset, aes(x = n_e, y = mean), inherit.aes = FALSE) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure", guide = guide_legend(order = 1)) +
#   scale_x_continuous(breaks = pretty) +
#   scale_linetype_manual(values = 2, name = NULL, guide = guide_legend(order = 2)) +
#   facet_wrap(~ trait, scales = "free_y") +
#   xlab("Number of training environments") +
#   ylab("Prediction accuracy") +
#   theme_presentation2() + theme(legend.position = "bottom", legend.direction = "horizontal")
# 
# # Save
# ggsave(filename = "cumulative_environment_prediction_STP16_1.jpg", plot = g_cumulative_pred_example,
#        height = 5, width = 10, dpi = 1000, path = fig_dir)
# 
# 
# 
# 
# ## For each trait and environment, calculate the difference between the accuracy
# ## using the nth training environment and the max accuracy
# ## Then summarize for each trait
# cumulative_pred_diff <- cumulative_pred_results %>% 
#   mutate(diff_accuracy = accuracy - max_accuracy_adding) %>%
#   group_by(trait, dist_method_abbr, n_e, iter) %>% 
#   summarize(diff_accuracy = mean(diff_accuracy), max_accuracy = mean(max_accuracy)) %>% ## Take the mean over all environments for a method/iteration
#   summarize_at(vars(diff_accuracy, max_accuracy), funs(mean, sd, n())) %>%
#   ungroup() %>%
#   select(trait:n_e, max_accuracy = max_accuracy_mean, mean = diff_accuracy_mean, sd = diff_accuracy_sd, n = diff_accuracy_n) %>%
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)
# 
# # Plot - Random
# g_diff_accuracy <- cumulative_pred_diff %>% 
#   ggplot(aes(x = n_e, y = mean, color = dist_method_abbr)) + 
#   geom_hline(aes(yintercept = 0, lty = "Accuracy using\nall environments")) +
#   geom_line() + 
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_linetype_manual(values = 2, name = NULL) +
#   scale_x_continuous(breaks = pretty) +
#   scale_y_continuous(breaks = pretty) +
#   facet_grid(trait ~ ., switch = "y", scales = "free_y") + 
#   xlab("Number of training environments") +
#   ylab("Prediction accuracy (relative to using all data)") +
#   theme_presentation() +
#   theme(strip.placement = "outside")
# 
# # Save
# ggsave(filename = "cumulative_environment_prediction_relative.jpg", plot = g_diff_accuracy,
#        height = 8, width = 7, dpi = 1000, path = fig_dir)
# 
# 
# 
# # Plot - examples
# g_diff_accuracy <- cumulative_pred_diff %>% 
#   mutate_at(vars(mean, lower, upper), funs(. + max_accuracy)) %>%
#   filter(trait != "HeadingDateAGDD", dist_method_abbr %in% c("GrCD", "PhnD", "TYEC", "Rndm")) %>% 
#   ggplot(aes(x = n_e, y = mean, color = dist_method_abbr)) + 
#   geom_hline(aes(yintercept = max_accuracy, lty = "Accuracy using\nall environments")) +
#   geom_point() +
#   geom_line() + 
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure", guide = guide_legend(order = 1)) +
#   scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_linetype_manual(values = 2, name = NULL, guide = guide_legend(order = 3)) +
#   scale_x_continuous(breaks = pretty) +
#   scale_y_continuous(breaks = pretty) +
#   facet_wrap(~ trait, scales = "free_y", nrow = 1) + 
#   xlab("Number of training environments") +
#   ylab("Prediction accuracy") +
#   theme_presentation2() +
#   theme(legend.position = "bottom")
# 
# # Save
# ggsave(filename = "cumulative_environment_prediction_relative_example.jpg", plot = g_diff_accuracy,
#        height = 5, width = 10, dpi = 1000, path = fig_dir)
# 
# 
# ## Choose an arbitrary number of environments and show the average accuracy advantage
# cumulative_pred_diff_example <- cumulative_pred_diff %>% 
#   filter(n_e %in% c(3, 5, 10)) %>% 
#   select(trait:n_e, mean) %>% 
#   mutate(mean = round(mean, 3)) %>% 
#   spread(dist_method_abbr, mean) %>% 
#   arrange(n_e)
# 
# # Save as an image
# t_cumulative_pred_diff <- grid.arrange(tableGrob(cumulative_pred_diff_example, rows = NULL, theme = ttheme_minimal()))
# 
# ggsave(filename = "cumulative_environment_prediction_example.jpg", plot = t_cumulative_pred_diff, path = fig_dir,
#        height = 5, width = 7, dpi = 1000)
# 
# 
# ## For each distance method, find the average number of environments in which the accuracy is maximized
# cumulative_pred_nE <- cumulative_pred_results %>% 
#   group_by(trait, dist_method_abbr, iter, environment) %>%
#   top_n(x = ., n = 1, wt = accuracy) %>% summarize(n_e = mean(n_e))
# 
# cumulative_pred_nE_summ <- cumulative_pred_nE %>%
#   summarize(n_e = mean(n_e)) %>%
#   summarize_at(vars(n_e), funs(mean, sd, n())) %>%
#   ungroup() %>%
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)
# 
# # Plot
# g_n_train <- cumulative_pred_nE_summ %>% 
#   ggplot(aes(x = trait, y = mean, color = dist_method_abbr)) +
#   geom_point(position = position_dodge(0.6), size = 3) +
#   geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6)) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_y_continuous(breaks = pretty) +
#   ylab("Peak accuracy training environments") +
#   theme_presentation() +
#   theme(axis.title.x = element_blank())
# 
# ggsave(filename = "cumulative_max_accuracy_nE.jpg", plot = g_n_train, path = fig_dir,
#        height = 5, width = 8, dpi = 1000)
# 
# 
# ## Using the phenotypic distance measure, calculate the average distance at which
# ## the maximum prediction accuracy is achieved
# 
# # First calculate the mean distance for each training set
# average_pheno_distance <- pred_env_dist_rank$tp %>% 
#   filter(model == "pheno_dist") %>% 
#   mutate(env_rank = map(env_rank, ~unlist(.) %>% data_frame(pred_environment = names(.), distance = .) %>% mutate(n_e = seq(nrow(.))))) %>% 
#   unnest() %>%
#   group_by(trait, environment) %>% 
#   mutate(distance = cummean(distance)) %>%
#   ungroup() %>%
#   select(-model)
# 
# 
# cumulative_pred_nE_summ_avg <- cumulative_pred_nE %>% 
#   left_join(., average_pheno_distance) %>%
#   summarize(distance = mean(distance)) %>%
#   summarize_at(vars(distance), funs(mean, sd, n())) %>%
#   ungroup() %>%
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)
# 
# 
# # Plot
# g_dist<- cumulative_pred_nE_summ_avg %>% 
#   ggplot(aes(x = trait, y = mean, color = dist_method_abbr)) +
#   geom_point(position = position_dodge(0.6), size = 3) +
#   geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6)) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_y_continuous(breaks = function(x) pretty(x, n = 3)) +
#   ylab("Peak accuracy training environments") +
#   facet_wrap(~trait, scales = "free", nrow = 1) +
#   theme_presentation() +
#   theme(axis.title.x = element_blank(), strip.text = element_blank())
# 
# ggsave(filename = "cumulative_max_accuracy_pheno_dist.jpg", plot = g_dist, path = fig_dir,
#        height = 5, width = 10, dpi = 1000)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## Window predictions
# 
# 
# # Rename the distance methods
# window_pred_results <- cluster_pred_out_window %>% 
#   mutate(out = map(out, ~mutate(., window = seq(nrow(.))))) %>%
#   unnest(out) %>%
#   rename(dist_method = model) %>%
#   mutate(iter = parse_number(dist_method),
#          iter = ifelse(str_detect(dist_method, "sample"), iter, 1),
#          dist_method = ifelse(str_detect(dist_method, "sample"), "sample", dist_method),
#          dist_method = str_replace_all(dist_method, dist_method_replace),
#          dist_method_abbr = factor(abbreviate(dist_method), levels = dist_method_abbr)) %>%
#   left_join(., select(loeo_accuracy, trait, environment, max_accuracy = base), by = c("environment", "trait")) 
# 
# 
# ## Summarize the accuracies and create CIs for the random method
# window_pred_results_summ <- window_pred_results %>%
#   group_by(environment, trait, dist_method, dist_method_abbr, window, max_accuracy) %>% 
#   summarize_at(vars(accuracy), funs(mean, sd, n())) %>% 
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat) %>%
#   ungroup()
# 
# 
# # Split by trait
# window_pred_results_split <- window_pred_results_summ %>% 
#   split(.$trait) %>% 
#   map(~mutate(., environment = factor(environment, levels = unique(.$environment[order(.$max_accuracy, decreasing = TRUE)]))))
# 
# # Plot
# g_plotlist <- window_pred_results_split %>%
#   map(~ggplot(data = ., aes(x = window, y = mean, color = dist_method_abbr)) + 
#         geom_hline(aes(yintercept = max_accuracy, lty = "Accuracy\nUsing\nAll Data")) +
#         # geom_point() + 
#         geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
#         geom_line() + 
#         scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#         scale_x_continuous(breaks = pretty) +
#         scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
#         facet_wrap(~ environment, ncol = 5) + 
#         xlab("Window number") +
#         ylab("Prediction accuracy") +
#         theme_minimal() )
# 
# # Save
# for (i in seq_along(g_plotlist)) {
#   ggsave(filename = paste0("window_environment_prediction_", names(g_plotlist)[i], ".jpg"), plot = g_plotlist[[i]],
#          height = 10, width = 8, dpi = 1000, path = fig_dir)
# }
# 
# 
# 
# 
# 
# 
# # For each trait and environment, calculate the difference between the accuracy
# ## using the nth training environment and the max accuracy
# ## Then summarize for each trait
# window_pred_diff <- window_pred_results %>% 
#   mutate(diff_accuracy = accuracy - max_accuracy) %>%
#   group_by(trait, dist_method_abbr, window, iter) %>% 
#   summarize(diff_accuracy = mean(diff_accuracy)) %>% ## Take the mean over all environments for a method/iteration
#   summarize_at(vars(diff_accuracy), funs(mean, sd, n())) %>%
#   ungroup() %>%
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)
# 
# # Plot - Random
# g_diff_accuracy <- window_pred_diff %>% 
#   ggplot(aes(x = window, y = mean, color = dist_method_abbr)) + 
#   geom_hline(aes(yintercept = 0, lty = "Accuracy\nUsing\nAll Data")) +
#   geom_line() + 
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill = dist_colors["Rndm"], alpha = 0.25) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
#   scale_x_continuous(breaks = pretty) +
#   scale_y_continuous(breaks = pretty) +
#   facet_grid(trait ~ ., switch = "y", scales = "free_y") + 
#   xlab("Window number") +
#   ylab("Prediction accuracy (relative to using all data)") +
#   theme_presentation() +
#   theme(strip.placement = "outside")
# 
# # Save
# ggsave(filename = "window_environment_prediction_relative.jpg", plot = g_diff_accuracy,
#        height = 8, width = 7, dpi = 1000, path = fig_dir)
# 
# 
# ## Fit smoothing lines
# g_diff_accuracy_smooth <- window_pred_diff %>% 
#   ggplot(aes(x = window, y = mean, color = dist_method_abbr)) + 
#   geom_hline(aes(yintercept = 0, lty = "Accuracy\nUsing\nAll Data")) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE) + 
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_linetype_manual(values = c("Accuracy\nUsing\nAll Data" = 2), name = NULL) +
#   scale_x_continuous(breaks = pretty) +
#   scale_y_continuous(breaks = pretty) +
#   facet_grid(trait ~ ., switch = "y", scales = "free_y") + 
#   xlab("Window number") +
#   ylab("Prediction accuracy (relative to using all data)") +
#   theme_presentation() +
#   theme(strip.placement = "outside")
# 
# # Save
# ggsave(filename = "window_environment_prediction_relative_smooth.jpg", plot = g_diff_accuracy_smooth,
#        height = 8, width = 7, dpi = 1000, path = fig_dir)
# 
# 
# 
# ## For each distance method, find the average window number in which the accuracy is maximized
# window_pred_wn <- window_pred_results %>% 
#   group_by(trait, dist_method_abbr, iter, environment) %>%
#   top_n(x = ., n = 1, wt = accuracy) %>% summarize(window = mean(window))
# 
# window_pred_wn_summ <- window_pred_wn %>%
#   summarize(window = mean(window)) %>%
#   summarize_at(vars(window), funs(mean, sd, n())) %>%
#   ungroup() %>%
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat)
# 
# # Plot
# g_window_train <- window_pred_wn_summ %>% 
#   ggplot(aes(x = trait, y = mean, color = dist_method_abbr)) +
#   geom_point(position = position_dodge(0.6), size = 3) +
#   geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6)) +
#   scale_color_manual(values = dist_colors, name = "Distance\nmeasure") +
#   scale_y_continuous(breaks = pretty) +
#   ylab("Peak accuracy training environments") +
#   theme_presentation() +
#   theme(axis.title.x = element_blank())
# 
# ggsave(filename = "window_max_accuracy_wn.jpg", plot = g_window_train, path = fig_dir,
#        height = 5, width = 8, dpi = 1000)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #### Environment covariance matrix prediction ####
# 
# # Load the results
# load(file.path(result_dir, "env_cov_mat_predictions.RData"))
# 
# 
# 
# ## MC predictions
# mc_pred_tidy <- environment_mc_predictions %>%
#   unnest() %>% unnest() %>% 
#   select(-trait1, -trait2) %>%
#   mutate_at(vars(trait, model, environment), as.factor)
# 
# ## Summarize the correlations across all environments for each iteration
# mc_pred_summ <- mc_pred_tidy %>% 
#   group_by(trait, pTrainEnv, model, environment, iter) %>% 
#   summarize(accuracy = cor(value, pgv))
# 
# # Now take the mean over iterations
# mc_pred_summ1 <- mc_pred_summ %>% 
#   summarize(accuracy = mean(accuracy))
# 
# # Plot for each model and trait
# g_model_acc <- mc_pred_summ1 %>%
#   ggplot(aes(x = trait, y = accuracy, fill = model)) +
#   geom_boxplot(position = "dodge", alpha = 0.5) +
#   xlab("Trait") + 
#   ylab("Prediction accuracy") +
#   scale_fill_discrete(name = "Model") +
#   facet_grid(~ pTrainEnv) +
#   theme_presentation2()
# 
# ggsave(filename = "environmental_cov_model_accuracy.jpg", plot = g_model_acc, path = fig_dir, width = 10, height = 4, dpi = 1000)
# 








