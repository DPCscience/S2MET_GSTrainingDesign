## S2MET Predictions
## 
## Script for analyzing cross-validation results
## 
## Author: Jeff Neyhart
## Last modified: January 28, 2019
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
library(broom)
library(car)

# Load the environmental distance df
load(file.path(result_dir, "distance_method_results.RData"))

# Load results
cv_results <- list.files(result_dir, full.names = TRUE, pattern = "cross_validation_")
lapply(cv_results, load, envir = .GlobalEnv)


## Significant level
alpha <- 0.05


## Create a factor for the distance methods
dist_method_replace <- c("pheno_dist" = "Phenotypic Distance", "pheno_loc_dist" = "Location Phenotypic Distance",  "great_circle_dist" = "Great Circle Distance", 
                         "OYEC_All" = "One Year All ECs", "OYEC_Mean" = "One Year Mean Cor EC", "OYEC_IPCA" = "One Year IPCA Cor EC", 
                         "MYEC_All" = "Multi Year All ECs", "MYEC_Mean" = "Multi Year Mean Cor EC", "MYEC_IPCA" = "Multi Year IPCA Cor EC",
                         "sample" = "Random")
dist_method_abbr <- abbreviate(dist_method_replace)

## Alternative abbreviations
dist_method_abbr <- setNames(c("PD", "LocPD", "GCD", "1Yr-All-EC", "1Yr-Mean-EC", "1Yr-IPCA-EC", "All-EC", "Mean-EC", "IPCA-EC", "Random"),
                             names(dist_method_replace))

colors <- umn_palette(3)
colors_use <- c(colors[2], "#FFDE7A", colors[c(1, 3, 8, 4, 9, 5, 10)], "grey75")
dist_colors <- setNames(colors_use, dist_method_abbr)






###
### Cross-validation results
### 
### 

# Regular cross-validation


## Calculate prediction accuracy and prepare for modelling
cv_prediction_accuracy <- cv_predictions %>% 
  unnest(prediction) %>% 
  group_by(trait, cv, model, environment, .id) %>% 
  summarize(accuracy = cor(value, pred_value)) %>%
  ungroup()
  

cv_prediction_accuracy_tomodel <- cv_prediction_accuracy %>%
  mutate(zscore = ztrans(accuracy),
         design = str_remove_all(cv, "[0-9]"),
         cv = str_extract(cv, "cv[0-9]{1,2}")) %>%
  mutate_at(vars(environment, model, cv), as.factor)

## Fit a model by trait and CV
## The model will estimate the fixed effect of model and random everything else
cv0_prediction_accuracy_fit <- cv_prediction_accuracy_tomodel %>%
  filter(cv == "cv0") %>%
  group_by(cv, trait, model, design) %>%
  summarize(mean_acc = mean(zscore))


## For each trait, fit a model then test for variance heterogeneity across CV schemes
cv_prediction_accuracy_tomodel %>%
  filter(cv != "cv0") %>% 
  group_by(trait) %>% 
  do(levene_test = leveneTest(y = .$zscore, group = .$cv)) %>% 
  pull()

## Results suggest heterogeneity of variance, so fit model for each CV scheme



# Other cv
other_cv_prediction_accuracy_fit <- cv_prediction_accuracy_tomodel %>%
  filter(cv != "cv0") %>% filter(cv != "cv00") %>%
  distinct(design, trait, cv) %>%
  # distinct(design, trait) %>%
  mutate(fit = list(NULL), effect = list(NULL), anova = list(NULL), ranova = list(NULL))

## Iterate over rows
for (i in seq(nrow(other_cv_prediction_accuracy_fit))) {
  row <- other_cv_prediction_accuracy_fit[i,]
  df <- cv_prediction_accuracy_tomodel %>% filter(trait == row$trait, cv == row$cv)

  
  # Fit the model
  fit <- lmer(zscore ~ 1 + model + (1|environment) + (1|model:environment), data = df)
  # fit <- lmer(zscore ~ 1 + model + cv + model:cv + (1|environment) + (1|model:environment) + (1|cv:environment) + 
  #               (1|cv:model:environment), data = df)
  
  ran <- ranova(fit)
  aov <- anova(fit)
  
  other_cv_prediction_accuracy_fit$fit[[i]] <- fit
  other_cv_prediction_accuracy_fit$anova[[i]] <- aov
  other_cv_prediction_accuracy_fit$ranova[[i]] <- ran
  other_cv_prediction_accuracy_fit$effect[[i]] <- as.data.frame(Effect("model", fit))
  # other_cv_prediction_accuracy_fit$effect[[i]] <- as.data.frame(Effect(c("model", "cv"), fit))
  
}



## Extract effects
cv_prediction_accuracy_effect <- bind_rows(
  unnest(other_cv_prediction_accuracy_fit, effect),
  rename(cv0_prediction_accuracy_fit, fit = mean_acc)) %>%
  mutate_at(vars(fit, lower, upper), funs(zexp(.)))

## A vector to rename models
model_replace <- c("M2" = "M1 (Main effect)", "M3" = "M2 (GxE)", "M4" = "M1 (Main effect)", "M5" = "M2 (GxE)")


## Plot
g_cv <- cv_prediction_accuracy_effect %>%
  mutate(model = str_replace_all(model, model_replace),
         cv = str_to_upper(cv)) %>%
  rename(Model = model) %>%
  ggplot(aes(x = cv, y = fit, color = Model, shape = Model)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.75), width = 0.5, color = "black") +
  geom_point(position = position_dodge(0.75), size = 2) + 
  scale_y_continuous(breaks = pretty, limits = c(0, 1), name = "Prediction accuracy") +
  xlab("Cross-validation scheme") +
  facet_grid(~ trait) +
  # theme_presentation2(base_size = 10) +
  theme_acs() +
  theme(legend.position = c(0.15, 0.80))

ggsave(filename = "cross_validation_accuracy.jpg", plot = g_cv, path = fig_dir, width = 5, height = 3, dpi = 1000)






## CV0 and CV00 - predicting future years

# Calculate accuracy
cv_zero_future_acc <- cv_zero_future_prediction %>% 
  unnest(prediction) %>% 
  group_by(trait, cv, model, .id, environment) %>% 
  summarize(accuracy = cor(value, pred_value)) %>%
  ungroup() %>%
  mutate(zscore = ztrans(accuracy)) %>%
  mutate_at(vars(model, environment), as.factor)

## Summarize for CV0
cv0_future_summ <- cv_zero_future_acc %>% 
  filter(cv == "cv0") %>% 
  group_by(cv, trait, model) %>% 
  summarize(fit = mean(accuracy)) %>%
  ungroup()

# trait       model accuracy
# 1 GrainYield  M4       0.486
# 2 GrainYield  M5_PD    0.469
# 3 HeadingDate M4       0.822
# 4 HeadingDate M5_PD    0.815
# 5 PlantHeight M4       0.497
# 6 PlantHeight M5_PD    0.489


## Fit a model for CV00
cv00_future_tomodel <- cv_zero_future_acc %>%
  filter(cv == "cv00") %>%
  group_by(cv, trait) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq(nrow(cv00_future_tomodel))) {

  df <- unnest(cv00_future_tomodel[i,], data)
  
  # Fit a model
  fit <- lmer(zscore ~ model + (1|environment) + (1|model:environment), data = df)
  
  # Significance testing
  anova_df <- tidy(anova(fit))
  ranova_df <- tidy(ranova(fit))
  effects_df <- as.data.frame(Effect("model", fit))
  
  cv00_future_tomodel$out[[i]] <- data_frame(anova = list(anova_df), ranova = list(ranova_df), effects = list(effects_df))
  
}

cv00_future_summ <- cv00_future_tomodel %>% 
  unnest(out) %>% 
  select(-data) %>%
  unnest(effects) %>%
  mutate_at(vars(fit, lower, upper), zexp)


## Combine and plot
bind_rows(cv0_future_summ, cv00_future_summ) %>%
  mutate(cv = toupper(cv)) %>%
  ggplot(aes(x = cv, y = fit, ymin = lower, ymax = upper, shape = model, color = model)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.5, color = "black") +
  geom_point(position = position_dodge(0.9), size = 2) +
  facet_grid(~ trait)







####
#### Parent-offspring cross-validation
#### 




## Temporary edits
pocv_predictions_tidy <- bind_rows(
  pocv_predictions %>% filter(cv %in% c("pocv1", "pocv2")) %>% select(-model, -prediction) %>% unnest(),
  pocv_predictions %>% filter(cv %in% c("pocv0", "pocv00")) %>%   select(-results) ) %>%
  unnest()

## Calculate prediction accuracy and prepare for modelling
pocv_prediction_accuracy <- pocv_predictions_tidy %>% 
  # unnest(prediction) %>% 
  group_by(trait, cv, model, environment, .id) %>% 
  summarize(accuracy = cor(value, pred_value)) %>%
  ungroup()


pocv_prediction_accuracy_tomodel <- pocv_prediction_accuracy %>%
  mutate(zscore = ztrans(accuracy),
         design = str_remove_all(cv, "[0-9]"),
         cv = str_extract(cv, "cv[0-9]{1,2}")) %>%
  mutate_at(vars(environment, model, cv), as.factor)

## Fit a model by trait and CV
## The model will estimate the fixed effect of model and random everything else
pocv0_prediction_accuracy_fit <- pocv_prediction_accuracy_tomodel %>%
  filter(cv == "cv0") %>%
  group_by(design, cv, trait, model) %>%
  summarize(mean_acc = mean(zscore))


## For each trait, fit a model then test for variance heterogeneity across CV schemes
pocv_prediction_accuracy_tomodel %>%
  filter(cv != "cv0") %>% 
  group_by(design, trait) %>% 
  do(levene_test = leveneTest(y = .$zscore, group = .$cv)) %>% 
  pull()

## Results suggest heterogeneity of variance, so fit models separately for all CV schema



# Other cv
other_pocv_prediction_accuracy_fit <- pocv_prediction_accuracy_tomodel %>%
  filter(cv != "cv0") %>% filter(cv != "cv00") %>%
  distinct(design, trait, cv) %>%
  # distinct(trait) %>%
  mutate(fit = list(NULL), effect = list(NULL), anova = list(NULL), ranova = list(NULL))

## Iterate over rows
for (i in seq(nrow(other_pocv_prediction_accuracy_fit))) {
  row <- other_pocv_prediction_accuracy_fit[i,]
  df <- pocv_prediction_accuracy_tomodel %>% filter(trait == row$trait, cv == row$cv)
  # df <- pocv_prediction_accuracy_tomodel %>% filter(trait == row$trait)
  
  
  # Fit the model
  fit <- lmer(zscore ~ 1 + model + (1|environment) + (1|model:environment), data = df)
  
  ran <- ranova(fit)
  aov <- anova(fit)
  
  other_pocv_prediction_accuracy_fit$fit[[i]] <- fit
  other_pocv_prediction_accuracy_fit$anova[[i]] <- aov
  other_pocv_prediction_accuracy_fit$ranova[[i]] <- ran
  other_pocv_prediction_accuracy_fit$effect[[i]] <- as.data.frame(Effect("model", fit))
  # other_pocv_prediction_accuracy_fit$effect[[i]] <- as.data.frame(Effect(c("model", "cv"), fit))
  
}



## Extract effects
pocv_prediction_accuracy_effect <- bind_rows(
  unnest(other_pocv_prediction_accuracy_fit, effect),
  rename(pocv0_prediction_accuracy_fit, fit = mean_acc)) %>%
  mutate_at(vars(fit, lower, upper), funs(zexp(.)))

## Plot
g_pocv <- pocv_prediction_accuracy_effect %>%
  mutate(model = str_replace_all(model, model_replace),
         cv = str_to_upper(cv)) %>%
  rename(Model = model) %>%
  ggplot(aes(x = cv, y = fit, color = Model, shape = Model)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.75), width = 0.5, color = "black") +
  geom_point(position = position_dodge(0.75), size = 2) + 
  scale_y_continuous(breaks = pretty, limits = c(0, 1), name = "Prediction accuracy") +
  xlab("Cross-validation scheme") +
  facet_grid(~ trait) +
  # theme_presentation2(base_size = 10) +
  theme_acs() +
  theme(legend.position = c(0.15, 0.80))

ggsave(filename = "parent_offspring_cross_validation_accuracy.jpg", plot = g_pocv, path = fig_dir, width = 5, height = 3, dpi = 1000)






#####
##### Parent-offspring validation
##### 

pov1_summ <- pov1_prediction %>% 
  unnest() %>%
  group_by(trait, model, environment) %>% 
  summarize(accuracy = cor(value, pred_value)) %>%
  summarize(fit = mean(accuracy)) %>%
  ungroup() %>%
  mutate(cv = "pov1")

## Summarize for POV0
pov0_future_summ <- pov0_future_prediction %>%
  unnest() %>%
  group_by(trait, model, environment) %>% 
  summarize(accuracy = cor(value, pred_value)) %>%
  summarize(fit = mean(accuracy)) %>%
  ungroup() %>%
  mutate(cv = "pov0_future")

pov0_loeo_summ <- pov0_loeo_prediction %>%
  unnest() %>%
  group_by(trait, model, environment) %>%
  summarize(accuracy = cor(value, pred_value)) %>%
  summarize(fit = mean(accuracy)) %>%
  ungroup() %>%
  mutate(cv = "pov0_loeo")

## Fit models to POV00
pov00_tomodel <- bind_rows(
  unnest(pov00_future_prediction) %>% mutate(cv = "pov00_future"),
  unnest(pov00_loeo_prediction) %>% mutate(cv = "pov00_loeo")
) %>%
  group_by(cv, trait, model, environment) %>%
  summarize(accuracy = cor(value, pred_value)) %>%
  ungroup() %>%
  mutate(model = as.factor(model),
         zscore = ztrans(accuracy)) %>%
  group_by(cv, trait) %>%
  nest() %>%
  mutate(out = list(NULL))


for (i in seq(nrow(pov00_tomodel))) {
  
  df <- unnest(pov00_tomodel[i,], data)
  
  # Fit a model
  fit <- lmer(zscore ~ model + (1|environment), data = df)
  
  # Significance testing
  anova_df <- tidy(anova(fit))
  ranova_df <- tidy(ranova(fit))
  effects_df <- as.data.frame(Effect("model", fit))
  
  pov00_tomodel$out[[i]] <- data_frame(anova = list(anova_df), ranova = list(ranova_df), effects = list(effects_df))
  
}

pov00_summ <- pov00_tomodel %>% 
  unnest(out) %>% 
  select(-data) %>%
  unnest(effects) %>%
  mutate_at(vars(fit, lower, upper), zexp)


## Combine and plot
bind_rows(pov1_summ, pov0_future_summ, pov0_loeo_summ, pov00_summ) %>%
  mutate(cv = toupper(cv)) %>%
  ggplot(aes(x = cv, y = fit, ymin = lower, ymax = upper, shape = model, color = model)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.5, color = "black") +
  geom_point(position = position_dodge(0.9), size = 2) +
  facet_grid(~ trait)

  
  
  



# Calculate accuracy
cv_zero_future_acc <- cv_zero_future_prediction %>% 
  unnest(prediction) %>% 
  group_by(trait, cv, model, .id, environment) %>% 
  summarize(accuracy = cor(value, pred_value)) %>%
  ungroup() %>%
  mutate(zscore = ztrans(accuracy)) %>%
  mutate_at(vars(model, environment), as.factor)

## Summarize for CV0
cv0_future_summ <- cv_zero_future_acc %>% 
  filter(cv == "cv0") %>% 
  group_by(cv, trait, model) %>% 
  summarize(fit = mean(accuracy)) %>%
  ungroup()












## Combine
all_pred_accuracy_effect <- bind_rows(cv_prediction_accuracy_effect, pocv_prediction_accuracy_effect)


## Plot
g_all_cv <- all_pred_accuracy_effect %>%
  mutate(model = str_replace_all(model, c("M2" = "M2 (Main effect)", "M3" = "M3 (GxE)")),
         cv = str_to_upper(cv), design = str_to_upper(design),
         group = paste0(design, "_", model)) %>%
  ggplot(aes(x = cv, y = fit, color = design, shape = model)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = group), position = position_dodge(0.75), width = 0.5, color = "black") +
  geom_point(position = position_dodge(0.75), size = 2) + 
  scale_color_discrete(name = "Design") +
  scale_shape_discrete(name = "Model") +
  scale_y_continuous(breaks = pretty, limits = c(0, 1), name = "Prediction accuracy") +
  xlab("Cross-validation scheme") +
  facet_grid(trait ~ .) +
  # theme_presentation2(base_size = 10) +
  theme_acs() +
  theme(legend.position = "bottom", legend.direction = "horizontal")


ggsave(filename = "all_cross_validation_accuracy.jpg", plot = g_all_cv, path = fig_dir, width = 4, height = 4, dpi = 1000)


#


































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








