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
library(broom)

# Load the environmental distance df
load(file.path(result_dir, "distance_method_results.RData"))

# Load results
load(file.path(result_dir, "distance_rank_predictions.RData"))
# Load the cluster results
load(file.path(result_dir, "cluster_predictions.RData"))


## Significant level
alpha <- 0.05

## Cutoff of prediction accuracy to remove an environment
env_cutoff <- 0.1
env_cutoff <- -Inf

environment_rank_pred_out <- bind_rows(environment_rank_pred_out)
environment_window_pred_out <- bind_rows(environment_window_pred_out)

## Unique models
unique(environment_rank_pred_out$model)


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










#### Leave-one-environment-out predictions ####


environment_loeo_predictions_analysis <- environment_loeo_predictions_geno_mean %>% 
  group_by(trait, set) %>% 
  summarize(accuracy = mean(accuracy)) %>%
  ungroup()


# trait       set       accuracy
# 1 GrainYield  complete     0.281
# 2 GrainYield  realistic    0.321
# 3 HeadingDate complete     0.421
# 4 HeadingDate realistic    0.418
# 5 PlantHeight complete     0.396
# 6 PlantHeight realistic    0.383


## Min, max, mean
environment_loeo_predictions_geno_mean %>% 
  group_by(trait, set) %>% 
  summarize_at(vars(accuracy), funs(min, max, mean)) 

# trait       set           min   max  mean
# 1 GrainYield  complete  -0.187  0.504 0.281
# 2 GrainYield  realistic  0.0456 0.494 0.321
# 3 HeadingDate complete   0.104  0.646 0.421
# 4 HeadingDate realistic  0.235  0.614 0.418
# 5 PlantHeight complete   0.0686 0.621 0.396
# 6 PlantHeight realistic  0.0362 0.589 0.383




#### Environmental distance predictions ####



## Predictions when adding one environment at a time


# Convert the accuracies to z scores
# Adjust the names of the distance models
environment_rank_pred_out1 <- environment_rank_pred_out %>%
  unnest(out) %>% 
  left_join(., distinct(environment_loeo_predictions_geno_mean, set, trait, nTrainEnv)) %>% # Filter out results for the max nEnv
  # filter(n_e < nTrainEnv) %>%
  mutate(zscore = ztrans(accuracy),
         nEnv = as.factor(n_e),
         environment = as.factor(validation_environment),
         model = ifelse(str_detect(model, "sample"), "sample", model),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr)) %>%
  filter(!mat_set %in% c("Jarquin", "MalosettiStand")) %>% # Filter out the relationship matrices
  # unite(model, model, mat_set, sep = "_") %>% 
  mutate(model = as.factor(str_remove_all(model, "_NA"))) %>%
  group_by(trait, set, environment, model, nEnv, n_e) %>% 
  summarize(zscore = mean(zscore)) %>%
  ungroup()



## Fit a model per trait

# Fit a model per trait
environment_rank_pred_fit <- environment_rank_pred_out1 %>%
  group_by(set, trait) %>%
  nest() %>%
  mutate(out = list(NULL))

## Iterate over sets and traits - the summaries of LMER don't work using do() or map(), so we have to use a loop
for (i in seq(nrow(environment_rank_pred_fit))) {
  df <- environment_rank_pred_fit$data[[i]] %>%
    droplevels()
  
  fit <- lmer(formula = zscore ~ 1 + model + nEnv + model:nEnv + (1|environment) + (1|environment:model) + (1|environment:nEnv), data = df)
  
  fit_summary <- data_frame(
    fitted = list(fit),
    nEnv_effects = list(Effect(focal.predictors = "nEnv", fit)),
    model_nEnv_effects = list(Effect(focal.predictors = c("model", "nEnv"), fit))) %>% 
    mutate_at(vars(contains("effects")), ~map(., ~as.data.frame(.)))
  
  # Add the summary to the DF
environment_rank_pred_fit$out[[i]] <- fit_summary
  
}

environment_rank_pred_fit <- unnest(environment_rank_pred_fit, out)

# Quick anova scan
environment_rank_pred_fit$fitted %>% map(anova)



## Collect the effects
environment_rank_pred_fit_effects <- environment_rank_pred_fit %>% 
  filter(! str_detect(trait, "AGDD")) %>%
  mutate_at(vars(contains("effects")), ~map(., as.data.frame))






## Plot the effect of model at an early number of training environments
## Highlight 1, 5, and 10 environments
nEnv_select <- c("1", "5", "10")
  
g_model_effect_nenv_list <-  environment_rank_pred_fit_effects %>% 
  split(.$set) %>%
  map(~{
    df <- unnest(., model_nEnv_effects)
    # Results for the maximum number of envs
    max_env <- mutate(df, nEnv = parse_number(nEnv)) %>% 
      group_by(trait) %>% 
      filter(nEnv == max(nEnv)) %>%
      slice(1) %>% 
      select(set, trait, fit) %>%
      mutate(mean_effects = zexp(fit))
    
    df %>%
      filter(nEnv %in% nEnv_select) %>%
      mutate_at(vars(fit, se, lower, upper), zexp) %>%
      mutate(model = factor(model, levels = dist_method_abbr),
             nEnv = factor(nEnv, levels = nEnv_select)) %>%
      ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) + 
      geom_hline(data = max_env, aes(yintercept = mean_effects, lty = "Accuracy using\nall data" )) +
      geom_point() + 
      geom_errorbar(width = 0.5) +
      facet_grid(trait ~ nEnv, scales = "free") +
      scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 3)) +
      scale_linetype_manual(values = 2, name = NULL) + 
      scale_y_continuous(breaks = pretty) +
      ylab("Prediction accuracy") +
      theme_presentation2() +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "bottom")
  })

## Plot this differently
g_model_effect_nenv_list <-  environment_rank_pred_fit_effects %>% 
  split(.$set) %>%
  map(~{
    df <- unnest(., model_nEnv_effects)
    # Results for the maximum number of envs
    max_env <- mutate(df, nEnv = parse_number(nEnv)) %>% 
      group_by(trait) %>% 
      filter(nEnv == max(nEnv)) %>%
      slice(1) %>% 
      select(set, trait, fit) %>%
      mutate(mean_effects = zexp(fit))
    
    df %>%
      filter(nEnv %in% nEnv_select) %>%
      mutate_at(vars(fit, se, lower, upper), zexp) %>%
      mutate(model = factor(model, levels = dist_method_abbr),
             nEnv = factor(nEnv, levels = nEnv_select)) %>%
      ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) + 
      geom_hline(data = max_env, aes(yintercept = mean_effects, lty = "Accuracy using\nall data" )) +
      geom_point() + 
      geom_errorbar(width = 0.5) +
      facet_grid(trait ~ nEnv, scales = "free") +
      scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 3)) +
      scale_linetype_manual(values = 2, name = NULL) + 
      scale_y_continuous(breaks = pretty) +
      ylab("Prediction accuracy") +
      theme_presentation2() +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "bottom")
  })




# Save
for (i in seq_along(g_model_effect_nenv_list)) {
  filename <- paste0("cumulative_pred_model_effect_envselect_", names(g_model_effect_nenv_list)[i], ".jpg")
  ggsave(filename = filename, plot = g_model_effect_nenv_list[[i]], path = fig_dir, height = 8, width = 8, dpi = 1000)
}



## Plot the interaction effect of model and number of environment
environment_rank_pred_fit_effects1 <- environment_rank_pred_fit_effects %>%
  unnest(model_nEnv_effects) %>%
  mutate_at(vars(fit, se, lower, upper), zexp) %>%
  mutate(nEnv = parse_number(nEnv), 
         model = factor(model, levels = dist_method_abbr), 
         set = str_to_title(set))


g_model_nenv_effect <- environment_rank_pred_fit_effects1 %>% 
  ggplot(aes(x = nEnv, y = fit, color = model, fill = model)) +
  geom_hline(data = group_by(environment_rank_pred_fit_effects1, set, trait) %>% filter(nEnv == max(nEnv)) %>% slice(1), 
             aes(yintercept = fit, lty = "Accuracy using\nall data")) +
  geom_point(size = 0.5) +
  geom_line(lwd = 0.5) +
  # geom_ribbon(aes(ymin = lower, ymax = upper),color = "grey75", alpha = 0.2) +
  facet_grid(trait ~ set, scales = "free", space = "free_x", switch = "y") +
  scale_color_manual(values = dist_colors, name = "Distance\nmeasure", guide = guide_legend(nrow = 3)) +
  scale_fill_manual(values = dist_colors, name = "Distance\nmeasure") +
  scale_y_continuous(breaks = pretty) +
  scale_x_continuous(breaks = function(x) seq(1, x[2], by = 5)) + 
  scale_linetype_manual(values = 2, name = NULL) +
  ylab("Prediction accuracy") +
  xlab("Number of training environments") +
  theme_presentation2(base_size = 10) +
  theme(legend.position = "bottom", legend.spacing.x = unit(x = 0.2, units = "line"), legend.key.height = unit(1, "line"))

# Save
ggsave(filename = "cumulative_pred_model_nEnv_effect.jpg", plot = g_model_nenv_effect, path = fig_dir, height = 5, width = 4, dpi = 1000)














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

## Re-organize, extract the random results
cluster_predictions_out <- cluster_predictions %>% 
  mutate(accuracy = map_dbl(out, "base"),
         random = map(out, "random") %>% map(~data_frame(rep = paste0("rep", seq_along(.)), accuracy = .)),
         model = ifelse(model == "pheno_location_dist", "pheno_loc_dist", model)) %>%
  select(-out)

  
cluster_predictions_base <- cluster_predictions_out %>%
  select(-random) %>%
  mutate(model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr),
         zscore = ztrans(accuracy)) %>% 
  mutate_at(vars(cluster, val_environment), as.factor)



## Model formula
## Models fitted individually for each set and trait
## Fixed effect of model and number of training environments
## zscore ~ model + nTrainEnv + (1|cluster:model) + (1|val_environment:cluster)
## 


# Fit a model per trait
cluster_predictions_base_analysis <- cluster_predictions_base %>%
  group_by(set, trait) %>%
  nest() %>%
  mutate(out = list(NULL))
  
## Iterate over sets and traits - the summaries of LMER don't work using do() or map(), so we have to use a loop
for (i in seq(nrow(cluster_predictions_base_analysis))) {
  df <- cluster_predictions_base_analysis$data[[i]] %>%
    droplevels()
  
  fit <- lmer(formula = zscore ~ 1 + model + nTrainEnv + (1|cluster:model) + (1|val_environment:cluster), data = df,
              contrasts = list(model = "contr.sum"))
  
  fit_summary <- data_frame(
        fitted = list(fit),
        mean_accuracy = fixef(fit)[1],
        model_effects = list(Effect("model", fit)),
        nTrainEnv_effects = list(Effect("nTrainEnv", fit)),
        fixef_test = list(tidy(anova(fit))),
        ranef_test = list(tidy(ranova(fit)))
      ) %>% mutate_at(vars(contains("effects")), ~map(., ~as.data.frame(.)))

  # Add the summary to the DF
  cluster_predictions_base_analysis$out[[i]] <- fit_summary
  
}
  
cluster_predictions_base_analysis <- cluster_predictions_base_analysis %>%
  select(-data) %>% 
  unnest(out)

## ANOVAs
unnest(cluster_predictions_base_analysis, fixef_test)

## Model was significant only for heading date and realistic


## RANOVAs
unnest(cluster_predictions_base_analysis, ranef_test) %>% filter(!str_detect(term, "none"))

## Environment %in% cluster and cluster %in% model was usually significant
## Not for heading date realistic





## Plot the effect of model
g_cluster_model_effect <- cluster_predictions_base_analysis %>%
  unnest(model_effects) %>%
  mutate(model = factor(model, levels = dist_method_abbr), set = str_to_title(set)) %>%
  mutate_at(vars(fit, lower, upper), zexp) %>%
  ## Plot
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) +
  geom_hline(data = group_by(environment_rank_pred_fit_effects1, set, trait) %>% filter(nEnv == max(nEnv)) %>% slice(1),
             aes(yintercept = fit, lty = "Accuracy using\nall data"), color = "black") +
  geom_point() +
  geom_errorbar(width = 0.5) + 
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 3)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty) + 
  facet_grid(trait ~ set, scales = "free_x", space = "free_x", switch = "y") +
  ylab("Prediction accuracy") + 
  theme_presentation2(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.box.margin = margin(), legend.spacing.x = unit(0.2, "lines"))
  
## Save
ggsave(filename = "cluster_model_effect.jpg", plot = g_cluster_model_effect, path = fig_dir, width = 4.5, height = 5, dpi = 1000)






## Subset for poster
dist_method_poster <- dist_method_abbr[c(1:3,7:9)]
# Rename the prediction scenario
set_replace <- c("complete" = "Leave-One-Out", "realistic" = "Time-Forward" )

environment_rank_pred_fit_effects1_use <- environment_rank_pred_fit_effects1 %>%
  mutate(set = str_to_lower(set),
         set = str_replace_all(set, set_replace)) %>%
  group_by(set, trait) %>% 
  filter(nEnv == max(nEnv)) %>% 
  slice(1)


g_cluster_model_effect_poster <- cluster_predictions_base_analysis %>%
  unnest(model_effects) %>%
  mutate(model = factor(model, levels = dist_method_abbr),
         set = str_replace_all(set, set_replace)) %>%
  filter(model %in% dist_method_poster) %>%
  mutate_at(vars(fit, lower, upper), zexp) %>%
  ## Plot
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) +
  geom_hline(data = environment_rank_pred_fit_effects1_use, aes(yintercept = fit, lty = "Accuracy using\nall data"), color = "black") +
  geom_point() +
  geom_errorbar(width = 0.5) + 
  scale_color_manual(values = dist_colors, name = "Distance\nmeasure", guide = guide_legend(ncol = 1)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty) + 
  facet_grid(trait ~ set, scales = "free_x", space = "free_x", switch = "y", labeller = labeller(trait = str_add_space)) +
  ylab("Prediction accuracy") + 
  labs(title = "Genomewide prediction accuracy") +
  theme_presentation2(base_size = 16) +
  theme(legend.position = "left", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.box.margin = margin(), legend.spacing.x = unit(0.2, "lines"))

## Save
ggsave(filename = "cluster_model_effect_poster.jpg", plot = g_cluster_model_effect_poster, path = fig_dir, width = 8.5, height = 6, dpi = 1000)






## Output a table
cluster_model_table <- cluster_predictions_base_analysis %>%
  unnest(model_effects) %>% 
  mutate_at(vars(fit, lower, upper), funs(formatC(., digits = 2))) %>% 
  mutate(annotation = paste0(fit, " (", lower, ", ", upper, ")")) %>% 
  bind_rows(., mutate(environment_loeo_predictions_analysis, annotation = formatC(accuracy, digits = 2), model = "AllData")) %>%
  select(trait, set, model, annotation) %>% 
  spread(model, annotation) %>%
  arrange(set)

write_csv(x = cluster_model_table, path = file.path(fig_dir, "cluster_model_predictions_table.csv"))




## Effect of n of training environments
cluster_predictions_base_analysis %>% 
  unnest(nTrainEnv_effects) %>%
  ggplot(aes(x = nTrainEnv, y = fit, ymin = lower, ymax = upper)) +
  geom_line() +
  geom_ribbon(color = "grey", alpha = 0.2) +
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 3)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty) + 
  facet_grid(trait ~ set, scales = "free_x", space = "free_x") +
  ylab("Prediction accuracy") + 
  theme_presentation2(base_size = 16) +
  theme(legend.position = "bottom", legend.box.margin = margin(), legend.spacing.x = unit(0.5, "lines"))






# ## Model the random iterations
# ## This analysis will model the difference between the base accuracy and the random accuracy
# cluster_predictions_random <- cluster_predictions_out %>%
#   unnest(random) %>%
#   rename(base_accuracy = accuracy, accuracy = accuracy1) %>%
#   mutate(dist_method = str_replace_all(model, dist_method_replace),
#          model = factor(abbreviate(dist_method), levels = dist_method_abbr)) %>%
#   mutate_at(vars(contains("accuracy")), funs(zscore = ztrans)) %>%
#   mutate(zscore_difference = base_accuracy_zscore - accuracy_zscore,
#          zscore_difference = ifelse(is.na(zscore_difference), 0, zscore_difference)) %>%
#   mutate_at(vars(cluster, val_environment), as.factor)
# 
# ## There will be some NAs from these conditions:
# ## 
# # 1 complete  MYICE PlantHeight
# # 2 complete  OYICE PlantHeight
# # 3 realistic LcPD  HeadingDate
# # 4 realistic MYICE GrainYield 
# # 
# # This is due to there being only one cluster - fill with zero
# 
# 
# ## Model formula
# ## Models fitted individually for each set and trait
# ## zscore ~ model + nTrainEnv + (1|cluster:model) + (1|val_environment:cluster)
# ## 
# 
# 
# # Fit a model per trait
# cluster_predictions_random_analysis <- cluster_predictions_random %>%
#   group_by(set, trait) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# ## Iterate over sets and traits - the summaries of LMER don't work using do() or map(), so we have to use a loop
# for (i in seq(nrow(cluster_predictions_random_analysis))) {
#   df <- cluster_predictions_random_analysis$data[[i]] %>%
#     droplevels()
#   
#   fit <- lmer(formula = zscore_difference ~ 1 + model + nTrainEnv + (1|cluster:model) + (1|val_environment:cluster) + (1|model:val_environment:cluster), 
#               data = df, contrasts = list(model = "contr.sum"))
#   
#   fit_summary <- data_frame(
#     fitted = list(fit),
#     mean_accuracy = fixef(fit)[1],
#     model_effects = list(Effect("model", fit)),
#     nTrainEnv_effects = list(Effect("nTrainEnv", fit)),
#     fixef_test = list(tidy(anova(fit))),
#     ranef_test = list(tidy(ranova(fit)))
#   ) %>% mutate_at(vars(contains("effects")), ~map(., ~as.data.frame(.)))
#   
#   # Add the summary to the DF
#   cluster_predictions_random_analysis$out[[i]] <- fit_summary
#   
# }
# 
# cluster_predictions_random_analysis <- cluster_predictions_random_analysis %>%
#   select(-data) %>% 
#   unnest(out)
# 
# 
# ## ANOVAs
# unnest(cluster_predictions_random_analysis, fixef_test)
# 
# ## Model was sometimes significant
# 
# 
# ## RANOVAs
# unnest(cluster_predictions_random_analysis, ranef_test) %>% filter(!str_detect(term, "none"))
# 
# ## Environment %in% cluster and cluster %in% model was more often not significant,
# ## which should be expected under random conditions



## Alternatively, analyze the accuracy of only the 2017 environments
# Fit a model per trait
cluster_predictions_analysis2017 <- cluster_predictions_out1 %>%
  filter(str_detect(val_environment, "17"), set != "realistic") %>%
  mutate(set = "complete2017") %>%
  droplevels() %>%
  group_by(set, trait) %>%
  nest() %>%
  mutate(out = list(NULL))

## Iterate over sets and traits - the summaries of LMER don't work using do() or map(), so we have to use a loop
for (i in seq(nrow(cluster_predictions_analysis2017))) {
  df <- cluster_predictions_analysis2017$data[[i]] %>%
    droplevels()
  
  fit <- lmer(formula = zscore ~ 1 + model + (1|cluster:model) + (1|val_environment:cluster), data = df,
              contrasts = list(model = "contr.sum"))
  
  fit_summary <- data_frame(
    fitted = list(fit),
    mean_accuracy = fixef(fit)[1],
    model_effects = list(Effect("model", fit)),
    fixef_test = list(tidy(anova(fit))),
    ranef_test = list(tidy(ranova(fit)))
  ) %>% mutate_at(vars(contains("effects")), ~map(., ~as.data.frame(.) %>% mutate(model = as.character(model))))
  
  # Add the summary to the DF
  cluster_predictions_analysis2017$out[[i]] <- fit_summary
  
}

cluster_predictions_analysis2017 <- cluster_predictions_analysis2017 %>%
  select(-data) %>% unnest(out)

## ANOVAs
unnest(cluster_predictions_analysis2017, fixef_test)

## Model was never significant

## RANOVAs
unnest(cluster_predictions_analysis2017, ranef_test) %>% filter(!str_detect(term, "none"))

## Environment %in% cluster and cluster %in% model was usually significant



## Plot the effect of model
g_cluster_model_effect1 <- bind_rows(cluster_predictions_analysis, cluster_predictions_analysis2017) %>%
  left_join(., bind_rows(environment_loeo_predictions_analysis, filter(environment_loeo_predictions_analysis2017, set == "complete2017"))) %>% # Add the all-data results
  rename(mean_effects = accuracy) %>%
  unnest(model_effects) %>%
  mutate(model = factor(model, levels = dist_method_abbr), set = str_to_title(set)) %>%
  mutate_at(vars(fit, lower, upper), zexp) %>%
  ## Plot
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) +
  geom_hline(aes(yintercept = mean_effects, lty = "Accuracy using\nall data"), color = "grey75") +
  geom_point() +
  geom_errorbar(width = 0.5) + 
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 3)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty) + 
  facet_grid(trait ~ set, scales = "free_x", space = "free_x") +
  ylab("Prediction accuracy") + 
  theme_presentation2(base_size = 16) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.box.margin = margin(), legend.spacing.x = unit(0.5, "lines"))

## Save
ggsave(filename = "cluster_model_effect_2017.jpg", plot = g_cluster_model_effect1, path = fig_dir, width = 8, height = 8, dpi = 1000)







## Cross validation cluster predictions
cluster_cv_predictions_out <- cluster_cv_predictions %>%
  mutate(cv_accuracy = map_dbl(out, 1),
         vp_accuracy = map_dbl(out, "base_vp")) %>%
  gather(design, accuracy, cv_accuracy, vp_accuracy) %>%
  mutate(zscore = ztrans(accuracy),
         model = ifelse(model == "pheno_location_dist", "pheno_loc_dist", model),
         dist_method = str_replace_all(model, dist_method_replace),
         model = factor(abbreviate(dist_method), levels = dist_method_abbr),
         design = as.factor(design)) %>%
  select(-out)
  


## Model formula
## Models fitted individually for each set and trait
## zscore ~ model + design + model:design + (1|cluster:model) + (1|val_environment:cluster)
## 


# Fit a model per trait
cluster_cv_predictions_analysis <- cluster_cv_predictions_out %>%
  group_by(set, trait) %>%
  nest() %>%
  mutate(out = list(NULL))

## Iterate over sets and traits - the summaries of LMER don't work using do() or map(), so we have to use a loop
for (i in seq(nrow(cluster_cv_predictions_analysis))) {
  df <- cluster_cv_predictions_analysis$data[[i]] %>%
    droplevels()
  
  fit <- lmer(formula = zscore ~ model + design + model:design + nTrainEnv + (1|cluster:model) + (1|val_environment:cluster) + 
                (1|cv_rep:val_environment), 
              data = df, contrasts = list(model = "contr.sum", design = "contr.sum"))

  fit_summary <- data_frame(
    fitted = list(fit),
    mean_accuracy = fixef(fit)[1],
    model_design_effects = list(Effect(c("model", "design"), fit)),
    design_effects = list(Effect("design", fit)),
    nTrainEnv_effects = list(Effect("nTrainEnv", fit)),
    fixef_test = list(tidy(anova(fit))),
    ranef_test = list(tidy(ranova(fit)))
  ) %>% mutate_at(vars(contains("effects")), ~map(., ~as.data.frame(.)))
  
  # Add the summary to the DF
  cluster_cv_predictions_analysis$out[[i]] <- fit_summary
  
}

cluster_cv_predictions_analysis <- cluster_cv_predictions_analysis %>%
  select(-data) %>% 
  unnest(out)


## ANOVAs
unnest(cluster_cv_predictions_analysis, fixef_test)

## Model was never significant
## Design was always signficiant
## Model x design was sometimes significant


## RANOVAs
unnest(cluster_cv_predictions_analysis, ranef_test) %>% filter(!str_detect(term, "none"))

## Environment %in% cluster and cluster %in% model was more often not significant,
## which should be expected under random conditions

## Plot the effect of model and prediction design
g_cluster_cv_model_effect_complete <- cluster_cv_predictions_analysis %>%
  filter(set == "complete") %>%
  unnest(model_design_effects) %>%
  mutate(model = factor(model, levels = dist_method_abbr), set = str_to_title(set),
         design = str_replace_all(design, c("vp_accuracy" = "POV", "cv_accuracy" = "CV"))) %>%
  mutate_at(vars(fit, lower, upper), zexp) %>%
  ## Plot
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) +
  geom_point() +
  geom_errorbar(width = 0.5) + 
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 2)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty) + 
  facet_grid(trait ~ design, scales = "free", space = "free_x") +
  ylab("Prediction accuracy") + 
  theme_presentation2(base_size = 16) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.box.margin = margin(), legend.spacing.x = unit(0.5, "lines"))

g_cluster_cv_model_effect_realistic <- cluster_cv_predictions_analysis %>%
  filter(set == "realistic") %>%
  unnest(model_design_effects) %>%
  mutate(model = factor(model, levels = dist_method_abbr), set = str_to_title(set),
         design = str_replace_all(design, c("vp_accuracy" = "POV", "cv_accuracy" = "CV"))) %>%
  mutate_at(vars(fit, lower, upper), zexp) %>%
  ## Plot
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper, color = model)) +
  geom_point() +
  geom_errorbar(width = 0.5) + 
  scale_color_manual(values = dist_colors, name = "Model", guide = guide_legend(nrow = 1)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_y_continuous(breaks = pretty) + 
  facet_grid(trait ~ design, scales = "free", space = "free_x") +
  ylab("Prediction accuracy") + 
  theme_presentation2(base_size = 16) +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.box.margin = margin(), legend.spacing.x = unit(0.5, "lines"))



## Save
ggsave(filename = "cluster_cv_model_design_effect_complete.jpg", plot = g_cluster_cv_model_effect_complete, path = fig_dir, 
       width = 6, height = 8, dpi = 1000)

## Save
ggsave(filename = "cluster_cv_model_design_effect_realistic.jpg", plot = g_cluster_cv_model_effect_realistic, path = fig_dir, 
       width = 6, height = 8, dpi = 1000)



## Plot the effect of design

g_cluster_cv_design_effect <- cluster_cv_predictions_analysis %>%
  unnest(design_effects) %>%
  mutate(set = str_to_title(set),
         design = str_replace_all(design, c("vp_accuracy" = "POV", "cv_accuracy" = "CV"))) %>%
  mutate_at(vars(fit, lower, upper), zexp) %>%
  ## Plot
  ggplot(aes(x = design, y = fit, ymin = lower, ymax = upper, color = design)) +
  geom_point() +
  geom_errorbar(width = 0.5) + 
  scale_color_brewer(name = "Prediction design", palette = "Set1", guide = FALSE) +
  scale_y_continuous(breaks = pretty) + 
  facet_grid(trait ~ set, scales = "free_y", space = "free_x") +
  ylab("Prediction accuracy") + 
  xlab("Design") + 
  theme_presentation2(base_size = 16) +
  theme(legend.position = "bottom", legend.box.margin = margin(), legend.spacing.x = unit(0.5, "lines"))


## Save
ggsave(filename = "cluster_cv_design_effect.jpg", plot = g_cluster_cv_design_effect, path = fig_dir, 
       width = 6, height = 8, dpi = 1000)







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








