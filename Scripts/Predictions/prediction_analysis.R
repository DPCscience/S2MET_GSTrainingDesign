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
library(gridExtra)
library(broom)
library(car)
library(cowplot)
library(grid)

# Load the environmental distance df
load(file.path(result_dir, "distance_method_results.RData"))

# Load results
load(file.path(result_dir, "all_predictions_00.RData"))
# load(file.path(result_dir, "distance_rank_predictions_pov.RData"))
load(file.path(result_dir, "distance_rank_predictions.RData"))





## Shorten the heritability results
env_herit <- stage_one_data %>% 
  filter(!str_detect(trial, "S2C1")) %>% 
  distinct(trait, val_environment = environment, heritability)






#### Cross-validation / POV using all data #####


## Combine the all-data prediction results
all_data_preds <- ls(pattern = "[0-9]{1,2}_predictions$") %>% 
  subset(., map_lgl(., ~inherits(get(.), "data.frame")))


cv_pov_results <- map(all_data_preds, get) %>%
  map(., ~rename_at(., vars(which(names(.) %in% "environment")), ~str_replace(., pattern = "environment", "val_environment"))) %>%
  map_df(., ~mutate_at(., vars(which(names(.) %in% "rep")), as.character)) %>%
  mutate(scheme = ifelse(scheme == "pcv00", "pocv0", scheme)) %>%
  ## Add CV/POV number
  mutate(scheme_number = as.character(str_extract(scheme, "[0-9]{1,2}")))



## Analyze together?
all_cv_pov_predictions_out <- cv_pov_results %>%
  left_join(., env_herit) %>%
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability),
         set = ifelse(is.na(set), "complete", set)) %>%
  mutate_at(vars(scheme, val_environment), as.factor)

## Test for homogeneity of variance across all schemes
homo_var_tests_schemes <- all_cv_pov_predictions_out %>% 
  group_by(set, trait) %>% 
  do(test = tidy(leveneTest(y = .$accuracy, group = .$scheme))) %>% 
  ungroup() %>%
  mutate(p_value = map_dbl(test, ~subset(., term == "group", p.value, drop = T)))

## Test for homogeneity of variance across scheme numbers
homo_var_tests_scheme_num <- all_cv_pov_predictions_out %>% 
  group_by(set, trait, scheme_number) %>% 
  do(test = tidy(leveneTest(y = .$accuracy, group = .$scheme))) %>% 
  ungroup() %>%
  mutate(p_value = map_dbl(test, ~subset(., term == "group", p.value, drop = T)))

### Evidence suggests testing each individually


# all_cv_pov_predictions_analysis <- all_cv_pov_predictions_out %>%
#   group_by(set, trait) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# 
# 
# ## Loop for each row
# for (i in seq(nrow(all_cv_pov_predictions_analysis))) {
#   
#   df <- all_cv_pov_predictions_analysis$data[[i]]
#   
#   ## Model
#   fit <- lmer(accuracy ~ scheme + (1|val_environment) + (1|val_environment:scheme), data = df)
#   
#   effects_keep <- list(as_data_frame(Effect(focal.predictors = "scheme", fit)))
#   
#   fit_summary <- data_frame(
#     fitted = list(fit),
#     scheme_effects = effects_keep,
#     anova = list(tidy(anova(fit))),
#     ranova = list(tidy(ranova(fit)))
#   )
#   
#   all_cv_pov_predictions_analysis$out[[i]] <- fit_summary
#   
#   
# }



## Separately
separate_cv_pov_predictions_analysis <- all_cv_pov_predictions_out %>%
  # group_by(set, trait, scheme) %>%
  group_by(set, trait, scheme_number) %>%
  nest() %>%
  mutate(model_rep = map(data, "rep") %>% map_lgl(~!all(is.na(.)))) %>%
  mutate(out = list(NULL))

## Loop for each row
for (i in seq(nrow(separate_cv_pov_predictions_analysis))) {
   
  df <- separate_cv_pov_predictions_analysis$data[[i]]
  
  ## If model_rep is FALSE, just take the mean
  if (!separate_cv_pov_predictions_analysis$model_rep[i]) {
    
    # boot_mean <- replicate(1000, mean(sample(df$accuracy, replace = TRUE)))
    # effects_keep <- summarize(df, accuracy = mean(accuracy)) %>%
    #   cbind(., matrix(quantile(x = boot_mean, c(alpha / 2, 1 - (alpha / 2))), ncol = 2)) %>% 
    #   rename_at(vars(-accuracy), ~c("lower", "upper"))
    # 
    # fit_summary <- data_frame(
    #   fitted = list(NULL),
    #   scheme_effects = list(effects_keep),
    #   anova = list(NULL),
    #   ranova = list(NULL)
    # )
    
    
    fit <- lmer(accuracy ~ 1 + scheme + (1|val_environment), data = df)
    
    effects_keep <- as.data.frame(Effect(focal.predictors = "scheme", fit))
    
    fit_summary <- data_frame(
      fitted = list(fit),
      scheme_effects = list(effects_keep),
      anova = list(tidy(anova(fit))),
      ranova = list(tidy(ranova(fit)))
    )
    
  } else {
    
    # ## Model
    # fit <- lmer(accuracy ~ 1 + (1|val_environment), data = df)
    # 
    # effects_keep <- cbind(data_frame(accuracy = fixef(fit)[1]), confint(fit, parm = "(Intercept)", method = "Wald")) %>%
    #   rename_at(vars(-accuracy), ~c("lower", "upper")) %>%
    #   as_data_frame()
    
    
    ## Model
    fit <- lmer(accuracy ~ 1 + scheme + (1|val_environment) + (1|val_environment:scheme), data = df)
    
    effects_keep <- as.data.frame(Effect(focal.predictors = "scheme", fit))
    
    fit_summary <- data_frame(
      fitted = list(fit),
      scheme_effects = list(effects_keep),
      anova = list(tidy(anova(fit))),
      ranova = list(tidy(ranova(fit)))
    )
    
  }
  
  separate_cv_pov_predictions_analysis$out[[i]] <- fit_summary
  
  
}


## Plot all
separate_cv_pov_predictions_analysis1 <- separate_cv_pov_predictions_analysis %>%
  mutate(effects = map(out, ~.$scheme_effects[[1]])) %>%
  unnest(effects) %>%
  mutate(scheme = factor(toupper(scheme), levels = cv_replace),
         set = factor(str_replace_all(set, set_replace), levels = set_replace)) %>%
  rename(accuracy = fit)

g_cv_pov_all_data <- separate_cv_pov_predictions_analysis1 %>%
  ggplot(aes(x = scheme, y = accuracy, ymin = lower, ymax = upper)) +
  geom_point() +
  geom_errorbar(width = 0.5) +
  scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
  xlab("Validation scheme") +
  facet_grid(trait ~ set, scales = "free_x", space = "free_x", labeller = labeller(set = str_to_title, trait = str_add_space), switch = "y") +
  theme_presentation2(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "all_cv_pov_predictions.jpg", plot = g_cv_pov_all_data, path = fig_dir, width = 4, height = 5, dpi = 1000)







#### Random environment subset predictions #####

random_pred <- ls(pattern = "[0-9]{1,2}_predictions_random") %>% 
  subset(., map_lgl(., ~inherits(get(.), "data.frame")))


cv_pov_random_results <- map(random_pred, get) %>%
  map(., ~rename_at(., vars(which(names(.) %in% "environment")), ~str_replace(., pattern = "environment", "val_environment"))) %>%
  map_df(., ~mutate_at(., vars(which(names(.) %in% "rep")), as.character)) %>%
  ## Add CV/POV number
  mutate(scheme_number = as.character(str_extract(scheme, "[0-9]{1,2}")),
         nTrainEnv = as.factor(nTrainEnv))



## Analyze together?
all_cv_pov_random_predictions_out <- cv_pov_random_results %>%
  left_join(., env_herit) %>%
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability),
         set = ifelse(is.na(set), "complete", set)) %>%
  mutate_at(vars(scheme, val_environment), as.factor)


## Separately
separate_cv_pov_random_predictions_analysis <- all_cv_pov_random_predictions_out %>%
  # group_by(set, trait, scheme) %>%
  group_by(set, trait, scheme_number) %>%
  nest() %>%
  mutate(model_rep = map(data, "rep") %>% map_lgl(~!all(is.na(.)))) %>%
  mutate(out = list(NULL))

## Loop for each row
for (i in seq(nrow(separate_cv_pov_random_predictions_analysis))) {
  
  df <- separate_cv_pov_random_predictions_analysis$data[[i]]
    
  fit <- lmer(accuracy ~ 1 + nTrainEnv + scheme + nTrainEnv:scheme + (1|val_environment) + (1|val_environment:scheme) + (1|val_environment:nTrainEnv), data = df)
  
  effects_keep <- as.data.frame(Effect(focal.predictors = c("scheme", "nTrainEnv"), fit))
  
  fit_summary <- data_frame(
    fitted = list(fit),
    scheme_effects = list(effects_keep),
    anova = list(tidy(anova(fit))),
    ranova = list(tidy(ranova(fit)))
  )
    
  separate_cv_pov_random_predictions_analysis$out[[i]] <- fit_summary
  
  
}


## Plot all
separate_cv_pov_random_predictions_analysis1 <- separate_cv_pov_random_predictions_analysis %>%
  mutate(effects = map(out, ~.$scheme_effects[[1]])) %>%
  unnest(effects) %>%
  mutate(scheme = ifelse(scheme == "pov0", "pocv0", scheme),
         scheme = factor(toupper(scheme), levels = cv_replace)) %>%
  rename(accuracy = fit)

g_cv_pov_random_data <- separate_cv_pov_random_predictions_analysis1 %>%
  ggplot(aes(x = nTrainEnv, y = accuracy, ymin = lower, ymax = upper)) +
  geom_point() +
  geom_errorbar(width = 0.5) +
  scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
  xlab("Number of training environments") +
  facet_grid(trait ~ set + scheme, scales = "free_x", space = "free_x", labeller = labeller(set = str_to_title, trait = str_add_space)) +
  theme_presentation2(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "all_cv_pov_random_predictions.jpg", plot = g_cv_pov_random_data, path = fig_dir, width = 10, height = 5, dpi = 1000)


## Add the accuracy when using all environments
g_cv_pov_random_data1 <- separate_cv_pov_random_predictions_analysis1 %>%
  left_join(., select(separate_cv_pov_predictions_analysis1, set, trait, scheme, fit = accuracy)) %>%
  ggplot(aes(x = nTrainEnv, y = accuracy, ymin = lower, ymax = upper)) +
  geom_hline(aes(yintercept = fit, lty = "All data")) +
  geom_point() +
  geom_errorbar(width = 0.5) +
  scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
  scale_linetype_manual(values = 2, name = NULL) +
  xlab("Number of training environments") +
  facet_grid(trait ~ set + scheme, scales = "free_x", space = "free_x", labeller = labeller(set = str_to_title, trait = str_add_space)) +
  theme_presentation2(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.05, 0.05), legend.margin = margin())

ggsave(filename = "all_cv_pov_random_predictions1.jpg", plot = g_cv_pov_random_data1, path = fig_dir, width = 10, height = 5, dpi = 1000)











#### Environmental distance predictions ####

## Make sure all observations were recorded
environment_rank_df <- pred_env_dist_rank %>%
  rename(val_environment = validation_environment) %>%
  filter(!mat_set %in% c("Jarquin", "MalosettiStand")) %>%
  filter(model %in% names(dist_method_abbr_use)) %>%
  select(-mat_set) %>%
  mutate(data = list(NULL)) %>%
  mutate(nTrainEnv = map(rank, seq)) %>%
  unnest(nTrainEnv)

## POV
# Regular
environment_rank_df %>% 
  left_join(., bind_rows(pov00_environment_rank_predictions)) %>%
  filter(is.na(accuracy))
# Random
environment_rank_df %>% 
  filter(model == "great_circle_dist") %>% 
  select(-model) %>% 
  crossing(., model = paste0("sample", seq(25))) %>%
  left_join(., bind_rows(pov00_environment_rank_random_predictions)) %>%
  filter(is.na(accuracy))


## CV
# Regular
environment_rank_df %>%
  crossing(rep = seq(25)) %>%
  left_join(., bind_rows(cv00_environment_rank_predictions)) %>%
  filter(is.na(cv00))

# Random
environment_rank_df %>% 
  filter(model == "great_circle_dist") %>% 
  select(-model) %>% 
  crossing(., rep = seq(25)) %>%
  mutate(model = paste0("sample", rep)) %>%
  left_join(., bind_rows(cv00_environment_rank_random_predictions)) %>%
  filter(is.na(cv00))





## Analyze pov00 first
pov00_predictions_out <- bind_rows(pov00_environment_rank_predictions, pov00_environment_rank_random_predictions) %>%
  mutate(model = ifelse(str_detect(model, "sample"), "sample", model)) %>%
  left_join(., env_herit) %>%
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, dist_method_abbr)) %>%
  mutate_at(vars(scheme, val_environment, nTrainEnv), as.factor) %>%
  filter(model %in% dist_method_abbr_use)


 
# Fit a model per trait
pov00_predictions_analysis <- pov00_predictions_out %>%
  group_by(set, trait, scheme) %>%
  nest() %>%
  mutate(out = list(NULL))

# ## Iterate over sets and traits - the summaries of LMER don't work using do() or map(), so we have to use a loop
# for (i in seq(nrow(pov00_predictions_analysis))) {
# 
#   df <- pov00_predictions_analysis$data[[i]] %>% droplevels()
#   df$nTrainEnv <- ordered(df$nTrainEnv)
#   
#   ## Fit a mixed model, where nTrainEnv is given a correlation covariance structure
#   fm1 <- lme(fixed = accuracy ~ 1 + model, data = df, random = ~ 1|val_environment + 1|nTrainEnv,
#              correlation = corAR1(form = ~ 1 |nTrainEnv))
#   
#   ## Check autocorrelation of residuals
#   acf(residuals(fm1, type = "normalized"))
#   
#   
#   
#   
#   fit <- lmer(formula = accuracy ~ 1 + model + nTrainEnv + model:nTrainEnv + (1|val_environment) + (1|val_environment:model) +
#                 (1|val_environment:nTrainEnv), data = df)
# 
#   fit_summary <- data_frame(
#     fitted = list(fit),
#     model_nEnv_effects = list(Effect(focal.predictors = c("model", "nTrainEnv"), fit))) %>%
#     mutate_at(vars(contains("effects")), ~map(., ~as.data.frame(.)))
# 
#   # Add the summary to the DF
#   pov00_predictions_analysis$out[[i]] <- fit_summary
# 
# }
# 
# pov00_predictions_analysis <- unnest(pov00_predictions_analysis, out)









## Just take averages, since everything is balanced
pov00_predictions_analysis <- pov00_predictions_out %>% 
  group_by(trait, set, model, nTrainEnv, scheme) %>%
  summarize(fit = mean(accuracy)) %>% 
  ungroup()



## Plot
pov00_predictions_analysis_toplot <- pov00_predictions_analysis %>%
  # unnest(model_nEnv_effects) %>%
  mutate(set = factor(str_replace_all(set, set_replace), levels = set_replace),
         model = factor(model, levels = rev(dist_method_abbr_use)),
         nTrainEnv = parse_number(nTrainEnv),
         scheme = toupper(scheme),
         size = model == "Random")

## Subset the final training environment accuracy
pov00_rank_pred_all_data <- pov00_predictions_analysis_toplot %>%
  group_by(set, trait, scheme) %>% 
  filter(model == model[1], nTrainEnv == max(nTrainEnv)) %>%
  select(trait, set, scheme, accuracy = fit)


pov00_predictions_analysis_toplot1 <- left_join(pov00_predictions_analysis_toplot, pov00_rank_pred_all_data)


  
g_pov00_rank_pred <- pov00_predictions_analysis_toplot1 %>%
  ggplot(aes(x = nTrainEnv, y = fit, color = model, group = model)) +
  geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
  geom_line(aes(size = size)) +
  facet_grid(trait ~ set, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
  scale_color_manual(values = dist_colors_use, name = "Distance\nmeasure",
                     guide = guide_legend(nrow = 3, keyheight = unit(0.5, "line"), reverse = TRUE)) +
  scale_linetype_manual(values = 2, name = NULL, guide = FALSE) +
  scale_size_manual(values = c(0.5, 1), guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
  scale_x_continuous(breaks = pretty, name = "Number of training set environments") +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom", strip.placement = "outside")

g_pov00_rank_pred1 <- g_pov00_rank_pred +
  scale_linetype_manual(values = 2, name = NULL, guide = guide_legend()) +
  scale_color_manual(values = dist_colors_use, guide = FALSE) +
  theme(legend.position = c(0.38, 0.05))

g_pov00_rank_pred2 <- plot_grid(g_pov00_rank_pred1, get_legend(g_pov00_rank_pred), ncol = 1, rel_heights = c(1, 0.07))


ggsave(filename = "cumulative_env_pred_pov00.jpg", plot = g_pov00_rank_pred2, path = fig_dir, width = 3.5, height = 6, dpi = 1000)




## What is the best/worst accuracy achievable for each model?
pov00_predictions_analysis_toplot1 %>%
  mutate(advantage = fit - accuracy, perc_advantage = advantage / accuracy) %>%
  group_by(set, trait, model) %>%
  mutate_at(vars(advantage), funs(min, max)) %>%
  filter(advantage == min | advantage == max) %>%
  ungroup() %>%
  arrange(set, trait, model)









## Cross-validation prediction
cv00_predictions_out <- bind_rows(cv00_environment_rank_predictions, cv00_environment_rank_random_predictions) %>%
  gather(scheme, accuracy, cv00, pocv00) %>%
  mutate(model = ifelse(str_detect(model, "sample"), "sample", model)) %>%
  left_join(., env_herit) %>%
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, dist_method_abbr),
         scheme_number = str_extract(scheme, "[0-9]{1,2}")) %>%
  mutate_at(vars(scheme, val_environment, nTrainEnv), as.factor) %>%
  filter(model %in% dist_method_abbr_use)

# ## Plot
# cv00_predictions_out %>% 
#   group_by(trait, set, nTrainEnv, scheme, model) %>% 
#   summarize(fit = mean(accuracy)) %>% 
#   ungroup() %>%
#   mutate(nTrainEnv = parse_number(nTrainEnv)) %>%
#   ggplot(aes(x = nTrainEnv, y = fit, color = model)) +
#   geom_line() +
#   facet_grid(trait ~ set + scheme)

# 
# ## Fit a model per trait
# cv00_predictions_analysis <- cv00_predictions_out %>%
#   bind_rows(., mutate(pov00_predictions_out, scheme_number = "00")) %>%
#   filter(scheme != "pocv00") %>%
#   mutate(scheme = as.factor(scheme)) %>%
#   group_by(set, trait, scheme_number) %>%
#   nest() %>%
#   mutate(out = list(NULL))
# 
# ## Iterate over sets and traits - the summaries of LMER don't work using do() or map(), so we have to use a loop
# for (i in seq(nrow(cv00_predictions_analysis))) {
#   
#   df <- cv00_predictions_analysis$data[[i]] %>% # mutate(nTrainEnv = parse_number(nTrainEnv)) %>%
#     droplevels()
#   
#   
#   if (any(summarize_at(df, vars(model, nTrainEnv), n_distinct) <= 1)) next
# 
#   fit <- lmer(formula = accuracy ~ 1 + scheme + model + nTrainEnv + scheme:model + nTrainEnv:model + scheme:model:nTrainEnv + 
#                 (1|val_environment) + (1|val_environment:model) + (1|val_environment:scheme) + (1|val_environment:model:scheme), 
#               data = filter(df, model %in% c("AMMI", "PD")))
#   
#   fit <- lmer(formula = accuracy ~ 1 + scheme + model + nTrainEnv + scheme:model + nTrainEnv:model + scheme:model:nTrainEnv + 
#                 (1|val_environment) + (1|val_environment:model) + (1|val_environment:scheme) + (1|val_environment:model:scheme) +
#                 (1|val_environment:nTrainEnv), 
#               data = filter(df, model %in% c("AMMI", "PD")))
#   
#   
#   fit_summary <- data_frame(
#     fitted = list(fit),
#     model_nEnv_effects = list(Effect(focal.predictors = c("model", "nTrainEnv", "scheme"), fit))) %>% 
#     mutate_at(vars(contains("effects")), ~map(., ~as.data.frame(.)))
#   
#   # Add the summary to the DF
#   cv00_predictions_analysis$out[[i]] <- fit_summary
#   
# }
# 
# cv00_predictions_analysis <- cv00_predictions_analysis %>%
#   filter(!map_lgl(out, is.null)) %>%
#   unnest(out)


## Just take averages, since everything is balanced
cv00_predictions_analysis <- cv00_predictions_out %>% 
  group_by(trait, set, model, nTrainEnv, scheme) %>%
  summarize(fit = mean(accuracy)) %>% 
  ungroup()


## Edit some columns for plotting
cv00_predictions_analysis_toplot <- cv00_predictions_analysis %>%
  # unnest(model_nEnv_effects) %>%
  mutate(set = factor(str_replace_all(set, set_replace), levels = set_replace),
         model = factor(model, levels = rev(dist_method_abbr_use)),
         nTrainEnv = parse_number(nTrainEnv),
         scheme = toupper(scheme),
         size = model == "Random")


## Subset the final training environment accuracy
cv00_rank_pred_all_data <- cv00_predictions_analysis_toplot %>%
  group_by(set, trait, scheme) %>% 
  filter(model == model[1], nTrainEnv == max(nTrainEnv)) %>%
  ungroup() %>%
  select(set, trait, scheme, accuracy = fit)


cv00_predictions_analysis_toplot1 <- left_join(cv00_predictions_analysis_toplot, cv00_rank_pred_all_data)



g_cv00_rank_pred <- cv00_predictions_analysis_toplot1 %>%
  ggplot(aes(x = nTrainEnv, y = fit, color = model, group = model)) +
  geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
  geom_line(aes(size = size)) +
  facet_grid(trait ~ set + scheme, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
  scale_color_manual(values = dist_colors_use, name = "Distance\nmeasure",
                     guide = guide_legend(nrow = 3, keyheight = unit(0.5, "line"), reverse = TRUE)) +
  scale_linetype_manual(values = 2, name = NULL, guide = FALSE) +
  scale_size_manual(values = c(0.5, 1), guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
  scale_x_continuous(breaks = pretty, name = "Number of training set environments") +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom", strip.placement = "outside")

g_cv00_rank_pred1 <- g_cv00_rank_pred +
  scale_linetype_manual(values = 2, name = NULL, guide = guide_legend()) +
  scale_color_manual(values = dist_colors_use, guide = FALSE) +
  theme(legend.position = c(0.20, 0.05))

g_cv00_rank_pred2 <- plot_grid(g_cv00_rank_pred1, get_legend(g_cv00_rank_pred), ncol = 1, rel_heights = c(1, 0.07))


ggsave(filename = "cumulative_env_pred_cv00.jpg", plot = g_cv00_rank_pred2, path = fig_dir, width = 7, height = 6, dpi = 1000)





## Combine CV and POV results
cv00_pov00_environment_rank_predictions_analysis_toplot <- bind_rows(cv00_predictions_analysis_toplot1, pov00_predictions_analysis_toplot1) %>%
  filter(scheme != "POCV00") %>%
  # Create limits
  group_by(trait, scheme) %>%
  mutate_at(vars(fit), funs(min, max)) %>%
  ungroup()
  

g_rank_pred <- cv00_pov00_environment_rank_predictions_analysis_toplot %>%
  mutate(group = paste0(set, " - ", scheme)) %>%
  ggplot(aes(x = nTrainEnv, y = fit, color = model, group = model)) +
  geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
  geom_line(aes(size = size)) +
  facet_grid(trait ~ group, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
  scale_color_manual(values = dist_colors_use, name = "Distance\nmeasure",
                     guide = guide_legend(nrow = 3, keyheight = unit(0.5, "line"), reverse = TRUE)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_size_manual(values = c(0.5, 1), guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
  scale_x_continuous(breaks = pretty, name = expression("Number of training set environments ("*italic(N[TE])*")")) +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom", strip.placement = "outside")

## Save
ggsave(filename = "cumulative_env_pred.jpg", plot = g_rank_pred, path = fig_dir, width = 7, height = 6, dpi = 1000)



## Plot editor
gg_vp <- function(x, let) {
  ggplot(data = x, aes(x = nTrainEnv, y = fit, color = model, group = model)) +
    geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
    geom_line(aes(size = size)) +
    geom_text(aes(x = Inf, y = -Inf, label = let), hjust = 2, vjust = -1, color = "black") +
    scale_color_manual(values = dist_colors_use, name = "Distance\nmeasure", guide = FALSE) +
    scale_linetype_manual(values = 2, name = NULL, guide = FALSE) +
    scale_size_manual(values = c(0.5, 1), guide = FALSE) +
    scale_y_continuous(breaks = function(x) pretty(x, 6)[c(3,6)], name = NULL, limits = c(unique(x$min), unique(x$max))) +
    # scale_y_continuous(breaks = pretty, name = NULL, limits = c(unique(x$min), unique(x$max))) +
    scale_x_continuous(breaks = pretty, name = NULL, labels = NULL) +
    theme_presentation2(base_size = 8) +
    theme(panel.border = element_rect(color = "black"))
}


## Add viewport of heading data/cv00/leave-one-out
g_rank_pred_vp1 <- cv00_pov00_environment_rank_predictions_analysis_toplot %>%
  filter(trait == "HeadingDate", set == "Leave-one-out", scheme == "CV00") %>%
  gg_vp(let = "A")

## Add viewport of heading data/pov00/leave-one-out
g_rank_pred_vp2 <- cv00_pov00_environment_rank_predictions_analysis_toplot %>%
  filter(trait == "HeadingDate", set == "Leave-one-out", scheme == "POV00")  %>%
  gg_vp(let = "B")

## Add viewport of heading data/cv00/time-forward
g_rank_pred_vp3 <- cv00_pov00_environment_rank_predictions_analysis_toplot %>%
  filter(trait == "HeadingDate", set == "Time-forward", scheme == "CV00")  %>%
  gg_vp(let = "C")

## Add viewport of heading data/pov00/time-forward
g_rank_pred_vp4 <- cv00_pov00_environment_rank_predictions_analysis_toplot %>%
  filter(trait == "HeadingDate", set == "Time-forward", scheme == "POV00")  %>%
  gg_vp(let = "D")


## Coordinates for viewports
vp1 <- viewport(x = 0.22, y = 0.52, width = 0.28, height = 0.2)
vp2 <- viewport(x = 0.51, y = 0.59, width = 0.28, height = 0.19)
vp3 <- viewport(x = 0.745, y = 0.52, width = 0.16, height = 0.2)
vp4 <- viewport(x = 0.91, y = 0.59, width = 0.16, height = 0.19)


jpeg(filename = file.path(fig_dir, "cumulative_env_pred_paper.jpg"), width = 7, height = 6, units = "in", res = 1000)
print(g_rank_pred)
print(g_rank_pred_vp1, vp = vp1)
print(g_rank_pred_vp2, vp = vp2)
print(g_rank_pred_vp3, vp = vp3)
print(g_rank_pred_vp4, vp = vp4)
dev.off()












## CV0 and POV0

cv0_pov0_predictions_out <- bind_rows(cv0_pocv0_environment_rank_predictions, cv0_pocv0_environment_rank_random_predictions) %>%
  mutate(model = ifelse(!is.na(rep), "sample", model),
         model = ifelse(str_detect(model, "sample"), "sample", model)) %>%
  left_join(., env_herit) %>%
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, dist_method_abbr),
         scheme_number = str_extract(scheme, "[0-9]{1,2}")) %>%
  mutate_at(vars(scheme, val_environment, nTrainEnv), as.factor) %>%
  filter(model %in% dist_method_abbr_use)

## Just take averages, since everything is balanced
cv0_pov0_predictions_analysis <- cv0_pov0_predictions_out %>% 
  group_by(trait, set, model, nTrainEnv, scheme) %>%
  summarize(fit = mean(accuracy)) %>% 
  ungroup()


## Edit some columns for plotting
cv0_pov0_predictions_analysis_toplot <- cv0_pov0_predictions_analysis %>%
  # unnest(model_nEnv_effects) %>%
  mutate(set = factor(str_replace_all(set, set_replace), levels = set_replace),
         model = factor(model, levels = rev(dist_method_abbr_use)),
         nTrainEnv = parse_number(nTrainEnv),
         scheme = toupper(scheme),
         size = model == "Random")


## Subset the final training environment accuracy
cv0_pov0_rank_pred_all_data <- cv0_pov0_predictions_analysis_toplot %>%
  group_by(set, trait, scheme) %>% 
  filter(model == model[1], nTrainEnv == max(nTrainEnv)) %>%
  ungroup() %>%
  select(set, trait, scheme, accuracy = fit)

cv0_pov0_predictions_analysis_toplot1 <- left_join(cv0_pov0_predictions_analysis_toplot, cv0_pov0_rank_pred_all_data)




g_rank_pred_0 <- cv0_pov0_predictions_analysis_toplot1 %>%
  mutate(group = paste0(set, " - ", scheme)) %>%
  ggplot(aes(x = nTrainEnv, y = fit, color = model, group = model)) +
  geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
  geom_line(aes(size = size)) +
  facet_grid(trait ~ group, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
  scale_color_manual(values = dist_colors_use, name = "Distance\nmeasure",
                     guide = guide_legend(nrow = 3, keyheight = unit(0.5, "line"), reverse = TRUE)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_size_manual(values = c(0.5, 1), guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
  scale_x_continuous(breaks = pretty, name = "Number of training set environments") +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom", strip.placement = "outside")

## Save
ggsave(filename = "cumulative_env_pred_cv0_pov0.jpg", plot = g_rank_pred_0, path = fig_dir, width = 7, height = 6, dpi = 1000)





## Reorganize by splitting on LOEO versus TF
cv_pov_rank_predictions1 <- bind_rows(cv0_pov0_predictions_analysis_toplot1, cv00_pov00_environment_rank_predictions_analysis_toplot) %>%
  mutate(scheme = ifelse(scheme == "POCV0", "POV0", scheme),
         scheme = factor(scheme, levels = cv_replace)) %>%
  droplevels()



### Stats

## What is the greatest decline in accuracy for CV and AMMI/PD?
# First find the peak for each CV/distance measure, trait
cv_pov_rank_predictions1 %>%
  filter(set == "Leave-one-environment-out", scheme %in% c("CV0", "CV00"), model %in% c("AMMI", "PD")) %>% 
  group_by(trait, set, model, scheme) %>% 
  mutate(max = max(fit)) %>% 
  filter(nTrainEnv == max(nTrainEnv)) %>% 
  mutate(decline = max - fit, per_decline = decline / max) %>%
  arrange(decline)

# trait       set                       model nTrainEnv scheme   fit size  accuracy    min   max  decline
# 1 HeadingDate Leave-one-environment-out AMMI         25 CV00   0.799 FALSE    0.799  0.761 0.799 0.000157
# 2 HeadingDate Leave-one-environment-out AMMI         25 CV0    0.882 FALSE    0.882 NA     0.885 0.00252 
# 3 HeadingDate Leave-one-environment-out PD           25 CV00   0.799 FALSE    0.799  0.761 0.805 0.00572 
# 4 HeadingDate Leave-one-environment-out PD           25 CV0    0.882 FALSE    0.882 NA     0.897 0.0147  
# 5 PlantHeight Leave-one-environment-out AMMI         26 CV0    0.749 FALSE    0.749 NA     0.810 0.0606  
# 6 PlantHeight Leave-one-environment-out AMMI         26 CV00   0.411 FALSE    0.406  0.244 0.496 0.0848  
# 7 PlantHeight Leave-one-environment-out PD           26 CV0    0.749 FALSE    0.749 NA     0.838 0.0883  
# 8 GrainYield  Leave-one-environment-out AMMI         22 CV0    0.609 FALSE    0.609 NA     0.700 0.0908  
# 9 PlantHeight Leave-one-environment-out PD           26 CV00   0.406 FALSE    0.406  0.244 0.507 0.101   
# 10 GrainYield  Leave-one-environment-out AMMI         22 CV00   0.393 FALSE    0.392  0.225 0.516 0.123   
# 11 GrainYield  Leave-one-environment-out PD           22 CV0    0.609 FALSE    0.609 NA     0.743 0.134   
# 12 GrainYield  Leave-one-environment-out PD           22 CV00   0.392 FALSE    0.392  0.225 0.549 0.157
# 











  

g_rank_pred_list <- cv_pov_rank_predictions1 %>%
  split(.$set) %>%
  map(~ggplot(data = ., aes(x = nTrainEnv, y = fit, color = model, group = model)) +
        geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
        geom_line(aes(size = size)) +
        facet_grid(trait ~ scheme, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
        scale_color_manual(values = dist_colors_use, name = "Distance\nmeasure",
                           guide = guide_legend(nrow = 3, keyheight = unit(0.5, "line"), reverse = TRUE, order = 1)) +
        scale_linetype_manual(values = 2, name = NULL) +
        scale_size_manual(values = c(0.5, 1), guide = FALSE) +
        scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
        scale_x_continuous(breaks = pretty, name = expression("Number of training set environments ("*italic(N[TE])*")")) +
        labs(subtitle = unique(.$set))  +
        theme_presentation2(base_size = 10) +
        theme(legend.position = "bottom") )


## Combine plots
# g_rank_pred_all <- plot_grid(plotlist = map(g_rank_pred_list, ~. + theme(legend.position = "none")), ncol = 1, labels = LETTERS[seq_along(g_rank_pred_list)])
g_rank_pred_all <- plot_grid(plotlist = map(g_rank_pred_list, ~. + theme(legend.position = "none")), ncol = 2, labels = LETTERS[seq_along(g_rank_pred_list)],
                             rel_widths = c(1, 0.75))

g_rank_pred_all1 <- plot_grid(g_rank_pred_all, get_legend(g_rank_pred_list[[1]]), ncol = 1, rel_heights = c(1, 0.07))
# ggsave(filename = "cumulative_env_pred_all.jpg", plot = g_rank_pred_all1, path = fig_dir, height = 8, width = 7, dpi = 1000)
ggsave(filename = "cumulative_env_pred_all.jpg", plot = g_rank_pred_all1, path = fig_dir, height = 5, width = 10, dpi = 1000)



## Two figures
ggsave(filename = "cumulative_env_pred_all_LOEO.jpg", plot = g_rank_pred_list$`Leave-one-environment-out` , path = fig_dir, height = 5.5, width = 7, dpi = 1000)
ggsave(filename = "cumulative_env_pred_all_TF.jpg", plot = g_rank_pred_list$`Time-forward`, path = fig_dir, height = 5.5, width = 5, dpi = 1000)


### Viewports for LOEO
g_vp_LOEO_list <- cv_pov_rank_predictions1 %>%
  filter(trait == "HeadingDate", set == "Leave-one-environment-out") %>%
  split(.$scheme) %>%
  map2(.x = ., .y = LETTERS[seq_along(.)], ~gg_vp(.x, let = " "))
  

## Coordinates for viewports
vp1 <- viewport(x = 0.20, y = 0.52, width = 0.215, height = 0.17)
vp2 <- viewport(x = 0.45, y = 0.59, width = 0.215, height = 0.12)
vp3 <- viewport(x = 0.65, y = 0.52, width = 0.215, height = 0.17)
vp4 <- viewport(x = 0.88, y = 0.56, width = 0.215, height = 0.17)

vp_list <- list(vp1, vp2, vp3, vp4)


jpeg(filename = file.path(fig_dir, "cumulative_env_pred_all_LOEO_paper.jpg"), width = 7, height = 5.5, units = "in", res = 1000)
print(g_rank_pred_list$`Leave-one-environment-out`)
for(i in seq_along(g_vp_LOEO_list)[-2]) print(g_vp_LOEO_list[[i]], vp = vp_list[[i]])
dev.off()



### Viewports for TF
g_vp_TF_list <- cv_pov_rank_predictions1 %>%
  filter(trait == "HeadingDate", set == "Time-forward") %>%
  split(.$scheme) %>%
  map2(.x = ., .y = LETTERS[seq_along(.)], ~gg_vp(.x, let = " "))


## Coordinates for viewports
vp1 <- viewport(x = 0.235, y = 0.52, width = 0.21, height = 0.17)
vp2 <- viewport(x = 0.45, y = 0.59, width = 0.21, height = 0.12)
vp3 <- viewport(x = 0.665, y = 0.52, width = 0.21, height = 0.17)
vp4 <- viewport(x = 0.88, y = 0.56, width = 0.21, height = 0.17)

vp_list <- list(vp1, vp2, vp3, vp4)


jpeg(filename = file.path(fig_dir, "cumulative_env_pred_all_TF_paper.jpg"), width = 5, height = 5.5, units = "in", res = 1000)
print(g_rank_pred_list$`Time-forward`)
for(i in seq_along(g_vp_TF_list)[-2]) print(g_vp_TF_list[[i]], vp = vp_list[[i]])
dev.off()










#### Environmental cluster predictions ####

cluster_pred <- ls(pattern = "[0-9]{1,2}_cluster_predictions") %>% 
  subset(., map_lgl(., ~inherits(get(.), "data.frame")))

cv00_cluster_predictions <- cv00_cluster_predictions %>%
  gather(scheme, accuracy, cv00, pocv00)

## Get and tidy
cluster_pred_df <- map(cluster_pred, get) %>%
  map(., ~mutate_at(., vars(which(names(.) %in% "rep")), as.character)) %>%
  map(., ~mutate_at(., vars(which(names(.) %in% "environment")), as.character))

cluster_pred_df[map_lgl(cluster_pred_df, ~"out" %in% names(.))] <- cluster_pred_df[map_lgl(cluster_pred_df, ~"out" %in% names(.))] %>%
  map(~mutate(., accuracy = map_dbl(out, "base")) %>% select(., -out))

## Add CV/POV number
cluster_pred_df <- bind_rows(cluster_pred_df) %>% 
  mutate(scheme_number = as.character(str_extract(scheme, "[0-9]{1,2}")),
         model = ifelse(model == "pheno_location_dist", "pheno_loc_dist", model))


  
cluster_predictions_base <- cluster_pred_df %>%
  left_join(., env_herit) %>% 
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability)) %>%
  mutate(model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr),
         zscore = accuracy) %>%
         # zscore = ztrans(accuracy)) %>% 
  mutate_at(vars(cluster, val_environment, scheme), as.factor) %>%
  filter(model %in% dist_method_abbr_use)


## Model formula
## Models fitted individually for each set and trait
## Fixed effect of model and number of training environments
## zscore ~ model + nTrainEnv + (1|cluster:model) + (1|val_environment:cluster)
## 



## Analyze separately
separate_cv_pov_cluster_predictions_analysis <- cluster_predictions_base %>%
  # group_by(set, trait, scheme) %>%
  group_by(set, trait, scheme_number) %>%
  nest() %>%
  mutate(model_rep = map_lgl(data, ~!all(is.na(.$rep)))) %>%
  mutate(out = list(NULL))

## Loop for each row
for (i in seq(nrow(separate_cv_pov_cluster_predictions_analysis))) {
  
  df <- separate_cv_pov_cluster_predictions_analysis$data[[i]]
  
  ## One model
  ## Model
  fit <- lmer(accuracy ~ 1 + model + scheme + model:scheme + (1|cluster:model) + (1|val_environment:cluster:model) + (1|nTrainEnv:cluster:model), data = df)
  
  effects_keep <- as.data.frame(Effect(c("model", "scheme"), fit))
  
  fit_summary <- data_frame(
    fitted = list(fit),
    scheme_effects = list(effects_keep),
    anova = list(tidy(anova(fit))),
    ranova = list(tidy(ranova(fit)))
  )
  
  separate_cv_pov_cluster_predictions_analysis$out[[i]] <- fit_summary
  
  
}


## What random effects are significant?
separate_cv_pov_cluster_predictions_analysis %>% 
  mutate(ranova = map(out, "ranova")) %>% 
  unnest(ranova) %>% 
  unnest(ranova) %>%
  select(set, trait, scheme_number, term, LRT, p.value) %>%
  filter(term != "<none>") %>% 
  as.data.frame() %>%
  filter(p.value <= alpha)





## Plot all
separate_cv_pov_cluster_predictions_analysis1 <- separate_cv_pov_cluster_predictions_analysis %>%
  mutate(effects = map(out, ~.$scheme_effects[[1]])) %>%
  unnest(effects) %>%
  mutate(scheme = factor(toupper(scheme), levels = cv_replace),
         model = factor(model, levels = dist_method_abbr_use),
         set = factor(str_replace_all(set, set_replace), levels = set_replace)) %>%
  left_join(., select(separate_cv_pov_predictions_analysis1, -lower, -upper, -se)) %>%
  mutate(group = paste0(set, " - ", scheme))
           

g_cv_pov_cluster <- separate_cv_pov_cluster_predictions_analysis1 %>%
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper)) +
  geom_hline(aes(yintercept = accuracy, lty = "All data"), color = "black") +
  geom_errorbar(width = 0.5) +
  geom_point(aes(color = model)) +
  scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
  scale_color_manual(values = dist_colors_use, guide = FALSE) +
  scale_linetype_manual(values = 2, name = NULL) +
  xlab("Distance method") +
  facet_grid(trait ~ group, scales = "free_x", space = "free_x", labeller = labeller(set = str_to_title, trait = str_add_space),
             switch = "y") +
  theme_presentation2(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside")

ggsave(filename = "all_cv_pov_cluster_predictions.jpg", plot = g_cv_pov_cluster, path = fig_dir, width = 10, height = 7, dpi = 1000)



# ### Add predictions with random environmental subsets
# separate_cv_pov_random_predictions_analysis2 <- separate_cv_pov_random_predictions_analysis1 %>%
#   select(set, trait, scheme, nTrainEnv, accuracy) %>% 
#   mutate(nTrainEnv = paste0("accuracy", nTrainEnv),
#          set = str_to_title(set)) %>%
#   spread(nTrainEnv, accuracy)
# 
# 
# dist_colors_random <- wesanderson::wes_palette(name = "Darjeeling1", 3)
# dist_colors_use_random <- c(dist_colors_use, "2 env" = dist_colors_random[1], "5 env" = dist_colors_random[2], "8 env" = dist_colors_random[3])
# 
# g_cv_pov_cluster1 <- separate_cv_pov_cluster_predictions_analysis1 %>%
#   left_join(., separate_cv_pov_random_predictions_analysis2) %>%
#   ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper)) +
#   geom_hline(aes(yintercept = accuracy, lty = "All data"), color = "black") +
#   geom_hline(aes(yintercept = accuracy2, lty = "2 env"), color = "black") +
#   geom_hline(aes(yintercept = accuracy5, lty = "5 env"), color = "black") +
#   geom_hline(aes(yintercept = accuracy8, lty = "8 env"), color = "black") +
#   geom_errorbar(width = 0.5) +
#   geom_point(aes(color = model)) +
#   scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
#   scale_color_manual(values = dist_colors_use_random, guide = FALSE) +
#   scale_linetype_manual(values = seq(2, 5), name = NULL) +
#   xlab("Distance method") +
#   facet_grid(trait ~ set + scheme, scales = "free_x", space = "free_x", labeller = labeller(set = str_to_title, trait = str_add_space)) +
#   theme_presentation2(base_size = 10) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggsave(filename = "all_cv_pov_cluster_predictions_with_random.jpg", plot = g_cv_pov_cluster1, path = fig_dir, width = 10, height = 7, dpi = 1000)
# 


## Paper version
g_cv_pov_cluster_list <- separate_cv_pov_cluster_predictions_analysis1 %>%
  mutate(scheme = ifelse(scheme == "POCV0", "POV0", as.character(scheme)),
         scheme = factor(scheme, levels = cv_replace)) %>%
  filter(!str_detect(scheme, "POCV")) %>%
  split(.$set) %>%
  map(~ggplot(data = ., aes(x = model, y = fit, ymin = lower, ymax = upper)) +
        geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data"), color = "grey") +
        geom_errorbar(width = 0.5, lwd = 0.5) +
        geom_point(aes(color = model), size = 1.25) +
        scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*r[MG]*")")) +
        scale_color_manual(values = dist_colors_use, guide = FALSE) +
        scale_linetype_manual(values = 1, name = NULL) +
        xlab("Distance measure") +
        facet_grid(trait ~ scheme, scales = "free", space = "free_x", labeller = labeller(set = str_to_title, trait = str_add_space), switch = "y") +
        labs(subtitle = unique(.$set))  +
        theme_presentation2(base_size = 8) +
        theme(axis.text.x = element_text(angle = 50, hjust = 1), legend.position = "bottom") )

g_cv_pov_cluster2 <- plot_grid(plotlist = map(g_cv_pov_cluster_list, ~. + theme(legend.position = "none")),
                               nrow = 1, labels = LETTERS[seq_along(g_cv_pov_cluster_list)], align = "hv")
g_cv_pov_cluster21 <- plot_grid(g_cv_pov_cluster2, get_legend(g_cv_pov_cluster_list[[1]]), ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "all_cv_pov_cluster_predictions_paper.jpg", plot = g_cv_pov_cluster21, path = fig_dir, width = 7, height = 3.5, dpi = 1000)










### Select CV00 and PV00

# Create a legend key
legend_key <- paste0(names(dist_colors_use), ": ", c("", "Phenotypic distance", "Location PD", "Geographic distance", "All covariates", "Mean-correlated ECs", "IPCA-correlated ECs"))


g_cv00_pov00_cluster <- separate_cv_pov_cluster_predictions_analysis1 %>%
  filter(scheme %in% c("CV00", "POV00")) %>%
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper)) +
  geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data"), color = "grey") +
  geom_errorbar(width = 0.5) +
  geom_point(aes(color = model), size = 1.5) +
  scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
  # scale_color_manual(values = dist_colors_use, name = "Cluster method", labels = legend_key) +
  scale_color_manual(values = dist_colors_use, name = "Cluster method", guide = FALSE) +
  scale_linetype_manual(values = 1, name = NULL) +
  xlab("Distance measure") +
  facet_grid(trait ~ set + scheme, scales = "free_x", space = "free_x", switch = "y", labeller = labeller(set = str_to_title, trait = str_add_space)) +
  # facet_grid(trait ~ group, scales = "free_x", space = "free_x", switch = "y", labeller = labeller(set = str_to_title, trait = str_add_space)) +
  theme_presentation2(base_size = 8) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right", legend.margin = margin())
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", legend.margin = margin())



ggsave(filename = "cv00_pov00_cluster_predictions.jpg", plot = g_cv00_pov00_cluster, path = fig_dir, width = 3.5, height = 5, dpi = 1000)



##












##### Random cluster prediction #####
##### 

# This will be analyzed by calculating the difference between the cluster prediction accuracy 
# and each replicate of the reandom cluster prediction accuracy, then modeling that difference

## Edit, then combine df
cv00_cluster_random_predictions <- cv00_cluster_random_predictions %>%
  gather(scheme, accuracy, cv00, pocv00)


random_cluster_pred <- ls(pattern = "[0-9]{1,2}_cluster_random_predictions") %>% 
  subset(., map_lgl(., ~inherits(get(.), "data.frame")))

## Get and tidy
random_cluster_pred_df <- map(random_cluster_pred, get) %>%
  map(., ~mutate_at(., vars(which(names(.) %in% "rep")), as.character)) %>%
  map(., ~mutate_at(., vars(which(names(.) %in% "environment")), as.character)) %>%
  bind_rows() %>% 
  mutate(scheme_number = as.character(str_extract(scheme, "[0-9]{1,2}")),
         model = ifelse(model == "pheno_location_dist", "pheno_loc_dist", model))

## Calculate prediction accuracy
random_cluster_predictions_base <- random_cluster_pred_df %>%
  left_join(., env_herit) %>% 
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability)) %>%
  mutate(model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr),
         zscore = accuracy) %>%
  # zscore = ztrans(accuracy)) %>% 
  mutate_at(vars(cluster, val_environment, scheme), as.factor) %>%
  filter(model %in% dist_method_abbr_use)

## Combine with the original clustering predictions, then calculate the advantage in accuracy
## First calculate averages over reps for the original predictions
random_base_cluster_predictions_base <- cluster_predictions_base %>% 
  group_by(set, trait, model, val_environment, nTrainEnv, cluster, scheme, scheme_number) %>% 
  summarize(base_accuracy = mean(accuracy)) %>%
  ungroup() %>%
  left_join(random_cluster_predictions_base,  by = c("set", "trait", "model", "val_environment", "cluster", "scheme", "scheme_number")) %>%
  mutate(advantage = base_accuracy - accuracy,
         scheme = as.factor(scheme)) %>%
  rename(nTrainEnv = nTrainEnv.x)




## Visualize
random_base_cluster_predictions_base %>%
  ggplot(aes(x = advantage, fill = model)) +
  geom_density() +
  facet_grid(trait ~ set)






## Model formula
## Models fitted individually for each set and trait
## Fixed effect of model and number of training environments
## zscore ~ model + nTrainEnv + (1|cluster:model) + (1|val_environment:cluster)
## 



## Analyze separately
random_cluster_predictions_analysis <- random_base_cluster_predictions_base %>%
  # group_by(set, trait, scheme) %>%
  group_by(set, trait, scheme_number) %>%
  nest() %>%
  filter(!map_lgl(data, ~all(is.na(.$advantage)))) %>%
  mutate(out = list(NULL))

## Loop for each row
for (i in seq(nrow(random_cluster_predictions_analysis))) {
  
  df <- random_cluster_predictions_analysis$data[[i]]
  
  ## Model
  fit <- lmer(advantage ~ 1 + model + scheme + model:scheme + (1|cluster:model) + (1|val_environment:cluster:model), data = df)
  
  effects_keep <- as.data.frame(Effect(c("model", "scheme"), fit))
    
  
  fit_summary <- data_frame(
    fitted = list(fit),
    scheme_effects = list(effects_keep),
    anova = list(tidy(anova(fit))),
    ranova = list(tidy(ranova(fit)))
  )
  
  random_cluster_predictions_analysis$out[[i]] <- fit_summary
  
  
}




## Plot all
random_cluster_predictions_analysis1 <- random_cluster_predictions_analysis %>%
  mutate(effects = map(out, ~.$scheme_effects[[1]])) %>%
  unnest(effects) %>%
  mutate(scheme = factor(toupper(scheme), levels = cv_replace),
         model = factor(model, levels = dist_method_abbr_use),
         scheme = ifelse(scheme == "POCV0", "POV0", as.character(scheme)),
         scheme = factor(scheme, levels = cv_replace),
         set = factor(str_replace_all(set, set_replace), levels = set_replace))


g_random_cv_pov_cluster <- random_cluster_predictions_analysis1 %>%
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(width = 0.5) +
  geom_point(aes(color = model)) +
  scale_y_continuous(breaks = pretty, name = "Advantage of clustering over random environments") +
  scale_color_manual(values = dist_colors_use, guide = FALSE) +
  scale_linetype_manual(values = 2, name = NULL) +
  xlab("Distance method") +
  facet_grid(trait ~ set + scheme, scales = "free", space = "free_x", labeller = labeller(set = str_to_title, trait = str_add_space)) +
  theme_presentation2(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "random_cv_pov_cluster_predictions.jpg", plot = g_random_cv_pov_cluster, path = fig_dir, width = 10, height = 7, dpi = 1000)



## Publication version is just CV and POV
g_random_cv_pov_cluster_list <- random_cluster_predictions_analysis1 %>%
  filter(!str_detect(scheme, "POCV")) %>%
  split(.$set) %>%
  map(~ggplot(data = ., aes(x = model, y = fit, ymin = lower, ymax = upper)) +
        geom_hline(yintercept = 0, color = "grey", lty = 5) +
        geom_errorbar(width = 0.5, lwd = 0.5) +
        geom_point(aes(color = model), size = 1.25) +
        scale_y_continuous(breaks = pretty, name = "Prediction accuracy versus random clusters") +
        scale_color_manual(values = dist_colors_use, guide = FALSE) +
        xlab("Distance measure") +
        labs(subtitle = unique(.$set)) +
        facet_grid(trait ~ scheme, scales = "free", space = "free_x", switch = "y", labeller = labeller(trait = str_add_space)) +
        theme_presentation2(base_size = 8) +
        theme(axis.text.x = element_text(angle = 50, hjust = 1)) )

g_random_cv_pov_cluster1 <- plot_grid(plotlist = g_random_cv_pov_cluster_list, nrow = 1, labels = LETTERS[seq_along(g_random_cv_pov_cluster_list)], 
                                      align = "hv")

ggsave(filename = "random_cv_pov_cluster_predictions_paper.jpg", plot = g_random_cv_pov_cluster1, path = fig_dir, width = 7, height = 3.5, dpi = 1000)




## Combine figure of cluster accuracy with random
g_cv_pov_cluster_LOEO <- plot_grid(g_cv_pov_cluster_list$`Leave-one-environment-out` + theme(legend.position = c(0.75, 1.12), legend.margin = margin()), 
                                   g_random_cv_pov_cluster_list$`Leave-one-environment-out` + labs(subtitle = NULL),
                                   nrow = 1, align = "hv", axis = "tblr", labels = LETTERS[1:2])

g_cv_pov_cluster_TF <- plot_grid(g_cv_pov_cluster_list$`Time-forward` + theme(legend.position = c(0.75, 1.12), legend.margin = margin()), 
                                 g_random_cv_pov_cluster_list$`Time-forward` + labs(subtitle = NULL),
                                 nrow = 1, align = "hv", axis = "tblr", labels = LETTERS[1:2])


ggsave(filename = "all_cv_pov_and_random_cluster_predictions_LOEO_paper.jpg", plot = g_cv_pov_cluster_LOEO, path = fig_dir, width = 7, height = 4, dpi = 1000)
ggsave(filename = "all_cv_pov_and_random_cluster_predictions_TF_paper.jpg", plot = g_cv_pov_cluster_TF, path = fig_dir, width = 7, height = 4, dpi = 1000)










##
