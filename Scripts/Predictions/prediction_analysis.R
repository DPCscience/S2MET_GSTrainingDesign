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
library(car)
library(cowplot)

# Load the environmental distance df
load(file.path(result_dir, "distance_method_results.RData"))

# Load results
# load(file.path(result_dir, "all_predictions.RData"))
load(file.path(result_dir, "all_predictions_00.RData"))
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
  facet_grid(trait ~ set, scales = "free_x", space = "free_x", labeller = labeller(set = str_to_title, trait = str_add_space)) +
  theme_presentation2(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "all_cv_pov_predictions.jpg", plot = g_cv_pov_all_data, path = fig_dir, width = 5, height = 7, dpi = 1000)







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

## Analyze pov00 first
pov00_predictions_out <- bind_rows(pov00_environment_rank_predictions, pov00_environment_rank_random_predictions) %>%
  mutate(model = ifelse(str_detect(model, "sample"), "sample", model)) %>%
  select(-train, -test) %>%
  left_join(., env_herit) %>%
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, dist_method_abbr)) %>%
  mutate_at(vars(scheme, val_environment, nTrainEnv), as.factor) %>%
  filter(model %in% dist_method_abbr_use)


## Fit a model per trait

# Fit a model per trait
pov00_predictions_analysis <- pov00_predictions_out %>%
  group_by(set, trait, scheme) %>%
  nest() %>%
  mutate(out = list(NULL))

## Iterate over sets and traits - the summaries of LMER don't work using do() or map(), so we have to use a loop
for (i in seq(nrow(pov00_predictions_analysis))) {
  
  df <- pov00_predictions_analysis$data[[i]] %>% droplevels()
  
  fit <- lmer(formula = accuracy ~ 1 + model + nTrainEnv + model:nTrainEnv + (1|val_environment) + (1|val_environment:model) + 
                (1|val_environment:nTrainEnv), data = df)
  
  fit_summary <- data_frame(
    fitted = list(fit),
    model_nEnv_effects = list(Effect(focal.predictors = c("model", "nTrainEnv"), fit))) %>% 
    mutate_at(vars(contains("effects")), ~map(., ~as.data.frame(.)))
  
  # Add the summary to the DF
  pov00_predictions_analysis$out[[i]] <- fit_summary
  
}

pov00_predictions_analysis <- unnest(pov00_predictions_analysis, out)


## Plot
pov00_predictions_analysis_toplot <- pov00_predictions_analysis %>%
  unnest(model_nEnv_effects) %>%
  mutate(set = factor(str_replace_all(set, set_replace), levels = set_replace),
         model = factor(model, levels = dist_method_abbr_use),
         nTrainEnv = parse_number(nTrainEnv),
         scheme = toupper(scheme),
         size = model == "Random") %>%
  ## Add results from all environment predictions
  left_join(., select(separate_cv_pov_predictions_analysis1, set, trait, scheme, accuracy))
  
g_pov00_rank_pred <- pov00_predictions_analysis_toplot %>%
  ggplot(aes(x = nTrainEnv, y = fit, color = model, group = model)) +
  geom_hline(aes(yintercept = accuracy, lty = "Accuracy\nusing\nall data")) + 
  geom_line(aes(size = size)) +
  facet_grid(trait ~ set, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
  scale_color_manual(values = dist_colors_use, name = "Distance\nmeasure", guide = guide_legend(nrow = 3)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_size_manual(values = c(0.5, 1), guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
  scale_x_continuous(breaks = pretty, name = "Number of training set environments") +
  theme_presentation2(base_size = 10) +
  theme(legend.position = "bottom", strip.placement = "outside")

ggsave(filename = "cumulative_env_pred_pov00.jpg", plot = g_pov00_rank_pred, path = fig_dir, width = 4.5, height = 6, dpi = 1000)




## What is the best/worst accuracy achievable for each model?
pov00_predictions_analysis_toplot %>%
  mutate(advantage = fit - accuracy, perc_advantage = advantage / accuracy) %>%
  group_by(set, trait, model) %>%
  mutate_at(vars(advantage), funs(min, max)) %>%
  filter(advantage == min | advantage == max) %>%
  ungroup() %>%
  arrange(set, trait, model)









## Cross-validation prediction
cv00_predictions_out <- bind_rows(cv00_environment_rank_predictions, cv00_environment_rank_random_predictions) %>%
  mutate(model = ifelse(str_detect(model, "sample"), "sample", model)) %>%
  left_join(., env_herit) %>%
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, dist_method_abbr)) %>%
  mutate_at(vars(scheme, val_environment, nTrainEnv), as.factor) %>%
  filter(model %in% dist_method_abbr_use)


## Fit a model per trait

# Fit a model per trait
cv00_predictions_analysis <- cv00_predictions_out %>%
  group_by(set, trait) %>%
  nest() %>%
  mutate(out = list(NULL))

## Iterate over sets and traits - the summaries of LMER don't work using do() or map(), so we have to use a loop
for (i in seq(nrow(cv00_predictions_analysis))) {
  
  df <- cv00_predictions_analysis$data[[i]] %>% # mutate(nTrainEnv = parse_number(nTrainEnv)) %>%
    droplevels()
  
  fit <- lmer(formula = accuracy ~ 1 + model + nTrainEnv + model:nTrainEnv + (1|val_environment) + (1|val_environment:model) + 
                (1|val_environment:nTrainEnv) + (1|val_environment:nTrainEnv:model), data = df)
  
  fit_summary <- data_frame(
    fitted = list(fit),
    model_nEnv_effects = list(Effect(focal.predictors = c("model", "nTrainEnv"), fit))) %>% 
    mutate_at(vars(contains("effects")), ~map(., ~as.data.frame(.)))
  
  # Add the summary to the DF
  cv00_predictions_analysis$out[[i]] <- fit_summary
  
}

cv00_predictions_analysis <- unnest(cv00_predictions_analysis, out)






















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
  left_join(., select(separate_cv_pov_predictions_analysis1, -lower, -upper, -se))
           

g_cv_pov_cluster <- separate_cv_pov_cluster_predictions_analysis1 %>%
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper)) +
  geom_hline(aes(yintercept = accuracy, lty = "All data"), color = "black") +
  geom_errorbar(width = 0.5) +
  geom_point(aes(color = model)) +
  scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
  scale_color_manual(values = dist_colors_use, guide = FALSE) +
  scale_linetype_manual(values = 2, name = NULL) +
  xlab("Distance method") +
  facet_grid(trait ~ set + scheme, scales = "free_x", space = "free_x", labeller = labeller(set = str_to_title, trait = str_add_space),
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
        geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data"), color = "black") +
        geom_errorbar(width = 0.5) +
        geom_point(aes(color = model)) +
        scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
        scale_color_manual(values = dist_colors_use, guide = FALSE) +
        scale_linetype_manual(values = 2, name = NULL) +
        xlab("Distance measure") +
        facet_grid(trait ~ scheme, scales = "free", space = "free_x", labeller = labeller(set = str_to_title, trait = str_add_space), switch = "y") +
        labs(subtitle = unique(.$set))  +
        theme_presentation2(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", strip.placement = "outside") )

g_cv_pov_cluster2 <- plot_grid(plotlist = map(g_cv_pov_cluster_list, ~. + theme(legend.position = "none")),
                               nrow = 2, labels = LETTERS[seq_along(g_cv_pov_cluster_list)], align = "hv")
g_cv_pov_cluster21 <- plot_grid(g_cv_pov_cluster2, get_legend(g_cv_pov_cluster_list[[1]]), ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "all_cv_pov_cluster_predictions_paper.jpg", plot = g_cv_pov_cluster21, path = fig_dir, width = 4.5, height = 8, dpi = 1000)










### Select CV00 and PV00

# Create a legend key
legend_key <- paste0(names(dist_colors_use), ": ", c("", "Phenotypic distance", "Location PD", "Geographic distance", "All covariates", "Mean-correlated ECs", "IPCA-correlated ECs"))



g_cv00_pov00_cluster <- separate_cv_pov_cluster_predictions_analysis1 %>%
  filter(scheme %in% c("CV00", "POV00"), set == "Complete") %>%
  mutate(scheme = str_replace_all(scheme, c("CV00" = "Cross-validation", "POV00" = "Parent-offspring validation"))) %>%
  ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper)) +
  geom_hline(aes(yintercept = accuracy, lty = "All data"), color = "black") +
  geom_errorbar(width = 0.5) +
  geom_point(aes(color = model), size = 2) +
  scale_y_continuous(breaks = pretty, name = "Prediction accuracy") +
  scale_color_manual(values = dist_colors_use, name = "Custer method", labels = legend_key) +
  scale_linetype_manual(values = 2, name = NULL) +
  xlab("Cluster method") +
  # facet_grid(trait ~ set + scheme, scales = "free_x", space = "free_x", switch = "y", labeller = labeller(set = str_to_title, trait = str_add_space)) +
  facet_grid(trait ~ scheme, scales = "free_x", space = "free_x", switch = "y", labeller = labeller(set = str_to_title, trait = str_add_space)) +
  theme_presentation2(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right", legend.margin = margin())


ggsave(filename = "cv00_pov00_cluster_predictions.jpg", plot = g_cv00_pov00_cluster, path = fig_dir, width = 7, height = 5, dpi = 1000)



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
         model = factor(model, levels = dist_method_abbr_use), set = str_to_title(set))


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
  mutate(scheme = ifelse(scheme == "POCV0", "POV0", as.character(scheme)),
         scheme = factor(scheme, levels = cv_replace)) %>%
  filter(!str_detect(scheme, "POCV")) %>%
  split(.$set) %>%
  map(~ggplot(data = ., aes(x = model, y = fit, ymin = lower, ymax = upper)) +
        geom_hline(yintercept = 0) +
        geom_errorbar(width = 0.5) +
        geom_point(aes(color = model)) +
        scale_y_continuous(breaks = pretty, name = "Prediction accuracy advantage") +
        scale_color_manual(values = dist_colors_use, guide = FALSE) +
        scale_linetype_manual(values = 2, name = NULL) +
        xlab("Distance method") +
        facet_grid(trait ~ scheme, scales = "free", space = "free_x", switch = "y", labeller = labeller(set = str_to_title, trait = str_add_space)) +
        theme_presentation2(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside") )

g_random_cv_pov_cluster1 <- plot_grid(plotlist = g_random_cv_pov_cluster_list, nrow = 2, labels = LETTERS[seq_along(g_random_cv_pov_cluster_list)])

ggsave(filename = "random_cv_pov_cluster_predictions_paper.jpg", plot = g_random_cv_pov_cluster1, path = fig_dir, width = 5, height = 8, dpi = 1000)







##























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




## Parent-offspring cross-validation

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


##


































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








