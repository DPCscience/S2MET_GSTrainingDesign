## S2MET Predictions
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
library(modelr)
library(viridis)
library(ggalt)



## Other functions and plot modifiers

## Pretty editor
pretty1 <- function(x) {
  x1 <- pretty.default(x = x, n = 6, min.n = 5, high.u.bias = 0.1)
  x1[seq(1, length(x1), by = 3)]
}

## Plot editor
gg_vp <- function(x, let) {
  x1 <- mutate_at(x, vars(fit), list(~min, ~max))
  
  ggplot(data = x1, aes(x = nTrainEnv, y = fit, color = model, group = model)) +
    geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
    geom_line(aes(size = size)) +
    geom_text(aes(x = Inf, y = -Inf, label = let), hjust = 2, vjust = -1, color = "black") +
    scale_color_manual(values = dist_colors_use, name = "Similarity\nmeasure", guide = FALSE) +
    scale_linetype_manual(values = 2, name = NULL, guide = FALSE) +
    scale_size_manual(values = c(0.5, 1), guide = FALSE) +
    scale_y_continuous(breaks = pretty1, name = NULL, limits = c(unique(x1$min), unique(x1$max))) +
    scale_x_continuous(breaks = pretty, name = NULL, labels = NULL) +
    theme_presentation2(base_size = 8) +
    theme(panel.border = element_rect(color = "black"))
}




# Load the environmental distance df
load(file.path(result_dir, "distance_method_results.RData"))

# Load results
load(file.path(result_dir, "all_data_cluster_predictions.RData"))
# load(file.path(result_dir, "distance_rank_predictions_pov.RData"))
load(file.path(result_dir, "distance_rank_predictions.RData"))



## Shorten the heritability results
env_herit <- stage_one_data %>% 
  filter(!str_detect(trial, "S2C1")) %>% 
  distinct(trait, val_environment = environment, heritability)


## Color scheme for schemes - less challenging to most challenging
scheme_color <- set_names(rev(wesanderson::wes_palette("Zissou1")[-3]), c("CV0", "POV0", "CV00", "POV00"))










#### Cross-validation / POV using all data #####

## Adjust pov00 and cv00
pov00_predictions <- gather(pov00_predictions, scheme, accuracy, pov00)
cv00_predictions <- gather(cv00_predictions, scheme, accuracy, cv00, pocv00)


## Combine the all-data prediction results
all_data_preds <- ls(pattern = "[0-9]{1,2}_predictions$") %>% 
  subset(., map_lgl(., ~inherits(get(.), "data.frame")))


## Grab the dfs containing the all-data prediction results
## Rename some schemes
cv_pov_results <- map(all_data_preds, get) %>%
  map(., ~rename_at(., vars(which(names(.) %in% "environment")), ~str_replace(., pattern = "environment", "val_environment"))) %>%
  map_df(., ~mutate_at(., vars(which(names(.) %in% "rep")), as.character)) %>%
  mutate(scheme = ifelse(scheme == "pcv00", "pocv0", scheme)) %>%
  ## Add CV/POV number
  mutate(scheme_number = as.character(str_extract(scheme, "[0-9]{1,2}")))



## Add the environmental heritabilities to the predictive ability estimate
## calculate prediction accuracy
all_cv_pov_predictions_out <- cv_pov_results %>%
  filter(trait %in% traits) %>%
  left_join(., env_herit) %>%
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability),
         set = ifelse(is.na(set), "complete", set)) %>%
  mutate_at(vars(scheme, val_environment), as.factor)



## Calculate average prediction accuracy per set, trait, scheme,
## and environment
separate_cv_pov_predictions_analysis <- all_cv_pov_predictions_out %>% 
  group_by(set, trait, scheme, val_environment) %>% 
  summarize_at(vars(ability, accuracy), mean) %>%
  ungroup() %>%
  mutate(set = str_remove(set, "[0-9]{4}"),
         set = str_replace_all(set, set_replace),
         scheme = factor(toupper(scheme), levels = cv_replace))


## There should be the same number of environments for the LOEO and LOYO per
## trait and scheme
separate_cv_pov_predictions_analysis %>%
  group_by(trait, trait, scheme) %>%
  summarize(n_distinct(val_environment))

# Good


## Group all years for LOYO
## Summarize prediction accuracy across sets, traits, and schemes
separate_cv_pov_predictions_analysis1 <- separate_cv_pov_predictions_analysis %>%
  group_by(set, trait, scheme) %>%
  summarize_at(vars(ability, accuracy), mean) %>%
  ungroup() 



## Plot all-data predictions
# g_cv_pov_all_data <- separate_cv_pov_predictions_analysis1 %>%
g_cv_pov_all_data <- separate_cv_pov_predictions_analysis %>%
  ggplot(aes(x = scheme, y = accuracy)) +
  # geom_col() +
  geom_boxplot() +
  scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
  facet_grid(trait ~ set, scales = "free", space = "free_x", switch = "y", 
             labeller = labeller(set = str_to_title, trait = str_add_space)) +
  theme_presentation2(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "all_cv_pov_predictions.jpg", plot = g_cv_pov_all_data, path = fig_dir, width = 4, height = 5, dpi = 1000)



## Paper version of the plot
## Color code from likely more accurate (warm) to likely less accurate (cold)
## Only plot POV and CV 0 and 00
# g_cv_pov_all_data_paper <- separate_cv_pov_predictions_analysis1 %>%
g_cv_pov_all_data_paper <- separate_cv_pov_predictions_analysis %>%
  filter(scheme != "POCV00") %>%
  mutate(scheme = str_replace(scheme, "POCV", "POV"),
         scheme = factor(scheme, levels = c("CV0", "POV0", "CV00", "POV00"))) %>%
  ggplot(aes(x = scheme, y = accuracy, fill = scheme)) +
  # geom_col() +
  geom_boxplot(alpha = 0.5) +
  scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
  scale_fill_manual(values = scheme_color, name = "Validation\nscheme") +
  facet_grid(trait ~ set, scales = "free", space = "free_x", labeller = labeller(trait = str_add_space), switch = "y") +
  theme_presentation2(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_blank(),
        axis.title.x = element_blank(), legend.position = "bottom", panel.grid.major.y = element_line(color = "grey85"))

ggsave(filename = "all_cv_pov_predictions_paper.jpg", plot = g_cv_pov_all_data_paper, path = fig_dir, 
       width = 4, height = 5, dpi = 1000)





## Output a table
separate_cv_pov_predictions_analysis1 %>%
  filter(scheme != "POCV00") %>%
  mutate(scheme = str_replace(scheme, "POCV", "POV"),
         scheme = factor(scheme, levels = c("CV0", "POV0", "CV00", "POV00"))) %>%
  mutate(accuracy = round(accuracy, 2), annotation = accuracy) %>%
  # mutate_at(vars(accuracy, lower, upper), round, 2) %>%
  # mutate(annotation = paste0(accuracy, " (", lower, ", ", upper, ")")) %>%
  select(set, trait, scheme, annotation) %>%
  spread(scheme, annotation) %>%
  # Output a table
  write_csv(x = ., path = file.path(fig_dir, "all_cv_pov_predictions.csv"))














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
  left_join(., pov00_environment_rank_predictions) %>%
  filter(is.na(accuracy))
# Random
environment_rank_df %>% 
  filter(model == "great_circle_dist") %>% 
  select(-model) %>% 
  crossing(., model = paste0("sample", seq(max(pov00_environment_rank_random_predictions$rep)))) %>%
  left_join(., mutate(pov00_environment_rank_random_predictions, model = paste0("sample", rep))) %>%
  filter(is.na(accuracy))


## CV
# Regular
environment_rank_df %>%
  crossing(rep = seq(25)) %>%
  left_join(., cv00_environment_rank_predictions) %>%
  filter(is.na(cv00))

# Random
environment_rank_df %>% 
  filter(model == "great_circle_dist") %>% 
  select(-model) %>% 
  crossing(., rep = seq(25)) %>%
  mutate(model = paste0("sample", rep)) %>%
  left_join(., cv00_environment_rank_random_predictions) %>%
  filter(is.na(cv00))



##### 

## Analyze pov00 first
pov00_predictions_out <- pov00_environment_rank_random_predictions %>%
  mutate(model = paste0("sample", rep)) %>%
  bind_rows(pov00_environment_rank_predictions, .) %>%
  mutate(model = ifelse(str_detect(model, "sample"), "sample", model)) %>%
  ## Add environment heritability to correct predictive ability
  left_join(., env_herit) %>%
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, dist_method_abbr)) %>%
  mutate_at(vars(scheme, val_environment, nTrainEnv), as.factor) %>%
  filter(model %in% dist_method_abbr_use)


## Make sure that each environment within a set, trait, and model has the same
## number of training environments
pov00_predictions_out %>% 
  group_by(trait, set, model, val_environment) %>% 
  summarize_at(vars(nTrainEnv), ~max(as.numeric(.))) %>% 
  summarize(n = n_distinct(nTrainEnv)) %>%
  filter(n > 1)

# Good

 
## Just take averages, since everything is balanced
## Do not group LOYO by years, given different numbers of training envs
pov00_predictions_analysis <- pov00_predictions_out %>% 
  mutate(set = f_set_replace(set)) %>%
  group_by(trait, set, model, nTrainEnv, scheme) %>%
  summarize(fit = mean(accuracy)) %>% 
  ungroup()



## Plot
pov00_predictions_analysis_toplot <- pov00_predictions_analysis %>%
  # unnest(model_nEnv_effects) %>%
  mutate(model = factor(model, levels = rev(dist_method_abbr_use)),
         nTrainEnv = parse_number(as.character(nTrainEnv)),
         scheme = toupper(scheme),
         size = model == "Random")
  # Add accuracy using all data
  

## Subset the final training environment accuracy
pov00_rank_pred_all_data <- pov00_predictions_analysis_toplot %>%
  group_by(set, trait, scheme) %>% 
  filter(nTrainEnv == max(nTrainEnv)) %>%
  distinct(trait, set, scheme, accuracy = fit)


## Add this final accuracy back to the table - use this to set a 
## baseline
pov00_predictions_analysis_toplot1 <- left_join(pov00_predictions_analysis_toplot, pov00_rank_pred_all_data)


## Plot 
g_pov00_rank_pred <- pov00_predictions_analysis_toplot1 %>%
  filter(trait %in% traits) %>%
  ggplot(aes(x = nTrainEnv, y = fit, color = model, group = model)) +
  geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
  geom_line(aes(size = size)) +
  facet_grid(trait ~ set, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
  scale_color_manual(values = dist_colors_use, name = "Similarity\nmeasure",
                     guide = guide_legend(nrow = 2, keyheight = unit(0.5, "line"), reverse = TRUE)) +
  scale_linetype_manual(values = 2, name = NULL, guide = FALSE) +
  scale_size_manual(values = c(0.5, 1), guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
  scale_x_continuous(breaks = pretty, name = "Number of training set environments") +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom", strip.placement = "outside")

g_pov00_rank_pred1 <- g_pov00_rank_pred +
  scale_linetype_manual(values = 2, name = NULL, guide = guide_legend()) +
  scale_color_manual(values = dist_colors_use, guide = FALSE) +
  theme(legend.position = c(0.05, -0.1))

g_pov00_rank_pred2 <- plot_grid(g_pov00_rank_pred1, get_legend(g_pov00_rank_pred), ncol = 1, rel_heights = c(1, 0.07))


ggsave(filename = "cumulative_env_pred_pov00.jpg", plot = g_pov00_rank_pred2, path = fig_dir, width = 7, height = 5, dpi = 1000)





## What is the best/worst accuracy achievable for each model?
pov00_predictions_analysis_toplot1 %>%
  mutate(advantage = fit - accuracy, perc_advantage = advantage / accuracy) %>%
  group_by(set, trait, model) %>%
  mutate_at(vars(advantage), list(~min, ~max)) %>%
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



## Make sure that each environment within a set, trait, and model has the same
## number of training environments
cv00_predictions_out %>% 
  group_by(trait, set, model, val_environment) %>% 
  summarize_at(vars(nTrainEnv), ~max(as.numeric(.))) %>% 
  summarize(n = n_distinct(nTrainEnv)) %>%
  filter(n > 1)

# good


## Just take averages, since everything is balanced
cv00_predictions_analysis <- cv00_predictions_out %>% 
  group_by(trait, set, model, nTrainEnv, scheme) %>%
  summarize(fit = mean(accuracy)) %>% 
  ungroup()


## Edit some columns for plotting
cv00_predictions_analysis_toplot <- cv00_predictions_analysis %>%
  # unnest(model_nEnv_effects) %>%
  mutate(set = f_set_replace(set),
         model = factor(model, levels = rev(dist_method_abbr_use)),
         nTrainEnv = parse_number(as.character(nTrainEnv)),
         scheme = toupper(scheme),
         size = model == "Random")


## Subset the final training environment accuracy
cv00_rank_pred_all_data <- cv00_predictions_analysis_toplot %>%
  group_by(set, trait, scheme) %>% 
  filter(model == "GCD", nTrainEnv == max(nTrainEnv)) %>% # Filter the GCD model because it is shared
  ungroup() %>%
  distinct(set, trait, scheme, accuracy = fit)

## Add this final prediction accuracy to the original df
## Use this as a baseline
cv00_predictions_analysis_toplot1 <- left_join(cv00_predictions_analysis_toplot, cv00_rank_pred_all_data)



g_cv00_rank_pred <- cv00_predictions_analysis_toplot1 %>%
  filter(trait %in% traits) %>%
  ggplot(aes(x = nTrainEnv, y = fit, color = model, group = model)) +
  geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
  geom_line(aes(size = size)) +
  facet_grid(trait ~ set + scheme, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
  scale_color_manual(values = dist_colors_use, name = "Similarity\nmeasure",
                     guide = guide_legend(nrow = 2, keyheight = unit(0.5, "line"), reverse = TRUE)) +
  scale_linetype_manual(values = 2, name = NULL, guide = FALSE) +
  scale_size_manual(values = c(0.5, 1), guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
  scale_x_continuous(breaks = pretty, name = "Number of training set environments") +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom", strip.placement = "outside")

g_cv00_rank_pred1 <- g_cv00_rank_pred +
  scale_linetype_manual(values = 2, name = NULL, guide = guide_legend()) +
  scale_color_manual(values = dist_colors_use, guide = FALSE) +
  theme(legend.position = c(0.05, -0.1))

g_cv00_rank_pred2 <- plot_grid(g_cv00_rank_pred1, get_legend(g_cv00_rank_pred), ncol = 1, rel_heights = c(1, 0.07))


ggsave(filename = "cumulative_env_pred_cv00.jpg", plot = g_cv00_rank_pred2, path = fig_dir, width = 10, height = 5, dpi = 1000)


# ggsave(filename = "cumulative_env_pred_cv00_withAGDD.jpg", plot = g_cv00_rank_pred2, path = fig_dir, width = 7, height = 8, dpi = 1000)



## Combine CV and POV results
cv00_pov00_environment_rank_predictions_analysis_toplot <- bind_rows(cv00_predictions_analysis_toplot1, 
                                                                     pov00_predictions_analysis_toplot1) %>%
  filter(scheme != "POCV00") %>%
  # Create limits
  group_by(trait, scheme) %>%
  mutate_at(vars(fit), funs(min, max)) %>%
  ungroup()
  

g_rank_pred <- cv00_pov00_environment_rank_predictions_analysis_toplot %>%
  filter(trait %in% traits) %>%
  mutate(group = paste0(set, " - ", scheme)) %>%
  ggplot(aes(x = nTrainEnv, y = fit, color = model, group = model)) +
  geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
  geom_line(aes(size = size)) +
  facet_grid(trait + scheme ~ set, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
  scale_color_manual(values = dist_colors_use, name = "Similarity\nmeasure",
                     guide = guide_legend(nrow = 2, keyheight = unit(0.5, "line"), reverse = TRUE)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_size_manual(values = c(0.5, 1), guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
  scale_x_continuous(breaks = pretty, name = expression("Number of training set environments ("*italic(N[TE])*")")) +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom", strip.placement = "outside")

## Save
ggsave(filename = "cumulative_env_pred.jpg", plot = g_rank_pred, path = fig_dir, width = 7, height = 9, dpi = 1000)








# 
# 
# ### Example for presentation
# traits_present <- c("GrainYield", "HeadingDate")
# model_present <- c("Random", "IPCA-EC", "GCD", "PD")
# 
# ### Just leave-one-environment-out
# g_rank_pred_base <- cv00_pov00_environment_rank_predictions_analysis_toplot %>%
#   # filter(trait %in% traits) %>%
#   filter(trait %in% traits_present) %>%
#   filter(model %in% model_present) %>%
#   filter(set == "Leave-one-environment-out") %>%
#   mutate(scheme = ifelse(scheme == "CV00", "Cross-validation", "Parent-offspring validation")) %>%
#   ggplot(aes(x = nTrainEnv, y = fit, color = model, group = model)) +
#   geom_line(aes(size = size), color = "white") +
#   geom_hline(aes(yintercept = accuracy, lty = "Accuracy using\nall data")) +
#   facet_grid(trait ~ scheme, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
#   scale_color_manual(values = dist_colors_use, name = "Similarity\nmeasure",
#                      guide = guide_legend(keyheight = unit(1, "line"), reverse = TRUE)) +
#   scale_linetype_manual(values = 2, name = NULL) +
#   scale_size_manual(values = c(0.5, 1), guide = FALSE) +
#   scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
#   scale_x_continuous(breaks = pretty, name = expression("Number of training set environments ("*italic(N[TE])*")")) +
#   theme_presentation2(base_size = 12) +
#   theme(legend.position = "right", strip.placement = "outside")
# 
# ggsave(filename = "cumulative_env_pred_LOEO_presentation_base.jpg", path = fig_dir, plot = g_rank_pred_base, 
#        width = 6.5, height = 5, units = "in", dpi = 1000)
# 
# 
# g_rank_pred <- g_rank_pred_base + 
#   geom_line(aes(size = size))
# 
# ## Add viewport of heading data/cv00/leave-one-out
# g_rank_pred_vp_list <- cv00_pov00_environment_rank_predictions_analysis_toplot %>%
#   filter(trait == "HeadingDate") %>%
#   filter(model %in% model_present) %>%
#   split(list(.$set, .$scheme)) %>%
#   map2(.x = ., .y = letters[seq_along(.)], ~gg_vp(x = .x, let = ""))
# 
# 
# ## Coordinates for viewports
# vp1 <- viewport(x = 0.29, y = 0.27, width = 0.29, height = 0.33)
# vp2 <- viewport(x = 0.60, y = 0.34, width = 0.29, height = 0.33)
# 
# 
# jpeg(filename = file.path(fig_dir, "cumulative_env_pred_LOEO_presentation.jpg"), width = 6.5, height = 5, units = "in", res = 1000)
# print(g_rank_pred)
# print(g_rank_pred_vp_list$`Leave-one-environment-out.CV00`, vp = vp1)
# print(g_rank_pred_vp_list$`Leave-one-environment-out.POV00`, vp = vp2)
# dev.off()




# 
# 
# 
# ## Example graphs to demonstrate concept
# g_rank_pred_example <- cv00_pov00_environment_rank_predictions_analysis_toplot %>%
#   filter(trait == "GrainYield", set == "Leave-one-environment-out", scheme == "CV00") %>%
#   filter(model %in% c("LocPD", "Random", "GCD")) %>%
#   mutate(model = factor(model, levels = unique(model))) %>%
#   ggplot(aes(x = nTrainEnv, y = fit, color = model, group = model)) +
#   geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) +
#   geom_line(lwd = 2) +
#   scale_color_discrete(guide = FALSE) +
#   scale_linetype_manual(values = 2, name = NULL) +
#   scale_size_manual(values = c(0.5, 1), guide = FALSE) +
#   scale_y_continuous(breaks = pretty, limits = c(0.23, 0.47), name = expression("Prediction accuracy ("*italic(r[MG])*")"), labels = NULL) +
#   scale_x_continuous(breaks = pretty, name = expression("Number of training set environments ("*italic(N[TE])*")")) +
#   theme_presentation2(base_size = 18) +
#   theme(legend.position = c(0.75, 0.10), strip.placement = "outside")
# 
# 
# ggsave(filename = "rank_prediction_example.jpg", plot = g_rank_pred_example, width = 7, height = 5.5, dpi = 1000,
#        path = "C:/Users/jln54/GoogleDrive/BarleyLab/ForKevin/Presentations/DefenseSeminar/figures/")
# ggsave(filename = "rank_prediction_example1.jpg", plot = g_rank_pred_example, width = 7, height = 5.5, dpi = 1000,
#        path = "C:/Users/jln54/GoogleDrive/BarleyLab/ForKevin/Presentations/DefenseSeminar/figures/")
# ggsave(filename = "rank_prediction_example2.jpg", plot = g_rank_pred_example, width = 7, height = 5.5, dpi = 1000,
#        path = "C:/Users/jln54/GoogleDrive/BarleyLab/ForKevin/Presentations/DefenseSeminar/figures/")









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


## Make sure that each environment within a set, trait, and model has the same
## number of training environments
cv0_pov0_predictions_out %>% 
  group_by(trait, set, model, val_environment) %>% 
  summarize_at(vars(nTrainEnv), ~max(as.numeric(.))) %>% 
  summarize(n = n_distinct(nTrainEnv)) %>%
  filter(n > 1)

## good


## Just take averages, since everything is balanced
cv0_pov0_predictions_analysis <- cv0_pov0_predictions_out %>% 
  group_by(trait, set, model, nTrainEnv, scheme) %>%
  summarize(fit = mean(accuracy)) %>% 
  ungroup()


## Edit some columns for plotting
cv0_pov0_predictions_analysis_toplot <- cv0_pov0_predictions_analysis %>%
  # unnest(model_nEnv_effects) %>%
  mutate(set = f_set_replace(set),
         model = factor(model, levels = rev(dist_method_abbr_use)),
         nTrainEnv = parse_number(as.character(nTrainEnv)),
         scheme = toupper(scheme),
         size = model == "Random")


## Subset the final training environment accuracy
cv0_pov0_rank_pred_all_data <- cv0_pov0_predictions_analysis_toplot %>%
  group_by(set, trait, scheme) %>% 
  filter(model == "GCD", nTrainEnv == max(nTrainEnv)) %>%
  ungroup() %>%
  distinct(set, trait, scheme, accuracy = fit)

## Combine final accuracy to the original df
cv0_pov0_predictions_analysis_toplot1 <- left_join(cv0_pov0_predictions_analysis_toplot, cv0_pov0_rank_pred_all_data)




g_rank_pred_0 <- cv0_pov0_predictions_analysis_toplot1 %>%
  filter(trait %in% traits) %>%
  mutate(group = paste0(set, " - ", scheme)) %>%
  ggplot(aes(x = nTrainEnv, y = fit, color = model, group = model)) +
  geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
  geom_line(aes(size = size)) +
  facet_grid(trait + scheme ~ set, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
  scale_color_manual(values = dist_colors_use, name = "Similarity\nmeasure",
                     guide = guide_legend(nrow = 2, keyheight = unit(0.5, "line"), reverse = TRUE)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_size_manual(values = c(0.5, 1), guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
  scale_x_continuous(breaks = pretty, name = "Number of training set environments") +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom", strip.placement = "outside")

## Save
ggsave(filename = "cumulative_env_pred_cv0_pov0.jpg", plot = g_rank_pred_0, path = fig_dir, width = 7, height = 9, dpi = 1000)






### Plot everything!!

## Combine the to_plot dfs
## Reorganize by splitting on LOEO versus TF
cv_pov_rank_predictions1 <- bind_rows(cv0_pov0_predictions_analysis_toplot1, 
                                      cv00_pov00_environment_rank_predictions_analysis_toplot) %>%
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

# trait           set                       model nTrainEnv scheme   fit size  accuracy    min   max   decline per_decline
# 1 HeadingDateAGDD Leave-one-environment-out AMMI         25 CV00   0.801 FALSE    0.802  0.763 0.801 0           0        
# 2 HeadingDate     Leave-one-environment-out AMMI         25 CV00   0.799 FALSE    0.799  0.760 0.799 0.0000443   0.0000554
# 3 HeadingDateAGDD Leave-one-environment-out AMMI         25 CV0    0.885 FALSE    0.885 NA     0.886 0.00109     0.00124  
# 4 HeadingDate     Leave-one-environment-out AMMI         25 CV0    0.882 FALSE    0.882 NA     0.885 0.00252     0.00285  
# 5 HeadingDate     Leave-one-environment-out PD           25 CV00   0.798 FALSE    0.799  0.760 0.804 0.00569     0.00708  
# 6 HeadingDateAGDD Leave-one-environment-out PD           25 CV00   0.801 FALSE    0.802  0.763 0.807 0.00595     0.00737  
# 7 HeadingDateAGDD Leave-one-environment-out PD           25 CV0    0.885 FALSE    0.885 NA     0.899 0.0144      0.0160   
# 8 HeadingDate     Leave-one-environment-out PD           25 CV0    0.882 FALSE    0.882 NA     0.897 0.0147      0.0164   
# 9 PlantHeight     Leave-one-environment-out AMMI         26 CV0    0.749 FALSE    0.749 NA     0.810 0.0606      0.0748   
# 10 PlantHeight     Leave-one-environment-out AMMI         26 CV00   0.412 FALSE    0.410  0.239 0.495 0.0834      0.168    
# 11 PlantHeight     Leave-one-environment-out PD           26 CV0    0.749 FALSE    0.749 NA     0.838 0.0883      0.105    
# 12 GrainYield      Leave-one-environment-out AMMI         22 CV0    0.609 FALSE    0.609 NA     0.700 0.0908      0.130    
# 13 PlantHeight     Leave-one-environment-out PD           26 CV00   0.409 FALSE    0.410  0.239 0.508 0.0992      0.195    
# 14 GrainYield      Leave-one-environment-out AMMI         22 CV00   0.393 FALSE    0.393  0.195 0.514 0.121       0.236    
# 15 GrainYield      Leave-one-environment-out PD           22 CV0    0.609 FALSE    0.609 NA     0.743 0.134       0.181    
# 16 GrainYield      Leave-one-environment-out PD           22 CV00   0.392 FALSE    0.393  0.195 0.546 0.155       0.283







## Plot all
g_rank_pred_all <- cv_pov_rank_predictions1 %>%
  ggplot(data = ., aes(x = nTrainEnv, y = fit, color = model, group = model)) +
  geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
  geom_line(aes(size = size)) +
  facet_grid(trait + scheme ~ set, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
  scale_color_manual(values = dist_colors_use, name = "Similarity\nmeasure",
                     guide = guide_legend(nrow = 2, keyheight = unit(0.5, "line"), reverse = TRUE, order = 1)) +
  scale_linetype_manual(values = 2, name = NULL) +
  scale_size_manual(values = c(0.5, 1), guide = FALSE) +
  scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
  scale_x_continuous(breaks = pretty, name = expression("Number of training set environments ("*italic(N[TE])*")")) +
  theme_presentation2(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(filename = "cumulative_env_pred_all.jpg", plot = g_rank_pred_all, path = fig_dir, height = 15, width = 8, dpi = 1000)


## Plot by set
g_rank_pred_list <- cv_pov_rank_predictions1 %>%
  split(.$set) %>%
  map(~{
    # Create annotation df - use some heuristics to guess the xy coordinates
    ann_df <- .x %>% 
      filter(trait == trait[1]) %>% 
      mutate_at(vars(fit, accuracy), list(~min, ~max)) %>%
      filter(scheme == scheme[1]) %>% 
      group_by(trait, scheme) %>% 
      summarize_at(vars(nTrainEnv, fit_min, fit_max, accuracy_max), max) %>%
      mutate(x = 0.70 * nTrainEnv, y = 0.5 * fit_max, 
             seg_x = 20, seg_y = 1.05 * y, seg_xend = nTrainEnv + 2, seg_yend = 0.97 * accuracy_max)
      
    
    ggplot(data = .x, aes(x = nTrainEnv, y = fit, color = model, group = model)) +
      geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data")) + 
      geom_line(aes(size = size)) +
      # Add an annotation for accuracy using all data
      geom_text(data = ann_df, aes(x = x, y = y, label = "Accuracy using all data"), inherit.aes = FALSE, size = 2.5) +
      geom_segment(data = ann_df, aes(x = seg_x, y = seg_y, xend = seg_xend, yend = seg_yend), inherit.aes = FALSE) +
      facet_grid(trait ~ scheme, labeller = labeller(trait = str_add_space), space = "free_x", scales = "free", switch = "y") +
      scale_color_manual(values = dist_colors_use, name = "Similarity\nmeasure",
                         guide = guide_legend(ncol = 1, keyheight = unit(1, "line"), reverse = TRUE, order = 1)) +
      scale_linetype_manual(values = 2, name = NULL, guide = FALSE) +
      scale_size_manual(values = c(0.5, 1), guide = FALSE) +
      scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
      scale_x_continuous(breaks = pretty, name = expression("Number of training set environments ("*italic(N[TE])*")")) +
      theme_presentation2(base_size = 10) +
      theme(legend.position = "right", legend.margin = margin())
  })




### Viewports for LOEO
# Only heading date
g_vp_LOEO_list <- cv_pov_rank_predictions1 %>%
  filter(trait == "HeadingDate", set == "Leave-one-environment-out") %>%
  split(.$scheme) %>%
  map2(.x = ., .y = letters[seq_along(.)], ~gg_vp(.x, let = " "))
  

## Coordinates for viewports
vp1 <- viewport(x = 0.185, y = 0.48, width = 0.185, height = 0.20)
vp2 <- viewport(x = 0.45, y = 0.59, width = 0.185, height = 0.20)
vp3 <- viewport(x = 0.58, y = 0.48, width = 0.185, height = 0.20)
vp4 <- viewport(x = 0.775, y = 0.55, width = 0.185, height = 0.20)

vp_list <- list(vp1, vp2, vp3, vp4)


jpeg(filename = file.path(fig_dir, "cumulative_env_pred_all_LOEO_paper.jpg"), width = 7.5, height = 5.5, units = "in", res = 1000)
print(g_rank_pred_list$`Leave-one-environment-out`)
for(i in seq_along(g_vp_LOEO_list)[-2]) print(g_vp_LOEO_list[[i]], vp = vp_list[[i]])
dev.off()


# 
# ### Viewports for TF
# g_vp_TF_list <- cv_pov_rank_predictions1 %>%
#   filter(trait == "HeadingDate", set == "Time-forward") %>%
#   split(.$scheme) %>%
#   map2(.x = ., .y = letters[seq_along(.)], ~gg_vp(.x, let = " "))
# 
# 
# ## Coordinates for viewports
# vp1 <- viewport(x = 0.24, y = 0.52, width = 0.21, height = 0.17)
# vp2 <- viewport(x = 0.45, y = 0.59, width = 0.21, height = 0.17)
# vp3 <- viewport(x = 0.67, y = 0.52, width = 0.21, height = 0.17)
# vp4 <- viewport(x = 0.76, y = 0.56, width = 0.21, height = 0.17)
# 
# vp_list <- list(vp1, vp2, vp3, vp4)
# 
# 
# jpeg(filename = file.path(fig_dir, "cumulative_env_pred_all_TF_paper.jpg"), width = 5, height = 5.5, units = "in", res = 1000)
# print(g_rank_pred_list$`Time-forward`)
# for(i in seq_along(g_vp_TF_list)[-2]) print(g_vp_TF_list[[i]], vp = vp_list[[i]])
# dev.off()
























#### Environmental cluster predictions ####

cluster_pred <- ls(pattern = "[0-9]{1,2}_cluster_predictions") %>% 
  subset(., map_lgl(., ~inherits(get(.), "data.frame")))

cv00_cluster_predictions <- cv00_cluster_predictions %>%
  gather(scheme, accuracy, cv00, pocv00)
pov00_cluster_predictions <- gather(pov00_cluster_predictions, scheme, accuracy, pov00)

## Get and tidy
cluster_pred_df <- map(cluster_pred, get) %>%
  map(., ~mutate_at(., vars(which(names(.) %in% "rep")), as.character)) %>%
  map(., ~mutate_at(., vars(which(names(.) %in% "environment")), as.character))

cluster_pred_df[map_lgl(cluster_pred_df, ~"out" %in% names(.))] <- 
  cluster_pred_df[map_lgl(cluster_pred_df, ~"out" %in% names(.))] %>%
  map(~mutate(., accuracy = map_dbl(out, "base")) %>% select(., -out))

## Add CV/POV number
cluster_pred_df <- bind_rows(cluster_pred_df) %>% 
  mutate(scheme_number = as.character(str_extract(scheme, "[0-9]{1,2}")),
         model = ifelse(model == "pheno_location_dist", "pheno_loc_dist", model)) %>%
  ## For each trait, only use common val environments between LOYO and LOEO
  split(.$trait) %>% 
  map_df(~filter(.x, val_environment %in% unique(filter(.x, set != "complete")$val_environment))) %>%
  # Remove pocv00
  filter(scheme != "pocv00")



## Format some variables as factors for modeling  
cluster_predictions_base <- cluster_pred_df %>%
  left_join(., env_herit) %>% 
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability)) %>%
  mutate(model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr),
         zscore = accuracy,
         # zscore = ztrans(accuracy),
         scheme = toupper(scheme),
         # Group all LOYO environments together
         set = str_remove(set, "[0-9]{4}"),
         set = str_replace_all(set, set_replace),
         cluster = as.factor(cluster)) %>%
  filter(model %in% dist_method_abbr_use)
  

## Summarize average accuracy per trait, model, set, scheme, and val environmet
separate_cv_pov_cluster_predictions_analysis <- cluster_predictions_base %>%
  group_by(trait, set, model, scheme, val_environment) %>%
  summarize_at(vars(accuracy, ability), mean) %>%
  ungroup() %>%
  mutate(fit = accuracy) %>%
  left_join(., rename_at(separate_cv_pov_predictions_analysis1, vars(ability, accuracy), ~paste0("max_", .))) %>%
  # Add cluster and nTrainEnv back in
  left_join(., distinct(cluster_predictions_base, set, model, trait, val_environment, cluster, nTrainEnv)) %>%
  mutate(scheme = str_replace(scheme, "POCV", "POV"),
         scheme = factor(scheme, levels = c("CV0", "POV0", "CV00", "POV00")),
         model = factor(model, levels = dist_method_abbr_use)) %>%
  filter(scheme %in% names(scheme_color)) %>%
  mutate_at(vars(cluster, val_environment, scheme, model), as.factor)



## Analyze using a mixed model

## Model formula
## Models fitted individually for each set and trait
## Fixed effect of model and number of training environments
## zscore ~ model + nTrainEnv + (1|cluster:model) + (1|val_environment:cluster)
## 

## List of formulae
form_list <- formulas(~ accuracy,
                      fit1 = ~ 1 + model + scheme + model:scheme,
                      fit2 = add_predictors(fit1, ~ (1|cluster:model)),
                      fit3 = add_predictors(fit1, ~ (1 + nTrainEnv|cluster:model)),
                      fit4 = add_predictors(fit2, ~ (1|val_environment:cluster:model)),
                      fit5 = add_predictors(fit3, ~ (1|val_environment:cluster:model)),
                      fit6 = add_predictors(fit2, ~ (1 + nTrainEnv|val_environment:cluster:model)),
                      fit7 = add_predictors(fit3, ~ (1 + nTrainEnv|val_environment:cluster:model)))

## Model building
# First fit by trait and set
model_building1 <- separate_cv_pov_cluster_predictions_analysis %>%
  group_by(trait, set) %>%
  do({
    df <- .
    
    fit_list <- setNames(vector("list", length(form_list)), names(form_list))
    fit_list$fit1 <- lm(formula = form_list$fit1, data = df)
    for (f in names(form_list[-1])) fit_list[[f]] <- lmer(formula = form_list[[f]], data = df)
    
    ## Return the nested anova and model-specific anova
    tibble(
      nested_anova = list(
        tidy(eval(parse(text = paste0("anova(", paste0("fit_list$", rev(names(fit_list)), collapse = ", "), ")"))))),
      anova = list(map(fit_list, anova) %>% map(tidy) %>% imap_dfr(~mutate(.x, model = .y)))
    )
    
           
  }) %>% ungroup()

## Decide what model to use based on AIC
model_building1 %>% 
  unnest(nested_anova) %>%
  group_by(trait, set) %>% 
  top_n(x = ., n = 1, wt = -AIC)

## Model 4
form_use_name <- "fit4"
form_use <- form_list[[form_use_name]]

## Check fit-specific anovas
model_building1 %>%
  unnest(anova) %>%
  filter(model == form_use_name,
         p.value < alpha)

## No interaction effect, so model together


cv_pov_cluster_predictions_model_analysis <- separate_cv_pov_cluster_predictions_analysis %>%
  group_by(set, trait) %>%
  nest() %>%
  mutate(model_rep = map_lgl(data, ~!all(is.na(.$rep)))) %>%
  mutate(out = list(NULL))

## Loop for each row
for (i in seq(nrow(cv_pov_cluster_predictions_model_analysis))) {

  df <- cv_pov_cluster_predictions_model_analysis$data[[i]]

  ## fit the final model
  fit <- lmer(form_use, data = df)

  effects_keep <- as.data.frame(Effect(c("model", "scheme"), fit))

  fit_summary <- data_frame(
    fitted = list(fit),
    scheme_effects = list(effects_keep),
    anova = list(tidy(anova(fit))),
    ranova = list(tidy(ranova(fit)))
  )

  cv_pov_cluster_predictions_model_analysis$out[[i]] <- fit_summary

}


## Unnest effects
cv_pov_cluster_predictions_model_analysis1 <- cv_pov_cluster_predictions_model_analysis %>%
  unnest(out) %>%
  unnest(scheme_effects) %>%
  # Add all-data accuracies back in
  left_join(., mutate(separate_cv_pov_predictions_analysis1, scheme = ifelse(scheme == "POCV0", "POV0", as.character(scheme))) %>% 
              rename_at(vars(ability, accuracy), ~paste0("max_", .))) %>%
  mutate(model = factor(model, levels = dist_method_abbr_use),
         scheme = factor(scheme, levels = c("CV0", "POV0", "CV00", "POV00")))
         


## Summarize the mean accuracy and 1 sd
separate_cv_pov_cluster_predictions_analysis1 <- separate_cv_pov_cluster_predictions_analysis %>%
  group_by(trait, set, model, scheme) %>%
  summarize_at(vars(fit, max_ability, max_accuracy), list(~mean, ~min, ~max)) %>%
  ungroup() %>%
  select(-matches("accuracy_max|ability_max|accuracy_min|ability_min")) %>%
  rename_at(vars(contains("_mean")), ~str_remove(., "_mean"))
         
         

## Determine relative widths based on number of similarity measures
rel_widths <- separate_cv_pov_cluster_predictions_analysis %>%
  group_by(set) %>% 
  summarize(n = n_distinct(model)) %>% 
  mutate(n / max(n)) %>%
  pull()

## Plot boxplot
g_cv_pov_cluster_box_list <- separate_cv_pov_cluster_predictions_analysis %>%
  filter(trait %in% traits) %>%
  ## Set min/max based on trait
  group_by(trait) %>%
  mutate_at(vars(fit), list(~min, ~max)) %>%
  split(.$set) %>%
  map(~{
    df1 <- group_by(.x, trait, scheme, model) %>%
      summarize(mean = mean(fit))
    
    ggplot(data = .x, aes(x = model, y = fit)) +
      ## Add blank points for force the y axis scale
      geom_point(data = .x, aes(x = model, y = min), color = "white") + 
      geom_point(data = .x, aes(x = model, y = max), color = "white") + 
      ##
      geom_hline(aes(yintercept = max_accuracy, lty = "Accuracy using all data"), color = "black") +
      # geom_boxplot(aes(fill = model), width = 0.5, alpha = 0.4, position = position_dodge(0.9), outlier.alpha = 0.1) +
      geom_violin(aes(fill = model), width = 0.5, alpha = 0.4, position = position_dodge(0.9),
                  draw_quantiles = 0.5) +
      # Add points for mean
      geom_point(data = df1, aes(y = mean), size = 2) +
      scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
      scale_fill_manual(values = dist_colors_use, name = "Similarity measure", guide = guide_legend(nrow = 1)) +
      scale_linetype_manual(values = 2, name = NULL) +
      xlab("Similarity measure") +
      facet_grid(trait ~ scheme, scales = "free", space = "free_x", switch = "y",
                 labeller = labeller(set = str_to_title, trait = str_add_space)) +
      labs(subtitle = unique(.x$set)) +
      theme_presentation2(base_size = 10) +
      theme(axis.text.x = element_blank(), strip.placement = "outside", legend.position = "bottom",
            legend.margin = margin(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  })

## Combine plots
g_cv_pov_cluster_box <- plot_grid(plotlist = map(g_cv_pov_cluster_box_list, ~. + theme(legend.position = "none")), 
                                  labels = letters[seq_along(g_cv_pov_cluster_box_list)], nrow = 1, rel_widths = rel_widths)
g_cv_pov_cluster_box1 <- plot_grid(g_cv_pov_cluster_box, get_legend(g_cv_pov_cluster_box_list[[1]]), ncol = 1, 
                                   rel_heights = c(1, 0.05))

ggsave(filename = "all_cv_pov_cluster_predictions_boxplot.jpg", plot = g_cv_pov_cluster_box1, path = fig_dir, 
       width = 10, height = 5, dpi = 1000)



## Plot points
g_cv_pov_cluster_list <- separate_cv_pov_cluster_predictions_analysis1 %>%
  filter(trait %in% traits) %>%
  ## Set min/max based on trait
  group_by(trait) %>%
  mutate_at(vars(fit), list(~min, ~max)) %>%
  split(.$set) %>%
  map(~{
    ggplot(data = .x, aes(x = model, y = fit)) +
      ## Add blank points for force the y axis scale
      geom_point(data = .x, aes(x = model, y = min), color = "white") + 
      geom_point(data = .x, aes(x = model, y = max), color = "white") + 
      ##
      geom_hline(aes(yintercept = max_accuracy, lty = "Accuracy using all data"), color = "black") +
      # Range bar
      # geom_linerange(aes(ymin = fit_min, ymax = fit_max)) +
      geom_point(aes(color = model, fill = model), alpha = 0.9, position = position_dodge(0.9), size = 2) +
      scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
      scale_fill_manual(values = dist_colors_use, name = "Similarity measure", guide = guide_legend(nrow = 1)) +
      scale_color_manual(values = dist_colors_use, name = "Similarity measure", guide = guide_legend(nrow = 1)) +
      scale_linetype_manual(values = 2, name = NULL) +
      xlab("Similarity measure") +
      facet_grid(trait ~ scheme, scales = "free", space = "free_x", switch = "y",
                 labeller = labeller(set = str_to_title, trait = str_add_space)) +
      labs(subtitle = unique(.x$set)) +
      theme_presentation2(base_size = 10) +
      theme(axis.text.x = element_blank(), strip.placement = "outside", legend.position = "bottom",
            legend.margin = margin(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  })

## Combine plots
g_cv_pov_cluster <- plot_grid(plotlist = map(g_cv_pov_cluster_list, ~. + theme(legend.position = "none")), 
                              labels = letters[seq_along(g_cv_pov_cluster_list)], nrow = 1, rel_widths = rel_widths)
g_cv_pov_cluster1 <- plot_grid(g_cv_pov_cluster, get_legend(g_cv_pov_cluster_list[[1]]), ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "all_cv_pov_cluster_predictions_point.jpg", plot = g_cv_pov_cluster1, path = fig_dir, 
       width = 7.5, height = 4, dpi = 1000)



## Plot points
g_cv_pov_cluster_model_list <- cv_pov_cluster_predictions_model_analysis1 %>%
  filter(trait %in% traits) %>%
  ## Set min/max based on trait
  group_by(trait) %>%
  mutate_at(vars(fit), list(~min, ~max)) %>%
  split(.$set) %>%
  map(~{
    ggplot(data = .x, aes(x = model, y = fit)) +
      ## Add blank points for force the y axis scale
      geom_point(data = .x, aes(x = model, y = min), color = "white") + 
      geom_point(data = .x, aes(x = model, y = max), color = "white") + 
      ##
      geom_hline(aes(yintercept = max_accuracy, lty = "Accuracy using all data"), color = "black") +
      # Range bar
      # geom_linerange(aes(ymin = fit_min, ymax = fit_max)) +
      geom_point(aes(color = model, fill = model), alpha = 0.9, position = position_dodge(0.9), size = 2) +
      scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
      scale_fill_manual(values = dist_colors_use, name = "Similarity measure", guide = guide_legend(nrow = 1)) +
      scale_color_manual(values = dist_colors_use, name = "Similarity measure", guide = guide_legend(nrow = 1)) +
      scale_linetype_manual(values = 2, name = NULL) +
      xlab("Similarity measure") +
      facet_grid(trait ~ scheme, scales = "free", space = "free_x", switch = "y",
                 labeller = labeller(set = str_to_title, trait = str_add_space)) +
      labs(subtitle = unique(.x$set)) +
      theme_presentation2(base_size = 10) +
      theme(axis.text.x = element_blank(), strip.placement = "outside", legend.position = "bottom",
            legend.margin = margin(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  })

## Combine plots
g_cv_pov_cluster_model <- plot_grid(plotlist = map(g_cv_pov_cluster_model_list, ~. + theme(legend.position = "none")), 
                              labels = letters[seq_along(g_cv_pov_cluster_model_list)], nrow = 1, rel_widths = rel_widths)
g_cv_pov_cluster_model1 <- plot_grid(g_cv_pov_cluster_model, get_legend(g_cv_pov_cluster_model_list[[1]]), 
                                     ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "all_cv_pov_cluster_predictions_model_point.jpg", plot = g_cv_pov_cluster_model1, path = fig_dir, 
       width = 7.5, height = 4, dpi = 1000)











### Select CV00 and PV00
# 
# # Create a legend key
# legend_key <- paste0(names(dist_colors_use), ": ", c("", "Phenotypic distance", "Location PD", "Geographic distance", "All covariates", "Mean-correlated ECs", "IPCA-correlated ECs"))
# 
# ## Example for presentation
# 
# g_cv00_pov00_cluster <- separate_cv_pov_cluster_predictions_analysis1 %>%
#   filter(scheme %in% c("CV00", "POV00"), set == "Leave-one-environment-out") %>%
#   mutate(scheme = ifelse(scheme == "CV00", "Cross-validation", "Parent-offspring validation")) %>%
#   filter(trait %in% traits_present, model %in% model_present) %>%
#   ggplot(aes(x = model, y = fit, ymin = lower, ymax = upper)) +
#   geom_hline(aes(yintercept = accuracy, lty = "Accuracy using all data"), color = "grey") +
#   geom_errorbar(width = 0.5) +
#   geom_point(aes(color = model), size = 3) +
#   scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
#   scale_color_manual(values = dist_colors_use, name = "Cluster method", guide = FALSE) +
#   scale_linetype_manual(values = 1, name = NULL) +
#   xlab("Similarity measure") +
#   labs(caption = "") +
#   facet_grid(trait ~ scheme, scales = "free_x", space = "free_x", switch = "y", labeller = labeller(set = str_to_title, trait = str_add_space)) +
#   theme_presentation2(base_size = 12) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.70, -0.29), legend.margin = margin())
# 
# ggsave(filename = "cv00_pov00_cluster_predictions_presentation.jpg", plot = g_cv00_pov00_cluster, path = fig_dir, 
#        width = 4.5, height = 4.5, dpi = 1000)







##### 
##### Random cluster prediction #####
##### 

# This will be analyzed by calculating the difference between the cluster prediction accuracy 
# and each replicate of the reandom cluster prediction accuracy, then modeling that difference

random_cluster_pred <- ls(pattern = "[0-9]{1,2}_cluster_random_predictions") %>% 
  subset(., map_lgl(., ~inherits(get(.), "data.frame")))

## Get and tidy
random_cluster_pred_df <- map(random_cluster_pred, get) %>%
  # Tidy if column names match (cv0 or pov0)
  map(., ~gather(., scheme, accuracy, matches("cv|po"))) %>%
  map(., ~mutate_at(., vars(which(names(.) %in% "rep")), as.character)) %>%
  map(., ~mutate_at(., vars(which(names(.) %in% "environment")), as.character)) %>%
  bind_rows() %>% 
  mutate(scheme_number = as.character(str_extract(scheme, "[0-9]{1,2}")),
         model = ifelse(model == "pheno_location_dist", "pheno_loc_dist", model))

## Reformat some variables
random_cluster_predictions_base <- random_cluster_pred_df %>%
  left_join(., env_herit) %>% 
  mutate(ability = accuracy, accuracy = ability / sqrt(heritability)) %>%
  mutate(model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr),
         zscore = accuracy,
         # zscore = ztrans(accuracy),
         scheme = toupper(scheme),
         # Group all LOYO environments together
         set = str_remove(set, "[0-9]{4}"),
         set = str_replace_all(set, set_replace)) %>%
  mutate_at(vars(cluster, val_environment, scheme), as.factor) %>%
  filter(model %in% dist_method_abbr_use)


## Combine with the original clustering predictions, then calculate the advantage in accuracy
## First calculate averages over reps for the original predictions
random_and_cluster_predictions_base <- cluster_predictions_base %>% 
  group_by(set, trait, model, val_environment, nTrainEnv, cluster, scheme, scheme_number) %>% 
  summarize(base_accuracy = mean(accuracy)) %>%
  ungroup() %>%
  # Now add the random prediction accuracies back in and calculate the difference between 
  # informed accuracy and the random accuracy
  left_join(., random_cluster_predictions_base,  
            by = c("set", "trait", "model", "val_environment", "cluster", "scheme", "scheme_number")) %>%
  mutate(accuracy_advantage = base_accuracy - accuracy,
         scheme = str_replace(scheme, "POCV", "POV"),
         scheme = factor(scheme, levels = c("CV0", "POV0", "CV00", "POV00"))) %>%
  rename(nTrainEnv = nTrainEnv.x)




## Summarize average accuracy per trait, model, set, scheme, and val environmet
random_cluster_predictions_analysis <- random_and_cluster_predictions_base %>%
  # Calculate the accuracy differential
  group_by(trait, set, model, scheme, val_environment) %>%
  summarize_at(vars(contains("accuracy")), mean) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = dist_method_abbr_use)) %>%
  filter(scheme %in% names(scheme_color)) %>%
  rename(fit = accuracy_advantage)


## Analyze using a mixed model

## Model formula
## Models fitted individually for each set and trait
## Fixed effect of model and number of training environments
## zscore ~ model + nTrainEnv + (1|cluster:model) + (1|val_environment:cluster)
## 

## List of formulae
form_list <- formulas(~ accuracy_advantage,
                      fit1 = ~ 1 + model + scheme + model:scheme,
                      fit2 = add_predictors(fit1, ~ (1|cluster:model)),
                      fit3 = add_predictors(fit1, ~ (1 + nTrainEnv|cluster:model)),
                      fit4 = add_predictors(fit2, ~ (1|val_environment:cluster:model)),
                      fit5 = add_predictors(fit3, ~ (1|val_environment:cluster:model)),
                      fit6 = add_predictors(fit2, ~ (1 + nTrainEnv|val_environment:cluster:model)),
                      fit7 = add_predictors(fit3, ~ (1 + nTrainEnv|val_environment:cluster:model)))

## Model building
# First fit by trait and set
model_building1 <- random_and_cluster_predictions_base %>%
  group_by(trait, set) %>%
  do({
    df <- .
    
    fit_list <- setNames(vector("list", length(form_list)), names(form_list))
    fit_list$fit1 <- lm(formula = form_list$fit1, data = df)
    for (f in names(form_list[-1])) fit_list[[f]] <- lmer(formula = form_list[[f]], data = df)
    
    ## Return the nested anova and model-specific anova
    tibble(
      nested_anova = list(
        tidy(eval(parse(text = paste0("anova(", paste0("fit_list$", rev(names(fit_list)), collapse = ", "), ")"))))),
      anova = list(map(fit_list, anova) %>% map(tidy) %>% imap_dfr(~mutate(.x, model = .y)))
    )
    
    
  }) %>% ungroup()

## Decide what model to use based on AIC
model_building1 %>% 
  unnest(nested_anova) %>%
  group_by(trait, set) %>% 
  top_n(x = ., n = 1, wt = -AIC)

## Model 6
form_use_name <- "fit6"
form_use <- form_list[[form_use_name]]

## Check fit-specific anovas
model_building1 %>%
  unnest(anova) %>%
  filter(model == form_use_name,
         p.value < alpha)

## Significant interaction effect, so model scheme separately
## Remove scheme and scheme:model from the formula
form_use <- as.formula(paste0("accuracy_advantage ~ ", str_remove_all(as.character(form_use)[3], "\\+ scheme|\\+ model:scheme")))

## Analyze separately
random_cluster_predictions_analysis_model <- random_and_cluster_predictions_base %>%
  group_by(set, trait, scheme) %>%
  nest() %>%
  mutate(model_rep = map_lgl(data, ~!all(is.na(.$rep)))) %>%
  mutate(out = list(NULL))

## Loop for each row
for (i in seq(nrow(random_cluster_predictions_analysis_model))) {
  
  df <- random_cluster_predictions_analysis_model$data[[i]]
  
  ## fit the final model
  fit <- lmer(form_use, data = df)

  
  effects_keep <- as.data.frame(Effect(c("model"), fit))
  
  fit_summary <- data_frame(
    fitted = list(fit),
    scheme_effects = list(effects_keep),
    anova = list(tidy(anova(fit))),
    ranova = list(tidy(ranova(fit)))
  )
  
  random_cluster_predictions_analysis_model$out[[i]] <- fit_summary
  
}


## Unnest effects
random_cluster_predictions_analysis_model1 <- random_cluster_predictions_analysis_model %>%
  unnest(out) %>%
  unnest(scheme_effects) %>%
  mutate(model = factor(model, levels = dist_method_abbr_use),
         scheme = factor(scheme, levels = c("CV0", "POV0", "CV00", "POV00")))



## Summarize the mean accuracy and 1 sd
random_cluster_predictions_analysis1 <- random_cluster_predictions_analysis %>%
  group_by(trait, set, model, scheme) %>%
  summarize_at(vars(fit), list(~mean, ~min, ~max)) %>%
  ungroup() %>%
  rename(fit = mean)



## Determine relative widths based on number of similarity measures
rel_widths <- random_cluster_predictions_analysis %>%
  group_by(set) %>% 
  summarize(n = n_distinct(model)) %>% 
  mutate(n / max(n)) %>%
  pull()

## Plot boxplot
g_cv_pov_cluster_diff_box_list <- random_cluster_predictions_analysis %>%
  filter(trait %in% traits) %>%
  ## Set min/max based on trait
  group_by(trait) %>%
  mutate_at(vars(fit), list(~min, ~max)) %>%
  split(.$set) %>%
  map(~{
    df1 <- group_by(.x, trait, scheme, model) %>%
      summarize(mean = mean(fit))
    
    ggplot(data = .x, aes(x = model, y = fit)) +
      ## Add blank points for force the y axis scale
      geom_point(data = .x, aes(x = model, y = min), color = "white") + 
      geom_point(data = .x, aes(x = model, y = max), color = "white") + 
      ##
      geom_hline(yintercept = 0, color = "grey85") +
      # geom_boxplot(aes(fill = model), width = 0.5, alpha = 0.4, position = position_dodge(0.9), outlier.alpha = 0.1) +
      geom_violin(aes(fill = model), width = 0.5, alpha = 0.4, position = position_dodge(0.9),
                  draw_quantiles = 0.5) +
      # Add points for mean
      geom_point(data = df1, aes(y = mean), size = 2) +
      scale_y_continuous(breaks = pretty, name = expression("Relative prediction accuracy ("*italic(r[MG])*")")) +
      scale_fill_manual(values = dist_colors_use, name = "Similarity measure", guide = guide_legend(nrow = 1)) +
      scale_linetype_manual(values = 2, name = NULL) +
      xlab("Similarity measure") +
      facet_grid(trait ~ scheme, scales = "free", space = "free_x", switch = "y",
                 labeller = labeller(set = str_to_title, trait = str_add_space)) +
      labs(subtitle = unique(.x$set)) +
      theme_presentation2(base_size = 10) +
      theme(axis.text.x = element_blank(), strip.placement = "outside", legend.position = "bottom",
            legend.margin = margin(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  })

## Combine plots
g_cv_pov_cluster_diff_box <- plot_grid(plotlist = map(g_cv_pov_cluster_diff_box_list, ~. + theme(legend.position = "none")), 
                                  labels = letters[seq_along(g_cv_pov_cluster_diff_box_list)], nrow = 1, rel_widths = rel_widths)
g_cv_pov_cluster_diff_box1 <- plot_grid(g_cv_pov_cluster_diff_box, get_legend(g_cv_pov_cluster_diff_box_list[[1]]), ncol = 1, 
                                   rel_heights = c(1, 0.05))

ggsave(filename = "all_cv_pov_cluster_predictions_advantage_boxplot.jpg", plot = g_cv_pov_cluster_diff_box1, path = fig_dir, 
       width = 10, height = 5, dpi = 1000)



## Plot points
g_cv_pov_cluster_diff_list <- random_cluster_predictions_analysis1 %>%
  filter(trait %in% traits) %>%
  ## Set min/max based on trait
  group_by(trait) %>%
  mutate_at(vars(fit), list(~min, ~max)) %>%
  split(.$set) %>%
  map(~{
    ggplot(data = .x, aes(x = model, y = fit)) +
      ## Add blank points for force the y axis scale
      geom_point(data = .x, aes(x = model, y = min), color = "white") + 
      geom_point(data = .x, aes(x = model, y = max), color = "white") + 
      ##
      geom_hline(yintercept = 0, color = "grey85") +
      # Range bar
      # geom_linerange(aes(ymin = fit_min, ymax = fit_max)) +
      geom_point(aes(color = model, fill = model), alpha = 0.9, position = position_dodge(0.9), size = 2) +
      scale_y_continuous(breaks = pretty, name = expression("Relative prediction accuracy ("*italic(r[MG])*")")) +
      scale_fill_manual(values = dist_colors_use, name = "Similarity measure", guide = guide_legend(nrow = 1)) +
      scale_color_manual(values = dist_colors_use, name = "Similarity measure", guide = guide_legend(nrow = 1)) +
      scale_linetype_manual(values = 2, name = NULL) +
      xlab("Similarity measure") +
      facet_grid(trait ~ scheme, scales = "free", space = "free_x", switch = "y",
                 labeller = labeller(set = str_to_title, trait = str_add_space)) +
      labs(subtitle = unique(.x$set)) +
      theme_presentation2(base_size = 10) +
      theme(axis.text.x = element_blank(), strip.placement = "outside", legend.position = "bottom",
            legend.margin = margin(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  })

## Combine plots
g_cv_pov_cluster_diff <- plot_grid(plotlist = map(g_cv_pov_cluster_diff_list, ~. + theme(legend.position = "none")), 
                              labels = letters[seq_along(g_cv_pov_cluster_diff_list)], nrow = 1, rel_widths = rel_widths)
g_cv_pov_cluster_diff1 <- plot_grid(g_cv_pov_cluster_diff, get_legend(g_cv_pov_cluster_diff_list[[1]]), 
                                    ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "all_cv_pov_cluster_predictions_advantage_point.jpg", plot = g_cv_pov_cluster_diff1, path = fig_dir, 
       width = 7.5, height = 4, dpi = 1000)



## Plot model marginal means
g_cv_pov_cluster_diff_model_list <- random_cluster_predictions_analysis_model1 %>%
  filter(trait %in% traits) %>%
  ## Set min/max based on trait
  group_by(trait) %>%
  mutate_at(vars(fit), list(~min, ~max)) %>%
  split(.$set) %>%
  map(~{
    ggplot(data = .x, aes(x = model, y = fit)) +
      ## Add blank points for force the y axis scale
      geom_point(data = .x, aes(x = model, y = min), color = "white") + 
      geom_point(data = .x, aes(x = model, y = max), color = "white") + 
      ##
      geom_hline(yintercept = 0, color = "grey85") +
      # Range bar
      # geom_linerange(aes(ymin = fit_min, ymax = fit_max)) +
      geom_point(aes(color = model, fill = model), alpha = 0.9, position = position_dodge(0.9), size = 2) +
      scale_y_continuous(breaks = pretty, name = expression("Relative prediction accuracy ("*italic(r[MG])*")")) +
      scale_fill_manual(values = dist_colors_use, name = "Similarity measure", guide = guide_legend(nrow = 1)) +
      scale_color_manual(values = dist_colors_use, name = "Similarity measure", guide = guide_legend(nrow = 1)) +
      scale_linetype_manual(values = 2, name = NULL) +
      xlab("Similarity measure") +
      facet_grid(trait ~ scheme, scales = "free", space = "free_x", switch = "y",
                 labeller = labeller(set = str_to_title, trait = str_add_space)) +
      labs(subtitle = unique(.x$set)) +
      theme_presentation2(base_size = 10) +
      theme(axis.text.x = element_blank(), strip.placement = "outside", legend.position = "bottom",
            legend.margin = margin(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
  })

## Combine plots
g_cv_pov_cluster_diff_model <- plot_grid(plotlist = map(g_cv_pov_cluster_diff_model_list, ~. + theme(legend.position = "none")), 
                                    labels = letters[seq_along(g_cv_pov_cluster_diff_model_list)], 
                                    nrow = 1, rel_widths = rel_widths)
g_cv_pov_cluster_diff_model1 <- plot_grid(g_cv_pov_cluster_diff_model, get_legend(g_cv_pov_cluster_diff_model_list[[1]]), 
                                     ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "all_cv_pov_cluster_predictions_advantage_model_point.jpg", plot = g_cv_pov_cluster_diff_model1, path = fig_dir, 
       width = 7.5, height = 4, dpi = 1000)








##############################
## Combine random and nominal cluster predictions
##############################


## Subset some silimarity measures for publication
similarity_measure_subset <- c("AMMI", "LocPD", "GCD", "IPCA-EC")

# Create a new df for plotting
cluster_prediction_analysis_model <- bind_rows(
  mutate(cv_pov_cluster_predictions_model_analysis1, analysis = "base"),
  mutate(random_cluster_predictions_analysis_model1, analysis = "random")
) %>% filter(model %in% similarity_measure_subset)


## Plot modifier
gg_add <- list(
  geom_point(aes(color = model, fill = model), alpha = 0.9, position = position_dodge(0.9), size = 2) ,
  scale_fill_manual(values = dist_colors_use, name = "Similarity measure", guide = guide_legend(nrow = 1)) ,
  scale_color_manual(values = dist_colors_use, name = "Similarity measure", guide = guide_legend(nrow = 1)) ,
  scale_linetype_manual(values = 2, name = NULL) ,
  xlab("Similarity measure") ,
  facet_grid(trait ~ scheme, scales = "free", space = "free_x", switch = "y",
             labeller = labeller(set = str_to_title, trait = str_add_space)) ,
  theme_presentation2(base_size = 10) ,
  theme(axis.text.x = element_blank(), strip.placement = "outside", legend.position = "bottom",
        legend.margin = margin(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
)


## Plot model marginal means of nominal prediction accuracy and
## random cluster prediction accuracy
g_cv_pov_cluster_model_combined <- cluster_prediction_analysis_model %>%
  filter(trait %in% traits) %>%
  ## Set min/max based on trait
  group_by(trait) %>%
  mutate_at(vars(fit), list(~min, ~max)) %>%
  split(list(.$set, .$analysis)) %>%
  map(~{
    ## Decide y axis label
    if (unique(.x$analysis) == "base") {
      y_axis_lab <- expression("Prediction accuracy ("*italic(r[MG])*")")
      
      ggplot(data = .x, aes(x = model, y = fit)) +
        ## Add blank points for force the y axis scale
        geom_point(data = .x, aes(x = model, y = min), color = "white") +
        geom_point(data = .x, aes(x = model, y = max), color = "white") + 
        ##
        geom_hline(aes(yintercept = max_accuracy), color = "black") +
        scale_y_continuous(breaks = pretty, name = y_axis_lab) +
        labs(subtitle = unique(.x$set)) +
        gg_add
      
    } else {
      y_axis_lab <- expression("Relative prediction accuracy ("*italic(r[MG])*")")
      
      ggplot(data = .x, aes(x = model, y = fit)) +
        ## Add blank points for force the y axis scale
        geom_point(data = .x, aes(x = model, y = min), color = "white") +
        geom_point(data = .x, aes(x = model, y = max), color = "white") + 
        ##
        geom_hline(yintercept = 0, color = "grey85") +
        scale_y_continuous(breaks = pretty, name = y_axis_lab) +
        labs(subtitle = unique(.x$set)) +
        gg_add
    }

  })

## Combine plots
g_cv_pov_cluster_combine_model <- plot_grid(plotlist = map(g_cv_pov_cluster_model_combined, ~. + theme(legend.position = "none")), 
                                         labels = letters[seq_along(g_cv_pov_cluster_model_combined)], 
                                         nrow = 2, rel_widths = rel_widths)
g_cv_pov_cluster_combine_model1 <- plot_grid(g_cv_pov_cluster_combine_model, get_legend(g_cv_pov_cluster_model_combined[[1]]), 
                                          ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "all_cv_pov_cluster_predictions_combined_model_point.jpg", plot = g_cv_pov_cluster_combine_model1,
       path = fig_dir, width = 6, height = 6, dpi = 1000)







## Plot model marginal means of nominal prediction accuracy and
## random cluster prediction accuracy
##  
## This time use dumbell graph
## 
g_cv_pov_cluster_model_combined_alt <- cluster_prediction_analysis_model %>%
  ## Subtract random from based
  select(set, trait, model, scheme, analysis, fit) %>% 
  spread(analysis, fit) %>% 
  mutate(random = base - random) %>% 
  gather(analysis, fit, base, random) %>% 
  left_join(., select(cluster_prediction_analysis_model, -fit)) %>%
  ##
  filter(trait %in% traits) %>%
  ## Set min/max based on trait
  group_by(trait) %>%
  mutate_at(vars(fit), list(~min, ~max)) %>%
  split(.$set) %>%
  map(~{
    ## Spread fit by analysis
    df1 <- select(.x, set, trait, model, scheme, analysis, fit) %>%
      spread(analysis, fit) %>%
      left_join(., select(.x, set, trait, model, scheme, max_accuracy))
    
    ggplot(data = df1, aes(x = model)) +
      ## Add blank points for force the y axis scale
      geom_point(data = .x, aes(y = min), color = "white") + 
      geom_point(data = .x, aes(y = max), color = "white") + 
      ##
      geom_hline(aes(yintercept = max_accuracy, lty = "Accuracy using\nall data"), color = "black") +
      # Segment between points
      geom_segment(aes(y = random, yend = base, x = model, xend = model), color = "grey85") +
      ## Points for base and random
      geom_point(aes(y = random, shape = "Random"), color = "grey85", size = 1) +
      geom_point(aes(y = base, shape = "Similarity", color = model), alpha = 0.5, size = 1.75) +
      scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
      scale_fill_manual(values = dist_colors_use, name = "Similarity\nmeasure", guide = guide_legend(nrow = 2)) +
      scale_color_manual(values = dist_colors_use, name = "Similarity\nmeasure", guide = guide_legend(nrow = 2)) +
      scale_shape_manual(name = "Cluster\ndesign", guide = guide_legend(nrow = 2), values = c("Similarity" = 16, "Random" = 15)) +
      scale_linetype_manual(values = 2, name = NULL, guide = guide_legend(order = 1)) +
      xlab("Similarity measure") +
      facet_grid(trait ~ scheme, scales = "free", space = "free_x", switch = "y",
                 labeller = labeller(set = str_to_title, trait = str_add_space)) +
      labs(subtitle = unique(.x$set)) +
      theme_presentation2(base_size = 8) +
      theme(axis.text.x = element_blank(), strip.placement = "outside", legend.position = "bottom",
            legend.margin = margin(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
            legend.key.height = unit(0.5, "lines"))
    
  })
    
 

## Combine plots
# Custom y_axis
g_cv_pov_cluster_combine_model_alt <- plot_grid(
  plotlist = map(g_cv_pov_cluster_model_combined_alt, ~. + theme(legend.position = "none", axis.title.y = element_blank())), 
  labels = letters[seq_along(g_cv_pov_cluster_model_combined_alt)], nrow = 1, rel_widths = rel_widths
)
# Add y axis
g_cv_pov_cluster_combine_model_alt1 <- plot_grid(
  grid::textGrob(label = expression("Prediction accuracy ("*italic(r[MG])*")"), rot = 90, gp = gpar(fontsize = 8)),
  g_cv_pov_cluster_combine_model_alt, rel_widths = c(0.05, 1), nrow = 1
)
g_cv_pov_cluster_combine_model_alt2 <- plot_grid(
  g_cv_pov_cluster_combine_model_alt1, get_legend(g_cv_pov_cluster_model_combined_alt[[1]]), 
  ncol = 1, rel_heights = c(1, 0.1)
)

ggsave(filename = "all_cv_pov_cluster_predictions_combined_model_paper.jpg", plot = g_cv_pov_cluster_combine_model_alt2,
       path = fig_dir, width = 6, height = 3.5, dpi = 1000)









## Save plotting data
save_list <- c("separate_cv_pov_predictions_analysis", "separate_cv_pov_predictions_analysis1", 
               "cv_pov_rank_predictions1", "separate_cv_pov_cluster_predictions_analysis", 
               "cv_pov_cluster_predictions_model_analysis1", "separate_cv_pov_cluster_predictions_analysis1",
               "random_cluster_predictions_analysis", "random_cluster_predictions_analysis_model1", 
               "random_cluster_predictions_analysis1")

save(list = save_list, file = file.path(result_dir, "prediction_plotting_data.RData"))









########
########
## Appendix
########
########


# #### Random environment subset predictions #####
# 
# random_pred <- ls(pattern = "[0-9]{1,2}_predictions_random") %>% 
#   subset(., map_lgl(., ~inherits(get(.), "data.frame")))
# 
# 
# cv_pov_random_results <- map(random_pred, get) %>%
#   map(., ~rename_at(., vars(which(names(.) %in% "environment")), ~str_replace(., pattern = "environment", "val_environment"))) %>%
#   map_df(., ~mutate_at(., vars(which(names(.) %in% "rep")), as.character)) %>%
#   ## Add CV/POV number
#   mutate(scheme_number = as.character(str_extract(scheme, "[0-9]{1,2}")),
#          nTrainEnv = as.factor(nTrainEnv))
# 
# 
# 
# ## Analyze together?
# all_cv_pov_random_predictions_out <- cv_pov_random_results %>%
#   left_join(., env_herit) %>%
#   mutate(ability = accuracy, accuracy = ability / sqrt(heritability),
#          set = ifelse(is.na(set), "complete", set)) %>%
#   mutate_at(vars(scheme, val_environment), as.factor)
# 
# 
# ## Separately
# separate_cv_pov_random_predictions_analysis <- all_cv_pov_random_predictions_out %>%
#   # group_by(set, trait, scheme) %>%
#   group_by(set, trait, scheme_number) %>%
#   nest() %>%
#   mutate(model_rep = map(data, "rep") %>% map_lgl(~!all(is.na(.)))) %>%
#   mutate(out = list(NULL))
# 
# ## Loop for each row
# for (i in seq(nrow(separate_cv_pov_random_predictions_analysis))) {
#   
#   df <- separate_cv_pov_random_predictions_analysis$data[[i]]
#   
#   fit <- lmer(accuracy ~ 1 + nTrainEnv + scheme + nTrainEnv:scheme + (1|val_environment) + (1|val_environment:scheme) + (1|val_environment:nTrainEnv), data = df)
#   
#   effects_keep <- as.data.frame(Effect(focal.predictors = c("scheme", "nTrainEnv"), fit))
#   
#   fit_summary <- data_frame(
#     fitted = list(fit),
#     scheme_effects = list(effects_keep),
#     anova = list(tidy(anova(fit))),
#     ranova = list(tidy(ranova(fit)))
#   )
#   
#   separate_cv_pov_random_predictions_analysis$out[[i]] <- fit_summary
#   
#   
# }
# 
# 
# ## Plot all
# separate_cv_pov_random_predictions_analysis1 <- separate_cv_pov_random_predictions_analysis %>%
#   mutate(effects = map(out, ~.$scheme_effects[[1]])) %>%
#   unnest(effects) %>%
#   mutate(scheme = ifelse(scheme == "pov0", "pocv0", scheme),
#          scheme = factor(toupper(scheme), levels = cv_replace)) %>%
#   rename(accuracy = fit)
# 
# g_cv_pov_random_data <- separate_cv_pov_random_predictions_analysis1 %>%
#   ggplot(aes(x = nTrainEnv, y = accuracy, ymin = lower, ymax = upper)) +
#   geom_point() +
#   geom_errorbar(width = 0.5) +
#   scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
#   xlab("Number of training environments") +
#   facet_grid(trait ~ set + scheme, scales = "free_x", space = "free_x", labeller = labeller(set = str_to_title, trait = str_add_space)) +
#   theme_presentation2(base_size = 10) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggsave(filename = "all_cv_pov_random_predictions.jpg", plot = g_cv_pov_random_data, path = fig_dir, width = 10, height = 5, dpi = 1000)
# 
# 
# ## Add the accuracy when using all environments
# g_cv_pov_random_data1 <- separate_cv_pov_random_predictions_analysis1 %>%
#   left_join(., select(separate_cv_pov_predictions_analysis1, set, trait, scheme, fit = accuracy)) %>%
#   ggplot(aes(x = nTrainEnv, y = accuracy, ymin = lower, ymax = upper)) +
#   geom_hline(aes(yintercept = fit, lty = "All data")) +
#   geom_point() +
#   geom_errorbar(width = 0.5) +
#   scale_y_continuous(breaks = pretty, name = expression("Prediction accuracy ("*italic(r[MG])*")")) +
#   scale_linetype_manual(values = 2, name = NULL) +
#   xlab("Number of training environments") +
#   facet_grid(trait ~ set + scheme, scales = "free_x", space = "free_x", labeller = labeller(set = str_to_title, trait = str_add_space)) +
#   theme_presentation2(base_size = 10) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.05, 0.05), legend.margin = margin())
# 
# ggsave(filename = "all_cv_pov_random_predictions1.jpg", plot = g_cv_pov_random_data1, path = fig_dir, width = 10, height = 5, dpi = 1000)

