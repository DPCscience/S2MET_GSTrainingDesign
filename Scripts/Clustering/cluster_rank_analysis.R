## S2MET Predictions
## 
## Analyze clustering and distance rank results
## 
## Last modified: 31 October 2019
## Author: Jeff Neyhart
## 
## 
## 

## Directories, packages, and such
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))


# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))
# Load environmental correlations
load(file.path(result_dir, "environmental_genetic_correlations.RData"))

## Filter out undesirable rankings
environment_rank_df <- bind_rows(pred_env_dist_rank, pred_env_rank_random) %>%
  rename(val_environment = validation_environment) %>%
  filter(!mat_set %in% c("Jarquin", "MalosettiStand")) %>%
  # filter(model %in% names(dist_method_abbr_use)) %>%
  filter(set %in% c("complete", "realistic2017")) %>%
  select(-mat_set)




## Compare environment rankings with the correlation matrix
environment_rank_df1 <- environment_rank_df %>%
  split(.$trait) %>%
  map2_df(.x = ., .y = env_cor_tp[traits], ~{
    mutate(.x, rank = map(rank, function(e) subset(e, e %in% colnames(.y)))) %>%
      mutate(rank2 = map2(val_environment, rank, function(ve, te) .y[ve, te]))
  }) %>%
  mutate(rank2_summ = map(rank2, cummean)) %>%
  # mutate(rank2 = map(rank2, ~tibble(step = seq_along(.), mean_cor = .)))
  mutate_at(vars(rank2, rank2_summ), ~map(., ~tibble(step = seq_along(.), mean_cor = .)))


## Summarize the correlation of the next environment at each step
environment_rank_next_cor <- environment_rank_df1 %>%
  unnest(rank2) %>%
  mutate(model = ifelse(str_detect(model, "sample"), "sample", model),
         model = str_replace_all(model, dist_method_abbr_use)) %>%
  group_by(trait, set, model, step) %>%
  summarize_at(vars(mean_cor), list(~mean, ~sd)) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = dist_method_abbr_use))

## Plot
environment_rank_next_cor %>%
  ggplot(aes(x = step, y = mean, color = model)) +
  # geom_line(data = filter(environment_rank_avg_cor, model == "Random"), lwd = 1.5) +
  geom_line(lwd = ifelse(environment_rank_next_cor$model == "Random", 1.25, 0.5)) +
  scale_color_manual(values = dist_colors, guide = guide_legend(override.aes = list(lwd = 1))) +
  facet_grid(trait ~ set, switch = "y", scales = "free_x", space = "free_x", 
             labeller = labeller(trait = str_add_space, set = ~str_replace_all(., set_replace))) +
  scale_x_continuous(name = "Number of training environments", breaks = pretty) +
  scale_y_continuous(name = "Correlation with testing environment", breaks = pretty) +
  labs(subtitle = "Average correlation of next-ranked training\nenvironment with testing environment") +
  theme_presentation2(12) +
  theme(legend.position = "bottom")
# Save
ggsave(filename = "environment_correlation_next.jpg", path = fig_dir, width = 5, height = 8, dpi = 1000)



## Summarize the average correlation at each step
environment_rank_avg_cor <- environment_rank_df1 %>%
  unnest(rank2_summ) %>%
  mutate(model = ifelse(str_detect(model, "sample"), "sample", model),
         model = str_replace_all(model, dist_method_abbr_use)) %>%
  group_by(trait, set, model, step) %>%
  summarize_at(vars(mean_cor), list(~mean, ~sd)) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = dist_method_abbr_use))


## Plot
environment_rank_avg_cor %>%
  ggplot(aes(x = step, y = mean, color = model)) +
  # geom_line(data = filter(environment_rank_avg_cor, model == "Random"), lwd = 1.5) +
  geom_line(lwd = ifelse(environment_rank_avg_cor$model == "Random", 1.25, 0.5)) +
  scale_color_manual(values = dist_colors, guide = guide_legend(override.aes = list(lwd = 1))) +
  facet_grid(trait ~ set, switch = "y", scales = "free_x", space = "free_x", 
             labeller = labeller(trait = str_add_space, set = ~str_replace_all(., set_replace))) +
  scale_x_continuous(name = "Number of training environments", breaks = pretty) +
  scale_y_continuous(name = "Correlation with testing environment", breaks = pretty) +
  labs(subtitle = "Average correlation of training environments\nwith testing environment") +
  theme_presentation2(12) +
  theme(legend.position = "bottom")
# Save
ggsave(filename = "environment_correlation_all.jpg", path = fig_dir, width = 5, height = 8, dpi = 1000)






