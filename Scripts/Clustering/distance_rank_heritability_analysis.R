## Heritability based on environmental distance - analysis
## 
## Author: Jeff Neyhart
## Last Updated: October 9, 2018
## 
## This script will look at prediction accuracies from adding environments after
## ranking based on specific distance metrics
## 
## 



## Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Significance level
alpha <- 0.05


# Read in the results
load(file.path(result_dir, "cluster_heritability_tp_TESTING.RData"))

# Summarize the heritability
cluster_herit1 <- cluster_herit_out %>% 
  mutate(out = map(out, ~{data_frame(n_e = seq_along(.) + 1, heritability = map_dbl(., "heritability"))})) %>% 
  unnest()

## Summarize the random samples as a baseline
cluster_herit_random <- cluster_herit1 %>%
  filter(str_detect(model, "sample")) %>% 
  group_by(environment, trait, n_e) %>% 
  summarize_at(vars(heritability), funs(mean, sd, n())) %>% 
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), 
         lower = mean - stat, upper = mean + stat) %>% 
  select(-se, -stat) %>%
  ungroup() %>%
  mutate(model = "random")

# Subset the non-random models
cluster_herit_model <- cluster_herit1 %>% 
  filter(!str_detect(model, "sample")) %>%
  rename(mean = heritability)






# Plot
cluster_herit_model %>% 
  ggplot(aes(x = n_e, y = mean, color = model, fill = model)) +
  geom_ribbon(data = cluster_herit_random, aes(ymin = lower, ymax = upper), alpha = 0.25, fill = "grey75", color = "grey75") + 
  geom_line(data = cluster_herit_random, color = "grey75") + 
  geom_line() + 
  facet_grid(trait ~ environment, scales = "free_y") +
  # scale_color_manual(values = c("Random" = "grey75"), name = "Model") + 
  # scale_fill_manual(values = c("Random" = "grey75"), name = "Model") + 
  theme_acs()



### Calculate heritability as relative to using all environments
cluster_herit_stand <- cluster_herit1 %>% 
  group_by(environment, trait, model) %>% 
  mutate(rel_heritability = heritability - heritability[n()])

# First average over environments
cluster_herit_stand_summ <- cluster_herit_stand %>% 
  group_by(trait, model, n_e) %>% 
  summarize(herit = mean(rel_heritability)) %>%
  # Summarize over replicates
  ungroup() %>%
  mutate(model = ifelse(str_detect(model, "sample"), "random", model)) %>%
  group_by(trait, model, n_e) %>%
  summarize_at(vars(herit), funs(mean, sd, n())) %>% 
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), 
         lower = mean - stat, upper = mean + stat) %>% 
  select(-se, -stat) %>%
  ungroup()


# Vector of colors
colors <- setNames(c(RColorBrewer::brewer.pal(length(unique(cluster_herit_model$model)), name = "Set3"), "grey75"), 
                     nm = c(unique(cluster_herit_model$model), "random"))


# Plot
g_cluster_herit_stand <- cluster_herit_stand_summ %>% 
  ggplot(aes(x = n_e, y = mean, ymin = lower, ymax = upper, color = model, fill = model)) +
  geom_ribbon(alpha = 0.25) + 
  geom_line() + 
  facet_grid(trait ~ ., scales = "free_y", switch = "y") +
  scale_color_manual(values = colors, name = "Distance\nmeasure") +
  scale_fill_manual(values = colors, name = "Distance\nmeasure") +
  scale_x_continuous(breaks = pretty) +
  ylab("Relative heritability") +
  xlab("Number of environments") +
  theme_acs()
 
ggsave(filename = "distance_rank_heritability.jpg", plot = g_cluster_herit_stand, path = fig_dir, width = 5, height = 7, dpi = 1000)



## What models outperform random?
cluster_herit_stand_compare <- cluster_herit_stand_summ %>%
  filter(model != "random") %>% 
  left_join(., filter(cluster_herit_stand_summ, model == "random"), c("trait", "n_e")) %>% 
  select(trait, model = model.x, n_e, mean = mean.x, random_mean = mean.y, lower = lower.y, upper = upper.y)

# For each model and trait, what is the lower n_e value that is significantly different than random?
cluster_herit_stand_compare_model <- cluster_herit_stand_compare %>% 
  filter(mean != 0, mean >= upper) %>% 
  group_by(trait, model) %>% 
  filter(n_e == min(n_e))

cluster_herit_stand_compare_model %>%
  arrange(trait, n_e) %>%
  as.data.frame()

# Of the models that are better than random early-on, which are common to all traits?
common_models <- cluster_herit_stand_compare_model %>% 
  ungroup() %>% 
  filter(n_e <= 5) %>% 
  split(.$trait) %>% 
  map("model") %>%
  reduce(intersect) %>%
  # Make sure pheno dist and GCD are included
  union(., c("great_circle_dist", "pheno_dist", "random"))


## Plot only those models
g_cluster_herit_stand_subset <- cluster_herit_stand_summ %>% 
  filter(model %in% common_models) %>%
  ggplot(aes(x = n_e, y = mean, ymin = lower, ymax = upper, color = model, fill = model)) +
  geom_ribbon(alpha = 0.25) + 
  geom_line() + 
  facet_grid(trait ~ ., scales = "free_y", switch = "y") +
  scale_color_manual(values = colors, name = "Distance\nmeasure") +
  scale_fill_manual(values = colors, name = "Distance\nmeasure") +
  scale_x_continuous(breaks = pretty) +
  ylab("Relative heritability") +
  xlab("Number of environments") +
  theme_acs()

ggsave(filename = "distance_rank_heritability_common.jpg", plot = g_cluster_herit_stand_subset, path = fig_dir, width = 5, height = 7, dpi = 1000)










