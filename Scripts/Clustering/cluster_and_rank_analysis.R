## Analyze the clustering and heritability results
## 
## Author: Jeff Neyhart
## Last modified: 21 November 2019
## 
## This script will perform different clustering procedures on the S2MET data
## and calculate the within and across-cluster heritabilities.
## 


# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load packages and the source script
library(lme4)
library(broom)
library(lubridate)
library(cowplot)
library(patchwork)

# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))
# Load the heritability results
load(file.path(result_dir, "cluster_heritability_results.RData"))
load(file.path(result_dir, "environmental_genetic_correlations.RData"))





## Assign traits to whether a higher value is favorable
higher_favorable_traits <- "GrainYield"

# Edit the traits to make higher value more favorable
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  mutate(value = if_else(trait %in% higher_favorable_traits, value, value * -1))

## Plot heatmaps
heat_colors <- wesanderson::wes_palette("Zissou1")



#############################
## Distance rank analysis
#############################

## Edit model names
pred_env_dist_rank_use <- pred_env_dist_rank %>%
  mutate(model = ifelse(model == "pheno_location_dist", "pheno_loc_dist", model),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr_use)) %>%
  filter(model %in% dist_method_abbr_use,
         set == "complete")


## For every environment, take the first 5 most similar environments.
## Look for patterns

pred_env_dist_rank_first_five <- pred_env_dist_rank_use %>%
  mutate(first_five = map(rank, ~.[1:5])) %>%
  # How many of these match the year of the validation environment?
  mutate(match_year = map2_dbl(validation_environment, first_five, 
                               ~sum(str_count(string = .y, pattern = str_extract(.x, "[0-9]{2}")))),
         # How many match location?
         match_loc = map2_dbl(validation_environment, first_five, 
                              ~sum(str_count(string = .y, pattern = str_extract(.x, "[A-Z]{3}")))))

## Find the number of shared first-five environments between ranks
pred_env_dist_rank_first_five_avg <- pred_env_dist_rank_first_five %>%
  mutate(nE = map_dbl(rank, length)) %>% # Number of available environments
  select(-rank, -contains("match")) %>%
  group_by(set, trait, validation_environment) %>%
  do(crossing(., select(., model, first_five))) %>%
  ungroup() %>%
  mutate(num_overlap = map2_dbl(.x = first_five, .y = first_five1, ~length(intersect(.x, .y))),
         prop_overlap = num_overlap / 5,
         prop_overlap_chance = 5 / nE) %>%
  select(-contains("first_five"))
  


## Find the pairwise spearman rank correlation between rankings for a particular test environment
pred_env_dist_rank_cor <- pred_env_dist_rank_use %>%
  group_by(trait, set, validation_environment) %>%
  do({
    df <- .
    
    # Crossing
    crossing(df, select(df, model, rank)) %>%
      # filter(model != model1) %>%
      # Use xtfrm to convert a factor into a numeric
      mutate(rank_cor = map2_dbl(.x = rank, .y = rank1, ~cor(x = xtfrm(.x), y = xtfrm(.y), method = "spearman"))) %>%
      select(-rank, -rank1)
    
  }) %>% ungroup()


## Add the proportion of overalp of first five
## Determine the average correlation between rankings for each pair of similarity
## measures for each trait
pred_env_dist_rank_cor_avg <- pred_env_dist_rank_cor %>%
  left_join(., pred_env_dist_rank_first_five_avg) %>%
  group_by(set, trait, model, model1) %>%
  summarize_at(vars(rank_cor, num_overlap, prop_overlap, prop_overlap_chance), mean) %>%
  ungroup() %>%
  split(list(.$set, .$trait)) %>%
  # Assign the lower triangle
  map_df(~mutate(., lower_triangle = duplicated(
    map2_chr(.x = model, .y = model1, ~paste(sort(c(as.character(.x), as.character(.y))), collapse = "_")))) )

## What is the range in correlation coefficients?
pred_env_dist_rank_cor_avg %>% 
  filter(model != model1) %>%
  group_by(trait) %>% 
  summarize_at(vars(rank_cor), list(~mean, ~min, ~max))
  
# trait         mean      min   max
# 1 GrainYield  0.128   0.00675 0.299
# 2 HeadingDate 0.0622 -0.0599  0.197
# 3 PlantHeight 0.0248 -0.0659  0.297


## Compare AMMI with others
pred_env_dist_rank_cor_avg %>%
  filter_at(vars(contains("model")), any_vars(. == "AMMI")) %>%
  filter(model != model1, !lower_triangle) %>%
  as.data.frame() %>%
  arrange(trait, set, desc(rank_cor)) %>%
  split(.$trait)


## AMMI is most correlated with PD or LocPD,
## but then the next model is often IPCA-EC

## Compare PD with others
pred_env_dist_rank_cor_avg %>%
  filter_at(vars(contains("model")), any_vars(. == "PD")) %>%
  filter(model != model1, !lower_triangle) %>%
  as.data.frame() %>%
  arrange(trait, set, desc(rank_cor)) %>%
  split(.$trait)


## PD is most correlated with PD or LocPD, but it is consistently
## correlated with IPCA-EC (0.073 - 0.15)


## Compare IPCA-EC with others
pred_env_dist_rank_cor_avg %>%
  filter_at(vars(contains("model")), any_vars(. == "IPCA-EC")) %>%
  filter(model != model1, !lower_triangle) %>%
  as.data.frame() %>%
  arrange(trait, set, desc(rank_cor)) %>%
  split(.$trait)

# IPCA-EC is not consistently correlated with anything,
# but it is often correlated with LocPD (0.16 - 0.30)



# Glyph modifier
draw_key_text1 <- function(data, params, size) {
  
  grobTree(
    textGrob(label = "3.3", x = 0.5, y = 0.5),
    gp = gpar(
      col = data$colour %||% "grey20",
      fill = alpha(data$fill %||% "grey85", data$alpha),
      lwd = (data$size %||% 0.5) * .pt,
      lty = data$linetype %||% 1
    )
  )
  
  
}

## Plot using a heatmap
# The upper triangle will be a heat color,
# the lower triangle will be the proportion of overapping first-five environments
g_dist_rank_cor_heatmap <- pred_env_dist_rank_cor_avg %>% 
  filter(set == "complete") %>%
  mutate(set = f_set_replace(set, abbr = T),
         # If lower triangle, make rank_cor blank
         rank_cor = case_when(lower_triangle ~ as.numeric(NA), model == model1 ~ as.numeric(NA), TRUE ~ rank_cor), 
         overlap_use = num_overlap,
         label = ifelse(lower_triangle, round(overlap_use, 1), NA)) %>% # If lower triangle, make the label
  ## Create plot
  ggplot(aes(x = model, y = model1, fill = rank_cor)) + 
  geom_tile() + 
  # Add alternative tile for the diagonal
  geom_tile(data = subset(pred_env_dist_rank_cor_avg, model == model1), fill = "white") +
  geom_text(aes(label = label, size = overlap_use), key_glyph = draw_key_text1) +
  scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[3], high = heat_colors[5], 
                       na.value = "grey95", name = "Rank\ncorrelation", midpoint = median(pred_env_dist_rank_cor_avg$rank_cor),
                       guide = guide_colorbar(title.position = "top")) +
  facet_grid(trait ~ ., labeller = labeller(trait = str_add_space), switch = "y", scales = "free", space = "free") +
  scale_x_discrete(position = "bottom", name = "Similarity measure") +
  scale_size_continuous(name = "Overlapping\nfirst five\nenvironments", labels = NULL, breaks = 3.3,
                        guide = guide_legend(override.aes = list(fill = "grey85"))) +
  # labs(subtitle = unique(.$set)) +
  theme_presentation2(base_size = 8) +
  theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1), strip.placement = "outside",
        legend.position = "right", panel.border = element_blank() )

## Save
ggsave(filename = "dist_rank_cor_heatmap.jpg", plot = g_dist_rank_cor_heatmap, path = fig_dir,
       height = 6, width = 3.5, dpi = 1000)

  
  
  
  
  
#### Compare distance rankings to phenotypic correlations #####


## Filter out undesirable rankings
environment_rank_df <- bind_rows(pred_env_dist_rank, pred_env_rank_random) %>%
  mutate(model = ifelse(model == "pheno_location_dist", "pheno_loc_dist", model),
         model = ifelse(str_detect(model, "sample"), "sample", model),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr_use)) %>%
  filter(model %in% dist_method_abbr_use,
         set == "complete") %>%
  rename(val_environment = validation_environment)


## Compare environment rankings with the correlation matrix
environment_rank_df1 <- environment_rank_df %>%
  split(.$trait) %>%
  map2_df(.x = ., .y = env_cor_tp[traits], ~{
    mutate(.x, rank = map(rank, function(e) subset(e, e %in% colnames(.y)))) %>%
      mutate(rank_cor_window = map2(val_environment, rank, function(ve, te) .y[ve, te]))
  }) %>%
  mutate(rank_cor_cummean = map(rank_cor_window, cummean)) %>%
  mutate_at(vars(rank_cor_window, rank_cor_cummean), ~map(., ~tibble(step = seq_along(.), mean_cor = .)))


## Summarize the correlation of the next environment at each step
environment_rank_cor_window <- environment_rank_df1 %>%
  unnest(rank_cor_window) %>%
  group_by(trait, set, model, step) %>%
  summarize_at(vars(mean_cor), list(~mean, ~sd, ~n())) %>%
  ungroup() %>%
  mutate(se = ifelse(model == "Random", sd / sqrt(n), NA))

## Plot
g_environment_rank_cor_window <- environment_rank_cor_window %>%
  ggplot(aes(x = step, y = mean, color = model)) +
  # geom_line(data = filter(environment_rank_avg_cor, model == "Random"), lwd = 1.5) +
  geom_line(lwd = ifelse(environment_rank_cor_window$model == "Random", 1.25, 0.5)) +
  scale_color_manual(values = dist_colors, guide = guide_legend(override.aes = list(lwd = 1)),
                     name = "Similarity measure") +
  facet_grid(trait ~ ., switch = "y", scales = "free", 
             labeller = labeller(trait = str_add_space)) +
  scale_x_continuous(name = "Number of training environments", breaks = pretty) +
  scale_y_continuous(name = "Phenotypic correlation with testing environment", breaks = pretty) +
  labs(subtitle = "Moving phenotypic correlation of training and test environments") +
  theme_presentation2(12) +
  theme(legend.position = "bottom")
# Save
ggsave(filename = "environment_correlation_window.jpg", plot = g_environment_rank_cor_window, 
       path = fig_dir, width = 5, height = 8, dpi = 1000)



## Summarize the average correlation at each step
environment_rank_avg_cor <- environment_rank_df1 %>%
  unnest(rank_cor_cummean) %>%
  group_by(trait, set, model, step) %>%
  summarize_at(vars(mean_cor), list(~mean, ~sd, lower = ~quantile(., probs = c(alpha / 2)),
                                    upper = ~quantile(., probs = 1 - (alpha / 2)))) %>%
  ungroup()


## Plot
g_environment_rank_avg_cor <- environment_rank_avg_cor %>%
  ggplot(aes(x = step, y = mean, color = model)) +
  # geom_ribbon(data = filter(environment_rank_avg_cor, model == "Random"), aes(ymin = lower, ymax = upper), fill = "grey85",
  #             color = NA, alpha = 0.5) +
  # geom_line(data = filter(environment_rank_avg_cor, model == "Random"), lwd = 1.5) +
  geom_line(lwd = ifelse(environment_rank_avg_cor$model == "Random", 1.25, 0.5)) +
  scale_color_manual(values = dist_colors, guide = guide_legend(override.aes = list(lwd = 1)),
                     name = "Similarity\nmeasure") +
  facet_grid(trait ~ ., switch = "y", scales = "free", space = "free_x", 
             labeller = labeller(trait = str_add_space)) +
  scale_x_continuous(name = "Number of training environments", breaks = pretty) +
  scale_y_continuous(name = "Phenotypic correlation with testing environment", breaks = pretty) +
  labs(subtitle = "Cumulative phenotypic correlation of\ntraining and test environments") +
  theme_presentation2(12) +
  theme(legend.position = "bottom")
# Save
ggsave(filename = "environment_correlation_cummean.jpg", plot = g_environment_rank_avg_cor, 
       path = fig_dir, width = 5, height = 8, dpi = 1000)














#############################
## cluster analysis
#############################




### Are clusters defined more by location or environment?
## Compare the proportion of each year/location in a cluster with the expectation of random
## environment assignments and the cluster assignments 
## What is the proportion of environments assigned to the same cluster across distance method
cluster_grouping <- cluster_df %>% 
  select(-data) %>% 
  group_by(set, trait, model) %>% 
  do({
    row <- .
    
    clus <- row$cluster[[1]] %>%
      mutate(location = str_extract(environment, "[A-Z]{3}"),
             year = str_extract(environment, "[0-9]{2}"))
    
    # List of all locations/years
    locs <- clus$location
    yrs <- clus$year
    
    ## Calculate proportions of locations and year
    clus_props <- clus %>% 
      group_by(cluster) %>% 
      summarize_at(vars(location, year), funs(list(as.data.frame(table(.))))) %>%
      mutate(nEnv = map(location, "Freq") %>% map_dbl(sum),
             nLoc = map_dbl(location, ~n_distinct(.$location)),
             nYr = map_dbl(year, ~n_distinct(.$year))) %>%
      # Calculate expectations
      mutate(expLoc = map_dbl(nEnv, ~mean(replicate(1000, n_distinct(sample(locs, .))))),
             expYr = map_dbl(nEnv, ~mean(replicate(1000, n_distinct(sample(yrs, .)))))) %>%
      # Calculate p value from exact binomial test
      mutate(LocTest = pmap_dbl(list(nLoc, nEnv, expLoc), ~binom.test(x = ..1, n = ..2, p = ..3 / ..2)$p.value),
             YrTest = pmap_dbl(list(nYr, nEnv, expYr), ~binom.test(x = ..1, n = ..2, p = ..3 / ..2)$p.value))
    
    
    ## Calculate some summary
    summarize_at(clus_props, vars(nEnv:expYr), mean) %>%
      mutate(LocTest = prop.test(x = nLoc, n = nEnv, p = expLoc / nEnv )$p.value,
             YrTest = prop.test(x = nYr, n = nEnv, p = expYr / nEnv)$p.value)
    
    
  }) %>% ungroup()



## Analyze tests
cluster_grouping %>%
  select(set:model, LocTest, YrTest) %>%
  gather(test, p_value, LocTest, YrTest) %>%
  # filter(p_value < alpha)
  arrange(p_value)


## Table
cluster_grouping %>%
  mutate(LocDev = nLoc - expLoc, YrDev = nYr - expYr) %>% 
  select(set:model, LocDev, YrDev) %>% 
  gather(term, deviation, LocDev, YrDev) %>% 
  split(list(.$set, .$trait, .$term)) %>%
  map_df(~arrange(., desc(abs(deviation)))) %>%
  write_csv(x = ., path = file.path(fig_dir, "cluster_location_year_bias.csv"))











## How many clusters per trait, set, and model?
cluster_summary <- cluster_df %>% 
  unnest(cluster) %>% 
  mutate(model = ifelse(model == "pheno_location_dist", "pheno_loc_dist", model),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr_use)) %>%
  filter(model %in% dist_method_abbr_use) %>%
  group_by(set, trait, model) %>% 
  summarize(nCluster = n_distinct(cluster)) %>% 
  spread(model, nCluster)




## What is the proportion of environments assigned to the same cluster across distance method
cluster_compare <- cluster_df %>%
  select(-data) %>% 
  group_by(set, trait) %>% 
  do(compare_clusters(x = .)) %>% 
  ungroup()



## Filter out
cluster_compare1 <- cluster_compare %>%
  mutate_at(vars(model, model1), funs(ifelse(. == "pheno_location_dist", "pheno_loc_dist", .))) %>%
  mutate_at(vars(model, model1), funs(str_replace_all(., dist_method_abbr))) %>%
  filter_at(vars(model1, model), all_vars(. %in% dist_method_abbr_use)) %>%
  mutate_at(vars(model, model1), funs(factor(., levels = dist_method_abbr_use))) %>%
  arrange(set, trait, model, model1) %>%
  # Designate upper or lower triangle
  split(list(.$set, .$trait)) %>%
  map_df(~mutate(., lower_triangle = duplicated(map2_chr(.x = model, .y = model1, ~paste(sort(c(as.character(.x), as.character(.y))), collapse = "_")))))

## Compare
cluster_compare1 %>% 
  group_by(set, model) %>% 
  summarize(mean = mean(prop_common))


cluster_compare1 %>% 
  filter_at(vars(model, model1), all_vars(. %in% c("AMMI", "PD", "LocPD"))) %>%
  group_by(set, model, model1) %>% 
  summarize(mean = mean(prop_common))

cluster_compare1 %>% 
  filter_at(vars(model, model1), all_vars(. %in% c("All-EC", "Mean-EC", "IPCA-EC"))) %>%
  group_by(set, model) %>% 
  summarize(mean = mean(prop_common))


## Subset some silimarity measures for publication
similarity_measure_subset <- c("AMMI", "LocPD", "GCD", "IPCA-EC")



cluster_compare_hm_list <- cluster_compare1 %>% 
  filter(!(str_detect(set,"realistic") & model == "AMMI")) %>%
  # filter(set %in% c("complete", "realistic2017")) %>%
  filter(trait %in% traits) %>% 
  filter_at(vars(contains("model")), ~. %in% similarity_measure_subset) %>%
  mutate(set = str_replace_all(set, set_replace),
         # Put the year in parentheses
         set = str_replace(string = set, pattern = "([A-Za-z-]*)([0-9]{4})", replacement = "\\1 (\\2)"),
         prop_common = ifelse(model == model1, NA, prop_common),
         annotation_use = ifelse(model == model1, annotation1, round(prop_common, 2))) %>%
  filter(!lower_triangle) %>%
  split(.$set) %>%
  map(~ggplot(data = ., aes(x = model, y = model1, fill = prop_common)) + 
        geom_tile() +  
        geom_text(aes(label = annotation_use), size = 1.75) + 
        scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[3], high = heat_colors[5], 
                             midpoint = 0.5, limits = c(0, 1), na.value = "grey95", name = "Proportion of overlap  ") +
        facet_wrap(~ trait, labeller = labeller(trait = str_add_space)) +
        scale_x_discrete(position = "top") +
        labs(subtitle = unique(.$set)) +
        theme_presentation2(base_size = 8) +
        theme(axis.title = element_blank(), axis.text.x = element_text(angle = 325, hjust = 1), strip.placement = "outside",
              legend.position = "top", legend.key.height = unit(0.75, "line"), panel.border = element_blank() ) )


## Combine plots
g_cluster_compare <- plot_grid(plotlist = map(cluster_compare_hm_list, ~. + theme(legend.position = "none")), ncol = 2, 
                               labels = letters[seq_along(cluster_compare_hm_list)], label_size = 12, align = "hv")
g_cluster_compare1 <- plot_grid(g_cluster_compare, get_legend(cluster_compare_hm_list[[1]] + theme(legend.position = "right")), 
                                nrow = 1, rel_widths = c(1, 0.15))

ggsave(filename = "cluster_compare_heatmap.jpg", plot = g_cluster_compare1, path = fig_dir, width = 11, height = 5, dpi = 1000)



## Portrait orientation
# Remove strip from all plots except first
cluster_compare_hm_list1 <- cluster_compare_hm_list
g_cluster_compare <- plot_grid(plotlist = map(cluster_compare_hm_list1, ~. + theme(legend.position = "none")), ncol = 1, 
                               labels = letters[seq_along(cluster_compare_hm_list1)], label_size = 10, align = "hv")
g_cluster_compare1 <- plot_grid(g_cluster_compare, get_legend(cluster_compare_hm_list1[[1]]), ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "cluster_compare_heatmap_portrait.jpg", plot = g_cluster_compare1, 
       path = fig_dir, width = 4, height = 8, dpi = 1000)



### Plot all scenarios as one
cluster_compare_hm <- cluster_compare1 %>% 
  filter(!(str_detect(set,"realistic") & model == "AMMI")) %>%
  # filter(set %in% c("complete", "realistic2017")) %>%
  filter(trait %in% traits) %>% 
  filter_at(vars(contains("model")), ~. %in% similarity_measure_subset) %>%
  mutate(set = f_set_replace(set, abbr = T),
         prop_common = ifelse(model == model1, NA, prop_common),
         annotation_use = ifelse(model == model1, annotation1, round(prop_common, 2))) %>%
  filter(!lower_triangle) %>%
  
  ## Create plot
  ggplot(data = ., aes(x = model, y = model1, fill = prop_common)) + 
  geom_tile() +  
  geom_text(aes(label = annotation_use), size = 1.75) + 
  scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[3], high = heat_colors[5], 
                       midpoint = 0.5, limits = c(0, 1), na.value = "grey95", name = "Proportion of\nenvironmental overlap  ") +
  facet_grid(trait ~ set, labeller = labeller(trait = str_add_space), switch = "both", scales = "free", space = "free") +
  scale_x_discrete(position = "top", name = "Similarity measure") +
  # labs(subtitle = unique(.$set)) +
  theme_presentation2(base_size = 8) +
  theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 325, hjust = 1), strip.placement = "outside",
        legend.position = "bottom", legend.key.height = unit(0.75, "line"), panel.border = element_blank() )

ggsave(filename = "cluster_compare_heatmap_combined.jpg", plot = cluster_compare_hm, 
       path = fig_dir, width = 4, height = 4, dpi = 1000)



# #### Greyscale version
# cluster_compare_hm_list <- cluster_compare1 %>% 
#   filter(!(str_detect(set,"realistic") & model == "AMMI")) %>%
#   filter(set %in% c("complete", "realistic2017")) %>%
#   filter(trait %in% traits) %>% 
#   mutate(set = str_replace_all(set, set_replace),
#          # set = str_replace_all(set, "Time-forward", "Leave-one-year-out "),
#          set = str_replace_all(set, "2017", ""),
#          prop_common = ifelse(model == model1, NA, prop_common),
#          annotation_use = ifelse(model == model1, annotation1, round(prop_common, 2))) %>%
#   filter(!lower_triangle) %>%
#   split(.$set) %>%
#   map(~ggplot(data = ., aes(x = model, y = model1, fill = prop_common)) + 
#         geom_tile() +  
#         geom_text(aes(label = annotation_use), size = 2) + 
#         scale_fill_gradient2(low = grey.colors(11)[11], mid = grey.colors(11)[5], high = grey.colors(10)[1], 
#                              midpoint = 0.5, limits = c(0, 1), na.value = "grey95",
#                              name = "Proportion of\noverlap") +
#         facet_wrap(~ trait, labeller = labeller(trait = str_add_space)) +
#         labs(subtitle = unique(.$set)) +
#         theme_presentation2(base_size = 10) +
#         theme(axis.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom") )
# 
# g_cluster_compare <- plot_grid(plotlist = map(cluster_compare_hm_list, ~. + theme(legend.position = "none")), ncol = 1, 
#                                labels = LETTERS[seq_along(cluster_compare_hm_list)], align = "hv")
# # g_cluster_compare1 <- plot_grid(g_cluster_compare, get_legend(cluster_compare_hm_list[[1]]), ncol = 1, rel_heights = c(1, 0.1))
# g_cluster_compare1 <- plot_grid(g_cluster_compare, get_legend(cluster_compare_hm_list[[1]] + theme(legend.position = "right")), 
#                                 nrow = 1, rel_widths = c(1, 0.15))
# 
# 
# ggsave(filename = "cluster_compare_heatmap_greyscale.jpg", plot = g_cluster_compare1, path = fig_dir, width = 6.5, height = 5, dpi = 1000)









#### Clusters versus non-cluster heritability


## Extract heritability
cluster_herit <- cluster_varcomp %>%
  mutate_at(vars(undivided, subdivided), ~map_dbl(., "heritability")) %>%
  rename(clustered = subdivided, unclustered = undivided) %>%
  mutate(model = ifelse(model  == "pheno_location_dist", "pheno_loc_dist", model),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr_use)) %>%
  filter(model %in% dist_method_abbr_use)

## Calculate ratio
cluster_herit %>%
  mutate(effectiveness = clustered / unclustered) %>%
  select(-unclustered, -clustered) %>%
  spread(model, effectiveness)

## Plot heritability
g_cluster_herit <- cluster_herit %>%
  filter(!(model == "AMMI" & set == "Time-forward")) %>%
  filter(model %in% similarity_measure_subset) %>%
  gather(group, heritability, -set, -model, -trait) %>%
  mutate(group = ifelse(group == "unclustered", paste0("H[1]~(", group, ")"), paste0("H[2]~(", group, ")")),
         set = f_set_replace(set, abbr = T)) %>%
  # Plot
  ggplot(aes(x = model, y = heritability, color = model, shape = group)) +
  geom_point(position = position_dodge2(0.5), size = 2) +
  scale_color_manual(values = dist_colors_use, name = "Distance\nmethod", guide = FALSE) +
  scale_shape_discrete(name = NULL, labels = function(x) parse(text = x), guide = guide_legend(nrow = 2, label.position = "left")) +
  facet_grid(trait ~ set, scales = "free_x", space = "free_x", switch = "y",
             labeller = labeller(trait = str_add_space)) +
  xlab("Distance measure") +
  scale_y_continuous(breaks = pretty, name = "Heritability") +
  theme_presentation2(base_size = 10) +
  theme(legend.position = c(0.12, 0.07), legend.margin = margin(), legend.key.height = unit(0.75, "line"), 
        axis.text.x = element_text(angle = 45, hjust = 1), strip.placement = "outside")

ggsave(filename = "cluster_heritability_complete.jpg", plot = g_cluster_herit, path = fig_dir, 
       width = 6, height = 4.5, dpi = 1000)



cluster_varcomp1 <- cluster_varcomp %>%
  mutate_at(vars(undivided, subdivided), ~map(., "var_comp")) %>%
  mutate(model = ifelse(model  == "pheno_location_dist", "pheno_loc_dist", model),
         model = str_replace_all(model, dist_method_abbr),
         model = factor(model, levels = dist_method_abbr_use),
         set = str_to_title(set)) %>%
  filter(model %in% dist_method_abbr_use) %>%
  select(-subdivided) %>%
  unnest() %>%
  ## Calculate variance proportions
  group_by(set, model, trait) %>%
  mutate(varprop = variance / sum(variance)) 







## Plot
g_cluster_varcomp <- cluster_varcomp1 %>%
  ggplot(aes(x = model, y = varprop, fill = source)) +
  geom_col() +
  facet_grid(trait ~ set, scales = "free_x", space = "free_x") +
  theme_presentation2() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1))



# Save a table
cluster_varcomp_table <- cluster_varcomp1 %>%
  mutate(annotation = list(NULL)) %>%
  mutate(annotation = str_trim(paste0(formatC(x = varprop, digits = 2))),
         annotation = str_remove_all(annotation, "NA")) %>%
  select(set, trait, model, source, annotation) %>%
  spread(source, annotation)

write_csv(x = cluster_varcomp_table, path = file.path(fig_dir, "cluster_varcomp_table.csv"))


## Of all trait-model combinations, which result in significant G X cluster?
# model       trait              term     variance        LRT df      p.value
# 1    great_circle_dist  GrainYield line_name:cluster 3.242287e+04  61.206548  1 5.139032e-15
# 2    great_circle_dist HeadingDate line_name:cluster 2.764723e-01  10.605145  1 1.127734e-03
# 3             MYEC_All  GrainYield line_name:cluster 3.479472e+04  56.852976  1 4.696428e-14
# 4             MYEC_All HeadingDate line_name:cluster 8.207656e-01  76.461140  1 2.245892e-18
# 5             MYEC_All PlantHeight line_name:cluster 7.876960e-01   7.698212  1 5.527556e-03
# 6            MYEC_IPCA  GrainYield line_name:cluster 2.189498e+04  43.586845  1 4.055561e-11
# 7            MYEC_IPCA HeadingDate line_name:cluster 6.151390e-01  43.976911  1 3.322725e-11
# 8            MYEC_Mean  GrainYield line_name:cluster 2.477421e+04  39.265137  1 3.699888e-10
# 9            MYEC_Mean HeadingDate line_name:cluster 8.828733e-01  88.245794  1 5.780602e-21
# 10           MYEC_Mean PlantHeight line_name:cluster 7.727645e-01  10.012722  1 1.554626e-03
# 11            OYEC_All  GrainYield line_name:cluster 1.649557e+04  19.217038  1 1.166675e-05
# 12            OYEC_All HeadingDate line_name:cluster 4.445435e-01  19.614002  1 9.477208e-06
# 13            OYEC_All PlantHeight line_name:cluster 1.575744e+00  29.782854  1 4.832471e-08
# 14           OYEC_IPCA  GrainYield line_name:cluster 4.357908e+04  42.842169  1 5.933914e-11
# 15           OYEC_IPCA HeadingDate line_name:cluster 4.605247e-01  33.263482  1 8.047986e-09
# 16           OYEC_Mean  GrainYield line_name:cluster 2.195249e+04  38.632631  1 5.115646e-10
# 17           OYEC_Mean HeadingDate line_name:cluster 7.379050e-01  63.731684  1 1.425730e-15
# 18           OYEC_Mean PlantHeight line_name:cluster 5.870112e-01   9.630360  1 1.913871e-03
# 19          pheno_dist  GrainYield line_name:cluster 7.501530e+04 163.011091  1 2.487701e-37
# 20          pheno_dist HeadingDate line_name:cluster 1.165701e+00  87.328155  1 9.192991e-21
# 21          pheno_dist PlantHeight line_name:cluster 3.291683e+00 199.800223  1 2.309020e-45
# 22 pheno_location_dist  GrainYield line_name:cluster 2.789433e+04  40.461515  1 2.005272e-10
# 23 pheno_location_dist HeadingDate line_name:cluster 1.871909e+00 169.755848  1 8.365500e-39
# 24 pheno_location_dist PlantHeight line_name:cluster 9.410029e-01  17.598128  1 2.728568e-05

## Line_name:cluster:


























