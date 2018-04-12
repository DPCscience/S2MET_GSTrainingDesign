## Analysis of predictions based on environmental distance
## 
## Author: Jeff Neyhart
## Last Updated: April 11, 2018
## 
## This script will look at prediction accuracies from adding environments after
## ranking based on specific distance metrics
## 


# # Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load some packages
library(ggforce)

# Load the results
load(file.path(result_dir, "environmental_distance_predictions.RData"))
load(file.path(result_dir, "environmental_distance_window_predictions.RData"))
load(file.path(result_dir, "environmental_distance_heritability.RData"))




## Bind the list elements together and unnest
cumulative_pred_results <- env_dist_predictions_out %>% 
  bind_rows() %>% 
  unnest()


# Re-scale the distance measurements to unit variance
cumulative_pred_toplot <- cumulative_pred_results %>%
  # Remove the environmental correlation method
  filter(dist_method %in% names(dist_method_replace)) %>%
  group_by(environment, trait, dist_method) %>% 
  mutate(scaled_distance = scale(distance),
         scaled_distance = scaled_distance - min(scaled_distance)) %>%
  ungroup() %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace),
         n_train_env = map_dbl(train_envs, length))





## For each prediction environment, trait, and distance method, find the prediction
## accuracy after adding all environments (i.e. the last value) and the local
## maximum prediction accuracy (if present)
## 
cumulative_pred_toplot_maxima <- cumulative_pred_toplot %>% 
  group_by(environment, trait, dist_method) %>% 
  mutate(local_max_accuracy = max(accuracy),
         local_max_n_train = n_train_env[accuracy == local_max_accuracy],
         local_max_scaled_dist = scaled_distance[accuracy == local_max_accuracy],
         all_env_acccuracy = accuracy[n_train_env == max(n_train_env)]) %>%
  ungroup() %>%
  select(-train_envs) %>%
  # Is the local maximum greater than the final accuracy?
  # What is the difference between the local maximum and the final
  mutate(local_adv = local_max_accuracy - all_env_acccuracy,
         local_better = (local_adv > 0))

## Filter for the observations in which the accuracy is the local maximum
cumulative_pred_toplot_maxima_filter <- cumulative_pred_toplot_maxima %>%
  filter(accuracy == local_max_accuracy)
  



## Distribution of the number of training environments or the scaled distance (cumulative) 
## of the local maximum
g_local_max <- cumulative_pred_toplot_maxima_filter %>% 
  ggplot(aes(x = n_train_env, fill = dist_method)) +
  # ggplot(aes(x = scaled_distance, fill = dist_method)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~ trait, ncol = 2) +
  theme_bw()

## Distribution of local advantage
g_local_adv <- cumulative_pred_toplot_maxima_filter %>% 
  ggplot(aes(x = local_adv, fill = dist_method)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~ trait, ncol = 2) +
  theme_bw()


## Calculate the average environmental distance (scaled) where the local maxima is
## found
cumulative_pred_maxima_summary <- cumulative_pred_toplot_maxima_filter %>% 
  group_by(trait, dist_method) %>% 
  summarize_at(vars(scaled_distance, n_train_env, local_adv, local_better), mean)




  

# ## Sample environments that have all traits
# cumulative_pred_toplot_sample <- cumulative_pred_toplot %>% 
#   group_by(environment) %>% filter(n_distinct(trait) == 3) %>% ungroup() %>%
#   filter(environment %in% sample(unique(.$environment), 5))

## Common plot modifier
g_mod <- list(
  geom_point(aes(x = local_max_scaled_dist, y = local_max_accuracy)),
  geom_line(lwd = 0.5),
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)),
  ylab("Prediction Accuracy"),
  xlab("Cumulative Distance of Environments in Training Set"),
  theme_bw(),
  theme(panel.grid = element_blank(),
        legend.key.height = unit(2, units = "line"),
        text = element_text(size = 8)),
  labs(title = "Cumulative Environmental Clusters")
)

# List over traits
g_cumulative_pred <- setNames(traits, traits) %>%
  map(function(tr) {
    # Subset the trait
    cumulative_pred_toplot_tr <- cumulative_pred_toplot_maxima %>%
      filter(trait == tr) %>%
      mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
    
    # Create one plot using facet_wrap
    g_total <- cumulative_pred_toplot_tr %>% 
      ggplot(aes(x = scaled_distance, y = accuracy, color = dist_method)) +
      g_mod +
      facet_wrap( ~ environment + trait) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = dist_method, color = NULL), alpha = 0.1)
    
    # Save
    save_file <- file.path(fig_dir, str_c("cumulative_env_dist_pred_", tr, ".jpg"))
    ggsave(filename = save_file, plot = g_total, height = 15, width = 12, dpi = 1000)
    
    # Return the plot
    # return(g_pages)
    return(g_total)
    
  })


### Plot a single distance method over many environments
# List over traits
g_cumulative_pred1 <- dist_method_replace %>%
  map(function(dm) {
    
    # Subset the distance method
    cumulative_pred_toplot_dm <- cumulative_pred_toplot_maxima %>%
      filter(dist_method == dm)
    
    # Iterate over the trait
    g_list <- setNames(traits, traits) %>%
      map(function(tr) {
        # Subset the trait data
        cumulative_pred_toplot_tr <- cumulative_pred_toplot_dm %>%
          filter(trait == tr) %>%
          mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
        
        # Extract the ylim values from the single plots created above
        ylim_tr <- ggplot_build(g_cumulative_pred[[tr]])$layout$panel_ranges[[1]]$y.range
        
        # Create the plot and return
        cumulative_pred_toplot_tr %>% 
          ggplot(aes(x = scaled_distance, y = accuracy)) +
          g_mod + 
          ylim(ylim_tr) +
          facet_wrap( ~ environment + trait + dist_method)  +
          geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.1)
        
      })
      
    # Iterate over the plots and save
    for (i in seq_along(g_list)) {
      # Save
      save_file <- file.path(
        fig_dir, str_c("cumulative_env_dist_pred_", names(dist_method_replace[dist_method_replace == dm]), 
                       "_", names(g_list)[i], ".jpg"))
      ggsave(filename = save_file, plot = g_list[[i]], height = 15, width = 12, dpi = 1000)
    }
    
    # Return the plot
    return(g_list)
    
  })











## Prediction accuracy using a sliding window
## Bind the list elements together and unnest
window_pred_results <- env_dist_window_predictions_out %>% 
  bind_rows() %>% 
  unnest()


# Re-scale the distance measurements to unit variance
window_pred_toplot <- window_pred_results %>%
  # Remove the environmental correlation method
  filter(dist_method %in% names(dist_method_replace)) %>%
  group_by(environment, trait, dist_method) %>% 
  mutate(scaled_distance = scale(distance),
         scaled_distance = scaled_distance - min(scaled_distance)) %>%
  ungroup() %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace),
         n_train_env = map_dbl(train_envs, length))


## Identify the largest prediction accuracy among the sliding windows
window_pred_toplot_maxima <- window_pred_toplot %>%
  group_by(environment, trait, dist_method) %>% 
  mutate(local_max_accuracy = max(accuracy),
         local_max_scaled_dist = scaled_distance[accuracy == local_max_accuracy]) %>%
  ungroup() %>%
  select(-train_envs)

## Filter for the local accuracy maxima
window_pred_toplot_maxima_filter <- window_pred_toplot_maxima %>% 
  filter(accuracy == local_max_accuracy)
  


## Distribution of the number of training environments or the scaled distance (cumulative) 
## of the local maximum
g_local_max <- window_pred_toplot_maxima_filter %>% 
  ggplot(aes(x = scaled_distance, fill = dist_method)) +
  geom_density(alpha = 0.25) + 
  facet_wrap(~ trait, ncol = 2) +
  theme_bw() +
  theme(legend.key.height = unit(2, "lines"))



## Calculate the average environmental distance (scaled) where the local maxima is
## found
window_pred_maxima_summary <- window_pred_toplot_maxima_filter %>% 
  group_by(trait, dist_method) %>% 
  summarize_at(vars(scaled_distance), mean)




## Common plot modifier
g_mod <- list(
  geom_point(aes(x = local_max_scaled_dist, y = local_max_accuracy)),
  geom_line(lwd = 0.5),
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)),
  ylab("Prediction Accuracy"),
  xlab("Mean Distance of Environments in the Training Set"),
  theme_bw(),
  theme(panel.grid = element_blank(),
        legend.key.height = unit(2, units = "line"),
        text = element_text(size = 8)),
  labs(title = "Sliding Window Environmental Clusters")
)



# List over traits
g_window_pred <- setNames(traits, traits) %>%
  map(function(tr) {
    # Subset the trait
    window_pred_toplot_tr <- window_pred_toplot_maxima %>%
      filter(trait == tr) %>%
      mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
    
    # Create one plot using facet_wrap
    g_total <- window_pred_toplot_tr %>% 
      ggplot(aes(x = scaled_distance, y = accuracy, color = dist_method)) +
      g_mod +
      facet_wrap( ~ environment + trait) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = dist_method, color = NULL), alpha = 0.1)
    
    # Save
    save_file <- file.path(fig_dir, str_c("window_env_dist_pred_", tr, ".jpg"))
    ggsave(filename = save_file, plot = g_total, height = 15, width = 12, dpi = 1000)
    
    # Return the plot
    # return(g_pages)
    return(g_total)
    
  })


### Plot a single distance method over many environments
# List over traits
g_window_pred1 <- dist_method_replace %>%
  map(function(dm) {
    
    # Subset the distance method
    window_pred_toplot_dm <- window_pred_toplot_maxima %>%
      filter(dist_method == dm)
    
    # Iterate over the trait
    g_list <- setNames(traits, traits) %>%
      map(function(tr) {
        # Subset the trait data
        window_pred_toplot_tr <- window_pred_toplot_dm %>%
          filter(trait == tr) %>%
          mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
        
        # Extract the ylim values from the single plots created above
        ylim_tr <- ggplot_build(g_window_pred[[tr]])$layout$panel_ranges[[1]]$y.range
        
        # Create the plot and return
        window_pred_toplot_tr %>% 
          ggplot(aes(x = scaled_distance, y = accuracy)) +
          g_mod + 
          ylim(ylim_tr) +
          facet_wrap( ~ environment + trait + dist_method)  +
          geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.1)
        
      })
    
    # Iterate over the plots and save
    for (i in seq_along(g_list)) {
      # Save
      save_file <- file.path(
        fig_dir, str_c("window_env_dist_pred_", names(dist_method_replace[dist_method_replace == dm]), 
                       "_", names(g_list)[i], ".jpg"))
      ggsave(filename = save_file, plot = g_list[[i]], height = 15, width = 12, dpi = 1000)
    }
    
    # Return the plot
    return(g_list)
    
  })




### Heritability across different ranks of environments

## Cumulative
# Bind rows and unnest
cumulative_heritability <- env_dist_heritability_out %>%
  bind_rows() %>% 
  unnest()

# Scale the distance measurements by environment and distance method
# Then subtract the minimum to have the values start at 0
cumulative_heritability_toplot <- cumulative_heritability %>%
  group_by(environment, trait, dist_method) %>%
  mutate(scaled_distance = scale(distance),
         scaled_distance = scaled_distance - min(scaled_distance)) %>%
  ungroup() %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace))


# ## Sample the same environments as those in the cumulative example above
# cumulative_heritability_toplot_sample <- cumulative_heritability_toplot %>%
#   filter(environment %in% unique(cumulative_pred_toplot_sample$environment))


# Plot
g_mod <- list(
  geom_point(),
  geom_line(),
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)), 
  ylab("Heritability"),
  xlab("Mean Distance of Environments in the Training Set"),
  ylim(c(0,1)),
  theme_bw(),
  theme(panel.grid = element_blank(),
        legend.key.height = unit(2, units = "line")),
  labs(title = "Cumulative Environmental Clusters")
)


# List over traits
g_cumulative_herit <- traits %>%
  map(function(tr) {
    # Subset the trait
    cumulative_heritability_toplot_tr <- cumulative_heritability_toplot %>%
      filter(trait == tr)
    
    # Calculate the number of pages
    ngrp <- n_distinct(cumulative_heritability_toplot_tr$environment)
    npgs <- ngrp / nenv_pg
    
    # Round up
    npgs_round <- ifelse(sign(round(npgs) - npgs) == 1, round(npgs), round(npgs + 1))
    
    # Iterate over the number of pages
    g_pages <- map(seq(npgs_round), function(pg) {
      cumulative_heritability_toplot_tr %>% 
        ggplot(aes(x = scaled_distance, y = heritability_out, color = dist_method)) +
        g_mod +
        facet_grid_paginate(environment + trait ~ ., scales = "free_y", ncol = 1, 
                            nrow = nenv_pg, page = pg)
    })
    
    # Iterate over pages and save
    for (i in seq_along(g_pages)) {
      # Save
      save_file <- file.path(fig_dir, str_c("cumulative_env_dist_herit_", tr, "_page", i, ".jpg"))
      ggsave(filename = save_file, plot = g_pages[[i]], height = 10, width = 8, dpi = 1000)
    }
    
    # Return the plot
    return(g_pages)
    
  })




## Sliding Window
# Bind rows and unnest
window_heritability <- env_dist_window_heritability_out %>%
  bind_rows() %>% 
  unnest()

# Scale the distance measurements by environment and distance method
# Then subtract the minimum to have the values start at 0
window_heritability_toplot <- window_heritability %>%
  group_by(environment, trait, dist_method) %>%
  mutate(scaled_distance = scale(distance),
         scaled_distance = scaled_distance - min(scaled_distance)) %>%
  ungroup() %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace))

# ## Sample the same environments as those in the cumulative example above
# window_heritability_toplot_sample <- window_heritability_toplot %>%
#   filter(environment %in% unique(cumulative_pred_toplot_sample$environment))

# Plot
g_mod <- list(
  geom_point(),
  geom_line(),
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)), 
  ylab("Heritability"),
  xlab("Mean Distance of Environments in the Training Set"),
  ylim(c(0,1)),
  theme_bw(),
  theme(panel.grid = element_blank(),
        legend.key.height = unit(2, units = "line")),
  labs(title = "Window Environmental Clusters")
)


# List over traits
g_window_herit <- traits %>%
  map(function(tr) {
    # Subset the trait
    window_heritability_toplot_tr <- window_heritability_toplot %>%
      filter(trait == tr)
    
    # Calculate the number of pages
    ngrp <- n_distinct(window_heritability_toplot_tr$environment)
    npgs <- ngrp / nenv_pg
    
    # Round up
    npgs_round <- ifelse(sign(round(npgs) - npgs) == 1, round(npgs), round(npgs + 1))
    
    # Iterate over the number of pages
    g_pages <- map(seq(npgs_round), function(pg) {
      window_heritability_toplot_tr %>% 
        ggplot(aes(x = scaled_distance, y = heritability_out, color = dist_method)) +
        g_mod +
        facet_grid_paginate(environment + trait ~ ., scales = "free_y", ncol = 1, 
                            nrow = nenv_pg, page = pg)
    })
    
    # Iterate over pages and save
    for (i in seq_along(g_pages)) {
      # Save
      save_file <- file.path(fig_dir, str_c("window_env_dist_herit_", tr, "_page", i, ".jpg"))
      ggsave(filename = save_file, plot = g_pages[[i]], height = 10, width = 8, dpi = 1000)
    }
    
    # Return the plot
    return(g_pages)
    
  })








