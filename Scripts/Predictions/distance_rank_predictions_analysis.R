## Analysis of predictions based on environmental distance
## 
## Author: Jeff Neyhart
## Last Updated: April 9, 2018
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

# Number of environments per page
nenv_pg <- 8



## Bind the list elements together and unnest
cumulative_pred_results <- env_dist_predictions_out %>% 
  bind_rows() %>% 
  unnest()


# Re-scale the distance measurements to unit variance
cumulative_pred_toplot <- cumulative_pred_results %>%
  group_by(environment, trait, dist_method) %>% 
  mutate(scaled_distance = scale(distance),
         scaled_distance = scaled_distance - min(scaled_distance)) %>%
  ungroup() %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace))
  

# ## Sample environments that have all traits
# cumulative_pred_toplot_sample <- cumulative_pred_toplot %>% 
#   group_by(environment) %>% filter(n_distinct(trait) == 3) %>% ungroup() %>%
#   filter(environment %in% sample(unique(.$environment), 5))

## Common plot modifier
g_mod <- list(
  geom_point(),
  geom_line(),
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)),
  ylab("Prediction Accuracy"),
  xlab("Cumulative Distance of Environments in Training Set"),
  theme_bw(),
  theme(panel.grid = element_blank(),
        legend.key.height = unit(2, units = "line")),
  labs(title = "Cumulative Environmental Clusters")
)

# List over traits
g_cumulative_pred <- traits %>%
  map(function(tr) {
    # Subset the trait
    cumulative_pred_toplot_tr <- cumulative_pred_toplot %>%
      filter(trait == tr) %>%
      mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
    
    # Calculate the number of pages
    ngrp <- n_distinct(cumulative_pred_toplot_tr$environment)
    npgs <- ngrp / nenv_pg
    
    # Round up
    npgs_round <- ifelse(sign(round(npgs) - npgs) == 1, round(npgs), round(npgs + 1))
    
    # Iterate over the number of pages
    g_pages <- map(seq(npgs_round), function(pg) {
      cumulative_pred_toplot_tr %>% 
        ggplot(aes(x = scaled_distance, y = accuracy, color = dist_method)) +
        g_mod +
        facet_grid_paginate(environment + trait ~ ., scales = "free_y", ncol = 1, 
                            nrow = nenv_pg, page = pg)
    })
    
    # Iterate over pages and save
    for (i in seq_along(g_pages)) {
      # Save
      save_file <- file.path(fig_dir, str_c("cumulative_env_dist_pred_", tr, "_page", i, ".jpg"))
      ggsave(filename = save_file, plot = g_pages[[i]], height = 10, width = 8, dpi = 1000)
    }
    
    # Return the plot
    return(g_pages)
    
  })


### Plot a single distance method over many environments
# List over traits
g_cumulative_pred1 <- dist_method_replace %>%
  map(function(dm) {
    
    # Subset the trait and the distance method
    g_list <- cumulative_pred_toplot %>% 
      filter(dist_method == dm) %>% 
      split(.$trait) %>% 
      list(., env_herit_rank) %>% 
      pmap(~mutate(.x, environment = factor(environment, levels = levels(.y$environment)))) %>%
      map(~ggplot(., aes(x = scaled_distance, y = accuracy, color = environment)) +
            g_mod +
            facet_grid(trait + dist_method ~ .) )
      
    # Iterate over the plots and save
    for (i in seq_along(g_list)) {
      # Save
      save_file <- file.path(
        fig_dir, str_c("cumulative_env_dist_pred_", names(dist_method_replace[dist_method_replace == dm]), 
                       "_", names(g_list)[i], ".jpg"))
      ggsave(filename = save_file, plot = g_list[[i]], height = 8, width = 8, dpi = 1000)
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
  group_by(environment, trait, dist_method) %>% 
  mutate(scaled_distance = scale(distance),
         scaled_distance = scaled_distance - min(scaled_distance)) %>%
  ungroup() %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace))


# ## Sample the same environments as those in the cumulative example above
# window_pred_toplot_sample <- window_pred_toplot %>%
#   filter(environment %in% unique(cumulative_pred_toplot_sample$environment))


# Plot
g_mod <- list(
  geom_point(),
  geom_line(),
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)), 
  ylab("Prediction Accuracy"),
  xlab("Mean Distance of Environments in the Training Set"),
  theme_bw(),
  theme(panel.grid = element_blank(),
        legend.key.height = unit(2, units = "line")),
  labs(title = "Sliding Window Environmental Clusters")
)


# List over traits
g_window_pred <- traits %>%
  map(function(tr) {
    # Subset the trait
    window_pred_toplot_tr <- window_pred_toplot %>%
      filter(trait == tr)
    
    # Calculate the number of pages
    ngrp <- n_distinct(window_pred_toplot_tr$environment)
    npgs <- ngrp / nenv_pg
    
    # Round up
    npgs_round <- ifelse(sign(round(npgs) - npgs) == 1, round(npgs), round(npgs + 1))
    
    # Iterate over the number of pages
    g_pages <- map(seq(npgs_round), function(pg) {
      window_pred_toplot_tr %>% 
        ggplot(aes(x = scaled_distance, y = accuracy, color = dist_method)) +
        g_mod +
        facet_grid_paginate(environment + trait ~ ., scales = "free_y", ncol = 1, 
                            nrow = nenv_pg, page = pg)
    })
    
    # Iterate over pages and save
    for (i in seq_along(g_pages)) {
      # Save
      save_file <- file.path(fig_dir, str_c("window_env_dist_pred_", tr, "_page", i, ".jpg"))
      ggsave(filename = save_file, plot = g_pages[[i]], height = 10, width = 8, dpi = 1000)
    }
    
    # Return the plot
    return(g_pages)
    
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








