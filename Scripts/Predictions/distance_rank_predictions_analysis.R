## Analysis of predictions based on environmental distance
## 
## Author: Jeff Neyhart
## Last Updated: April 23, 2018
## 
## This script will look at prediction accuracies from adding environments after
## ranking based on specific distance metrics
## 


# # Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# # Load some packages
library(lubridate)
library(ggforce)
library(ggridges)

# Load the results
# load(file.path(result_dir, "environmental_distance_predictions.RData"))
# load(file.path(result_dir, "environmental_distance_window_predictions.RData"))
# load(file.path(result_dir, "environmental_distance_heritability.RData"))

load(file.path(result_dir, "environmental_distance_predictions_TESTING.RData"))


## Bind the list elements together and unnest
cumulative_pred_results <- env_dist_predictions_out %>% 
  bind_rows() %>% 
  unnest() %>%
  rename(dist_method = model)

## Create a color scheme for the distance methods
dist_method_unique <- cumulative_pred_results$dist_method %>%
  unique() %>% 
  .[!str_detect(., "sample")] 

dist_method_replace <- dist_method_unique %>%
  str_replace_all("_", " ") %>% 
  str_to_title() %>% 
  str_replace_all(" ", "") %>%
  set_names(., dist_method_unique)
  

dist_colors <- c(setNames(umn_palette(3, length(dist_method_replace)), dist_method_replace), "Random" = "grey75")

## Significant level
alpha <- 0.05


# Re-scale the distance measurements to unit variance
cumulative_pred_adj <- cumulative_pred_results %>%
  group_by(environment, trait, dist_method) %>% 
  mutate(scaled_distance = scale(distance),
         scaled_distance = scaled_distance - min(scaled_distance)) %>%
  ungroup() %>%
  mutate(n_train_env = map_dbl(train_envs, length))

### First, for each distance method, find the rank of each environment
train_env_rank <- cumulative_pred_adj %>% 
  select(environment, trait, dist_method, n_train_env, train_envs) %>%
  group_by(environment, trait, dist_method) %>% 
  filter(n_train_env == max(n_train_env)) %>%
  unnest(train_envs) %>%
  # Add the rank
  mutate(train_env_rank = seq(n())) %>%
  ungroup()






# ## Bind the elements of the LRT results list
# env_dist_lrt_results <- env_dist_lrt_predictions_out %>%
#   bind_rows() %>%
#   select(environment:dist_method, results_out) %>%
#   unnest()
# 
# ## Select when it is first significant ("first")
# env_dist_lrt_results_first <- env_dist_lrt_results %>% 
#   group_by(environment, trait, dist_method) %>% 
#   mutate(is_sig = p_value <= alpha | n_env == max(n_env)) %>% 
#   filter(is_sig) %>%
#   mutate(which_sig = which(is_sig)) %>% 
#   filter(which_sig == min(which_sig)) %>%
#   ungroup()
# 





## For each environment and trait, order the environments according to the 
## prediction accuracy obtained after using all the data
## 
cumulative_env_pred_order <- cumulative_pred_adj %>% 
  filter(dist_method %in% names(dist_method_replace)) %>% 
  group_by(environment, trait) %>% 
  filter(n_train_env == max(n_train_env)) %>%
  slice(1) %>% 
  select(trait, environment, accuracy) %>% 
  ungroup() %>% split(.$trait) %>% 
  map(~mutate(., environment = factor(environment, levels = environment[order(accuracy, decreasing = TRUE)])))





## Separate the original samples from the random samples
## Then summarize the mean and quantiles of the prediction accuracies
## for each additional training population environment
cumulative_pred_random <- cumulative_pred_adj %>% 
  filter(str_detect(dist_method, "sample")) %>%
  group_by(environment, trait, n_train_env) %>% 
  summarize_at(vars(accuracy), funs(mean = mean(.), lower = quantile(., probs = 0.025),
                                    upper = quantile(., probs = 0.975))) %>%
  ungroup() %>%
  rename(accuracy = mean) %>%
  mutate(dist_method = "Random")

# Now grab the original runs
cumulative_pred_orig <- cumulative_pred_adj %>% 
  filter(dist_method %in% names(dist_method_replace)) %>%
  select(environment, trait, dist_method, scaled_distance, n_train_env, accuracy)

## Combine and add heritability information
# Remove environment-trait combinations with < h2 = 0.10
cumulative_pred_toplot <- bind_rows(cumulative_pred_orig, cumulative_pred_random) %>%
  mutate(n_train_env = as.integer(n_train_env)) %>%
  left_join(., select(stage_one_data, environment, trait, heritability)) %>%
  filter(heritability >= 0.10) %>%
  mutate(pred_ability = accuracy,
         accuracy = accuracy / sqrt(heritability),
         acc_lower = lower / sqrt(heritability),
         acc_upper = upper / sqrt(heritability))

# ## Alternatively, center by the minimum and scale by the heritability
# # This should center at 0
# cumulative_pred_toplot <- bind_rows(cumulative_pred_orig, cumulative_pred_random) %>%
#   mutate(n_train_env = as.integer(n_train_env)) %>%
#   left_join(., select(stage_one_data, environment, trait, heritability)) %>%
#   filter(heritability >= 0.10) %>%
#   group_by(trait, environment, dist_method) %>%
#   mutate(pred_ability = accuracy,
#          accuracy = scale(accuracy, center = min(accuracy),  scale = sqrt(unique(heritability))),
#          lower = scale(lower, center = min(accuracy), scale = sqrt(unique(heritability))),
#          upper = scale(upper, center = min(accuracy), scale = sqrt(unique(heritability)))) %>%
#   ungroup()




## Find the local maximum and create a data.frame
cumulative_pred_toplot_maximum <- cumulative_pred_toplot %>%  
  group_by(environment, trait, dist_method) %>% 
  filter(accuracy == max(accuracy)) %>%
  ungroup() %>%
  mutate(annotation = "local_maximum")

Ntrain_envs <- c(1, 5, 10)

## Find the accuray at n = 1, 5, 10 training environments
cumulative_pred_toplot_Nenv <- cumulative_pred_toplot %>%  
  filter(n_train_env %in% Ntrain_envs) %>%
  mutate(annotation = str_c(n_train_env, "envs"))





# ## Use the LRT results to select the accuracy when the LRT is deemed significant.
# ## First use the random results to identify the points along the 'add-one-environment' line
# 
# 
# ## Separate the original samples from the random samples
# cumulative_pred_toplot_lrt_random <- env_dist_lrt_results_first %>% 
#   filter(str_detect(dist_method, "sample")) %>%
#   # filter(dist_method %in% names(dist_method_replace)) %>%
#   left_join(., filter(cumulative_pred_toplot, dist_method == "Random"),
#             by = c("environment", "trait", "n_env" = "n_train_env")) %>%
#   group_by(environment, trait) %>% 
#   summarize_at(vars(n_env, scaled_distance, accuracy), mean) %>% 
#   mutate(dist_method = "Random", annotation = "lrt") %>%
#   rename(n_train_env = n_env)
#   
# # Original
# cumulative_pred_toplot_lrt_orig <- cumulative_pred_toplot %>%
#   inner_join(., env_dist_lrt_results_first, by = c("environment", "trait", "dist_method", "n_train_env" = "n_env")) %>%
#   select(environment:upper) %>% 
#   mutate(annotation = "lrt")
# 
# # Combine
# cumulative_pred_toplot_lrt <- bind_rows(cumulative_pred_toplot_lrt_orig, cumulative_pred_toplot_lrt_random)
#   

# 
# 
# ## A vector to replace the annotation names
# annotation_replace <- c(
#   "local_maximum" = "Maximum", 
#   setNames(str_replace(string = str_c(Ntrain_envs, "envs"), pattern = "e", replacement = " E"), 
#            str_c(Ntrain_envs, "envs")), 
#   "lrt" = "Likelihood Ratio Test")
# 
# 
# ## Combine the data together
# cumulative_pred_toplot_annotate <- bind_rows(cumulative_pred_toplot_maximum, 
#                                              cumulative_pred_toplot_Nenv,
#                                              cumulative_pred_toplot_lrt) %>%
#   # Calculate the difference between the annotated accuracy and the "final" accuracy
#   mutate(annotation = as_replaced_factor(x = annotation, replacement = annotation_replace),
#          dist_method = as_replaced_factor(x = dist_method, c(dist_method_replace, "Random" = "Random")))
# 
# 
# # Change the distance methods to a factor
# cumulative_pred_toplot_final <- cumulative_pred_toplot %>%
#   mutate(dist_method = as_replaced_factor(x = dist_method, replacement = c(dist_method_replace, "Random" = "Random")),
#          # Remove the upper/lower for non random dist_methods
#          lower = if_else(dist_method == "Random", lower, as.numeric(NA)),
#          upper = if_else(dist_method == "Random", upper, as.numeric(NA)))
# 





##### Plotting ######






## Common plot modifier
g_mod <- list(
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = dist_method, color = NULL), alpha = 0.25),
  # geom_point(aes(x = local_max_scaled_dist, y = local_max_accuracy)),
  geom_line(lwd = 0.5),
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)),
  ylab("Prediction Accuracy"),
  xlab("Number of Environments in Training Set"),
  scale_color_manual(values = dist_colors, name = NULL),
  scale_fill_manual(values = dist_colors, name = NULL),
  theme_acs(),
  theme(legend.key.height = unit(2, units = "line"))
)


## Testing
g_cumupred_test <- cumulative_pred_toplot %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace)) %>%
  ggplot(aes(x = n_train_env, y = accuracy, color = dist_method)) + 
  g_mod +
  facet_grid(trait ~ environment, scale = "free_y")

ggsave(filename = "cumulative_env_dist_pred_testing.jpg", plot = g_cumupred_test, path = fig_dir,
       height = 6, width = 10, dpi = 1000)




# # List over traits
# g_cumulative_pred <- setNames(traits, traits) %>%
#   map(function(tr) {
#     # Subset the trait
#     cumulative_pred_toplot_tr <- cumulative_pred_toplot_final %>%
#       filter(trait == tr) %>%
#       mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
#     
#     cumulative_pred_toplot_annotate_tr <- cumulative_pred_toplot_annotate %>%
#       filter(trait == tr) %>%
#       mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
#       
#     # Create one plot using facet_wrap
#     g_total <- cumulative_pred_toplot_tr %>% 
#       # ggplot(aes(x = scaled_distance, y = accuracy, color = dist_method)) +
#       ggplot(aes(x = n_train_env, y = accuracy, color = dist_method)) +
#       geom_point(data = cumulative_pred_toplot_annotate_tr, aes(shape = annotation)) +
#       g_mod +
#       facet_wrap( ~ environment + trait)
#       # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = dist_method, color = NULL), alpha = 0.1)
#     
#     # Save
#     save_file <- file.path(fig_dir, str_c("cumulative_env_dist_pred_", tr, ".jpg"))
#     ggsave(filename = save_file, plot = g_total, height = 15, width = 12, dpi = 1000)
#     
#     # Return the plot
#     # return(g_pages)
#     return(g_total)
#     
#   })
# 
# 
# cumumlative_blank <- list()
# 
# ### Plot a single distance method over many environments
# # List over traits
# g_cumulative_pred1 <- dist_method_replace %>%
#   map(function(dm) {
#     
#     # Subset the distance method
#     cumulative_pred_toplot_dm <- cumulative_pred_toplot_final %>%
#       filter(dist_method %in% c(dm, "Random"))
#     
#     cumulative_pred_toplot_annotate_dm <- cumulative_pred_toplot_annotate %>%
#       filter(dist_method %in% c(dm, "Random"))
#     
#     # Iterate over the trait
#     g_list <- setNames(traits, traits) %>%
#       map(function(tr) {
#         # Subset the trait data
#         cumulative_pred_toplot_tr <- cumulative_pred_toplot_dm %>%
#           filter(trait == tr) %>%
#           # mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
#           mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
#         
#         cumulative_pred_toplot_annotate_tr <- cumulative_pred_toplot_annotate_dm %>%
#           filter(trait == tr) %>%
#           mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
#         
#           
#         # Extract the ylim values from the single plots created above
#         ylim_tr <- ggplot_build(g_cumulative_pred[[tr]])$layout$panel_ranges[[1]]$y.range
#         
#         # Create the plot and return
#         cumulative_pred_toplot_tr %>% 
#           # ggplot(aes(x = scaled_distance, y = accuracy)) +
#           ggplot(aes(x = n_train_env, y = accuracy, color = dist_method)) +
#           geom_point(data = cumulative_pred_toplot_annotate_tr, aes(shape = annotation)) +
#           g_mod + 
#           ylim(ylim_tr) +
#           facet_wrap( ~ environment + trait)
#           # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.1)
#         
#       })
#       
#     # Iterate over the plots and save
#     for (i in seq_along(g_list)) {
#       # Save
#       save_file <- file.path(
#         fig_dir, str_c("cumulative_env_dist_pred_", names(dist_method_replace[dist_method_replace == dm]), 
#                        "_", names(g_list)[i], ".jpg"))
#       ggsave(filename = save_file, plot = g_list[[i]], height = 15, width = 12, dpi = 1000)
#     }
#     
#     # Return the plot
#     return(g_list)
#     
#   })





## Prediction accuracy using a sliding window
## Bind the list elements together and unnest
window_pred_results <- env_dist_window_predictions_out %>% 
  unnest() %>% unnest() %>%
  rename(dist_method = model)

## Create a color scheme for the distance methods
dist_method_unique <- window_pred_results$dist_method %>%
  unique() %>% 
  .[!str_detect(., "sample")] 

dist_method_replace <- dist_method_unique %>%
  str_replace_all("_", " ") %>% 
  str_to_title() %>% 
  str_replace_all(" ", "") %>%
  set_names(., dist_method_unique)

dist_colors <- c(setNames(umn_palette(3, length(dist_method_replace)), dist_method_replace), "Random" = "grey75")


# Re-scale the distance measurements to unit variance
window_pred_adj <- window_pred_results %>%
  group_by(window_size, environment, trait, dist_method) %>% 
  mutate(scaled_distance = scale(distance),
         scaled_distance = scaled_distance - min(scaled_distance),
         ## Add the rank of the training set
         train_rank = seq(n())) %>%
  ungroup()


## For each environment and trait, order the environments according to the 
## prediction accuracy obtained after using all data (i.e. use the same order as above)
## 
window_env_pred_order <- cumulative_env_pred_order





## Separate the original samples from the random samples
## Then summarize the mean and quantiles of the prediction accuracies
## for each additional training population environment
window_pred_random <- window_pred_adj %>% 
  filter(str_detect(dist_method, "sample")) %>%
  group_by(window_size, environment, trait, train_rank) %>% 
  summarize_at(vars(accuracy), funs(mean = mean(.), lower = quantile(., probs = 0.025),
                                    upper = quantile(., probs = 0.975))) %>%
  ungroup() %>%
  rename(accuracy = mean) %>%
  mutate(dist_method = "Random")

# Now grab the original runs
window_pred_orig <- window_pred_adj %>% 
  filter(dist_method %in% names(dist_method_replace)) %>%
  select(window_size, environment, trait, dist_method, scaled_distance, train_rank, accuracy)

## Combine and add heritability information
# Remove environment-trait combinations with < h2 = 0.10
window_pred_toplot <- bind_rows(window_pred_orig, window_pred_random) %>%
  mutate(n_train_env = as.integer(train_rank),
         window_size = as.factor(window_size)) %>%
  left_join(., select(stage_one_data, environment, trait, heritability)) %>%
  filter(heritability >= 0.10) %>%
  mutate(pred_ability = accuracy,
         accuracy = accuracy / sqrt(heritability),
         acc_lower = lower / sqrt(heritability),
         acc_upper = upper / sqrt(heritability))




# ## Identify the largest prediction accuracy among the sliding windows
# window_pred_toplot_maxima <- window_pred_toplot %>%
#   group_by(environment, trait, dist_method, n_train_env) %>% 
#   mutate(local_max_accuracy = max(accuracy),
#          local_max_window_rank = window_rank[accuracy == local_max_accuracy]) %>%
#   ungroup() %>%
#   mutate(dist_method = as_replaced_factor(dist_method, c(dist_method_replace, "Random" = "Random")),
#          window_rank = as.integer(window_rank))
# 
# ## Filter for the local accuracy maxima
# window_pred_toplot_maxima_filter <- window_pred_toplot_maxima %>% 
#   filter(accuracy == local_max_accuracy)
  



# cumumlative_blank <- cumulative_pred_toplot %>% 
#   mutate(annotation = "cumulative") %>%
#   filter(dist_method != "Random") %>%
#   split(.$trait) %>%
#   map(~mutate(., environment = factor(environment, levels(cumulative_env_pred_order[[unique(.$trait)]]$environment)))) %>% 
#   map(~select(., annotation, trait, dist_method, environment) %>% 
#         distinct() %>% 
#         arrange(trait, dist_method, environment))
# 
# window_blank <- window_pred_toplot %>%
#   mutate(annotation = "window") %>%
#   filter(dist_method != "Random") %>%
#   split(.$trait) %>%
#   map(~mutate(., environment = factor(environment, levels(cumulative_env_pred_order[[unique(.$trait)]]$environment)))) %>% 
#   map(~select(., annotation, trait, dist_method, environment) %>% 
#         distinct() %>% 
#         arrange(trait, dist_method, environment))






# ## Distribution of the number of training environments or the scaled distance (cumulative) 
# ## of the local maximum
# g_local_max <- window_pred_toplot_maxima_filter %>% 
#   # filter(dist_method != "Random") %>%
#   ggplot(aes(x = scaled_distance, fill = dist_method)) +
#   # ggplot(aes(x = window_rank, fill = dist_method)) +
#   # geom_density(alpha = 0.25) + 
#   geom_density_ridges(aes(y = dist_method), alpha = 0.25) +
#   scale_fill_manual(values = dist_colors, guide = FALSE) +
#   xlab("Scaled Average Environmental Distance at Maximum Accuracy") +
#   labs(title = "Average Environmental Distance at Maximum Accuracy") +
#   facet_wrap(~ trait, ncol = 2) +
#   theme_bw() +
#   theme(legend.key.height = unit(2, "lines"))
# 
# save_file <- file.path(fig_dir, "window_pred_max_dist.jpg")
# ggsave(filename = save_file, plot = g_local_max, height = 6, width = 8, dpi = 1000)
# 
# 
# g_local_max <- window_pred_toplot_maxima_filter %>% 
#   # filter(dist_method != "Random") %>%
#   # ggplot(aes(x = scaled_distance, fill = dist_method)) +
#   ggplot(aes(x = window_rank, fill = dist_method)) +
#   # geom_density(alpha = 0.25) + 
#   geom_density_ridges(aes(y = dist_method), alpha = 0.25) +
#   scale_fill_manual(values = dist_colors, guide = FALSE) +
#   xlab("Rank of Sliding Window at Maximum Accuracy") +
#   labs(title = "Rank of Sliding Window at Maximum Accuracy") +
#   facet_wrap(~ trait, ncol = 2) +
#   theme_bw() +
#   theme(legend.key.height = unit(2, "lines"))
# 
# save_file <- file.path(fig_dir, "window_pred_max_rank.jpg")
# ggsave(filename = save_file, plot = g_local_max, height = 6, width = 8, dpi = 1000)
# 
# 
# 
# ## Calculate the average environmental distance (scaled) where the local maxima is
# ## found and calculate the average rank of the window where the maximum is found
# window_pred_maxima_summary <- window_pred_toplot_maxima_filter %>% 
#   group_by(trait, dist_method, n_train_env) %>% 
#   summarize_at(vars(window_rank, scaled_distance), mean)


##### Plotting ######






## Common plot modifier
g_mod <- list(
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = dist_method, color = NULL), alpha = 0.25),
  # geom_point(aes(x = local_max_scaled_dist, y = local_max_accuracy)),
  geom_line(lwd = 0.5),
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)),
  ylab("Prediction Accuracy"),
  xlab("Rank of training set (by distance)"),
  scale_color_manual(values = dist_colors, name = NULL),
  scale_fill_manual(values = dist_colors, name = NULL),
  theme_acs(),
  theme(legend.key.height = unit(2, units = "line"))
)


## Testing
g_window_test <- window_pred_toplot %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace)) %>%
  ggplot(aes(x = train_rank, y = pred_ability, color = dist_method, lty = window_size)) + 
  g_mod +
  facet_grid(trait ~ environment + window_size, scale = "free")

ggsave(filename = "window_env_dist_pred_testing.jpg", plot = g_window_test, path = fig_dir,
       height = 6, width = 12, dpi = 1000)



# ## Common plot modifier
# g_mod <- list(
#   geom_ribbon(aes(ymin = lower, ymax = upper, fill = dist_method, color = NULL), alpha = 0.25),
#   # geom_point(aes(x = local_max_scaled_dist, y = local_max_accuracy)),
#   geom_point(aes(x = local_max_window_rank, y = local_max_accuracy)),
#   geom_line(lwd = 0.5),
#   # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)),
#   scale_color_manual(values = dist_colors),
#   scale_fill_manual(values = dist_colors),
#   ylab("Prediction Accuracy"),
#   xlab("Distance Rank of Sliding Window"),
#   theme_bw(),
#   theme(panel.grid = element_blank(),
#         legend.key.height = unit(2, units = "line"),
#         text = element_text(size = 8)),
#   labs(title = "Sliding Window Environmental Clusters")
# )






# # List over traits
# g_window_pred <- setNames(traits, traits) %>%
#   map(function(tr) {
#     # Subset the trait
#     window_pred_toplot_tr <- window_pred_toplot_maxima %>%
#       filter(trait == tr) %>%
#       # mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
#       mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
#       
#     # Create one plot using facet_wrap
#     g_total <- window_pred_toplot_tr %>% 
#       # ggplot(aes(x = scaled_distance, y = accuracy, color = dist_method)) +
#       ggplot(aes(x = window_rank, y = accuracy, color = dist_method)) +
#       g_mod +
#       facet_wrap( ~ environment + trait)
#       # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = dist_method, color = NULL), alpha = 0.1)
#     
#     # Save
#     save_file <- file.path(fig_dir, str_c("window_env_dist_pred_", tr, ".jpg"))
#     ggsave(filename = save_file, plot = g_total, height = 15, width = 12, dpi = 1000)
#     
#     # Return the plot
#     # return(g_pages)
#     return(g_total)
#     
#   })
# 
# 
# ### Plot a single distance method over many environments
# # List over traits
# g_window_pred1 <- dist_method_replace %>%
#   map(function(dm) {
#     
#     # Subset the distance method
#     window_pred_toplot_dm <- window_pred_toplot_maxima %>%
#       filter(dist_method %in% c(dm, "Random"))
#     
#     # Iterate over the trait
#     g_list <- setNames(traits, traits) %>%
#       map(function(tr) {
#         # Subset the trait data
#         window_pred_toplot_tr <- window_pred_toplot_dm %>%
#           filter(trait == tr) %>%
#           # mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
#           mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
#           
#         # Extract the ylim values from the single plots created above
#         ylim_tr <- ggplot_build(g_window_pred[[tr]])$layout$panel_ranges[[1]]$y.range
#         
#         # Create the plot and return
#         window_pred_toplot_tr %>% 
#           # ggplot(aes(x = scaled_distance, y = accuracy)) +
#           ggplot(aes(x = window_rank, y = accuracy, color = dist_method)) +
#           g_mod + 
#           ylim(ylim_tr) +
#           facet_wrap( ~ environment + trait)
#           # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.1)
#         
#       })
#     
#     # Iterate over the plots and save
#     for (i in seq_along(g_list)) {
#       # Save
#       save_file <- file.path(
#         fig_dir, str_c("window_env_dist_pred_", names(dist_method_replace[dist_method_replace == dm]), 
#                        "_", names(g_list)[i], ".jpg"))
#       g_print <- g_list[[i]] + geom_smooth(method = "lm", se = FALSE)
#       
#       ggsave(filename = save_file, plot = , height = 15, width = 12, dpi = 1000)
#     }
#     
#     # Return the plot
#     return(g_list)
#     
#   })

















####### Statistical Testing #########

## For each prediction environment, trait, and distance method, find the prediction
## accuracy after adding all environments (i.e. the last value) and the local
## maximum prediction accuracy (if present)
## 
## This is to be used for analysis only
## 

cumulative_pred_analysis_terminal <- cumulative_pred_toplot %>% 
  group_by(environment, trait) %>% 
  filter(n_train_env == max(n_train_env)) %>% 
  slice(1) %>%
  ungroup() %>%
  select(-dist_method, -lower, -upper, -acc_lower, -acc_upper) %>%
  rename_at(vars(accuracy, pred_ability, n_train_env, scaled_distance), ~str_c(., "_final"))

# Now subset the dataframe for use as a base
cumulative_pred_analysis_final <- cumulative_pred_toplot %>% 
  # select(-scaled_distance, -lower:-upper) %>% 
  left_join(., cumulative_pred_analysis_terminal, by = c("environment", "trait", "heritability")) %>%
  mutate(accuracy_advantage = accuracy - accuracy_final,
         pred_ability_advantage = pred_ability - pred_ability_final)

## Calculate the average advantage for each distance method, then compare
## with random
cumulative_pred_analysis_final_summ <- cumulative_pred_analysis_final %>% 
  group_by(trait, dist_method, n_train_env) %>% 
  summarize_at(vars(contains("advantage")), funs(mean, adv_lower = quantile(., probs = 0.025), 
                                     adv_upper = quantile(., probs = 0.975))) %>%
  ungroup()

cumulative_pred_analysis_final_summ1 <- left_join(
  filter(cumulative_pred_analysis_final_summ, dist_method != "Random"),
  filter(cumulative_pred_analysis_final_summ, dist_method == "Random") %>% 
    select(-dist_method) %>% 
    rename_at(vars(-trait, -n_train_env), ~str_c(., "_random"))
)

## For each n_train_env, plot the average advantage across all environments using a
## boxplot
g_advantage <- cumulative_pred_analysis_final %>% 
  ggplot(aes(x = n_train_env, y = pred_ability_advantage, group = n_train_env)) +
  geom_hline(yintercept = 0) + 
  geom_boxplot() + 
  facet_grid(trait ~ dist_method) +
  theme_bw()


## Common plot modifier
g_mod <- list(
  geom_hline(yintercept = 0), 
  geom_ribbon(aes(ymin = pred_ability_advantage_adv_lower_random, ymax = pred_ability_advantage_adv_upper_random,
                  color = "Random", fill = "Random"), alpha = 0.5), 
  # geom_ribbon(aes(ymin = adv_lower, ymax = adv_upper), alpha = 0.5),
  geom_line(aes(y = pred_ability_advantage_mean_random, color = "Random"), lwd = 1),
  geom_line(lwd = 1), 
  scale_fill_manual(values = dist_colors, name = NULL), 
  scale_color_manual(values = dist_colors, name = NULL), 
  # facet_grid(trait ~ dist_method),
  theme_acs(),
  theme(legend.key.height = unit(x = 2, units = "lines"))
)


## Plot the mean advantage for each distance method, with a ribbon to illustrate
## the confidence interval of the random
g_mean_advantage <- cumulative_pred_analysis_final_summ1 %>% 
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace)) %>%
  ggplot(aes(x = n_train_env, y = pred_ability_advantage_mean, col = dist_method, fill = dist_method)) + 
  g_mod +
  facet_grid(trait ~ dist_method) +
  ylab("Prediction Advantage Over Using All Data") +
  xlab("Number of Training Environments")

## Save
ggsave(filename = "cumulative_mean_pred_ability_advantage_testing.jpg", plot = g_mean_advantage,
       path = fig_dir, height = 5, width = 8, dpi = 1000)








## Window results
window_pred_analysis_terminal <- window_pred_toplot %>% 
  group_by(window_size, environment, trait) %>% 
  filter(train_rank == max(train_rank)) %>% 
  slice(1) %>%
  ungroup() %>%
  select(-dist_method, -lower, -upper, -acc_lower, -acc_upper) %>%
  rename_at(vars(accuracy, pred_ability, train_rank, scaled_distance), ~str_c(., "_final"))

# Now subset the dataframe for use as a base
window_pred_analysis_final <- window_pred_toplot %>% 
  # select(-scaled_distance, -lower:-upper) %>% 
  left_join(., window_pred_analysis_terminal, by = c("window_size", "environment", "trait", "heritability")) %>%
  mutate(accuracy_advantage = accuracy - accuracy_final,
         pred_ability_advantage = pred_ability - pred_ability_final)

## Calculate the average advantage for each distance method, then compare
## with random
window_pred_analysis_final_summ <- window_pred_analysis_final %>% 
  group_by(window_size, trait, dist_method, train_rank) %>% 
  summarize_at(vars(contains("advantage")), funs(mean, adv_lower = quantile(., probs = 0.025), 
                                                 adv_upper = quantile(., probs = 0.975))) %>%
  ungroup()

window_pred_analysis_final_summ1 <- left_join(
  filter(window_pred_analysis_final_summ, dist_method != "Random"),
  filter(window_pred_analysis_final_summ, dist_method == "Random") %>% 
    select(-dist_method) %>% 
    rename_at(vars(-window_size, -trait, -train_rank), ~str_c(., "_random"))
)


## For each n_train_env, plot the average advantage across all environments using a
## boxplot
g_advantage_window <- window_pred_analysis_final %>% 
  ggplot(aes(x = train_rank, y = pred_ability_advantage, group = train_rank)) +
  geom_hline(yintercept = 0) + 
  geom_boxplot() + 
  facet_grid(trait ~ dist_method) +
  theme_bw()


## Plot the mean advantage for each distance method, with a ribbon to illustrate
## the confidence interval of the random
g_mean_advantage_window <- window_pred_analysis_final_summ1 %>% 
  filter(window_size == 5) %>% # Filter for this window size because it shows the most apparent differences
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace)) %>%
  ggplot(aes(x = train_rank, y = pred_ability_advantage_mean, col = dist_method, fill = dist_method)) + 
  g_mod +
  facet_grid(trait ~ dist_method + window_size) +
  ylab("Prediction Advantage Over Using All Data") +
  xlab("Distance Rank of Sliding Window")

## Save
ggsave(filename = "window_mean_pred_ability_advantage_testing.jpg", plot = g_mean_advantage_window,
       path = fig_dir, height = 5, width = 8, dpi = 1000)



## What is going on with the uptick in accuracy at the end of the sliding windows?
## These are the C1R-only environments, so they have one extra window (I trimmed these
## off)











## Plot
## Distribution of the number of training environments or the scaled distance (cumulative) 
## of the local maximum
g_local_max <- cumulative_pred_toplot_annotate %>% 
  filter(annotation %in% c("Maximum")) %>%
  ggplot(aes(x = n_train_env, fill = dist_method)) +
  # ggplot(aes(x = scaled_distance, fill = dist_method)) + 
  # geom_density(alpha = 0.25) + 
  geom_density_ridges(aes(y = dist_method), alpha = 0.25) +
  scale_fill_manual(values = dist_colors, guide = FALSE) +
  xlab("Number of Training Environments to Reach Maximum") +
  labs(title = "Training Environments to Reach Maximum Accuracy") +
  facet_wrap(~ trait, ncol = 2) +
  theme_bw() +
  theme(legend.key.height = unit(2, "lines"))

save_file <- file.path(fig_dir, "cumulative_pred_max_ntrain.jpg")
ggsave(filename = save_file, plot = g_local_max, height = 6, width = 8, dpi = 1000)



# ## Distribution of local advantage
# g_local_adv <- cumulative_pred_toplot_annotate %>% 
#   ggplot(aes(x = advantage, fill = dist_method)) + 
#   # geom_density(alpha = 0.25) + 
#   geom_density_ridges(aes(y = dist_method), alpha = 0.25) +
#   scale_fill_manual(values = dist_colors, guide = FALSE) +
#   xlab("Prediction Accuracy Advantage") +
#   labs(title = "Advantage of Selected Environments Over All Environments") +
#   facet_grid(annotation ~ trait) +
#   theme_bw() +
#   theme(legend.key.height = unit(2, "lines"))

## Distribution of local advantage
g_local_adv <- cumulative_pred_toplot_annotate %>% 
  ggplot(aes(x = dist_method, y = advantage, fill = annotation)) + 
  geom_boxplot(alpha = 0.5) + 
  # scale_fill_manual(values = dist_colors, guide = FALSE) +
  ylab("Accuracy Advantage Over Using All Data") +
  labs(title = "Advantage of Selected Environments Over All Environments") +
  facet_grid( ~ trait) +
  theme_bw() +
  theme(legend.key.height = unit(2, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

save_file <- file.path(fig_dir, "cumulative_pred_max_advantage.jpg")
ggsave(filename = save_file, plot = g_local_adv, height = 6, width = 8, dpi = 1000)








## Determine whether the local maximum, or the optimal accuracy based on subsets 
## of environments is significantly different from the terminal accuracy
cumulative_pred_significant <- cumulative_pred_toplot_annotate %>% 
  mutate(sig_better = accuracy > terminal_upper,
         sig_worse = accuracy < terminal_lower) %>%
  group_by(trait, dist_method, annotation) %>% 
  summarize_at(vars(starts_with("sig")), mean)


# Calculate the proportion of significant accuracies across all observations
cumulative_pred_toplot_final %>%   
  mutate(sig_better = accuracy > terminal_upper,
         sig_worse = accuracy < terminal_lower) %>% 
  group_by(trait, dist_method) %>% 
  summarize_at(vars(starts_with("sig")), mean)



## For the prediction accuracy determined using the sliding window, fit a linear
## model for each of the responses. Then determine the number of coefficients
## that are significantly greater than or less than zero
window_pred_toplot_fit <- window_pred_toplot_maxima %>% 
  filter(dist_method != "Random") %>%
  group_by(environment, trait, dist_method) %>% 
  do(fit = lm(accuracy ~ scaled_distance, data = .)) %>%
  ungroup() %>% 
  mutate(coef = map(fit, ~broom::tidy(.)),
         df = map_dbl(fit, df.residual))

# Extract the coefficients and perform hypothesis tests
window_pred_toplot_fit_coef <- window_pred_toplot_fit %>% 
  unnest(coef) %>% 
  filter(term != "(Intercept)") %>%
  mutate(p_value_gt0 = pt(q = statistic, df = df, lower.tail = FALSE), 
         p_value_lt0 = pt(q = statistic, df = df))

## Determine the proportion of p_values that are less than 0, greater than 0
window_pred_toplot_fit_sig <- window_pred_toplot_fit_coef %>% 
  mutate_at(vars(starts_with("p_value")), ~. <= alpha) %>%
  group_by(trait, dist_method) %>%
  summarize_at(vars(starts_with("p_value")), funs(mean))

# Plot
window_pred_toplot_fit_sig %>% 
  gather(test, prop_sig, starts_with("p_value")) %>% 
  ggplot(aes(x = dist_method, y = prop_sig, fill = test)) + 
  geom_col(position = "dodge") + 
  facet_grid(~trait)







#########

## The distribution of the ranks should follow a uniform distribution under the
## null hypothesis of idependence of rank and training environment. I will use a 
## KS test to test the observed rank distribution versus a uniform distribution
## 

train_env_rank_test <- train_env_rank %>% 
  filter(dist_method %in% names(dist_method_replace)) %>%
  group_by(trait, dist_method, train_envs) %>% 
  do({data_frame(
    mean_rank = mean(.$train_env_rank),
    sd_rank = sd(.$train_env_rank),
    # The less than alternative tests if the distribution is stochastically larger (higher rank)
    # The greater than alternative tests if the distribution is stochastically less (lower rank)
    ks_test_low = list(ks.test(x = .$train_env_rank, y = "punif", min = 1, max = nrow(.), alternative = "l")),
    ks_test_high = list(ks.test(x = .$train_env_rank, y = "punif", min = 1, max = nrow(.), alternative = "g"))
  )}) %>%
  # Extract the p_value and adjust
  group_by(trait, dist_method) %>% 
  mutate(p_value_lt = map_dbl(ks_test_low, "p.value"), 
         p_value_gt = map_dbl(ks_test_high, "p.value"),
         n_test = n(),
         p_adj_lt = p.adjust(p = p_value_lt, method = "bonf"),
         p_adj_gt = p.adjust(p = p_value_gt, method = "bonf")) %>%
  ungroup()

## Examine the significance results
# Get the training environments with the lowest and highest rank
train_env_rank_test_sig <- train_env_rank_test %>% 
  select(trait:sd_rank, contains("p_adj")) %>%
  gather(test_type, p_adj, contains("p_adj")) %>%
  filter(p_adj <= 0.05) %>% 
  group_by(trait, dist_method, test_type) %>% 
  filter(mean_rank == min(mean_rank) | mean_rank == max(mean_rank))

## Count the number of times a training environment is significant
train_env_rank_test_count <- train_env_rank_test %>% 
  select(trait:sd_rank, contains("p_adj")) %>%
  gather(test_type, p_adj, contains("p_adj"))  %>% 
  group_by(trait, train_envs, test_type) %>% 
  summarize(n_times_sig = sum(p_adj <= 0.05)) %>%
  arrange(desc(n_times_sig))

## Add heritability information for comparison
train_env_rank_test %>% 
  select(trait, dist_method, train_envs, mean_rank) %>% 
  left_join(., bind_rows(env_herit_rank), by = c("trait", "train_envs" = "environment")) %>% 
  group_by(trait, dist_method) %>% 
  summarize(rank_herit_cor = cor(mean_rank, heritability))

## Heritability only seems to be correlated with the D phenotypic distance metric,
## which makes sense.
# Plot this - actually evidence is not strong enough to suggest a trend between
# heritability and the rank
train_env_rank_test %>% 
  select(trait, dist_method, train_envs, mean_rank) %>% 
  left_join(., bind_rows(env_herit_rank), by = c("trait", "train_envs" = "environment")) %>% 
  ggplot(aes(x = mean_rank, heritability)) +
  geom_point() + 
  geom_smooth(method = "lm", se = F) + 
  facet_grid(trait ~ dist_method, scales = "free_y")


## For each prediction environment, take the mean rank of environments from the same
## year and the mean rank of environments not from the same year
train_env_rank_year <- train_env_rank %>%
  filter(dist_method %in% names(dist_method_replace)) %>%
  mutate_at(vars(environment, train_envs), funs(year = str_extract(string = ., pattern = "[0-9]{2}") %>%
                                                  parse_date_time(orders = "y") %>% year())) %>%
  mutate(same_year = environment_year == train_envs_year) %>% 
  group_by(trait, dist_method, same_year, environment) %>% 
  summarize_at(vars(train_env_rank), funs(mean, sd)) %>% 
  mutate(overall_mean = mean(mean)) %>%
  ungroup()

## Plot
train_env_rank_year %>% 
  ggplot(aes(x = dist_method, y = mean, fill = same_year)) + 
  geom_boxplot(position = "dodge") + 
  facet_grid(~trait) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






### Fit a model to compare the effects of the different factors
cumulative_pred_models <- cumulative_pred_toplot_annotate %>% 
  mutate(environment = as.factor(environment)) %>% 
  group_by(trait) %>%
  do(fit = lm(advantage ~ environment + dist_method + annotation, data = .))

# dist_method:annotation is not signficant

cumulative_pred_effects <- cumulative_pred_models %>% 
  mutate(effects = list(effects::allEffects(fit)))

cumulative_pred_effects1 <- bind_cols(
  cumulative_pred_effects, 
  cumulative_pred_effects$effects %>% map(~map(., as.data.frame)) %>% transpose() %>% as_data_frame())

## Plot
cumulative_pred_effects1 %>% 
  unnest(dist_method) %>% 
  ggplot(aes(x = dist_method, y = fit, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_errorbar() +
  ylab("effect") + 
  facet_wrap(~trait) + 
  theme_bw()

cumulative_pred_effects1 %>% 
  unnest(annotation) %>% 
  ggplot(aes(x = annotation, y = fit, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_errorbar() +
  ylab("effect") + 
  facet_wrap(~trait) + 
  theme_bw()



















##### Heritability across different ranks of environments ###### 

## Cumulative
# Bind rows and unnest
cumulative_heritability <- env_dist_heritability_out %>%
  bind_rows() %>% 
  unnest()

# Re-scale the distance measurements to unit variance
cumulative_heritability_adj <- cumulative_heritability %>%
  group_by(environment, trait, dist_method) %>% 
  mutate(scaled_distance = scale(distance),
         scaled_distance = scaled_distance - min(scaled_distance)) %>%
  ungroup() %>%
  mutate(n_train_env = map_dbl(train_envs, length))


## Separate the original samples from the random samples
## Then summarize the mean and quantiles of the prediction accuracies
## for each additional training population environment
cumulative_heritability_random <- cumulative_heritability_adj %>% 
  filter(str_detect(dist_method, "sample")) %>%
  group_by(environment, trait, n_train_env) %>% 
  summarize_at(vars(heritability_out), funs(mean = mean(.), lower = quantile(., probs = 0.025),
                                            upper = quantile(., probs = 0.975))) %>%
  ungroup() %>%
  rename(heritability_out = mean) %>%
  mutate(dist_method = "Random")

# Now grab the original runs
cumulative_heritability_orig <- cumulative_heritability_adj %>% 
  filter(dist_method %in% names(dist_method_replace)) %>%
  select(environment, trait, dist_method, scaled_distance, n_train_env, heritability_out)

## Combine
cumulative_heritability_toplot <- bind_rows(cumulative_heritability_orig, cumulative_heritability_random) %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace),
         dist_method = factor(dist_method, levels = names(dist_colors)),
         n_train_env = as.integer(n_train_env))




## For each prediction environment, trait, and distance method, find the heritability
## after adding all environments (i.e. the last value) and the local
## maximum heritability (if present)
## 
cumulative_heritability_toplot_maxima <- cumulative_heritability_toplot %>% 
  group_by(environment, trait, dist_method) %>% 
  mutate(local_max_heritability = max(heritability_out),
         local_max_n_train = n_train_env[heritability_out == local_max_heritability],
         local_max_scaled_dist = scaled_distance[heritability_out == local_max_heritability],
         all_env_heritability = heritability_out[n_train_env == max(n_train_env)]) %>%
  ungroup() %>%
  # Is the local maximum greater than the final heritability_out?
  # What is the difference between the local maximum and the final
  mutate(local_adv = local_max_heritability - all_env_heritability,
         local_better = (local_adv > 0))

## Filter for the observations in which the accuracy is the local maximum
cumulative_heritability_toplot_maxima_filter <- cumulative_heritability_toplot_maxima %>%
  filter(heritability_out == local_max_heritability)


## Common plot modifier
g_mod <- list(
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = dist_method, color = NULL), alpha = 0.25),
  # geom_point(aes(x = local_max_scaled_dist, y = local_max_accuracy)),
  geom_point(aes(x = local_max_n_train, y = local_max_heritability)),
  geom_line(lwd = 0.5),
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)),
  ylab("Heritability"),
  xlab("Number of Environments in Training Set"),
  scale_color_manual(values = dist_colors),
  scale_fill_manual(values = dist_colors),
  theme_bw(),
  theme(panel.grid = element_blank(),
        legend.key.height = unit(2, units = "line"),
        text = element_text(size = 8)),
  labs(title = "Cumulative Environmental Clusters")
)



# List over traits
g_cumulative_herit <- setNames(traits, traits) %>%
  map(function(tr) {
    # Subset the trait
    cumulative_heritability_toplot_tr <- cumulative_heritability_toplot_maxima %>%
      filter(trait == tr) %>%
      # mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
      mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
    
    
    # Create one plot using facet_wrap
    g_total <- cumulative_heritability_toplot_tr %>% 
      # ggplot(aes(x = scaled_distance, y = heritability_out, color = dist_method)) +
      ggplot(aes(x = n_train_env, y = heritability_out, color = dist_method)) +
      g_mod +
      facet_wrap( ~ environment + trait)
    # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = dist_method, color = NULL), alpha = 0.1)
    
    # Save
    save_file <- file.path(fig_dir, str_c("cumulative_env_dist_herit_", tr, ".jpg"))
    ggsave(filename = save_file, plot = g_total, height = 15, width = 12, dpi = 1000)
    
    # Return the plot
    # return(g_pages)
    return(g_total)
    
  })


### Plot a single distance method over many environments
# List over traits
g_cumulative_herit1 <- dist_method_replace %>%
  map(function(dm) {
    
    # Subset the distance method
    cumulative_heritability_toplot_dm <- cumulative_heritability_toplot_maxima %>%
      filter(dist_method %in% c(dm, "Random"))
    
    # Iterate over the trait
    g_list <- setNames(traits, traits) %>%
      map(function(tr) {
        # Subset the trait data
        cumulative_heritability_toplot_tr <- cumulative_heritability_toplot_dm %>%
          filter(trait == tr) %>%
          # mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
          mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
        
        
        # Extract the ylim values from the single plots created above
        ylim_tr <- ggplot_build(g_cumulative_herit[[tr]])$layout$panel_ranges[[1]]$y.range
        
        # Create the plot and return
        cumulative_heritability_toplot_tr %>% 
          # ggplot(aes(x = scaled_distance, y = heritability_out)) +
          ggplot(aes(x = n_train_env, y = heritability_out, color = dist_method)) +
          g_mod + 
          ylim(ylim_tr) +
          facet_wrap( ~ environment + trait)
        # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.1)
        
      })
    
    # Iterate over the plots and save
    for (i in seq_along(g_list)) {
      # Save
      save_file <- file.path(
        fig_dir, str_c("cumulative_env_dist_herit_", names(dist_method_replace[dist_method_replace == dm]), 
                       "_", names(g_list)[i], ".jpg"))
      ggsave(filename = save_file, plot = g_list[[i]], height = 15, width = 12, dpi = 1000)
    }
    
    # Return the plot
    return(g_list)
    
  })











##### Heritability using a sliding window #####


## Bind the list elements together and unnest
window_heritability <- env_dist_window_heritability_out %>%
  bind_rows() %>% 
  unnest()

# Re-scale the distance measurements to unit variance
window_heritability_adj <- window_heritability %>%
  group_by(environment, trait, dist_method) %>% 
  mutate(scaled_distance = scale(distance),
         scaled_distance = scaled_distance - min(scaled_distance),
         window_rank = seq(n())) %>% # Create a vector specifying the order of the windows
  ungroup() %>%
  mutate(n_train_env = map_dbl(train_envs, length))


## Separate the original samples from the random samples
## Then summarize the mean and quantiles of the prediction accuracies
## for each additional training population environment
window_heritability_random <- window_heritability_adj %>% 
  filter(str_detect(dist_method, "sample")) %>%
  group_by(environment, trait, window_rank) %>% 
  summarize_at(vars(heritability_out), funs(mean = mean(.), lower = quantile(., probs = 0.025),
                                            upper = quantile(., probs = 0.975))) %>%
  ungroup() %>%
  rename(heritability_out = mean) %>%
  mutate(dist_method = "Random")

# Now grab the original runs
window_heritability_orig <- window_heritability_adj %>% 
  filter(dist_method %in% names(dist_method_replace)) %>%
  select(environment, trait, dist_method, scaled_distance, window_rank, heritability_out)


## Combine
window_heritability_toplot <- bind_rows(window_heritability_orig, window_heritability_random) %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace),
         dist_method = factor(dist_method, levels = names(dist_colors)),
         window_rank = as.integer(window_rank))


## For each prediction environment, trait, and distance method, find the heritability
## after adding all environments (i.e. the last value) and the local
## maximum heritability (if present)
## 
window_heritability_toplot_maxima <- window_heritability_toplot %>% 
  group_by(environment, trait, dist_method) %>% 
  mutate(local_max_heritability = max(heritability_out),
         local_max_window_rank = window_rank[heritability_out == local_max_heritability],
         local_max_scaled_dist = scaled_distance[heritability_out == local_max_heritability],
         all_env_heritability = heritability_out[window_rank == max(window_rank)]) %>%
  ungroup()

## Filter for the observations in which the accuracy is the local maximum
window_heritability_toplot_maxima_filter <- window_heritability_toplot_maxima %>%
  filter(heritability_out == local_max_heritability)




## Common plot modifier
g_mod <- list(
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = dist_method, color = NULL), alpha = 0.25),
  # geom_point(aes(x = local_max_scaled_dist, y = local_max_heritability)),
  geom_point(aes(x = local_max_window_rank, y = local_max_heritability)),
  geom_line(lwd = 0.5),
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)),
  scale_color_manual(values = dist_colors),
  scale_fill_manual(values = dist_colors),
  ylab("Heritability"),
  xlab("Distance Rank of Sliding Window"),
  theme_bw(),
  theme(panel.grid = element_blank(),
        legend.key.height = unit(2, units = "line"),
        text = element_text(size = 8)),
  labs(title = "Sliding Window Environmental Clusters")
)


# List over traits
g_window_herit <- setNames(traits, traits) %>%
  map(function(tr) {
    # Subset the trait
    window_heritability_toplot_tr <- window_heritability_toplot_maxima %>%
      filter(trait == tr) %>%
      # mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
      mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
    
    
    # Create one plot using facet_wrap
    g_total <- window_heritability_toplot_tr %>% 
      # ggplot(aes(x = scaled_distance, y = heritability_out, color = dist_method)) +
      ggplot(aes(x = window_rank, y = heritability_out, color = dist_method)) +
      g_mod +
      facet_wrap( ~ environment + trait)
    # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = dist_method, color = NULL), alpha = 0.1)
    
    # Save
    save_file <- file.path(fig_dir, str_c("window_env_dist_herit_", tr, ".jpg"))
    ggsave(filename = save_file, plot = g_total, height = 15, width = 12, dpi = 1000)
    
    # Return the plot
    # return(g_pages)
    return(g_total)
    
  })


### Plot a single distance method over many environments
# List over traits
g_window_herit1 <- dist_method_replace %>%
  map(function(dm) {
    
    # Subset the distance method
    window_heritability_toplot_dm <- window_heritability_toplot_maxima %>%
      filter(dist_method %in% c(dm, "Random"))
    
    # Iterate over the trait
    g_list <- setNames(traits, traits) %>%
      map(function(tr) {
        # Subset the trait data
        window_heritability_toplot_tr <- window_heritability_toplot_dm %>%
          filter(trait == tr) %>%
          # mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
          mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
        
        
        # Extract the ylim values from the single plots created above
        ylim_tr <- ggplot_build(g_window_herit[[tr]])$layout$panel_ranges[[1]]$y.range
        
        # Create the plot and return
        window_heritability_toplot_tr %>% 
          # ggplot(aes(x = scaled_distance, y = heritability_out)) +
          ggplot(aes(x = window_rank, y = heritability_out, color = dist_method)) +
          g_mod + 
          ylim(ylim_tr) +
          facet_wrap( ~ environment + trait)
        # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.1)
        
      })
    
    # Iterate over the plots and save
    for (i in seq_along(g_list)) {
      # Save
      save_file <- file.path(
        fig_dir, str_c("window_env_dist_herit_", names(dist_method_replace[dist_method_replace == dm]), 
                       "_", names(g_list)[i], ".jpg"))
      ggsave(filename = save_file, plot = g_list[[i]], height = 15, width = 12, dpi = 1000)
    }
    
    # Return the plot
    return(g_list)
    
  })



