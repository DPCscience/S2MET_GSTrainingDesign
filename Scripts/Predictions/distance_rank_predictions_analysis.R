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

# # Load some packages
library(lubridate)
library(ggforce)
library(ggridges)

# Load the results
load(file.path(result_dir, "environmental_distance_predictions.RData"))
load(file.path(result_dir, "environmental_distance_window_predictions.RData"))
load(file.path(result_dir, "environmental_distance_heritability.RData"))

# Load the results of the model-based clustering
load(file.path(result_dir, "environmental_distance_lrt.RData"))

## Create a color scheme for the distance methods
dist_colors <- c(setNames(umn_palette(3, length(dist_method_replace)), dist_method_replace), "Random" = "grey75")
## Significant level
alpha <- 0.05


## Bind the elements of the LRT results list
env_dist_lrt_results <- env_dist_lrt_predictions_out %>%
  bind_rows() %>%
  select(environment:dist_method, results_out) %>%
  unnest()

## Two different tests using the LRT results
## Either select when it is first significant ("first"), or select when it is MOST significant ("most")
## If no points are significant, take the last point.
env_dist_lrt_results_first <- env_dist_lrt_results %>% 
  group_by(environment, trait, dist_method) %>% 
  mutate(is_sig = p_value <= alpha | n_env == max(n_env)) %>% 
  filter(is_sig) %>%
  mutate(which_sig = which(is_sig)) %>% 
  filter(which_sig == min(which_sig)) %>%
  ungroup()


## Bind the list elements together and unnest
cumulative_pred_results <- env_dist_predictions_out %>% 
  bind_rows() %>% 
  unnest()

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


##########







### Plots of prediction accuracy.



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

## Combine
cumulative_pred_toplot <- bind_rows(cumulative_pred_orig, cumulative_pred_random) %>%
  mutate(n_train_env = as.integer(n_train_env))




## For each prediction environment, trait, and distance method, find the prediction
## accuracy after adding all environments (i.e. the last value) and the local
## maximum prediction accuracy (if present)
## 
cumulative_pred_toplot_final <- cumulative_pred_toplot %>% 
  group_by(environment, trait, dist_method) %>% 
  mutate(all_env_acccuracy = accuracy[n_train_env == max(n_train_env)]) %>%
  ungroup() %>% 
  mutate(accuracy_gt_final = accuracy > all_env_acccuracy)


## Find the local maximum and create a data.frame
cumulative_pred_toplot_maximum <- cumulative_pred_toplot_final %>%  
  group_by(environment, trait, dist_method) %>% 
  filter(accuracy == max(accuracy)) %>%
  ungroup() %>%
  mutate(annotation = "local_maximum")

Ntrain_envs <- c(1, 5, 10)

## Find the accuray at n = 1, 5, 10 training environments
cumulative_pred_toplot_Nenv <- cumulative_pred_toplot_final %>%  
  filter(n_train_env %in% Ntrain_envs) %>%
  mutate(annotation = str_c(n_train_env, "envs"))

## Use the LRT results to select the accuracy when the LRT is deemed significant.
cumulative_pred_toplot_lrt <- cumulative_pred_toplot_final %>%
  inner_join(., env_dist_lrt_results_first, by = c("environment", "trait", "dist_method", "n_train_env" = "n_env")) %>% 
  select(environment:all_env_acccuracy) %>% 
  mutate(annotation = "lrt")


## A vector to replace the annotation names
annotation_replace <- c("local_maximum" = "Maximum", 
                        setNames(str_replace(string = str_c(Ntrain_envs, "envs"), pattern = "e", replacement = " E"), 
                                 str_c(Ntrain_envs, "envs")), 
                        "lrt" = "Likelihood Ratio Test")


## Combine the data together
cumulative_pred_toplot_annotate <- bind_rows(cumulative_pred_toplot_maximum, 
                                             cumulative_pred_toplot_Nenv,
                                             cumulative_pred_toplot_lrt) %>%
  # Calculate the difference between the annotated accuracy and the "final" accuracy
  mutate(advantage = accuracy - all_env_acccuracy,
         pos_advantage = advantage > 0,
         annotation = as_replaced_factor(x = annotation, replacement = annotation_replace))


## Plot
## Distribution of the number of training environments or the scaled distance (cumulative) 
## of the local maximum
g_local_max <- cumulative_pred_toplot_annotate %>% 
  filter(annotation %in% c("local_maximum")) %>%
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
  xlab("Prediction Accuracy Advantage") +
  labs(title = "Advantage of Selected Environments Over All Environments") +
  facet_grid( ~ trait) +
  theme_bw() +
  theme(legend.key.height = unit(2, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1))

save_file <- file.path(fig_dir, "cumulative_pred_max_advantage.jpg")
ggsave(filename = save_file, plot = g_local_adv, height = 6, width = 8, dpi = 1000)


####### ANALYSIS ######
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



  

## Common plot modifier
g_mod <- list(
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = dist_method, color = NULL), alpha = 0.25),
  # geom_point(aes(x = local_max_scaled_dist, y = local_max_accuracy)),
  geom_line(lwd = 0.5),
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)),
  ylab("Prediction Accuracy"),
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
g_cumulative_pred <- setNames(traits, traits) %>%
  map(function(tr) {
    # Subset the trait
    cumulative_pred_toplot_tr <- cumulative_pred_toplot_final %>%
      filter(trait == tr) %>%
      # mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
      mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
    
    cumulative_pred_toplot_annotate_tr <- cumulative_pred_toplot_annotate %>%
      filter(trait == tr) %>%
      mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
      
    # Create one plot using facet_wrap
    g_total <- cumulative_pred_toplot_tr %>% 
      # ggplot(aes(x = scaled_distance, y = accuracy, color = dist_method)) +
      ggplot(aes(x = n_train_env, y = accuracy, color = dist_method)) +
      geom_point(data = cumulative_pred_toplot_annotate_tr, aes(shape = annotation)) +
      g_mod +
      facet_wrap( ~ environment + trait)
      # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = dist_method, color = NULL), alpha = 0.1)
    
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
    cumulative_pred_toplot_dm <- cumulative_pred_toplot_final %>%
      filter(dist_method %in% c(dm, "Random"))
    
    cumulative_pred_toplot_annotate_dm <- cumulative_pred_toplot_annotate %>%
      filter(dist_method %in% c(dm, "Random"))
    
    # Iterate over the trait
    g_list <- setNames(traits, traits) %>%
      map(function(tr) {
        # Subset the trait data
        cumulative_pred_toplot_tr <- cumulative_pred_toplot_dm %>%
          filter(trait == tr) %>%
          # mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
          mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
        
        cumulative_pred_toplot_annotate_tr <- cumulative_pred_toplot_annotate_dm %>%
          filter(trait == tr) %>%
          mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
        
          
        # Extract the ylim values from the single plots created above
        ylim_tr <- ggplot_build(g_cumulative_pred[[tr]])$layout$panel_ranges[[1]]$y.range
        
        # Create the plot and return
        cumulative_pred_toplot_tr %>% 
          # ggplot(aes(x = scaled_distance, y = accuracy)) +
          ggplot(aes(x = n_train_env, y = accuracy, color = dist_method)) +
          geom_point(data = cumulative_pred_toplot_annotate_tr, aes(shape = annotation)) +
          g_mod + 
          ylim(ylim_tr) +
          facet_wrap( ~ environment + trait)
          # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.1)
        
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
window_pred_adj <- window_pred_results %>%
  group_by(environment, trait, dist_method) %>% 
  mutate(scaled_distance = scale(distance),
         scaled_distance = scaled_distance - min(scaled_distance),
         window_rank = seq(n())) %>% # Create a vector specifying the order of the windows
  ungroup() %>%
  mutate(n_train_env = map_dbl(train_envs, length))


## Separate the original samples from the random samples
## Then summarize the mean and quantiles of the prediction accuracies
## for each additional training population environment
window_pred_random <- window_pred_adj %>% 
  filter(str_detect(dist_method, "sample")) %>%
  group_by(environment, trait, window_rank) %>% 
  summarize_at(vars(accuracy, scaled_distance), 
               funs(mean = mean(.), lower = quantile(., probs = 0.025), upper = quantile(., probs = 0.975))) %>%
  ungroup() %>%
  select(environment:window_rank, accuracy = accuracy_mean, scaled_distance = scaled_distance_mean,
         lower = accuracy_lower, upper = accuracy_upper) %>%
  mutate(dist_method = "Random")

# Now grab the original runs
window_pred_orig <- window_pred_adj %>% 
  filter(dist_method %in% names(dist_method_replace)) %>%
  select(environment, trait, dist_method, scaled_distance, window_rank, accuracy)

## Combine
window_pred_toplot <- bind_rows(window_pred_orig, window_pred_random) %>%
  mutate(dist_method = as_replaced_factor(dist_method, c(dist_method_replace, "Random" = "Random")),
         window_rank = as.integer(window_rank))




## Identify the largest prediction accuracy among the sliding windows
window_pred_toplot_maxima <- window_pred_toplot %>%
  group_by(environment, trait, dist_method) %>% 
  mutate(local_max_accuracy = max(accuracy),
         local_max_window_rank = window_rank[accuracy == local_max_accuracy]) %>%
  ungroup()

## Filter for the local accuracy maxima
window_pred_toplot_maxima_filter <- window_pred_toplot_maxima %>% 
  filter(accuracy == local_max_accuracy)
  



## Distribution of the number of training environments or the scaled distance (cumulative) 
## of the local maximum
g_local_max <- window_pred_toplot_maxima_filter %>% 
  # filter(dist_method != "Random") %>%
  ggplot(aes(x = scaled_distance, fill = dist_method)) +
  # ggplot(aes(x = window_rank, fill = dist_method)) +
  # geom_density(alpha = 0.25) + 
  geom_density_ridges(aes(y = dist_method), alpha = 0.25) +
  scale_fill_manual(values = dist_colors, guide = FALSE) +
  xlab("Scaled Average Environmental Distance at Maximum Accuracy") +
  labs(title = "Average Environmental Distance at Maximum Accuracy") +
  facet_wrap(~ trait, ncol = 2) +
  theme_bw() +
  theme(legend.key.height = unit(2, "lines"))

save_file <- file.path(fig_dir, "window_pred_max_dist.jpg")
ggsave(filename = save_file, plot = g_local_max, height = 6, width = 8, dpi = 1000)


g_local_max <- window_pred_toplot_maxima_filter %>% 
  # filter(dist_method != "Random") %>%
  # ggplot(aes(x = scaled_distance, fill = dist_method)) +
  ggplot(aes(x = window_rank, fill = dist_method)) +
  # geom_density(alpha = 0.25) + 
  geom_density_ridges(aes(y = dist_method), alpha = 0.25) +
  scale_fill_manual(values = dist_colors, guide = FALSE) +
  xlab("Rank of Sliding Window at Maximum Accuracy") +
  labs(title = "Rank of Sliding Window at Maximum Accuracy") +
  facet_wrap(~ trait, ncol = 2) +
  theme_bw() +
  theme(legend.key.height = unit(2, "lines"))

save_file <- file.path(fig_dir, "window_pred_max_rank.jpg")
ggsave(filename = save_file, plot = g_local_max, height = 6, width = 8, dpi = 1000)



## Calculate the average environmental distance (scaled) where the local maxima is
## found and calculate the average rank of the window where the maximum is found
window_pred_maxima_summary <- window_pred_toplot_maxima_filter %>% 
  group_by(trait, dist_method) %>% 
  summarize_at(vars(window_rank, scaled_distance), mean)




## Common plot modifier
g_mod <- list(
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = dist_method, color = NULL), alpha = 0.25),
  # geom_point(aes(x = local_max_scaled_dist, y = local_max_accuracy)),
  geom_point(aes(x = local_max_window_rank, y = local_max_accuracy)),
  geom_line(lwd = 0.5),
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)),
  scale_color_manual(values = dist_colors),
  scale_fill_manual(values = dist_colors),
  ylab("Prediction Accuracy"),
  xlab("Distance Rank of Sliding Window"),
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
      # mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
      mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
      
    # Create one plot using facet_wrap
    g_total <- window_pred_toplot_tr %>% 
      # ggplot(aes(x = scaled_distance, y = accuracy, color = dist_method)) +
      ggplot(aes(x = window_rank, y = accuracy, color = dist_method)) +
      g_mod +
      facet_wrap( ~ environment + trait)
      # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = dist_method, color = NULL), alpha = 0.1)
    
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
      filter(dist_method %in% c(dm, "Random"))
    
    # Iterate over the trait
    g_list <- setNames(traits, traits) %>%
      map(function(tr) {
        # Subset the trait data
        window_pred_toplot_tr <- window_pred_toplot_dm %>%
          filter(trait == tr) %>%
          # mutate(environment = factor(environment, levels = levels(env_herit_rank[[tr]]$environment)))
          mutate(environment = factor(environment, levels = levels(cumulative_env_pred_order[[tr]]$environment)))
          
        # Extract the ylim values from the single plots created above
        ylim_tr <- ggplot_build(g_window_pred[[tr]])$layout$panel_ranges[[1]]$y.range
        
        # Create the plot and return
        window_pred_toplot_tr %>% 
          # ggplot(aes(x = scaled_distance, y = accuracy)) +
          ggplot(aes(x = window_rank, y = accuracy, color = dist_method)) +
          g_mod + 
          ylim(ylim_tr) +
          facet_wrap( ~ environment + trait)
          # geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.1)
        
      })
    
    # Iterate over the plots and save
    for (i in seq_along(g_list)) {
      # Save
      save_file <- file.path(
        fig_dir, str_c("window_env_dist_pred_", names(dist_method_replace[dist_method_replace == dm]), 
                       "_", names(g_list)[i], ".jpg"))
      g_print <- g_list[[i]] + geom_smooth(method = "lm", se = FALSE)
      
      ggsave(filename = save_file, plot = , height = 15, width = 12, dpi = 1000)
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









## Prediction accuracy using a sliding window
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











