## Analysis of predictions based on a realistic breeding scenario
## 
## Author: Jeff Neyhart
## Last Updated: April 24, 2018
## 
## This script will look at prediction accuracies from adding environments after
## ranking based on specific distance metrics, in the context of a realistic
## breeding scenario. This means that one year (e.g. 2016) was predicted with 
## only data from previous years and only using distance metrics that would
## allow for new environments to be predicted.
## 

# # Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Prepare the phenotype BLUEs for analysis
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  filter(line_name %in% vp_geno)


# Load the true results
load(file.path(result_dir, "realistic_breeding_predictions.RData"))


## Create a character vector to replace the criteria for environmental selection
criteria_replace <- c("lrt" = "LRT", "first10" = "First 10 Envs", "first5" = "First 5 Envs", 
                      "first1" = "First 1 Envs", "all_env" = "All Environments")


## Treat the two years separately (because the number of training environments
## is so different)
breeding_prediction_results1 <- breeding_prediction_results %>% 
  list(., names(.)) %>%
  pmap_df(~mutate(.x, year = parse_number(.y)))
  

pred_results <- breeding_prediction_results1 %>% 
  mutate(., boot = map(results_out, "boot"),
         n_train_env = map_int(train_environment, length)) %>% 
  unnest(boot) %>%
  select(-population, -train_environment, -results_out) %>%
  rename(accuracy = cor) %>%
  ## For the 2016 year, drop the first 5 and first 10 criteria
  filter(!(criteria %in% c("first5", "first10") & year == 2016))




## For each prediction environment, calculate the difference between the accuracy
## after using all data and the accuracy after using a specific group of environments
## 
## This will require some re-arranging of data 
pred_results1 <- pred_results %>% 
  group_by(trait, val_environment, dist_method) %>% 
  mutate(terminal_accuracy = tail(accuracy, 1),
         advantage = accuracy - terminal_accuracy) %>%
  # Also add the confidence interval for the "all environment" criteria
  mutate(terminal_lower = tail(ci_lower, 1), terminal_upper = tail(ci_upper, 1)) %>%
  ungroup()
    


## Edit some variables
pred_results_toplot <- pred_results1 %>%
  mutate(., dist_method = str_replace_all(dist_method, dist_method_replace),
         criteria = as_replaced_factor(x = criteria, replacement = criteria_replace),
         criteria = factor(criteria, levels = rev(criteria_replace))) %>%
  group_by(year, trait, val_environment, dist_method) %>%
  mutate(max_train_env = max(n_train_env),
         prop_max_train_env = n_train_env / max_train_env) %>%
  ungroup()



# Plot the distribution of prediction accuracy advantages for each criteria and for 
# each distance method and trait
g_pred_advantage <- pred_results_toplot %>% 
  filter(criteria != "All Environments") %>% 
  ggplot(aes(x = dist_method, y = advantage, fill = criteria)) +
  geom_boxplot(alpha = 0.5) + 
  facet_grid(year ~ trait) +
  ylab("Accuracy Advantage") +
  xlab("Distance Method") +
  scale_fill_discrete(name = "Environment Selection\nCriteria") + 
  theme_bw()

# Save
ggsave(filename = "realistic_prediction_advantage_distribution.jpg", 
       plot = g_pred_advantage, path = fig_dir, width = 8, height = 6, dpi = 1000)


## Plot the proportion of environments in which the optimal training set resulted
## in a prediction accuracy that was not signficantly different from using all
## environments
pred_results_sig <- pred_results_toplot %>% 
  rowwise() %>% 
  mutate(not_different_than_all = between(x = accuracy, left = terminal_lower, right = terminal_upper), 
         sig_less_than_all = accuracy < terminal_lower, 
         sig_greater_than_all = accuracy > terminal_upper) %>% ungroup()


## Note that no accuracies were significantly greater than the terminal accuracy
g_pred_sig <- pred_results_sig %>% 
  filter(criteria != "All Environments") %>%
  group_by(trait, dist_method, criteria, year) %>% 
  summarize_at(vars(not_different_than_all:sig_greater_than_all), mean) %>% 
  # gather(test, prop, not_different_than_all:sig_greater_than_all) %>% 
  ggplot(aes(x = dist_method, y = not_different_than_all, col = criteria)) + 
  geom_point(position = position_dodge(0.5)) +
  ylab("Proportion of Accuracies Not Different from Terminal") +
  xlab("Distance Method") +
  scale_color_discrete(name = "Environment Selection\nCriteria") + 
  facet_grid(year ~ trait) + 
  theme_bw()

ggsave(filename = "realistic_prediction_significant_accuracy.jpg", 
       plot = g_pred_sig, path = fig_dir, width = 8, height = 6, dpi = 1000)



## Plot the proportion of the maximum number of training environments that were
## used to create the optimal trainin environment
pred_results_toplot %>% 
  filter(criteria != "All Environments") %>% 
  ggplot(aes(x = dist_method, y = prop_max_train_env, fill = criteria)) +
  geom_boxplot(alpha = 0.5) + 
  facet_grid(year ~ trait)




## 
## Fit a model to see the average effect of each criteria
pred_fit <- pred_results1 %>% 
  filter(criteria != "all_env") %>%
  mutate(year = as.factor(str_extract(string = val_environment, pattern = "[0-9]{2}")),
         val_environment = as.factor(val_environment),
         criteria = factor(criteria, levels = names(rev(criteria_replace)))) %>%
  group_by(trait, dist_method) %>%
  do(fit = lm(advantage ~ criteria + val_environment + n_train_env, data = .))











### Examine the results of accumulating environments to predict a new environment

## This may need to be changed
breeding_prediction_iterative_results <- breeding_prediction_iterative_results_use

# Bind rows
breeding_prediction_iterative_results1 <- breeding_prediction_iterative_results %>% 
  bind_rows() %>% 
  unnest(results_out) %>% 
  mutate(n_train_env = map_int(train_envs, length))





### Examine the results of randomly adding environments
# First add a designator for the sample number
breeding_prediction_iterative_results_random1 <- breeding_prediction_iterative_results_random %>% 
  map(~group_by(., trait) %>% mutate(dist_method = str_c("sample", seq(n()))) %>% ungroup) %>%
  # Add the year designator
  list(., names(.)) %>%
  pmap_df(~mutate(.x, year = parse_number(.y)))

# Combine with the observed phenotypic values
breeding_prediction_iterative_results_random2 <- breeding_prediction_iterative_results_random1 %>% 
  unnest(results_out) %>% 
  mutate(n_train_env = map_int(train_env, length)) %>%
  unnest(predict_out) %>% 
  inner_join(S2_MET_BLUEs_use, .) %>% 
  select(-trial, -location, -std_error)

## For the random predicitions accuracies, we first need to calculate the prediction
## accuracy by correlating the predicted values with the observed
## 
## Summarize the results of the random environmental additions by calculating 
## the mean accuracy and 95% confidence interval
breeding_prediction_iterative_results_random3 <- breeding_prediction_iterative_results_random2 %>%
  group_by(environment, trait, n_train_env, dist_method) %>% 
  summarize(accuracy = cor(value, pred_value)) %>%
  summarize_at(vars(accuracy), funs(cor = mean(.), lower = quantile(., probs = 0.025),
                                    upper = quantile(., probs = 0.975)))
  

## Create a df to plot
breeding_prediction_iterative_results_random_toplot <- breeding_prediction_iterative_results_random3 %>%
  mutate(dist_method = "Random") %>%
  rename(val_environment = environment) %>%
  ungroup()



## Split the datasets and iterate in parallel
g_breeding_iterative <- list(breeding_prediction_iterative_results1 %>% split(.$trait),
                             breeding_prediction_iterative_results_random_toplot %>% split(.$trait)) %>%
  pmap(~{
    
    # Sort environments by final prediction accuracy
    sorted_envs <- .x %>% 
      group_by(val_environment) %>% 
      filter(n_train_env == max(n_train_env),
             dist_method == dist_method[1]) %>% 
      ungroup() %>% 
      select(val_environment, dist_method, cor) %>% 
      arrange(desc(cor)) %>% 
      mutate(val_environment = factor(val_environment, levels = .$val_environment))
    
    ## Change the order of the environments
    dfx <- mutate(.x, val_environment = factor(val_environment, levels = levels(sorted_envs$val_environment)))
    dfy <- mutate(.y, val_environment = factor(val_environment, levels = levels(sorted_envs$val_environment)))
    
    # Plot
    g <- dfx %>% 
      mutate(val_environment = factor(val_environment, levels = levels(sorted_envs$val_environment))) %>%
      ggplot(aes(x = n_train_env, y = cor, col = dist_method, fill = dist_method)) + 
      geom_ribbon(data = dfy, aes(ymin = lower, ymax = upper), alpha = 0.25) +
      geom_line(lwd = 1) + 
      geom_line(data = dfy, lwd = 1) +
      facet_wrap(~ val_environment + trait) + 
      ylab("Prediction Accuracy") +
      xlab("Number of Training Environments") +
      theme_bw()
    
    ## Save
    save_file <- file.path(fig_dir, str_c("realistic_cumulative_predictions_", unique(dfx$trait), ".jpg"))
    ggsave(filename = save_file, plot = g, height = 12, width = 12, dpi = 1000)
    
    return(g)
    
  })
    












