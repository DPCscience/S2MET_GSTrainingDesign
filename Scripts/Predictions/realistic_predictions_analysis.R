## Analysis of predictions based on a realistic breeding scenario
## 
## Author: Jeff Neyhart
## Last Updated: April 22, 2018
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

# Load the results
load(file.path(result_dir, "realistic_breeding_predictions.RData"))

## Bind the rows of the list and unnest the results
breeding_prediction_results1 <- breeding_prediction_results %>%
  bind_rows()

pred_results <- breeding_prediction_results1 %>% 
  mutate(accuracy = map_dbl(results_out, "accuracy"),
         n_train_env = map_int(train_environment, length)) %>% 
  select(-train_environment, -results_out)




## For each prediction environment, calculate the difference between the accuracy
## after using all data and the accuracy after using a specific group of environments
## 
## This will require some re-arranging of data 
pred_results1 <- pred_results %>% 
  left_join(., pred_results %>% filter(dist_method == "all_env") %>% 
              spread(dist_method, accuracy) %>% select(-criteria, -n_train_env),
            by = c("trait", "val_environment", "population")) %>%
  mutate(advantage = accuracy - all_env)

## Edit some variables
pred_results_toplot <- pred_results1 %>%
  mutate(dist_method = str_replace_all(dist_method, dist_method_replace),
         criteria = str_replace_all(criteria, c("lrt" = "LRT", "first10" = "First 10 Envs",
                                                "first5" = "First 5 Envs", "first1" = "First 1 Envs")),
         criteria = factor(criteria, levels = c("LRT", "First 1 Envs", "First 5 Envs", "First 10 Envs")))



## Plot
pred_results_toplot %>% filter(dist_method != "all_env") %>% ggplot(aes(x = dist_method, y = advantage, fill = criteria)) + geom_boxplot() + facet_grid(~ trait)


## Summarize
pred_results_summ <- pred_results_toplot %>% 
  filter(dist_method != "all_env") %>% 
  group_by(trait, dist_method, criteria) %>% 
  summarize(mean_adv = mean(advantage),
            sd_adv = sd(advantage))


pred_results_summ %>% ggplot(aes(x = dist_method, y = mean_adv, fill = criteria)) + geom_col(position = "dodge") + geom_errorbar(aes(ymin = mean_adv - sd_adv, ymax = mean_adv + sd_adv), width = 0.5, position = position_dodge(0.9)) + facet_grid(~ trait)


## 
## Fit a model to see the average effect of each criteria
pred_fit <- pred_results1 %>% 
  filter(dist_method != "all_env") %>%
  mutate(year = as.factor(str_extract(string = val_environment, pattern = "[0-9]{2}")),
         val_environment = as.factor(val_environment)) %>%
  group_by(trait, dist_method) %>%
  do(fit = lm(advantage ~ criteria + val_environment, data = .))















