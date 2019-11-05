## s2MET Predictions
## Growth stage environmental average
## 
## Author: Jeff Neyhart
## Last modified: 29 July 2019
## 
## This R script will predict the environmental average growth stage
## based on GDD.
## 

# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load other libraries
library(lubridate)

## Load the weather dataset
load(file.path(data_dir, "Climate_Data/NOAA_Data/noaa_stations_trial_data.RData"))

## Number of years for the multiyear window
multiyear_window <- 10




### One-year environmental weather data



oneyear_env_data_unnest <- noaa_trial_data_oneyear_complete %>%
  unnest(data)

oneyear_trial_env_data <- trial_info %>% 
  distinct(environment, location, year, planting_date) %>%
  # mutate_all(parse_guess) %>%
  left_join(., oneyear_env_data_unnest) %>%
  mutate_at(vars(planting_date, date), ymd) %>%
  # Convert tenth of value to value
  mutate(value = value / 10)

## Select TMIN and TMAX
oneyear_temp_data <- oneyear_trial_env_data %>%
  filter(datatype %in% c("TMIN", "TMAX")) %>%
  filter(date >= planting_date) %>% 
  mutate(dap = as.numeric(date - planting_date))


# Calculate the mean daily temperature
# Use a loess smooth to impute missing temperature data
one_year_daily_temp_impute <- oneyear_temp_data %>%
  arrange(environment, datatype, dap) %>%
  group_by(environment, datatype) %>%
  do({
    df <- .
    
    # cat(unique(df$environment), unique(as.character(df$datatype)), unique(df$year), "\n")
    
    n_missing_start <- n_missing <- sum(is.na(df$value))
    
    # Reassign df to df1 (overwrite if missing data is present)
    df1 <- df
    
    # Iterate over window sizes until the missing data is imputed
    wind <- 2
    while(n_missing > 0) {
      
      df1 <- df %>%
        mutate(pred_value = window_mean(x = dap, y = value, window = wind)) %>%
        mutate(value = ifelse(is.na(value), pred_value, value)) %>%
        select(dap, date, value)
      
      n_missing <- sum(is.na(df1$value))
      wind <- wind + 1
      
    }
    
    # Remove grouping variables
    df2 <- df1[,c("dap", "date", "value")]
    
    # Return the data.frame
    tibble(imputed_data = list(df2), window = wind - 1, n_missing = n_missing_start)
    
  }) %>% ungroup()



## Calculate growing degree days (GDD) based on: https://ndawn.ndsu.nodak.edu/help-barley-growing-degree-days.html

## Calculate GDD, AGDD, haun stage, and growth stage for all environments using
## one-year data
one_year_growth_staging <- one_year_daily_temp_impute %>% 
  unnest() %>% 
  # Remove dap = 0
  filter(dap > 0) %>%
  select(-window, -n_missing) %>%
  spread(datatype, value) %>%
  split(.$environment) %>%
  map_df(~mutate(., gdd = gdd(tmin = TMIN, tmax = TMAX), agdd = cumsum(gdd),
                 growth_stage1 = stage_growth(x = agdd), growth_stage2 = stage_growth(x = agdd, method = "montana")) )


## Plot
one_year_growth_staging %>%
  # filter(environment %in% sample(unique(environment), 3)) %>%
  # ggplot(aes(x = dap, y = environment, color = growth_stage1)) +
  ggplot(aes(x = dap, y = environment, color = growth_stage2)) +
  geom_point()


## Range in length of each stage
one_year_growth_staging %>% 
  # group_by(growth_stage1, environment) %>%
  group_by(growth_stage2, environment) %>%
  summarize(length = n()) %>% 
  summarize_at(vars(length), list(min = min, max = max))
  







### 10-year average


multiyear_env_data_unnest <- noaa_trial_data_multiyear_complete %>%
  unnest(data)

## Find the average planting date per location
location_planting_date <- trial_info %>% 
  filter(!is.na(planting_date)) %>%
  mutate(mean_planting_day = yday(ymd(planting_date))) %>% 
  group_by(location) %>% 
  summarize_at(vars(mean_planting_day), ~floor(mean(.)))


multiyear_trial_env_data <- trial_info %>% 
  filter(!is.na(planting_date)) %>%
  distinct(environment, location, year) %>%
  left_join(., location_planting_date) %>%
  rename(env_year = year) %>%
  left_join(., multiyear_env_data_unnest) %>%
  mutate(day = yday(ymd(date)),
         value = value / 10)


## Select TMIN and TMAX
multiyear_temp_data <- multiyear_trial_env_data %>%
  filter(datatype %in% c("TMIN", "TMAX")) %>%
  filter(day >= mean_planting_day) %>% 
  mutate(dap = as.numeric(day - mean_planting_day))


# Calculate the mean daily temperature
# Use a loess smooth to impute missing temperature data
multiyear_daily_temp_impute <- multiyear_temp_data %>%
  arrange(environment, datatype, year, dap) %>%
  group_by(environment, year, datatype) %>%
  do({
    df <- .
    
    # cat(unique(df$environment), unique(as.character(df$datatype)), unique(df$year), "\n")
    
    n_missing_start <- n_missing <- sum(is.na(df$value))
    
    # Reassign df to df1 (overwrite if missing data is present)
    df1 <- df
    
    # Iterate over window sizes until the missing data is imputed
    wind <- 2
    while(n_missing > 0) {
      
      df1 <- df %>%
        mutate(pred_value = window_mean(x = dap, y = value, window = wind)) %>%
        mutate(value = ifelse(is.na(value), pred_value, value)) %>%
        select(dap, date, value)
      
      n_missing <- sum(is.na(df1$value))
      wind <- wind + 1
      
    }
    
    # Remove grouping variables
    df2 <- df1[,c("dap", "date", "value")]
    
    # Return the data.frame
    tibble(imputed_data = list(df2), window = wind - 1, n_missing = n_missing_start)
    
  }) %>% ungroup()





## Calculate growing degree days (GDD) based on: https://ndawn.ndsu.nodak.edu/help-barley-growing-degree-days.html

## Calculate GDD, AGDD, haun stage, and growth stage for all environments using
## multi-year data
## Here, we will average the GDD for each day after the average planting, then calculate
## AGDD from those GDD values

multiyear_mean_gdd <- multiyear_daily_temp_impute %>% 
  ## Add env_year back in
  left_join(., distinct(multiyear_trial_env_data, environment, env_year, mean_planting_day)) %>%
  unnest() %>% 
  ## Filter observations that are 10 years preceding the trial year (not including trial year)
  filter(year >= env_year - multiyear_window, year <= env_year - 1) %>% 
  # Filter dap > 0
  filter(dap > 0) %>%
  select(-window, -n_missing) %>%
  spread(datatype, value) %>%
  group_by(environment, year) %>%
  do( mutate(., gdd = gdd(tmin = TMIN, tmax = TMAX)) %>% select(mean_planting_day, dap:gdd) ) %>%
  # Group by environment and dap and calculate average GDD
  group_by(environment, dap) %>%
  summarize_at(vars(gdd, mean_planting_day), mean) %>%
  ungroup() %>%
  rename(avg_gdd = gdd)



## Plot to verify
qplot(x = dap, y = avg_gdd, data = multiyear_mean_gdd, color = environment) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE)

## Calculate AGDD per environment, then growth stage
multiyear_growth_staging <- multiyear_mean_gdd %>%
  split(.$environment) %>%
  map_df(~mutate(., avg_agdd = cumsum(avg_gdd), growth_stage1 = stage_growth(x = avg_agdd), 
                 growth_stage2 = stage_growth(x = avg_agdd, method = "montana")) )


## Plot
multiyear_growth_staging %>%
  # filter(environment %in% sample(unique(environment), 3)) %>%
  # ggplot(aes(x = dap, y = environment, color = growth_stage1)) +
  ggplot(aes(x = dap, y = environment, color = growth_stage2)) +
  geom_point()







### Validate by comparing the predicted flowering time via AGDD with the average
### heading date in an environment

## Load phenotypic data
load(file.path(result_dir, "genotype_environment_phenotypic_analysis.RData"))


data_to_model <- training_sets_twoway %>%
  select(-varE) %>% 
  unnest() %>%
  distinct(set, trait, environment, line_name) %>%
  left_join(., S2_MET_BLUEs)


# Fit the model per trait, extract the environmental means
full_model_fits <- data_to_model %>%
  mutate_at(., vars(line_name, environment), as.factor) %>%
  group_by(trait, set) %>%
  do(model = {
    
    df <- droplevels(.)
    
    # Control and weights
    control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", calc.derivs = FALSE)
    wts <- df$std_error^2
    
    # Formula
    form <- value ~ (1|line_name) + environment
    # Fit the model
    lmer(formula = form, data = df, control = control, weights = wts, contrasts = list(environment = "contr.sum"))
    
  })

# Extract the environmental means
env_means_all <- full_model_fits %>%
  group_by(trait, set) %>%
  do({
    
    df <- .
    mod <- df$model[[1]]
    
    mf <- model.frame(mod)
    
    # Get the fixed coefficients
    env_eff <- fixef(mod)
    # Tidy
    tibble(environment = levels(mf$environment), 
           h = c(env_eff[-1], -sum(env_eff[-1])),
           mean = h + env_eff[1])  
    
  }) %>% ungroup()


## Use the growth staging to determine the dap when flowering begins
flowering_oneyear <- one_year_growth_staging %>% 
  # filter(growth_stage1 == "flowering") %>% 
  filter(growth_stage2 == "flowering") %>% 
  group_by(environment) %>% 
  # filter(dap == min(dap)) %>%
  # filter(dap == max(dap)) %>%
  summarize(dap = mean(dap)) %>%
  ungroup()

## Compare average flowering with predicted
compare_flowering_oneyear <- env_means_all %>% 
  filter(trait == "HeadingDate") %>%
  left_join(., flowering_oneyear) %>%
  # Calculate bias
  mutate(bias = dap - mean) %>%
  # Prepare for plotting
  filter(set %in% c("complete", "realistic2017")) %>%
  mutate(set = str_remove_all(set, "2017"),
         set = str_replace_all(set, set_replace))

aggregate(bias ~ set, compare_flowering_oneyear, mean)


## Summarize
compare_flowering_oneyear %>%
  group_by(set) %>%
  do(test = cor.test(.$mean, .$dap)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"))

# Plot
limits <- select(compare_flowering_oneyear, mean, dap) %>% 
  {c(min(.), max(.))}

g_compare_flowering_oneyear <- compare_flowering_oneyear %>%
  ggplot(aes(x = dap, y = mean)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  ggrepel::geom_text_repel(data = subset(compare_flowering_oneyear, abs(bias) - 5 > 5), aes(label = environment), size = 3) +
  facet_wrap(~ set, ncol = 2) +
  scale_y_continuous(name = "Average observed heading date", breaks = pretty, limits = limits) +
  scale_x_continuous(name = "GDD predicted flowering date", breaks = pretty, limits = limits) +
  theme_presentation2()

## Save
ggsave(filename = "one_year_gdd_predicted_FT.jpg", plot = g_compare_flowering_oneyear,
       path = fig_dir, width = 6, height = 4, dpi = 1000)


## Repeat for multiyear growth staging
## Compare average flowering with predicted
flowering_multiyear <- multiyear_growth_staging %>% 
  # filter(growth_stage1 == "flowering") %>% 
  filter(growth_stage2 == "flowering") %>% 
  group_by(environment) %>% 
  # filter(dap == min(dap)) %>%
  # filter(dap == max(dap)) %>%
  summarize(dap = mean(dap)) %>%
  ungroup()


compare_flowering_multiyear <- env_means_all %>% 
  filter(trait == "HeadingDate") %>%
  left_join(., flowering_multiyear) %>%
  # Calculate bias
  mutate(bias = dap - mean) %>%
  # Prepare for plotting
  filter(set %in% c("complete", "realistic2017")) %>%
  mutate(set = str_remove_all(set, "2017"),
         set = str_replace_all(set, set_replace))

## Average bias
aggregate(bias ~ set, compare_flowering_multiyear, mean)


## Summarize
compare_flowering_multiyear %>%
  group_by(set) %>%
  summarize(cor = cor(mean, dap))
# Plot
limits <- select(compare_flowering_multiyear, mean, dap) %>% 
  {c(min(.), max(.))}

g_compare_flowering_multiyear <- compare_flowering_multiyear %>%
  ggplot(aes(x = dap, y = mean)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  ggrepel::geom_text_repel(data = subset(compare_flowering_multiyear, abs(bias) - 5 > 5), 
                           aes(label = environment), size = 3) +
  facet_wrap(~ set, ncol = 2) +
  scale_y_continuous(name = "Average observed heading date", breaks = pretty, limits = limits) +
  scale_x_continuous(name = "GDD predicted flowering date", breaks = pretty, limits = limits) +
  theme_presentation2()

## Save
ggsave(filename = "multiyear_gdd_predicted_FT.jpg", plot = g_compare_flowering_multiyear,
       path = fig_dir, width = 6, height = 4, dpi = 1000)




# Save the data
save_file <- file.path(data_dir, "agdd_growth_staging.RData")
save("one_year_growth_staging", "multiyear_growth_staging", file = save_file)

