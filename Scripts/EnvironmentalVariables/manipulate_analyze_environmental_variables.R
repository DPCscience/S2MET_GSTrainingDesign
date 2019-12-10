## Manipulate and Analyze Environmental Covariates
## 
## Author: Jeff Neyhart
## Last modified: June 11, 2018
## 
## This R Script outlines the methods to manipulate and analyze some of the collected
## environmental covariables
## 


# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load other libraries
library(lubridate)


### Transformation and Compilation 

# Load both climate and soil datasets
load(file.path(data_dir, "Climate_Data/NOAA_Data/noaa_stations_trial_data.RData"))
load(file.path(data_dir, "Soil_Data/complete_trial_soil_data.RData"))
load(file.path(data_dir, "Climate_Data/USNO_daylength_data.RData"))

## Load the growth stage predictions
load(file.path(data_dir, "agdd_growth_staging.RData"))


## Soil data
# Edit the soil data (define topsoil and subsoil, then scale and center). Note in cm
depth_cutoff <- 20

# Create a function that calculates the the proportion of overlap given two ranges
range_overlap <- function(x, y) {
  is_engulfed <- between(x = y[1], x[1], x[2]) & between(x = y[2], x[1], x[2])
  if (is_engulfed) {
    return(1)
  } else if (between(x = y[1], x[1], x[2])) {
    (x[2] - y[1]) / diff(y)
  } else {
    return(0)
  }}



# Calculate the proportion of each layer within the topsoil
# Use a weighted mean of the proportion of each horizon within topsoil vs subsoil
# to calculate the value in each layer
soil_data_complete1 <- soil_data_complete %>% 
  mutate(prop_subsoil = 1 - pmap_dbl(select(., hzdept_r, hzdepb_r),
                                     ~range_overlap(x = c(0, depth_cutoff), y = c(.x, .y)))) %>% 
  group_by(environment, location, variable) %>% 
  # Weight the mean value in each soil layer by the proportion of the horizon in
  # the layer
  summarize(subsoil = weighted.mean(x = value, w = prop_subsoil, na.rm = T), 
            topsoil = weighted.mean(x = value, w = 1 - prop_subsoil, na.rm = T)) %>%
  gather(layer, value, -environment:-variable) %>%
  # Center and scale
  unite(variable, variable, layer, sep = "_") %>%
  ungroup() %>%
  select(environment, variable, value)



## Create summary environmental variables. We will use those outlined in @Anderson2016. 
## Essentially, we will calculate the monthly min, mean, and max for temperature and 
## precipitation. We will also calculate the min mean and max of one variables during
## the wettest/driest or the warmest/coldest month.

## Use the USGS Bioclimatic Predictors for Supporting Ecological Applications in the Conterminous United States

### One-year environmental covariates
oneyear_env_data_unnest <- noaa_trial_data_oneyear_complete %>%
  unnest(data)

oneyear_trial_env_data <- trial_info %>% 
  distinct(environment, location, year, planting_date) %>%
  left_join(., oneyear_env_data_unnest) %>%
  mutate_at(vars(planting_date, date), ymd) %>%
  # Convert tenth of value to value
  mutate(value = value / 10)

# Summarize longitude/latitude
coord_clim_summ <- oneyear_trial_env_data %>% 
  select(environment, latitude, longitude, elevation) %>% 
  distinct() %>%
  gather(variable, value, -environment) %>%
  filter(!is.na(value))

## Subset data only after planting, then add 30-day intervals
oneyear_trial_env_data_interval <- oneyear_trial_env_data %>% 
  select(environment, planting_date, date, month, datatype, value) %>% 
  filter(date >= planting_date) %>% 
  mutate(dap = as.numeric(date - planting_date)) %>% 
  group_by(environment, datatype) %>% 
  mutate(interval = rep(1:12, each = 30, length.out = n())) %>%
  ungroup() %>%
  ## Add growth stage information
  left_join(one_year_growth_staging) %>%
  select(environment:interval, growth_stage2)



## Calculate summary statistics for temperature and precip

## First calculate some basic summary statistics:
# monthly mean of daily max temperature
# monthly mean of daily min temperature
# average monthly temperature ([Tmax + Tmin] / 2)
# total monthly precip

daily_stats_temp <- oneyear_trial_env_data_interval %>% 
  select(environment, dap, interval, growth_stage2, datatype, value) %>% 
  filter(datatype != "PRCP") %>% 
  spread(datatype, value) %>% 
  mutate(TAVG = (TMAX + TMIN) / 2,
         TRANGE = TMAX - TMIN)


## For each  statistic, summarize by interval or growth stage; then combine

summ_stats_temp_interval <- daily_stats_temp %>%
  group_by(environment, interval) %>% 
  summarize_at(vars(TMAX, TMIN, TAVG, TRANGE), list(mean = mean, sd = sd), na.rm = TRUE) %>% 
  rename_at(vars(contains("mean")), ~str_remove(., "_mean")) %>%
  select(environment:TRANGE, TSEASON = TAVG_sd) %>%
  gather(variable, value, -environment, -interval)

summ_stats_temp_stage <- daily_stats_temp %>%
  filter(!is.na(growth_stage2)) %>%
  group_by(environment, growth_stage2) %>% 
  summarize_at(vars(TMAX, TMIN, TAVG, TRANGE), list(mean = mean, sd = sd), na.rm = TRUE) %>% 
  rename_at(vars(contains("mean")), ~str_remove(., "_mean")) %>%
  select(environment:TRANGE, TSEASON = TAVG_sd) %>%
  gather(variable, value, -environment, -growth_stage2)



summ_stats_prcp_interval <- oneyear_trial_env_data_interval %>% 
  select(environment, dap, interval, datatype, value) %>% 
  filter(datatype == "PRCP") %>% 
  group_by(environment, interval) %>% 
  summarize(PPT = sum(value, na.rm = TRUE)) %>% 
  gather(variable, value, -environment, -interval)

summ_stats_prcp_stage <- oneyear_trial_env_data_interval %>% 
  select(environment, dap, growth_stage2, datatype, value) %>% 
  filter(!is.na(growth_stage2)) %>%
  filter(datatype == "PRCP") %>% 
  group_by(environment, growth_stage2) %>% 
  summarize(PPT = sum(value, na.rm = TRUE)) %>% 
  gather(variable, value, -environment, -growth_stage2)


# Combine
interval_summ_stats <- ls(pattern = "^summ[_a-z]*interval$") %>% 
  map_df(get) %>% 
  ungroup() %>%
  mutate(interval = str_c("interval_", interval)) %>% 
  unite(variable, interval, variable, sep = "_")

stage_summ_stats <- ls(pattern = "^summ[_a-z]*stage$") %>% 
  map_df(get) %>% 
  ungroup() %>%
  unite(variable, growth_stage2, variable, sep = "_")


## Calculate secondary summary statistics
# annual mean temperature (average of monthly temperatures; we will rename this growing annual mean temperature
# annual mean diurnal range (average of the difference between monthly mean min and max temperatures; we will rename this growing annual mean diurnal range)
# Max temperature of warmest month (max TMAX)
# Min temperature of coldest month (min TMIN)
# Annual temperature range (max(TMAX) - min(TMIN); renamed growing annual temperature range)
# isothermality (annual mean diurnal range / Annual temperature range) * 100
# temperature seasonality (standard deviation of TAVG)
# Annual precipitation (sum of PPT; renamed growing annual precipitation)
# Precipitation of wettest month (max PPT)
# Precipitation of driest month (min PPT)
            
## Again, calculate using intervals and then using growth stages

# Annual mean temperature
annual_TAVG_interval <- summ_stats_temp_interval %>% 
 filter(variable == "TAVG") %>% 
 group_by(environment) %>% 
 summarize(annual_TAVG = mean(value, na.rm = TRUE), 
           annual_TSEASON = sd(value, na.rm = TRUE)) %>%
 gather(variable, value, -environment)

annual_TAVG_stage <- summ_stats_temp_stage %>% 
  filter(variable == "TAVG") %>% 
  # Exclude pre-emergence and maturity
  filter(growth_stage2 %in% c("vegetative", "flowering", "grainfill")) %>%
  group_by(environment) %>% 
  summarize(annual_TAVG = mean(value, na.rm = TRUE), 
            annual_TSEASON = sd(value, na.rm = TRUE)) %>%
  gather(variable, value, -environment)

## Annual mean max, min, and range
annual_means_interval <- summ_stats_temp_interval %>% 
  filter(variable %in% c("TMAX", "TMIN", "TRANGE")) %>%
  spread(variable, value) %>%
  group_by(environment) %>% 
  summarize_at(vars(TMAX, TMIN, TRANGE), mean, na.rm = TRUE) %>%
  rename_at(vars(contains("T", ignore.case = FALSE)), ~str_c("annual_", .)) %>%
  gather(variable, value, -environment)

annual_means_stage <- summ_stats_temp_stage %>% 
  filter(variable %in% c("TMAX", "TMIN", "TRANGE")) %>%
  filter(growth_stage2 %in% c("vegetative", "flowering", "grainfill")) %>%
  spread(variable, value) %>%
  group_by(environment) %>% 
  summarize_at(vars(TMAX, TMIN, TRANGE), mean, na.rm = TRUE) %>%
  rename_at(vars(contains("T", ignore.case = FALSE)), ~str_c("annual_", .)) %>%
  gather(variable, value, -environment)

## Yearly max, yearly min, yearly range
annual_range_min_max_interval <- summ_stats_temp_interval %>% 
  filter(variable %in% c("TMAX", "TMIN")) %>%
  spread(variable, value) %>% 
  group_by(environment) %>% 
  summarize(max_TMAX = max(TMAX, na.rm = TRUE), min_TMIN = min(TMIN, na.rm = TRUE)) %>% 
  mutate(annual_TRANGE_max = max_TMAX - min_TMIN) %>% 
  gather(variable, value, -environment)

annual_range_min_max_stage <- summ_stats_temp_stage %>% 
  filter(variable %in% c("TMAX", "TMIN")) %>%
  spread(variable, value) %>% 
  group_by(environment) %>% 
  summarize(max_TMAX = max(TMAX, na.rm = TRUE), min_TMIN = min(TMIN, na.rm = TRUE)) %>% 
  mutate(annual_TRANGE_max = max_TMAX - min_TMIN) %>% 
  gather(variable, value, -environment)


## Isothermality
isotherm_interval <- bind_rows(filter(annual_means_interval, variable == "annual_TRANGE"), 
                               filter(annual_range_min_max_interval, variable == "annual_TRANGE_max")) %>% 
  spread(variable, value) %>% 
  group_by(environment) %>% 
  summarize(isothermality = (annual_TRANGE / annual_TRANGE_max) * 100) %>%
  gather(variable, value, -environment)

isotherm_stage <- bind_rows(filter(annual_means_stage, variable == "annual_TRANGE"), 
                            filter(annual_range_min_max_stage, variable == "annual_TRANGE_max")) %>% 
  spread(variable, value) %>% 
  group_by(environment) %>% 
  summarize(isothermality = (annual_TRANGE / annual_TRANGE_max) * 100) %>%
  gather(variable, value, -environment)


# Precipitation
annual_prcp_interval <- summ_stats_prcp_interval %>% 
  group_by(environment) %>% 
  summarize_at(vars(value), list(sum = sum, min = min, max = max), na.rm = TRUE) %>%
  select(environment, annual_PPT = sum, min_PPT = min, max_PPT = max) %>%
  gather(variable, value, -environment)

annual_prcp_stage <- summ_stats_prcp_stage %>% 
  group_by(environment) %>% 
  summarize_at(vars(value), list(sum = sum, min = min, max = max), na.rm = TRUE) %>%
  select(environment, annual_PPT = sum, min_PPT = min, max_PPT = max) %>%
  gather(variable, value, -environment)




## Create photothermal time stats - based on AGDD
phototherm_stats <- one_year_growth_staging %>% 
  mutate(day = yday(date)) %>%
  select(environment, agdd, day, dap) %>% 
  # Add location
  left_join(., distinct(trial_info, environment, location)) %>% 
  # Add daylength
  left_join(., daylength_loc_df) %>% 
  mutate(phototherm = agdd * hours) %>% 
  select(environment, dap, phototherm)



# Combine the oneyear data
oneyear_summary_df_interval <- map_df(c("coord_clim_summ", "interval_summ_stats", ls(pattern = "(annual|isotherm)[_a-zA-Z]*interval$")), get)

oneyear_summary_df_stage <- map_df(c("coord_clim_summ", "stage_summ_stats", ls(pattern = "(annual|isotherm)[_a-zA-Z]*stage$")), get)







### 10-year average

### One-year environmental covariates
multiyear_env_data_unnest <- noaa_trial_data_multiyear_complete %>%
  unnest(data)

## Find the average planting date per location
location_planting_date <- trial_info %>% 
  mutate(mean_planting_day = yday(ymd(planting_date))) %>% 
  group_by(location) %>% 
  summarize_at(vars(mean_planting_day), ~floor(mean(.)))


multiyear_env_data <- trial_info %>% 
  distinct(environment, location, year) %>%
  left_join(., location_planting_date) %>%
  rename(env_year = year) %>%
  left_join(., multiyear_env_data_unnest) %>%
  mutate(day = yday(ymd(date)),
         value = value / 10)


## Subset data only after planting, then add 30-day intervals
multiyear_env_data_interval <- multiyear_env_data %>% 
  select(environment, env_year, mean_planting_day, day, month, year, datatype, value) %>% 
  filter(day >= mean_planting_day) %>% 
  mutate(dap = as.numeric(day - mean_planting_day)) %>% 
  group_by(environment, year, datatype) %>% 
  mutate(interval = rep(1:12, each = 30, length.out = n())) %>%
  ungroup() %>%
  ## Add growth stage information
  left_join(multiyear_growth_staging) %>%
  select(environment:interval, growth_stage2)
  


## Proportion of missing data by year
prop_NA <- multiyear_env_data_interval %>% 
  group_by(environment, env_year, year, datatype) %>% 
  summarize(prop_NA = mean(is.na(value))) %>%
  arrange(desc(prop_NA))

## By environment, how many years are removed?
prop_NA %>% 
  group_by(environment, datatype) %>% 
  summarize(n_removed = sum(prop_NA > 0.20)) %>% 
  arrange(desc(n_removed))

obs_tokeep <- prop_NA %>% 
  filter(prop_NA <= 0.20) %>%
  ungroup() %>%
  select(-prop_NA)


## Remove location-datatype-year combinations with too much missing data. This should be no more than 3 years
multiyear_env_data_interval1 <- left_join(obs_tokeep, multiyear_env_data_interval)


## Calculate summary statistics for temperature and precip
daily_stats_temp <- multiyear_env_data_interval1 %>% 
  select(environment, env_year, year, dap, interval, growth_stage2, datatype, value) %>% 
  filter(datatype != "PRCP") %>% 
  spread(datatype, value) %>% 
  mutate(TAVG = (TMAX + TMIN) / 2,
         TRANGE = TMAX - TMIN)


  
  
summ_stats_temp_interval <- daily_stats_temp %>%
  group_by(environment, env_year, year, interval) %>% 
  summarize_at(vars(TMAX, TMIN, TAVG, TRANGE), list(mean = mean, sd = sd), na.rm = TRUE) %>% 
  rename_at(vars(contains("mean")), ~str_remove(., "_mean")) %>%
  select(environment:TRANGE, TSEASON = TAVG_sd) %>%
  gather(variable, value, -environment, -interval, -env_year, -year)

summ_stats_temp_stage <- daily_stats_temp %>%
  group_by(environment, env_year, year, growth_stage2) %>% 
  summarize_at(vars(TMAX, TMIN, TAVG, TRANGE), list(mean = mean, sd = sd), na.rm = TRUE) %>% 
  rename_at(vars(contains("mean")), ~str_remove(., "_mean")) %>%
  select(environment:TRANGE, TSEASON = TAVG_sd) %>%
  gather(variable, value, -environment, -growth_stage2, -env_year, -year)


## Precipitation summary
summ_stats_prcp_interval <- multiyear_env_data_interval1 %>% 
  select(environment, env_year, year, dap, interval, datatype, value) %>% 
  filter(datatype == "PRCP") %>% 
  group_by(environment, env_year, year, interval) %>% 
  summarize(PPT = sum(value, na.rm = TRUE)) %>% 
  gather(variable, value, -environment, -interval, -env_year, -year)

summ_stats_prcp_stage <- multiyear_env_data_interval1 %>% 
  select(environment, env_year, year, dap, growth_stage2, datatype, value) %>% 
  filter(datatype == "PRCP") %>% 
  group_by(environment, env_year, year, growth_stage2) %>% 
  summarize(PPT = sum(value, na.rm = TRUE)) %>% 
  gather(variable, value, -environment, -growth_stage2, -env_year, -year)

# Combine
interval_summ_stats <- ls(pattern = "summ_stats[_a-z]*_interval$") %>%
  map_df(get) %>%
  mutate(interval = str_c("interval_", interval)) %>% 
  unite(variable, interval, variable, sep = "_") %>%
  ungroup()

stage_summ_stats <- ls(pattern = "summ_stats[_a-z]*_stage$") %>%
  map_df(get) %>%
  unite(variable, growth_stage2, variable, sep = "_") %>%
  ungroup()




## Calculate secondary summary statistics
# annual mean temperature (average of monthly temperatures; we will rename this growing annual mean temperature
# annual mean diurnal range (average of the difference between monthly mean min and max temperatures; we will rename this growing annual mean diurnal range)
# Max temperature of warmest month (max TMAX)
# Min temperature of coldest month (min TMIN)
# Annual temperature range (max(TMAX) - min(TMIN); renamed growing annual temperature range)
# isothermality (annual mean diurnal range / Annual temperature range) * 100
# temperature seasonality (standard deviation of TAVG)
# Annual precipitation (sum of PPT; renamed growing annual precipitation)
# Precipitation of wettest month (max PPT)
# Precipitation of driest month (min PPT)

  

# Annual mean temperature
annual_TAVG_interval <- summ_stats_temp_interval %>% 
  filter(variable == "TAVG") %>% 
  group_by(environment, env_year, year) %>% 
  summarize(annual_TAVG = mean(value, na.rm = TRUE), 
            annual_TSEASON = sd(value, na.rm = TRUE)) %>%
  gather(variable, value, annual_TAVG, annual_TSEASON)

annual_TAVG_stage <- summ_stats_temp_stage %>% 
  filter(variable == "TAVG") %>% 
  filter(growth_stage2 %in% c("vegetative", "flowering", "grainfill")) %>%
  group_by(environment, env_year, year) %>% 
  summarize(annual_TAVG = mean(value, na.rm = TRUE), 
            annual_TSEASON = sd(value, na.rm = TRUE)) %>%
  gather(variable, value, annual_TAVG, annual_TSEASON)



## Annual mean max, min, and range
annual_means_interval <- summ_stats_temp_interval %>% 
  filter(variable %in% c("TMAX", "TMIN", "TRANGE")) %>%
  spread(variable, value) %>%
  group_by(environment, env_year, year) %>% 
  summarize_at(vars(TMAX, TMIN, TRANGE), mean, na.rm = TRUE) %>%
  rename_at(vars(contains("T", ignore.case = FALSE)), ~str_c("annual_", .)) %>%
  gather(variable, value, annual_TMAX, annual_TMIN, annual_TRANGE)

annual_means_stage <- summ_stats_temp_stage %>% 
  filter(variable %in% c("TMAX", "TMIN", "TRANGE")) %>%
  spread(variable, value) %>%
  group_by(environment, env_year, year) %>% 
  summarize_at(vars(TMAX, TMIN, TRANGE), mean, na.rm = TRUE) %>%
  rename_at(vars(contains("T", ignore.case = FALSE)), ~str_c("annual_", .)) %>%
  gather(variable, value, annual_TMAX, annual_TMIN, annual_TRANGE)


## Annual min, max, range
annual_range_min_max_interval <- summ_stats_temp_interval %>% 
  filter(variable %in% c("TMAX", "TMIN")) %>%
  spread(variable, value) %>% 
  group_by(environment, env_year, year) %>% 
  summarize(max_TMAX = max(TMAX, na.rm = TRUE), min_TMIN = min(TMIN, na.rm = TRUE)) %>% 
  filter_at(vars(max_TMAX, min_TMIN), all_vars(. != "-Inf")) %>%
  mutate(annual_TRANGE_max = max_TMAX - min_TMIN) %>% 
  gather(variable, value, max_TMAX, min_TMIN, annual_TRANGE_max)

annual_range_min_max_stage <- summ_stats_temp_stage %>% 
  filter(variable %in% c("TMAX", "TMIN")) %>%
  spread(variable, value) %>% 
  group_by(environment, env_year, year) %>% 
  summarize(max_TMAX = max(TMAX, na.rm = TRUE), min_TMIN = min(TMIN, na.rm = TRUE)) %>% 
  filter_at(vars(max_TMAX, min_TMIN), all_vars(. != "-Inf")) %>%
  mutate(annual_TRANGE_max = max_TMAX - min_TMIN) %>% 
  gather(variable, value, max_TMAX, min_TMIN, annual_TRANGE_max)


## Isothermality
isotherm_interval <- bind_rows(filter(annual_means_interval, variable == "annual_TRANGE"), 
                               filter(annual_range_min_max_interval, variable == "annual_TRANGE_max")) %>% 
  spread(variable, value) %>% 
  group_by(environment, env_year, year) %>% 
  summarize(isothermality = (annual_TRANGE / annual_TRANGE_max) * 100) %>%
  gather(variable, value, isothermality)

isotherm_stage <- bind_rows(filter(annual_means_interval, variable == "annual_TRANGE"), 
                            filter(annual_range_min_max_stage, variable == "annual_TRANGE_max")) %>% 
  spread(variable, value) %>% 
  group_by(environment, env_year, year) %>% 
  summarize(isothermality = (annual_TRANGE / annual_TRANGE_max) * 100) %>%
  gather(variable, value, isothermality)

# Precipitation
annual_prcp_interval <- summ_stats_prcp_interval %>% 
  group_by(environment, env_year, year) %>% 
  summarize_at(vars(value), list(sum = sum, min = min, max = max), na.rm = TRUE) %>%
  select(environment, env_year, year, annual_PPT = sum, min_PPT = min, max_PPT = max) %>%
  gather(variable, value, annual_PPT, min_PPT, max_PPT)

annual_prcp_stage <- summ_stats_prcp_stage %>% 
  group_by(environment, env_year, year) %>% 
  summarize_at(vars(value), list(sum = sum, min = min, max = max), na.rm = TRUE) %>%
  select(environment, env_year, year, annual_PPT = sum, min_PPT = min, max_PPT = max) %>%
  gather(variable, value, annual_PPT, min_PPT, max_PPT)



## Create photothermal time stats - based on AGDD
phototherm_stats <- multiyear_growth_staging %>% 
  select(environment, avg_agdd, mean_planting_day, dap) %>% 
  # Add location
  left_join(., distinct(trial_info, environment, location)) %>% 
  # Add daylength
  left_join(., daylength_loc_df, by= c("location", "mean_planting_day" = "day")) %>% 
  mutate(avg_phototherm = avg_agdd * hours) %>% 
  select(environment, dap, avg_phototherm)



# Combine the oneyear data
multiyear_summary_df_interval <- c("coord_clim_summ", "interval_summ_stats", ls(pattern = "(annual|isotherm)[_a-zA-Z]*interval$")) %>%
  map_df(get) %>%
  ## Take the average variable for each environment
  group_by(environment, variable) %>%
  summarize(value = mean(value, na.rm = TRUE)) %>%
  ungroup()
  

multiyear_summary_df_stage <- c("coord_clim_summ", "stage_summ_stats", ls(pattern = "(annual|isotherm)[_a-zA-Z]*stage$")) %>%
  map_df(get) %>%
  ## Take the average variable for each environment
  group_by(environment, variable) %>%
  summarize(value = mean(value, na.rm = TRUE)) %>%
  ungroup()


  

## Combine climate and soil data

# Create a df to store everything

environ_covariate_data <- tibble(
  timeframe = rep(c("interval", "growth_stage"), each = 2), 
  ec_group = rep(c("oneyear", "multiyear"), 2) ) %>%
  mutate(data = list(
    bind_rows(oneyear_summary_df_interval, soil_data_complete1),
    bind_rows(multiyear_summary_df_interval, soil_data_complete1),
    bind_rows(oneyear_summary_df_stage, soil_data_complete1),
    bind_rows(multiyear_summary_df_stage, soil_data_complete1)
  )) %>%
  # Replace NA with the variable mean across environments
  mutate(data = map(data, ~group_by(., variable) %>% mutate(value = ifelse(is.na(value), mean(value, na.rm = TRUE), value)) %>% 
                      ungroup() %>% arrange(environment, variable) ))


# Save the data
save_file <- file.path(data_dir, "environmental_data_compiled.RData")
save("environ_covariate_data", file = save_file)

