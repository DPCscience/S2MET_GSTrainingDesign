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

## Also calculate GDD according to the growth model here:  https://ndawn.ndsu.nodak.edu/help-barley-growing-degree-days.html

### One-year environmental covariates
oneyear_env_data_unnest <- noaa_trial_data_oneyear_complete %>%
  unnest(data)

oneyear_trial_env_data <- trial_info %>% 
  distinct(environment, location, year) %>%
  mutate_all(parse_guess) %>%
  left_join(., oneyear_env_data_unnest) %>%
  # Convert tenth of value to value
  mutate(value = value / 10)

# Summarize longitude/latitude
coord_clim_summ <- oneyear_trial_env_data %>% 
  select(environment, latitude, longitude, elevation) %>% 
  distinct() %>%
  gather(variable, value, -environment)


## Calculate summary statistics for temperature and precip

## First calculate some basic summary statistics:
# monthly mean of daily max temperature
# monthly mean of daily min temperature
# average monthly temperature ([Tmax + Tmin] / 2)
# total monthly precip

summ_stats_temp <- oneyear_trial_env_data %>% 
  select(environment, date, month, datatype, value) %>% 
  filter(datatype != "PRCP") %>% 
  spread(datatype, value) %>% 
  group_by(environment, month) %>% 
  summarize_at(vars(TMAX, TMIN), mean, na.rm = TRUE) %>% 
  mutate(TAVG = (TMAX + TMIN) / 2) %>%
  gather(variable, value, -environment, -month)

summ_stats_prcp <- oneyear_trial_env_data %>% 
  select(environment, date, month, datatype, value) %>% 
  filter(datatype == "PRCP") %>% 
  group_by(environment, month) %>% 
  summarize(PPT = sum(value, na.rm = TRUE)) %>% 
  gather(variable, value, -environment, -month)

# Combine
mon_summ_stats <- bind_rows(summ_stats_temp, summ_stats_prcp) %>%
  mutate(month = str_c("month_", month)) %>% 
  unite(variable, month, variable, sep = "_")


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
annual_TAVG <- summ_stats_temp %>% 
 filter(variable == "TAVG") %>% 
 group_by(environment) %>% 
 summarize(annual_TAVG = mean(value, na.rm = TRUE), temp_seasonality = sd(value, na.rm = TRUE)) %>%
 gather(variable, value, -environment)

annual_diurnal_range <- summ_stats_temp %>% 
 filter(variable != "TAVG") %>% 
 spread(variable, value) %>% 
 group_by(environment) %>% 
 mutate(diurnal_range = TMAX - TMIN) %>% 
 summarize(annual_diurnal_range = mean(diurnal_range, na.rm = TRUE)) %>%
 gather(variable, value, -environment)

range_min_max <- summ_stats_temp %>% 
 filter(variable != "TAVG") %>% 
 spread(variable, value) %>% 
 group_by(environment) %>% 
 summarize(max_TMAX = max(TMAX, na.rm = TRUE), min_TMIN = min(TMIN, na.rm = TRUE)) %>% 
 mutate(annual_temperature_range = max_TMAX - min_TMIN) %>% 
 gather(variable, value, -environment)

isotherm <- bind_rows(annual_diurnal_range, 
                     filter(range_min_max, variable == "annual_temperature_range")) %>% 
 spread(variable, value) %>% 
 group_by(environment) %>% 
 summarize(isothermality = (annual_diurnal_range / annual_temperature_range) * 100) %>%
 gather(variable, value, -environment)

# Precipitation
annual_prcp <- summ_stats_prcp %>% 
 group_by(environment) %>% 
 summarize(annual_precipitation = sum(value, na.rm = TRUE), max_PPT = max(value, na.rm = TRUE),
           min_PPT = min(value, na.rm = TRUE)) %>% 
 gather(variable, value, -environment)


## Calculate growing degree days (GDD) based on: https://ndawn.ndsu.nodak.edu/help-barley-growing-degree-days.html
# Calculate the mean daily temperature
# Use a loess smooth to impute missing temperature data
one_year_daily_temp <- oneyear_trial_env_data %>% 
  select(environment, date, month, datatype, value) %>% 
  filter(datatype != "PRCP") %>%
  mutate(day = yday(date))

one_year_daily_temp_impute <- one_year_daily_temp %>%
  arrange(environment, datatype, date) %>%
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
        mutate(pred_value = window_mean(x = day, y = value, window = wind)) %>%
        mutate(value = ifelse(is.na(value), pred_value, value)) %>%
        select(date, month, day, value)
      
      n_missing <- sum(is.na(df1$value))
      wind <- wind + 1
      
    }
    
    # Remove grouping variables
    df2 <- df1[,c("date", "month", "value", "day")]
    
    # Return the data.frame
    data_frame(imputed_data = list(df2), window = wind - 1, n_missing = n_missing_start)
    
  }) %>% ungroup()



gdd_stats <- one_year_daily_temp_impute %>% 
  unnest() %>% 
  select(-window, -n_missing) %>%
  spread(datatype, value) %>%
  # Convert from tenths of C to F
  mutate_at(vars(TMIN, TMAX), ~(. * (9/5)) + 32) %>%
  # If the daily min or max is < 32 F, it is set to 32
  mutate_at(vars(TMIN, TMAX), ~ifelse(. < 32, 32, .)) %>%
  mutate(TAVG = (TMIN + TMAX) / 2,
         GDD = TAVG - 32)

gdd_stats1 <- gdd_stats %>% 
  select(environment, day, GDD)



## Use planting date to calculate accumulated growing degree days
agdd_stats <- gdd_stats %>% 
  left_join(., distinct(trial_info, environment, planting_date)) %>%
  # Is the recorded date before the planting date?
  mutate(date = ymd(date), planting_date = ymd(planting_date), 
         planted = date >= planting_date) %>% 
  filter(planted) %>% 
  group_by(environment) %>% 
  mutate(AGDD = cumsum(GDD)) %>%
  ungroup()

# The AGDD to reach Haun stage 2.0
haun2 <- 384

## Adjust the accumulated GDD using the restrictions found in the above website
agdd_stats_adj <- agdd_stats %>%
  mutate(TMAX = ifelse(AGDD < haun2 & TMAX > 70, 70, TMAX),
         TMAX = ifelse(AGDD >= haun2 & TMAX > 95, 95, TMAX),
         TAVG = (TMIN + TMAX) / 2,
         GDD = TAVG - 32) %>%
  group_by(environment) %>% 
  mutate(AGDD = cumsum(GDD)) %>%
  ungroup()

# Tidy
agdd_stats_adj1 <- agdd_stats_adj %>%
  mutate(day = yday(date)) %>%
  select(environment, day, AGDD) %>%
  # Convert day to days from planting
  group_by(environment) %>%
  mutate(day_from_planting = seq(n_distinct(day))) %>%
  ungroup() %>%
  select(environment, day_from_planting, names(.), -day)


## Create photothermal time stats - based on AGDD
phototherm_stats <- agdd_stats_adj %>% 
  mutate(day = yday(date)) %>% 
  select(environment, AGDD, day) %>% 
  left_join(., distinct(trial_info, environment, location)) %>% 
  left_join(., daylength_loc_df) %>% 
  mutate(phototherm = AGDD * hours) %>% 
  select(environment, day, phototherm)

# Tidy
phototherm_stats1 <- phototherm_stats %>%
  select(environment, day, phototherm) %>%
  # Convert day to days from planting
  group_by(environment) %>%
  mutate(day_from_planting = seq(n_distinct(day))) %>%
  ungroup() %>%
  select(environment, day_from_planting, names(.), -day)


## Now create photothermal time information based on GDD
phototherm_stats_gdd <- gdd_stats %>%
  mutate(day = yday(date)) %>% 
  select(environment, GDD, day) %>%
  left_join(., distinct(trial_info, environment, location)) %>% 
  left_join(., daylength_loc_df) %>% 
  mutate(phototherm = GDD * hours) %>% 
  select(environment, day, phototherm)

# Tidy
phototherm_stats_gdd1 <- phototherm_stats_gdd %>%
  select(environment, day, phototherm_daily = phototherm)




# Combine the oneyear data
one_year_summary_df <- bind_rows(
  mon_summ_stats, annual_TAVG, 
  annual_diurnal_range, 
  range_min_max, 
  isotherm, 
  annual_prcp) %>% ungroup()

# Create a daily summary
one_year_daily_summary <- list(
  agdd = agdd_stats_adj1,
  phototherm = phototherm_stats1,
  gdd = gdd_stats1,
  phototherm_gdd = phototherm_stats_gdd1
)


# Set the maximum number of days in the interval length
max_days <- 90

## For each environment and daily covariable, build intervals of the average variable
one_year_daily_summary_interval <- one_year_daily_summary %>%
  map(~{
    df <- .
    
    # Find the minimum day among environments.
    min_day <- group_by(df, environment) %>% 
      summarize_at(vars(contains("day")), min) %>% 
      pull() %>% 
      max()
    
    # Build the intervals
    intervals_df <- crossing(begin = seq(min_day, min_day + max_days), end = seq(min_day + max_days)) %>% 
      filter(begin <= end)
    
    # For each environment, convert the day to starting from 1 (planting day)
    df1 <- df %>% 
      group_by(environment)
    
    # Day name
    dayname <- names(df1)[2]
    names(df1)[2] <- "day"
    
    # Iterate over the intervals
    intervals_summ <- pmap(intervals_df, ~filter(df1, between(day, .x, .y)) %>% 
                             summarize_at(vars(-environment, -day), mean))
    
    
    mutate(intervals_df, variable = tail(names(df), 1), out = intervals_summ)
    
  })










### 10-year average

### One-year environmental covariates
multiyear_env_data_unnest <- noaa_trial_data_multiyear_complete %>%
  unnest(data)

multiyear_trial_env_data <- trial_info %>% 
  distinct(environment, location, year) %>%
  mutate_all(parse_guess) %>%
  rename(env_year = year) %>%
  left_join(., multiyear_env_data_unnest) %>%
  # Convert tenth of value to value
  mutate(value = value / 10)

# Summarize longitude/latitude
coord_clim_summ <- multiyear_trial_env_data %>% 
  select(environment, latitude, longitude, elevation) %>% 
  distinct() %>%
  gather(variable, value, -environment)


## Proportion of missing data by year
prop_NA <- multiyear_trial_env_data %>% 
  group_by(environment, env_year, year, datatype) %>% 
  summarize(prop_NA = mean(is.na(value)))

obs_tokeep <- prop_NA %>% 
  filter(prop_NA <= 0.20) %>%
  ungroup() %>%
  select(-prop_NA)


## Remove location-datatype-year combinations with too much missing data. This should be no more than 3 years
multiyear_trial_env_data1 <- left_join(obs_tokeep, multiyear_trial_env_data)


## Calculate summary statistics for temperature and precip


summ_stats_temp <- multiyear_trial_env_data1 %>% 
  filter(datatype != "PRCP") %>% 
  spread(datatype, value) %>% 
  group_by(environment, env_year, year, month) %>% 
  summarize_at(vars(TMAX, TMIN), mean, na.rm = TRUE) %>% 
  mutate(TAVG = (TMAX + TMIN) / 2) %>%
  gather(variable, value, -environment, -env_year, -year, -month)

summ_stats_prcp <- multiyear_trial_env_data1 %>% 
  filter(datatype == "PRCP") %>% 
  group_by(environment, env_year, year, month) %>%  
  summarize(PPT = sum(value, na.rm = TRUE)) %>% 
  gather(variable, value, -environment, -env_year, -year, -month)

# Combine
mon_summ_stats <- bind_rows(summ_stats_temp, summ_stats_prcp) %>%
  mutate(month = str_c("month_", month)) %>% 
  unite(variable, month, variable, sep = "_")




# Calculate secondary summary statistics
  
# Annual mean temperature
annual_TAVG <- summ_stats_temp %>% 
  filter(variable == "TAVG") %>% 
  group_by(environment, env_year, year) %>% 
  summarize(annual_TAVG = mean(value, na.rm = TRUE), temp_seasonality = sd(value, na.rm = TRUE)) %>%
  gather(variable, value, -environment:-year)

annual_diurnal_range <- summ_stats_temp %>% 
  filter(variable != "TAVG") %>% 
  spread(variable, value) %>% 
  group_by(environment, env_year, year) %>% 
  mutate(diurnal_range = TMAX - TMIN) %>% 
  summarize(annual_diurnal_range = mean(diurnal_range, na.rm = TRUE)) %>%
  gather(variable, value, -environment:-year)

range_min_max <- summ_stats_temp %>% 
  filter(variable != "TAVG") %>% 
  spread(variable, value) %>% 
  group_by(environment, env_year, year) %>% 
  summarize(max_TMAX = max(TMAX, na.rm = TRUE), min_TMIN = min(TMIN, na.rm = TRUE)) %>% 
  mutate_at(vars(max_TMAX, min_TMIN), ~ifelse(is.infinite(.), NA, .)) %>%
  mutate(annual_temperature_range = max_TMAX - min_TMIN) %>% 
  gather(variable, value, -environment:-year)

isotherm <- bind_rows(annual_diurnal_range, 
                      filter(range_min_max, variable == "annual_temperature_range")) %>% 
  spread(variable, value) %>% 
  group_by(environment, env_year, year) %>% 
  summarize(isothermality = (annual_diurnal_range / annual_temperature_range) * 100) %>%
  gather(variable, value, -environment:-year)

# Precipitation
annual_prcp <- summ_stats_prcp %>% 
  group_by(environment, env_year, year) %>% 
  summarize(annual_precipitation = sum(value, na.rm = TRUE), 
            max_PPT = max(value, na.rm = TRUE), min_PPT = min(value, na.rm = TRUE)) %>% 
  gather(variable, value, -environment:-year)


## Calculate growing degree days (GDD) based on: https://ndawn.ndsu.nodak.edu/help-barley-growing-degree-days.html
# Calculate the mean daily temperature
# Use a loess smooth to impute missing temperature data
multi_year_daily_temp <- multiyear_trial_env_data1 %>% 
  select(environment, env_year, date, year, month, datatype, value) %>% 
  filter(datatype != "PRCP") %>%
  mutate(day = yday(date))


multi_year_daily_temp_impute <- multi_year_daily_temp %>%
  # filter(environment == "STP16") %>%
  arrange(environment, year, datatype, date) %>%
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
        mutate(pred_value = window_mean(x = day, y = value, window = wind)) %>%
        mutate(value = ifelse(is.na(value), pred_value, value)) %>%
        select(env_year, date, month, day, value)
        
      n_missing <- sum(is.na(df1$value))
      wind <- wind + 1
    
    }
    
    # Remove grouping variables
    df2 <- df1[,c("env_year", "date", "month", "value", "day")]
    
    # Return the data.frame
    data_frame(imputed_data = list(df2), window = wind - 1, n_missing = n_missing_start)
    
  }) %>% ungroup()



gdd_stats <- multi_year_daily_temp_impute %>%
  unnest(imputed_data) %>%
  select(-window, -n_missing) %>%
  spread(datatype, value) %>%
  # Convert from tenths of C to F
  mutate_at(vars(TMIN, TMAX), ~(. * (9/5)) + 32) %>%
  # If the daily min or max is < 32 F, it is set to 32
  mutate_at(vars(TMIN, TMAX), ~ifelse(. < 32, 32, .)) %>%
  mutate(TAVG = (TMIN + TMAX) / 2,
         GDD = TAVG - 32)

gdd_stats1 <- gdd_stats %>% 
  select(environment, env_year, year, day, GDD)




## Use planting date to calculate accumulated growing degree days
# Average these AGDDs over years
agdd_stats <- gdd_stats %>% 
  left_join(., distinct(trial_info, environment, planting_date)) %>%
  # Is the recorded yday before or after the planting yday?
  mutate(date = ymd(date), planting_date = ymd(planting_date), 
         planted = day >= yday(planting_date)) %>% 
  filter(planted) %>% 
  group_by(environment, year) %>% 
  mutate(AGDD = cumsum(GDD)) %>%
  ungroup()

# The AGDD to reach Haun stage 2.0
haun2 <- 384

## Adjust the accumulated GDD using the restrictions found in the above website
agdd_stats_adj <- agdd_stats %>%
  mutate(TMAX = ifelse(AGDD < haun2 & TMAX > 70, 70, TMAX),
         TMAX = ifelse(AGDD >= haun2 & TMAX > 95, 95, TMAX),
         TAVG = (TMIN + TMAX) / 2,
         GDD = TAVG - 32) %>%
  group_by(environment, year) %>% 
  mutate(AGDD = cumsum(GDD)) %>%
  ungroup()

# Tidy
agdd_stats_adj1 <- agdd_stats_adj %>%
  mutate(day = yday(date), plant_day = yday(planting_date)) %>%
  select(environment, env_year, year, day, plant_day, AGDD) %>%
  # Convert day to days from planting
  group_by(environment, year) %>%
  mutate(day_from_planting = day - plant_day) %>%
  ungroup() %>%
  select(environment, day_from_planting, names(.), -day, -plant_day)






## Now create photothermal time information based on GDD
phototherm_stats_gdd <- gdd_stats %>%
  select(environment, env_year, year, GDD, day) %>%
  left_join(., distinct(trial_info, environment, location)) %>% 
  left_join(., distinct(daylength_loc_df, location, day, hours)) %>% 
  filter(!is.na(hours)) %>%
  mutate(phototherm = GDD * hours) %>% 
  select(environment, env_year, year, day, phototherm)

# Tidy
phototherm_stats_gdd1 <- phototherm_stats_gdd %>%
  select(environment, env_year, year, day, phototherm_daily = phototherm)


## Now create photothermal time information based on AGDD
phototherm_stats_agdd <- agdd_stats_adj %>%
  select(environment, env_year, year, AGDD, planting_date, day) %>%
  left_join(., distinct(trial_info, environment, location), by = "environment") %>% 
  left_join(., distinct(daylength_loc_df, location, day, hours), by = c("day", "location")) %>% 
  filter(!is.na(hours)) %>%
  mutate(phototherm = AGDD * hours) %>% 
  select(environment, env_year, year, day, planting_date, phototherm)

# Tidy
phototherm_stats_agdd1 <- phototherm_stats_agdd %>%
  mutate(plant_day = yday(planting_date)) %>%
  # Convert day to days from planting
  group_by(environment, year) %>%
  mutate(day_from_planting = day - plant_day) %>%
  ungroup() %>%
  select(environment, day_from_planting, names(.), -day, -plant_day, -planting_date)



  
# Combine the monthly and annual data
multi_year_summary_df <- bind_rows(mon_summ_stats, annual_TAVG, annual_diurnal_range, 
                                   range_min_max, isotherm, annual_prcp) %>%
  # Summarize over the 10-year average
  ungroup() %>% 
  filter(year >= env_year - 10, year <= env_year - 1) %>% 
  group_by(environment, variable) %>%
  summarize(value = mean(value, na.rm = TRUE)) %>%
  ungroup()

## Average GDD and photothermality for each day over years
gdd_stats_avg <- gdd_stats1 %>% 
  filter(year >= env_year - 10, year <= env_year - 1) %>% 
  group_by(environment, day) %>% 
  summarize(GDD = mean(GDD, na.rm = TRUE)) %>%
  ungroup()

agdd_stats_avg <- agdd_stats_adj1 %>% 
  filter(year >= env_year - 10, year <= env_year - 1) %>% 
  group_by(environment, day_from_planting) %>% 
  summarize(AGDD = mean(AGDD, na.rm = TRUE)) %>%
  ungroup()

phototherm_daily_stats_avg <- phototherm_stats_gdd1 %>% 
  filter(year >= env_year - 10, year <= env_year - 1) %>% 
  group_by(environment, day) %>% 
  summarize(phototherm_daily = mean(phototherm_daily, na.rm = TRUE)) %>%
  ungroup()

phototherm_stats_avg <- phototherm_stats_agdd1 %>% 
  filter(year >= env_year - 10, year <= env_year - 1) %>% 
  group_by(environment, day_from_planting) %>% 
  summarize(phototherm = mean(phototherm, na.rm = TRUE)) %>%
  ungroup()


multi_year_daily_summary <- list(
  gdd = gdd_stats_avg,
  agdd = agdd_stats_avg,
  phototherm_gdd = phototherm_daily_stats_avg,
  phototherm = phototherm_stats_avg
)




# Set the maximum number of days in the interval length
max_days <- 90

## For each environment and daily covariable, build intervals of the average variable
multi_year_daily_summary_interval <- multi_year_daily_summary %>%
  map(~{
    df <- .
    
    # Find the minimum day among environments.
    min_day <- group_by(df, environment) %>% 
      summarize_at(vars(contains("day")), min) %>% 
      pull() %>% 
      max()
    
    # Build the intervals
    intervals_df <- crossing(begin = seq(min_day, min_day + max_days), end = seq(min_day + max_days)) %>% 
      filter(begin <= end)
    
    # For each environment, convert the day to starting from 1 (planting day)
    df1 <- df %>% 
      group_by(environment)
    
    # Day name
    dayname <- names(df1)[2]
    names(df1)[2] <- "day"
    
    # Iterate over the intervals
    intervals_summ <- pmap(intervals_df, ~filter(df1, between(day, .x, .y)) %>% 
                             summarize_at(vars(-environment, -day), mean))
    
    
    mutate(intervals_df, variable = tail(names(df), 1), out = intervals_summ)
    
  })













  
## Combine climate and soil data

# Replace NA with the variable mean across environments
one_year_env_df <- bind_rows(one_year_summary_df, soil_data_complete1) %>%
  group_by(variable) %>% 
  mutate(value = ifelse(is.na(value), mean(value, na.rm = TRUE), value)) %>% 
  ungroup()
  
multi_year_env_df <- bind_rows(multi_year_summary_df, soil_data_complete1) %>%
  group_by(variable) %>% 
  mutate(value = ifelse(is.na(value), mean(value, na.rm = TRUE), value)) %>% 
  ungroup()
  
# # Convert to matrices
# one_year_env_mat <- one_year_env_df %>% 
#   spread(variable, value) %>% 
#   as.data.frame() %>%
#   remove_rownames() %>%
#   column_to_rownames("environment") %>% 
#   as.matrix()
#   
# multi_year_env_mat <- multi_year_env_df %>% 
#   spread(variable, value) %>% 
#   as.data.frame() %>%
#   remove_rownames() %>%
#   column_to_rownames("environment") %>% 
#   as.matrix()

# Save the data
save_file <- file.path(data_dir, "environmental_data_compiled.RData")
save("one_year_env_df", "multi_year_env_df", "one_year_daily_summary", "one_year_daily_summary_interval",
     "multi_year_daily_summary", "multi_year_daily_summary_interval", file = save_file)

