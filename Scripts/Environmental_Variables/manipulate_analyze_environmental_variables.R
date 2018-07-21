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
  select(environment, latitude, longitude) %>% 
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
gdd_stats <- oneyear_trial_env_data %>% 
  select(environment, date, month, datatype, value) %>% 
  filter(datatype != "PRCP") %>% 
  spread(datatype, value) %>%
  # Convert from tenths of C to F
  mutate_at(vars(TMIN, TMAX), ~(. * (9/5)) + 32) %>%
  # If the daily min or max is < 32 F, it is set to 32
  mutate_at(vars(TMIN, TMAX), ~ifelse(. < 32, 32, .)) %>%
  mutate(TAVG = (TMIN + TMAX) / 2,
         GDD = TAVG - 32)

## Use planting date to calculate accumulated growing degree days
agdd_stats <- gdd_stats %>% 
  left_join(., distinct(trial_info, environment, planting_date)) %>%
  # Is the recorded date before the planting date?
  mutate(date = ymd(date), planting_date = ymd(planting_date), 
         planted = date >= planting_date) %>% 
  filter(planted) %>% 
  group_by(environment) %>% 
  mutate(GDD = ifelse(is.na(GDD), 0, GDD),
         AGDD = cumsum(GDD)) %>%
  ungroup()

# The AGDD to reach Haun stage 2.0
haun2 <- 384

## Adjust the accumulated GDD using the restrictions found in the above website
agdd_stats_adj <- agdd_stats %>%
  mutate(TMAX = ifelse(AGDD < haun2 & TMAX > 70, 70, TMAX),
         TMAX = ifelse(AGDD >= haun2 & TMAX > 95, 95, TMAX),
         TAVG = (TMIN + TMAX) / 2,
         GDD = TAVG - 32,
         GDD = ifelse(is.na(GDD), 0, GDD)) %>%
  group_by(environment) %>% 
  mutate(AGDD = cumsum(GDD)) %>%
  ungroup()

# Tidy
agdd_stats_adj1 <- agdd_stats_adj %>%
  mutate(day = yday(date),
         variable = str_c("AGDD_day", day)) %>%
  select(environment, variable, value = AGDD)


# Combine the oneyear data

one_year_summary_df <- bind_rows(mon_summ_stats, annual_TAVG, annual_diurnal_range, 
                                 range_min_max, isotherm, annual_prcp, agdd_stats_adj1) %>%
 ungroup()






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
  select(environment, latitude, longitude) %>% 
  distinct() %>%
  gather(variable, value, -environment)

## Calculate summary statistics for temperature and precip


summ_stats_temp <- multiyear_trial_env_data %>% 
  filter(datatype != "PRCP") %>% 
  spread(datatype, value) %>% 
  group_by(environment, env_year, year, month) %>% 
  summarize_at(vars(TMAX, TMIN), mean, na.rm = TRUE) %>% 
  mutate(TAVG = (TMAX + TMIN) / 2) %>%
  gather(variable, value, -environment, -env_year, -year, -month)

summ_stats_prcp <- trial_env_data %>% 
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



  
# Combine the monthly and annual data
multi_year_summary_df <- bind_rows(mon_summ_stats, annual_TAVG, annual_diurnal_range, 
                                   range_min_max, isotherm, annual_prcp) %>%
  # Summarize over the 10-year average
  ungroup() %>% 
  filter(year >= env_year - 10, year <= env_year - 1) %>% 
  group_by(environment, variable) %>%
  summarize(value = mean(value, na.rm = TRUE)) %>%
  ungroup()

  
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
save("one_year_env_df", "multi_year_env_df", file = save_file)

