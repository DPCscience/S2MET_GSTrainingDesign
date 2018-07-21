## Collect and Gather Environmental Covariates
## 
## Author: Jeff Neyhart
## Last modified: May 25, 2018
## 
## This R Script outlines the methods to gather environmental covariables on soil
## and weather data
## 


# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))


# Load packages and set directories
library(lubridate)
library(rnoaa)
library(geosphere)
library(soilDB)
library(foreign)
library(raster)




## Data Collection

### NOAA Weather Data

# Weather data will be collected by querying the NOAA API using the latitude/longitude coordinates
# from the trial metadata. These coordinates will be used to locate the nearest weather station
# and weather data will be drawn from these stations.

# Data will be collected on the basis of two approaches:
# 1. Long-term (10 years of data will be collected)
# 2. Short term (only data within the year of the environment will be used)

# Here is some code for downloading the GHCND stations and saving those station as an `.RData` file

# # Load the list of cooperative stations - ONLY RUN ONCE -
# station_list <- ghcnd_stations()
# 
# save_file <- file.path(env_var_dir, "NOAA_Weather_Stations/GHCND_stations.RData")
# save("station_list", file = save_file)



# API token
token <- 'ubmUkwrtxZOQLeyXZdtCMjkKcCTpPDcK'

# Start and end year
end_year <- max(trial_info$year)
start_year <- min(trial_info$year) - 10

# Collect data from March 1 to August 31 of each year
start_date <- "03-01"
end_date <- "08-31"

# Set a threshold for missing data on a station-datatype basis
max_missing <- 0.20


# Desired data types
desired_vars <- c("TMIN", "TMAX", "PRCP")

# Load the list of cooperative stations
load(file.path(data_dir, "Climate_Data/NOAA_Weather_Stations/GHCND_stations.RData"))


# Remove stations with no id and keep stations with the USC prefix
station_list1 <- station_list %>% 
  filter(str_detect(id, pattern = "USC|CA"),
         first_year <= start_year, 
         last_year >= end_year) %>%
  # Filter stations on whether all of the desired variables are present
  group_by(id) %>% 
  filter(all(desired_vars %in% element))

# Grab the unique stations
stations_unique <- station_list1 %>%
  distinct(.keep_all = TRUE) %>%
  ungroup()

# Create a data.frame of unique trial locations
locations_unique <- trial_info %>%
  filter(!is.na(latitude), !is.na(longitude),
         !location %in% c("Sidney_MT", "Driscoll_ND")) %>%
  distinct(location, .keep_all = T) %>% 
  dplyr::select(location, latitude, longitude)

# Create a list of distances from each trial location to all of the stations
# The distance is the great circle distance given the WGS84 ellipsoid
loc_dist <- locations_unique %>%
  dplyr::select(-location) %>%
  pmap(~{
    
    # Extract coordiates
    coord <- c(.y, .x)
    
    stations_unique %>% 
      group_by(id) %>% 
      summarize(dist_to_site = distGeo(p1 = coord, p2 = c(unique(longitude), unique(latitude))))
    
  })

trial_loc_dist <- locations_unique %>%
  mutate(station_dist = loc_dist)


# Find the 10 stations that are closest to each location
closest_station <- trial_loc_dist %>% 
  unnest() %>% 
  group_by(location) %>% 
  top_n(n = -10, wt = dist_to_site) %>% 
  arrange(location, dist_to_site) %>%
  ungroup()

# Group by environment and pull out data
trial_station <- closest_station %>%
  # Filter by the closest station. If tied, take the first
  group_by(location) %>% 
  top_n(n = -1, wt = dist_to_site) %>% 
  dplyr::slice(1) %>%
  ungroup()

trial_station_data <- map(trial_station$id, ~{
  
  # Change the station name
  station_id <- str_c("GHCND:", .)
  
  # Year sequences
  yrs <- seq(start_year, end_year)
  
  # Create a list of length number of years
  data_list <- vector("list", length(yrs)) %>%
    structure(names = yrs)
  
  # Iterate over the requested years
  for (yr in seq(start_year, end_year)) {
    
    # Format the date
    start_ymd <- ymd(str_c(yr, start_date))
    end_ymd <- ymd(str_c(yr, end_date))
    
    # Pull data
    ncdc_data <- ncdc(datasetid = 'GHCND', datatypeid = desired_vars, stationid = station_id, 
                      startdate = as.character(start_ymd), enddate = as.character(end_ymd), 
                      limit = 1000, token = token, )
    
    # Reformat and determine the number of days since the start date
    data_list[[as.character(yr)]] <- ncdc_data$data %>%
      as_data_frame()
    
  } # Close the loop
  
  # Print the completion of the location
  print(str_c("Data for collected for station ID: ", .))
  
  # Return the list
  bind_rows(data_list)
  
})

# Combine with the trial location data.frame
trial_station_df <- trial_station %>% 
  mutate(out = trial_station_data)
  

# Quality-control the data
trial_station_data_qc <- trial_station_df %>% 
  unnest() %>%
  mutate(date = ymd(ymd_hms(date))) %>%
  group_by(location, station) %>% 
  filter(fl_q == "") %>%
  # Reorder
  dplyr::select(location, id = station, date, datatype, value)


# Complete the data.frame for all possible observations
trial_station_data_complete <- trial_station_data_qc %>%
  ungroup() %>%
  mutate_at(vars(datatype, date), as.factor) %>%
  group_by(location, id) %>%
  complete(datatype, date) %>% 
  mutate(year = year(date), 
         month = month(date),
         day = yday(date) - yday(date[1])) %>% 
  dplyr::select(location:date, year:day, value)
  




# Deal with missing data by first assessing the level of missingness. If the level of missing 
# data is above the threshold provided above, reject the data from that environment. Then look 
# for data from the next closest station for each environment and extract data.


## Find the level of missing data for each environment
trial_station_missing <- trial_station_data_complete %>%
  group_by(location, datatype, year) %>%
  summarize(prop_miss = mean(is.na(value)))

## Plot this
g_station_missing <- trial_station_missing %>%
  ggplot(aes(x = year, y = prop_miss, color = datatype)) + 
  geom_line() + 
  facet_wrap(~location, ncol = 4) +
  theme_bw()




## Pull out environments that need re-gathering - total
redo_env_all <- trial_station_missing %>% 
  group_by(location, datatype) %>% 
  summarize(prop_miss = mean(prop_miss)) %>% 
  filter(any(prop_miss > max_missing)) %>%
  ungroup() %>% 
  distinct(location) %>%
  pull()

## Pull out environments that need re-gathering for 2015, 2016, or 2017
redo_env_rel <- trial_station_missing %>% 
  group_by(location, datatype) %>% 
  filter(year >= 2015) %>%
  filter(any(prop_miss > max_missing)) %>% 
  ungroup() %>% 
  distinct(location) %>%
  pull()

# Common environments
redo_env <- union(redo_env_all, redo_env_rel)

## Save the data from those environments that will not be re-done for 2015-2017
trial_station_data_oneyear_complete1 <- trial_station_data_complete %>%
  ungroup() %>%
  filter(year >= 2015, !location %in% redo_env_rel)

# Create a copy of the data
trial_station_data_multiyear_complete1 <- trial_station_data_complete %>%
  ungroup() %>%
  filter(!location %in% redo_env_all)

# Look for the second-closest station and proceed
dist_rank <- 2

# What is the max rank?
max_rank <- closest_station %>%
  group_by(location) %>% 
  summarize(n = n()) %>% 
  pull(n) %>% 
  unique()

# If the number of redo environments is greater than 0, proceed
while (length(redo_env) >= 1) {
  
  # IS the rank greater than the max rank? If so, stop
  if (dist_rank > max_rank) {
    print("Rank is above the max. Stopping.")
    break
  }
  
  # Go back to the closest_station data frame and pull out the next closest station
  trial_redo_station<- closest_station %>% 
    filter(location %in% redo_env) %>% 
    group_by(location) %>%
    top_n(n = -dist_rank, dist_to_site) %>% 
    top_n(n = 1, wt = dist_to_site) %>%
    # If tied, select the second
    top_n(n = 1, wt = id) 
  
  
  trial_redo_station_data <- map(trial_redo_station$id, ~{
    
    # Change the station name
    station_id <- str_c("GHCND:", .)
    
    # Year sequences
    yrs <- seq(start_year, end_year)
    
    # Create a list of length number of years
    data_list <- vector("list", length(yrs)) %>%
      structure(names = yrs)
    
    # Iterate over the requested years
    for (yr in seq(start_year, end_year)) {
      
      # Format the date
      start_ymd <- ymd(str_c(yr, start_date))
      end_ymd <- ymd(str_c(yr, end_date))
      
      # Pull data
      ncdc_data <- ncdc(datasetid = 'GHCND', datatypeid = desired_vars, stationid = station_id, 
                        startdate = as.character(start_ymd), enddate = as.character(end_ymd), 
                        limit = 1000, token = token)
      
      # Reformat and determine the number of days since the start date
      data_list[[as.character(yr)]] <- ncdc_data$data %>%
        as_data_frame()
      
    } # Close the loop
    
    # Print the completion of the location
    print(str_c("Data for collected for station ID: ", .))
    
    # Return the list
    bind_rows(data_list)
    
  })
  
  
  trial_redo_station_df <- trial_redo_station %>% 
    ungroup() %>%
    mutate(out = trial_redo_station_data)
  
  # Quality-control the data
  trial_redo_station_data_qc <- trial_redo_station_df %>% 
    unnest() %>%
    mutate(date = ymd(ymd_hms(date))) %>%
    group_by(location, station) %>% 
    filter(fl_q == "") %>%
    # Reorder
    dplyr::select(location, id = station, date, datatype, value)
  
  
  # Complete the data.frame for all possible observations
  trial_redo_station_data_complete <- trial_redo_station_data_qc %>%
    ungroup() %>%
    mutate_at(vars(datatype, date), as.factor) %>%
    group_by(location, id) %>%
    complete(datatype, date) %>% 
    mutate(year = year(date), 
           month = month(date),
           day = yday(date) - yday(date[1])) %>% 
    dplyr::select(location:date, year:day, value)
  
  ## Find the level of missing data for each environment
  trial_station_missing <- trial_redo_station_data_complete %>%
    group_by(location, datatype, year) %>%
    summarize(prop_miss = mean(is.na(value)))
  
  
  ## Pull out environments that need re-gathering - total
  redo_env_all <- trial_station_missing %>% 
    filter(any(mean(prop_miss) > max_missing)) %>% 
    ungroup() %>% 
    distinct(location) %>%
    pull()
  
  ## Pull out environments that need re-gathering for 2015, 2016, or 2017
  redo_env_rel <- trial_station_missing %>% 
    filter(year >= 2015) %>%
    filter(any(prop_miss > max_missing)) %>% 
    ungroup() %>% 
    distinct(location) %>%
    pull()
  
  
  ## Combine the redo data with the saved data
  # First do this for environments that were re-done for 2015-2017
  trial_station_data_oneyear_complete1 <- trial_station_data_oneyear_complete1 %>%
    bind_rows(., filter(trial_redo_station_data_complete, year >= 2015, 
                        !location %in% trial_station_data_oneyear_complete1$location,
                        !location %in% redo_env_rel)) %>%
    arrange(location)

  # Now do this for environments that were re-done for the whole period
  trial_station_data_multiyear_complete1 <- trial_station_data_multiyear_complete1 %>% 
    bind_rows(., filter(trial_redo_station_data_complete, 
                        !location %in% trial_station_data_multiyear_complete1$location,
                        !location %in% redo_env_all)) %>%
    arrange(location)

  # Create an object with the environments to be re-sampled
  redo_env <- union(redo_env_all, redo_env_rel) %>%
    unlist()
  
  # Increase the rank counter
  dist_rank <- dist_rank + 1
  
}


# Nest the data for each location and add station metadata
# Include the closest station and the station that was actually used
noaa_trial_data_oneyear_complete <- trial_station_data_oneyear_complete1 %>%
  nest(datatype:value) %>% 
  left_join(., mutate(closest_station, id = str_c("GHCND:", id)), by = c("location", "id")) %>%
  left_join(., by = c("location"), y = mutate(trial_station, id = str_c("GHCND:", id)) %>%
              dplyr::select(location, closest_station_id = id, dist_to_closest_station = dist_to_site)) %>%
  separate(id, c("database", "id"), sep = ":") %>%
  left_join(., distinct(station_list1, id, elevation, state), by = "id")


noaa_trial_data_multiyear_complete <- trial_station_data_multiyear_complete1 %>%
  nest(datatype:value) %>% 
  left_join(., mutate(closest_station, id = str_c("GHCND:", id)), by = c("location", "id")) %>%
  left_join(., by = c("location"), y = mutate(trial_station, id = str_c("GHCND:", id)) %>%
              dplyr::select(location, closest_station_id = id, dist_to_closest_station = dist_to_site)) %>%
  separate(id, c("database", "id"), sep = ":") %>%
  left_join(., distinct(station_list1, id, elevation, state), by = "id")

save_file <- file.path(data_dir, "Climate_Data/NOAA_Data/noaa_stations_trial_data.RData")
save("noaa_trial_data_oneyear_complete", "noaa_trial_data_multiyear_complete", file = save_file)






### Soil Data

# Soil data for sites in the U.S. will be queried from the NRCS Soil Survey, and soil data for 
# sites in Canada will be queried from the National Soil Database

# First filter the `trial_info` data.frame for U.S. versus Canadian locations

# Define the grid size around a latitude/longitude coordinate
grid_size <- 1e-6

# Filter out the Canadian sites and group by location
canadian <- c("EON", "CPE", "PQC")

trial_info_us <- trial_info %>%
  filter(!str_detect(environment, str_c(canadian, collapse = "|")))

trial_info_can <- trial_info %>%
  filter(str_detect(environment, str_c(canadian, collapse = "|")))




# Iterate over the locations
trial_info_us_mukey <- trial_info_us %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  distinct(environment, .keep_all = TRUE) %>% 
  split(.$environment) %>%
  purrr::map(function(loc_info) {
    
    # Create the coordinate grid
    b <- c(loc_info$longitude - grid_size,
           loc_info$latitude - grid_size,
           loc_info$longitude + grid_size,
           loc_info$latitude + grid_size)
    
    # Find information based on the coordinates and grid size
    mapunit <- mapunit_geom_by_ll_bbox(b)
    
    # Return the mukey
    loc_info %>%
      mutate(mukey = unique(mapunit$mukey)) }) %>%
  # Bind rows
  bind_rows()

# For each location and using the mukey, gather relevant soil variables
trial_info_us_vars <- trial_info_us_mukey %>%
  distinct(location, mukey) %>%
  group_by(location, mukey) %>%
  do({
    
    # Create the query string to get the cokey
    query <- str_c("SELECT * FROM component WHERE mukey = '", .$mukey, "'")
    
    # Query and get the cokey of the most abundant component
    cokey <- SDA_query(query) %>%
      top_n(n = 1, wt = comppct_r) %>%
      # If a tie, take the first one
      head(1) %>%
      dplyr::select(cokey, comppct_r)
    
    query <- str_c("SELECT * FROM chorizon WHERE cokey = '", cokey$cokey, "'")
    # Use the cokey to get characteristic data
    SDA_query(query) %>%
      gather(variable, value, -hzname, -hzdept_r, -hzdepb_r) %>%
      mutate(comppct_r = cokey$comppct_r) })

# Recombine data with the trial info
trial_info_us_vars1 <- left_join(trial_info_us_mukey, trial_info_us_vars, by = c("location", "mukey")) %>%
  

# Subset the variables
# Also parse the values of the variables to numbers
desired_vars <- c("hzdept_r", "hzdepb_r", "sandtotal_r", "silttotal_r", "claytotal_r", "om_r", "ph1to1h2o_r")

trial_info_us_vars_sub <- trial_info_us_vars1 %>%
  filter(variable %in% desired_vars) %>%
  mutate(value = parse_number(value))

# What is the level of missing data per variable?
trial_info_us_vars_sub %>%
  group_by(variable) %>%
  summarize(prop_miss = mean(is.na(value)))

# Per location?
trial_info_us_vars_sub %>%
  group_by(environment) %>%
  summarize(prop_miss = mean(is.na(value)))

# Merge data with the remaining trial info
nrcs_trial_data_complete <- trial_info_us_vars_sub

# Save the data
save_file <- file.path(env_var_dir, "Soil_Data/USDA_NRCS/nrcs_trial_soil_data.RData")
save("nrcs_trial_data_complete", file = save_file)



# Set the directory for pulling the Canadian soil data
nsdb_dir <- file.path(env_var_dir, "Soil_Data/NSDB/")

# Designate the shapefile with polygon information
shape_file <- file.path(nsdb_dir, "ca_all_slc_v3r2")
# Designate the dbf with soil component information
comp_file <- file.path(nsdb_dir, "ca_all_slc_v3r2_cmp.dbf")
# Designate the dbf file with soil layer information
layer_file <- file.path(nsdb_dir, "soil_layer_canada_v2r20170602.dbf")

# Read in files
# nsdb_poly <- readOGR(dsn = nsdb_dir, layer = "ca_all_slc_v3r2")
# Use raster instead
nsdb_poly <- shapefile(shape_file)

# Transform projection information
nsdb_poly <- spTransform(x = nsdb_poly, CRSobj = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")


# Read in the dbf files
nsdb_comp <- read.dbf(file = comp_file, as.is = TRUE)
nsdb_layer <- read.dbf(file = layer_file, as.is = TRUE)

## For each trial, find the soil information corresponding to the lat/long
trial_info_can_vars <- trial_info_can %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  group_by(environment) %>%
  do({
    
    # Create the box
    bbx <- c(.$longitude - 1e-6,
             .$longitude + 1e-6,
             .$latitude - 1e-6,
             .$latitude + 1e-6)
    
    # Crop the polygon
    nsdb_poly_crop <- crop(x = nsdb_poly, y = extent(bbx))
    
    # Get the POLY_ID
    poly_id <- nsdb_poly_crop$POLY_ID
    
    # Subset the components file for that POLY_ID and take the soil component that is most abudant
    comp_select <- nsdb_comp %>%
      filter(POLY_ID == poly_id) %>%
      top_n(n = 1, wt = PERCENT)
    
    # Now use the SOIL_ID from this selction to obtain information about that soil using the
    # nsdb_layer file
    soil_id <- as.character(comp_select$SOIL_ID)
    
    soil_select <- nsdb_layer %>% 
      filter(SOIL_ID == soil_id)
    
    # Combine the soil data back to the component data
    full_join(comp_select, soil_select, by = c("SOIL_ID", "PROVINCE", "SOIL_CODE", "MODIFIER", "PROFILE")) %>%
      gather(variable, value, -HZN_MAS, -UDEPTH, -LDEPTH) %>%
      select(HZN_MAS, UDEPTH, LDEPTH, variable, value) })

# Vector of desired variables
desired_vars <- c("TSAND", "TSILT", "TCLAY", "ORGCARB", "PH2")

# Subset the data for the desired variables
trial_info_can_vars_sub <- trial_info_can_vars %>%
  filter(variable %in% desired_vars) %>%
  mutate(value = parse_number(value),
         value = ifelse(value == -9, NA, value))

# What is the level of missing data per variable?
trial_info_can_vars_sub %>%
  group_by(variable) %>%
  summarize(prop_miss = mean(is.na(value)))

# Per location?
trial_info_can_vars_sub %>%
  group_by(environment) %>%
  summarize(prop_miss = mean(is.na(value)))

# Merge data with the remaining trial info
nsdb_trial_data_complete <- full_join(trial_info, trial_info_can_vars_sub, "environment")

# Save the data
save_file <- file.path(env_var_dir, "Soil_Data/NSDB/nsdb_trial_soil_data.RData")
save("nsdb_trial_data_complete", file = save_file)


# Combine the soil data and create common variable names


nrcs_file <- file.path(env_var_dir, "Soil_Data/USDA_NRCS/nrcs_trial_soil_data.RData")
nsdb_file <- file.path(env_var_dir, "Soil_Data/NSDB/nsdb_trial_soil_data.RData")

load(nrcs_file)
load(nsdb_file)

# Expand each df
nrcs_expand <- nrcs_trial_data_complete %>%
  filter(!is.na(variable)) %>%
  select(-comppct_r, -mukey, -notes) %>%
  spread(variable, value)

nsdb_expand <- nsdb_trial_data_complete %>%
  filter(!is.na(variable)) %>%
  select(-notes) %>%
  spread(variable, value)

# Combine
soil_data_complete <- full_join(nrcs_expand, nsdb_expand, by = c("trial", "environment", "location", "year", "latitude", "longitude", "planting_date", "hzname" = "HZN_MAS", "hzdept_r"= "UDEPTH", "hzdepb_r" = "LDEPTH", "claytotal_r" = "TCLAY", "om_r" = "ORGCARB", "ph1to1h2o_r" = "PH2", "sandtotal_r" = "TSAND", "silttotal_r" = "TSILT"))

# Re-tidy
soil_data_complete <- soil_data_complete %>%
  select(-trial) %>%
  gather(variable, value, -environment:-hzdepb_r)

# Export
save_file <- file.path(env_var_dir, "Soil_Data/complete_trial_soil_data.RData")
save("soil_data_complete", file = save_file)




























