## S2MET Phenotypic Adjustment
## 
## Author: Jeff Neyhart
## Last updated: March 26, 2018
## 
## This notebook outlines procedures for calculating adjusted phenotypic means of
## entries in trials that belong to the `S2MET` experiment.
## 

# Load libraries and directories

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(modelr)
library(broom)
library(pbr)
library(ggridges)

# Repository directory
repo_dir <- getwd()

# Source the main project script
source(file.path(repo_dir, "source.R"))




# Project directory
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET/"


# Load the tidy S2 data
load("C:/Users/jln54/GoogleDrive/BarleyLab/Breeding/PhenotypicData/Final/MasterPhenotypes/S2_tidy_pheno.RData")

# Subset the traits and relevant trials
S2_MET_tidy <- s2_tidy_pheno %>%
  filter(year %in% 2015:2017, !str_detect(trial, "FHB|PVV")) %>%
  mutate(location = ifelse(location == "CNY", "HNY", location),
         environment = ifelse(environment == "CNY15", "HNY15", environment)) %>%
  group_by(trait) %>%
  filter(n_distinct(environment) > 35) %>%
  # filter(dplyr::n_distinct(environment) >= 4, !trait %in% c("Moisture", "BarleyColor", "ThinGrains", "LodgingDegree", "MaturityDate", "KernelWeight")) %>%
  ungroup()




## Create a data.frame of data to model
# Set a dummy variable for check or test
# Also add correct replicate numbers
data_to_model <- S2_MET_tidy %>%
  group_by(trial, trait, line_name) %>% 
  ungroup() %>%
  mutate(line_name = as.character(line_name),
         line = ifelse(!line_name %in% checks, line_name, "00check"),
         check = ifelse(line_name %in% checks, line_name, "00line")) %>%
  mutate_at(vars(rep, line, check, row, column, blk), as.factor)


## Add heading date as AGDD
## Load the environmental covariables
load(file.path(data_dir, "environmental_data_compiled.RData"))

# Select the accumulated growing degree days
agdd <- one_year_daily_summary$agdd

# Pull out heading date, round to the nearest whole number, then assign the AGDD value
data_to_model_HD <- data_to_model %>%
  filter(trait == "HeadingDate") %>% 
  mutate(day_from_planting = round(value, 0)) %>% 
  left_join(., agdd) %>%
  mutate(trait = "HeadingDateAGDD") %>%
  select(-value, -day_from_planting) %>%
  rename(value = AGDD)


data_to_model1 <- bind_rows(data_to_model, data_to_model_HD)



## Phenotypic adjustment will be performed so that a two-stage analysis of 
## multiple environments is possible. This analysis pipeline is inspired by @Mohring2009.

## Fit models most appropriate for the trial.

## Stage-One

# Stage-one models are environment specific and estimate the standard error in each
# environment

# I will include row, column, and blk as random effects and use backwards elimination to reduce
# the models
rand_eff <- c("row", "column", "blk")

# The base formula will have line and check as fixed effects
base_form <- value ~ -1 + line + check

## Create a empty data frame of the distinct trials and traits
stage_one_results <- data_to_model1 %>% 
  distinct(trait, trial) %>%
  mutate(out = list(NULL))

# Iterate over rows
for (i in seq(nrow(stage_one_results))) {

  # Create a separate df object - this greatly improves the speed of the model fitting
  df <- data_to_model1 %>%
    filter(trait == stage_one_results$trait[i], trial == stage_one_results$trial[i]) %>%
    droplevels()
  
  # Print the trait
  cat(c(unique(df$trait), unique(df$trial)))
  
  # Determine the random effects that are available
  rand_eff_tomodel <- rand_eff[which(!is.na(df[1,rand_eff]))]
  
  # Create the full formula
  full_form <- base_form %>%
    add_predictors(as.formula(str_c("~ ", str_c("(1|", rand_eff_tomodel, ")", collapse = " + "))))
  
  # Fit the model
  fit <- lmer(formula = full_form, data = df, contrasts = list(line = "contr.sum"))
  # Backwards elimination
  back_elim <- step(fit, alpha.random = 0.05, reduce.fixed = FALSE)
  fit_final <- get_model(back_elim)
  
  # Tidy up
  if (class(fit_final) == "lm") {
    fit_final_tidy <- tidy(fit_final)
    
  } else {
    fit_final_tidy <- structure(fit_final, class = "merMod") %>%
      tidy() %>%
      filter(group == "fixed")
    
  }
  
  # Get the levels of the checks
  check_levels <- levels(df$check)
  
  fit_tidy1 <- fit_final_tidy %>% 
    mutate(term = if_else(term == "line00check", tail(check_levels, 1), 
                          str_replace_all(term, pattern = "line|check", replacement = "")), 
           estimate = if_else(term %in% head(check_levels, -1), estimate + estimate[1], estimate)) %>% 
    select(line = term, estimate, std.error)

  ## Fit a new model to calulate heritability
  ran_form <- str_replace_all(string = formula(fit_final), pattern = "line", replacement = "(1|line)") %>%
    str_replace_all("- 1", "") %>%
    tail(1) %>% 
    str_c("value ~ ", .) %>% 
    as.formula()
  
  fit_ran <- lmer(formula = ran_form, data = df)
  
  # Find the harmonic mean of the number of replicates
  n_r <- table(df$line_name) %>%
    harm_mean()
  
  # Calculate heritability
  h2 <- herit(object = fit_ran, exp = "line / (line + (Residual / n_r))", n_r = n_r)[[1]]
  
  # Return the coefficients and the variance of the adjusted means
  stage_one_results$out[[i]] <- data_frame(
    BLUE = list(fit_tidy1), 
    terms = list(formula(fit_final)),
    heritability = h2) 

}


# Unnest
stage_one_data <- stage_one_results %>%
  unnest()


# Extract the heritability estimates
stage_one_herit <- stage_one_data %>%
  select(trial, trait, heritability)



# Plot heritability
stage_one_herit %>%
  ggplot(aes(x = trial, y = heritability )) +
  geom_col() +
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(trait ~ .) +
  labs(title = "Broad-Sense Heritability of Each Trial") +
  ylab("Heritability") +
  xlab("Environment") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))


## Plot heritability of heading date (DAP) to heading date (AGDD)
stage_one_herit %>% 
  filter(str_detect(trait, "HeadingDate")) %>%
  spread(trait, heritability) %>% 
  qplot(x = HeadingDate, y = HeadingDateAGDD, data = .) + 
  geom_abline(slope = 1)


## How many trials/traits had co-factors in the final model?
stage_one_data <- stage_one_data %>%
  mutate(nCofactors = map_dbl(terms, ~str_count(string = ., pattern = paste0(rand_eff, collapse = "|")) %>% sum()))




# Combine trial information from the metadata (ie. change some trial/environment names)
stage_one_data <- stage_one_data %>%
  left_join(., select(trial_info, trial, environment, location, year), by = ("trial")) %>%
  select(trial, environment, location, year, names(.))

# Extract BLUEs
S2_MET_BLUEs_all <- stage_one_data %>% 
  select(trial:trait, BLUE) %>% 
  unnest() %>% 
  select(trial:year, line_name = line, trait, value = estimate, std_error = std.error) %>%
  # filter(line_name %in% entries) %>%
  ungroup() %>%
  droplevels() %>%
  arrange(trait, environment, line_name)

S2_MET_BLUEs <- S2_MET_BLUEs_all %>%
  filter(line_name %in% entries) %>%
  droplevels()

# Is the standard error of the BLUEs reflective of the heritability within an environment?
g_herit_se <- stage_one_data %>% 
  select(trial, environment, trait, BLUE, heritability) %>% 
  unnest(BLUE) %>% 
  group_by(environment, trait) %>%
  mutate(mean_se = mean(std.error)) %>%
  distinct(environment, trait, heritability, mean_se) %>% 
  ggplot(aes(x = heritability, y = mean_se, group = trait)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~ trait, scales = "free_y", ncol = 2) +
  ylab("Standard Error") +
  xlab("Heritability")



# Save
save_file <- file.path(data_dir, "S2_MET_BLUEs.RData")
save("S2_MET_BLUEs_all", "S2_MET_BLUEs", "stage_one_data", file = save_file)


# # Save
# save_file <- file.path(data_dir, "S2_MET_BLUEs_malt_quality.RData") 
# save("S2_MET_BLUEs_all", "S2_MET_BLUEs", "stage_one_data", file = save_file)
# 





