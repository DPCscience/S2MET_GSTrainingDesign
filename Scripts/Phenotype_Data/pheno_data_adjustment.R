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
library(stringr)
library(readxl)
library(lme4)
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

# Traits
traits <- c("GrainYield", "HeadingDate", "PlantHeight")

# Load the tidy S2 data
load("C:/Users/Jeff/Google Drive/Barley Lab/Breeding//Phenotypic Data/Final/Master Phenotypes/S2_MET_tidy.RData")


## Create a data.frame of data to model
# Set a dummy variable for check or test
# Also add correct replicate numbers
data_to_model <- S2_MET_tidy %>% 
  # Filter for relevant traits
  filter(trait %in% traits) %>%
  group_by(trial, trait, line_name) %>% 
  mutate(rep = seq(sum(duplicated(line_name)) + 1)) %>%
  ungroup() %>%
  mutate(line_name = as.character(line_name),
         line = ifelse(!line_name %in% checks, line_name, "00check"),
         check = ifelse(line_name %in% checks, line_name, "00line")) %>%
  mutate_at(vars(rep, line, check), as.factor)



## Phenotypic adjustment will be performed so that a two-stage analysis of 
## multiple environments is possible. This analysis pipeline is inspired by @Mohring2009.

## Fit models most appropriate for the trial.

## Stage-One

# Stage-one models are environment specific and estimate the standard error in each
# environment

# Designate two models, depending on the environment and the experimental design
forms <- formulas(~ value,
                  form1 = ~ -1 + line + check,
                  form2 = ~ -1 + line + check + (1|blk))

# Fit the models
stage_one <- data_to_model %>%
  group_by(trial, trait) %>%
  do({
    
    # Create a separate df object - this greatly improves the speed of the model fitting
    df <- droplevels(.)
    
    # Print the trait
    print(unique(df$trait))
    
    # Extract the function and the formula, then fit
    if (is.na(df$blk[1])) {
      f <- forms$form1
      f_mean <- str_replace(string = f, pattern = "-1 \\+", "") %>%
        tail(-1) %>% 
        str_c(collapse = " ~")
        
      fit <- lm(formula = f, data = df)
      fit_mean <- lm(formula = f_mean, data = df)
      
      # Tidy
      fit_tidy <- tidy(fit) %>% 
        mutate(group = "fixed")
      
      fit_mean_tidy <- tidy(fit_mean) %>% 
        mutate(group = "fixed")
      
    } else {
      f <- forms$form2
      f_mean <- str_replace(string = f, pattern = "-1 \\+", "") %>%
        tail(-1) %>% 
        str_c(collapse = " ~")
      
      fit <- lmer(formula = f, data = df)
      fit_mean <- lmer(formula = f_mean, data = df)
      
      fit_tidy <- tidy(fit)
      fit_mean_tidy <- tidy(fit_mean)
      
    }
    
    # Get the levels of the checks
    check_levels <- levels(df$check)
    
    fit_tidy1 <- fit_tidy %>% 
      filter(group == "fixed") %>% 
      mutate(term = if_else(term == "line00check", tail(check_levels, 1), 
                            str_replace_all(term, pattern = "line|check", replacement = "")), 
             estimate = if_else(term %in% head(check_levels, -1), estimate + estimate[1], estimate)) %>% 
      select(line = term, estimate, std.error)
    
    # Estimate z scores
    z_scores <- fit_mean_tidy %>% 
      filter(str_detect(term, "line")) %>% 
      transmute(line = str_replace_all(term, "line", ""), z = estimate / std.error)

    
    ## Fit a new model to calulate heritability
    f_ran <- str_replace_all(string = f, pattern = "line", replacement = "(1|line)") %>% 
      tail(1) %>% 
      str_c("value ~ ", .) %>% 
      as.formula()
    
    fit_ran <- lmer(formula = f_ran, data = df)
    
    # Find the harmonic mean of the number of replicates
    n_r <- table(df$line_name) %>%
      harm_mean()
    
    # Add the z scores to the fit_tidy
    fit_tidy2 <- left_join(x = fit_tidy1, y = z_scores, by = "line")
    
    # Return the coefficients and the variance of the adjusted means
    data_frame(fit = list(fit), fit_ran = list(fit_ran), BLUE = list(fit_tidy2), 
               n_r = n_r) })


## Calculate the heritability in each environment

# Calculate heritability
stage1_herit <- stage_one %>% 
  do(suppressWarnings(herit_boot(object = .$fit_ran[[1]], exp = "line / (line + (Residual / n_r))", 
                                 n_r = .$n_r, boot.reps = 100)))
  
# Plot
stage1_herit %>%
  ggplot(aes(x = trial, y = heritability - bias)) +
  geom_col() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(trait ~ .) +
  labs(title = "Broad-Sense Heritability of Each Trial") +
  ylab("Heritability") +
  xlab("Environment") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))



# Combine trial information from the metadata (ie. change some trial/environment names)
stage_one_data <- full_join(stage_one, stage1_herit) %>% 
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




