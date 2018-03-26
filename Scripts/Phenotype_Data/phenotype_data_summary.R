## S2MET Phenotypic Data Summary
## 
## Author: Jeff Neyhart
## Last updated: March 15, 2018
## 
## This notebook outlines procedures for calculating adjusted phenotypic means of
## entries in trials that belong to the `S2MET` experiment.
## 

# Load libraries and directories


## This notebook will provide some phenotyping data summaries for the S2MET project. 
## It will include:
## 1. Basic model for g + e + gxe
## 2. Heritability estimates
## 3. Trait correlations among environments

library(tidyverse)
library(broom)
library(stringr)
library(readxl)
library(modelr)
library(pbr)
library(rrBLUP)
library(ggridges)

# The head directory
repo_dir <- getwd()

source(file.path(repo_dir, "source.R"))


## Basic Summaries

## Number of lines per environment per trait

## Find the total number of possible line x environment combinations and find
## the proportion that are observed for each trait
## If a trait was not observed in an entire environment, that environment is not
## included in these calculations
(prob_observed <- S2_MET_BLUEs %>% 
    distinct(trait, environment, line_name) %>%
    group_by(trait) %>%
    mutate(observed = TRUE) %>% 
    complete(trait, environment, line_name, fill = list(observed = FALSE)) %>%
    group_by(trait) %>%
    summarize(prop_obs = mean(observed)))



## Find the proportion of unbalance in the dataset
## Do this for the TP, then for the VP, then for both

# Group by trait and complete the observations,
# then group by class and measure the completeness
S2_MET_BLUEs_completeness <- S2_MET_BLUEs %>%
  mutate(observed = TRUE) %>%
  select(line_name, environment, trait, observed) %>% 
  group_by(trait) %>% 
  complete(line_name, environment, fill = list(observed = FALSE)) %>%
  left_join(., select(entry_list, line_name = Line, class = Class)) %>%
  group_by(trait, class) %>%
  summarize(completeness = mean(observed))




## Visualization of distributions

# Sort on grain yield environmental mean
env_order <- S2_MET_BLUEs %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(environment, trait) %>% 
  mutate(env_mean = mean(value, na.rm = TRUE)) %>% 
  filter(trait == "GrainYield") %>% 
  complete(environment) %>%
  arrange(env_mean) %>%
  pull(environment) %>% 
  unique()

S2_MET_BLUEs_toplot <- S2_MET_BLUEs %>%
  mutate(environment = parse_factor(environment, levels = env_order),
         trait = str_replace_all(trait, trait_replace))
  


g_met_dist <- S2_MET_BLUEs_toplot %>%
  ggplot(aes(x = value, y = environment, fill = environment)) +
  geom_density_ridges() +
  facet_grid(. ~ trait, scales = "free_x") +
  scale_fill_discrete(guide = FALSE) +
  ylab("Environment") +
  xlab("Phenotypic Value") +
  labs(title = "Trait Distributions in All Environments") +
  theme_bw()

# Alternatively, plot using boxplot
g_met_boxplot <- S2_MET_BLUEs_toplot %>%
  ggplot(aes(x = environment, y = value, fill = environment)) +
  geom_boxplot(width = 1) +
  facet_grid(trait ~ ., scales = "free_y", switch = "y") +
  scale_fill_discrete(guide = FALSE) +
  ylab("Phenotypic Value") +
  xlab("Environment") +
  labs(title = "Trait Distributions in All Environments") + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1))


# Save it
save_file <- file.path(fig_dir, "met_trait_dist.jpg")
ggsave(filename = save_file, plot = g_met_dist, width = 5, height = 5, dpi = 1000)

# Save it
save_file <- file.path(fig_dir, "met_trait_dist_boxplot.jpg")
ggsave(filename = save_file, plot = g_met_boxplot, width = 9, height = 7, dpi = 1000)





## Stage-Two analysis

# Mixed model formula - everything is random!
random_terms <- c("(1|line_name)", "(1|location)", "(1|year)", "(1|line_name:location)",
                  "(1|line_name:year)", "(1|line_name:location:year)")

# Create formulas
reduced_forms <- combn(x = random_terms, m = length(random_terms) - 1, paste, simplify = FALSE, collapse = " + ") %>%
  map(~str_c("value ~ ", .) %>% as.formula) %>%
  set_names(rev(random_terms))

# Add the full formula
forms <- c(full = as.formula(str_c("value ~ ", str_c(random_terms, collapse = " + "))),
           reduced_forms)


# Load a new optimizer
library(optimx)

# Group by trait and fit the multi-environment model
stage_two_fits <- S2_MET_BLUEs %>%
  mutate_at(vars(location:line_name), as.factor) %>%
  group_by(trait) %>%
  do({
    
    df <- droplevels(.)
    
    # Table of lines by environments (i.e. plots)
    plot_table <- xtabs(formula = ~ line_name + location + year, data = df)
    
    ## Harmonic means
    # Locations
    harm_loc <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Year
    harm_year <- apply(X = plot_table, MARGIN = c(1,3), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Reps
    harm_rep <- apply(X = plot_table, MARGIN = c(1,2,3), sum) %>% 
      harm_mean()

    # Lmer control
    lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore",
                                calc.derivs = FALSE,
                                optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))
    
    # lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
        
    # Get the weights
    wts <- pull(df, std_error)^2
    
    # Fit the models and return
    fits <- fit_with(data = df, lmer, forms, weights = wts, control = lmer_control)
    
    ## Likelihood ratio tests
    lrt <- fits %>% 
      map_df(~logLik(.) %>% as.numeric) %>% gather(reduced_model_term, log_lik) %>% 
      mutate(statistic = (-2) * (log_lik - log_lik[1]), 
             pvalue = pchisq(q = statistic, df = 1, lower.tail = FALSE) / 2,
             variation_source = str_replace_all(string = reduced_model_term, pattern = "\\(1\\||\\)", "")) %>%
      select(variation_source, statistic, pvalue) %>%
      filter(variation_source != "full")
    
    # Calculate heritability
    h2 <- herit(object = fits$full, n_l = harm_loc, n_y = harm_year, n_r = harm_rep,
                exp = "line_name / (line_name + (line_name:year / n_y) + (line_name:location / n_l) + 
                (line_name:location:year / (n_l * n_y)) + (Residual / (n_l * n_y * n_r)))")
    
    # Return data_frame
    data_frame(fit = list(fits$full), lrt = list(lrt), h2 = h2, n_l = harm_loc, 
               n_y = harm_year, n_r = harm_rep) })



# Extract the full model and calculate heritability
stage_two_herit <- stage_two_fits %>% 
  unnest(out) %>%
  ungroup() %>% 
  mutate(full_fit = lapply(.$fits, "[[", "full")) %>%
  group_by(trait) %>%
  do({
    # Calculate heritability
    suppressWarnings(herit_boot(object = .$full_fit[[1]], 
                        exp = "line_name / (line_name + (line_name:environment / n_e) + (Residual / n_r))", 
                        boot.reps = 500, n_e = .$n_e, n_r = .$n_r)) })
    

# Plot
stage_two_herit %>% 
  mutate(H_corr = heritability - bias) %>%
  ggplot(aes(x = trait, y = H_corr, ymin = ci_lower, ymax = ci_upper)) +
  geom_col() +
  geom_errorbar(width = 0.5) +
  ylab("Heritability") +
  xlab("Trait") +
  labs(
    title = "Broad-Sense Heritability of Each Trait Across All Environments",
    subtitle = "Estimates have been corrected for bias and error bars reflect\na 95% confidence interval of 500 bootstrap replications.",
    caption = expression("Equation: "~frac(sigma[g]^2, sigma[g]^2 + frac(sigma[ge]^2, n[e]) + 
                                                  frac(sigma[epsilon]^2, n[e]*n[r]))))
  )


  

## What is the proportion of each variance component to the total phenotypic variance?
stage_two_varProp <- stage_two_fits %>% 
  ungroup() %>% 
  mutate(varcor = map(fit, ~as.data.frame(VarCorr(.)))) %>% 
  unnest(varcor) %>% 
  group_by(trait) %>%
  mutate(varProp = vcov / sum(vcov)) %>% 
  select(trait, h2, variation = grp, varProp)

# Plot this
# Plot
g_varprop <- stage_two_varProp %>%
  # Re-order the levels
  mutate(variation = factor(variation, levels = c("line_name", "location", "year", "line_name:location",
                                                  "line_name:year", "line_name:location:year", "Residual"))) %>%
  arrange(trait, desc(variation)) %>%
  ggplot(aes(x = trait, y = varProp, fill = variation)) +
  geom_bar(stat = "identity") +
  ylab("Proportion") +
  xlab("Trait") + 
  scale_fill_discrete(guide = guide_legend(title = "Variance Component")) +
  theme_bw()

# Save this plot
save_file <- file.path(fig_dir, "trait_varprop.jpg")
ggsave(filename = save_file, plot = g_varprop, height = 6, width = 8)


# Save the stage-two results
save_file <- file.path(result_dir, "stage_two_data.RData")
save("stage_two_fits", "stage_two_varProp", file = save_file)





## Environmental Correlations

## Estimate genetic correlations using the BLUEs from each environment

# Use all data
env_cor <- S2_MET_BLUEs %>% 
  split(.$trait) %>% 
  map(~select(., line_name, environment, value) %>% spread(environment, value) %>% 
        as.data.frame() %>% remove_rownames() %>% column_to_rownames("line_name") %>% 
        cor(., use = "pairwise.complete.obs"))

# Convert to a data.frame for visualization 
env_cor_df <- env_cor %>% 
  map(~as.data.frame(.) %>% rownames_to_column("environment1") %>% 
        gather(environment2, correlation, -environment1) %>% 
        filter(environment1 != environment2)) %>%
  list(., names(.)) %>%
  pmap_df(~mutate(.x, trait = .y))

# Plot
g_env_cor <- env_cor_df %>%
  ggplot(aes(x = environment1, y = environment2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "transparent", 
                       breaks = c(-1, 0, 1), limits = c(-1, 1), 
                       guide = guide_colorbar(title = "Genetic Correlation")) +
  ylab("Environment 2") +
  xlab("Environment 1") +
  facet_wrap(~ trait, nrow = 2, ncol = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = c(0.75, 0.25),
        legend.direction = "horizontal")

# Save this plot
save_file <- file.path(fig_dir, "environmental_correlation.jpg")
ggsave(filename = save_file, plot = g_env_cor, height = 8, width = 8)




## Save the correlation matrices for further use

save_file <- file.path(result_dir, "environmental_genetic_correlations.RData")
save("env_cor_df", "env_cor", file = save_file)

