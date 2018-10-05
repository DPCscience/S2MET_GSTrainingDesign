## S2MET Phenotypic Data Summary
## 
## Author: Jeff Neyhart
## Last updated: May 2, 2018
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
# Load a new optimizer
library(optimx)
library(lme4qtl)

# The head directory
repo_dir <- getwd()

source(file.path(repo_dir, "source.R"))

# Load the distance matrices
load(file.path(result_dir, "distance_method_results.RData"))


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


# trait           class completeness
# 1 GrainYield      S2C1R        0.964
# 2 GrainYield      S2TP         0.840
# 3 HeadingDate     S2C1R        0.932
# 4 HeadingDate     S2TP         0.893
# 5 HeadingDateAGDD S2C1R        0.932
# 6 HeadingDateAGDD S2TP         0.893
# 7 PlantHeight     S2C1R        0.964
# 8 PlantHeight     S2TP         0.895

# The VP data appears to be more balanced, but this is probably because there are more
# environments in which only the VP was grown than environments in which only the TP was grown



## Visualization of distributions

# Sort on grain yield environmental mean
env_order <- S2_MET_BLUEs %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(environment, trait) %>% 
  mutate(env_mean = mean(value, na.rm = TRUE)) %>% 
  filter(trait == "HeadingDate") %>% 
  complete(environment) %>%
  arrange(env_mean) %>%
  pull(environment) %>% 
  unique()

S2_MET_BLUEs_toplot <- S2_MET_BLUEs %>%
  mutate(environment = parse_factor(environment, levels = env_order))
  


g_met_dist <- S2_MET_BLUEs_toplot %>%
  ggplot(aes(x = value, y = environment, fill = environment)) +
  geom_density_ridges() +
  facet_grid(. ~ trait, scales = "free_x") +
  scale_fill_discrete(guide = FALSE) +
  ylab("Environment") +
  xlab("Phenotypic Value") +
  theme_acs()

# Save it
ggsave(filename = "met_trait_dist_agdd.jpg", plot = g_met_dist, path = fig_dir, width = 6, height = 5, dpi = 1000)


## Remove HD in AGDD
g_met_dist <- S2_MET_BLUEs_toplot %>%
  filter(trait != "HeadingDateAGDD") %>%
  ggplot(aes(x = value, y = environment, fill = environment)) +
  geom_density_ridges() +
  facet_grid(. ~ trait, scales = "free_x") +
  scale_fill_discrete(guide = FALSE) +
  ylab("Environment") +
  xlab("Phenotypic Value") +
  theme_acs()

# Save it
ggsave(filename = "met_trait_dist.jpg", plot = g_met_dist, path = fig_dir, width = 4.5, height = 5, dpi = 1000)






## Stage-Two analysis

# Fit a random model where environment is defined as location-year combinations
random_terms <- c("(1|line_name)", "(1|environment)", "(1|line_name:environment)")

# Create formulas
reduced_forms <- combn(x = random_terms, m = length(random_terms) - 1, paste, simplify = FALSE, collapse = " + ") %>%
  map(~str_c("value ~ ", .) %>% as.formula) %>%
  set_names(rev(random_terms))

# Add the full formula
forms <- c(full = as.formula(str_c("value ~ ", str_c(random_terms, collapse = " + "))),
           reduced_forms)


# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP
stage_two_fits_GE <- bind_rows(mutate(S2_MET_BLUEs, population = "all"), mutate(filter(S2_MET_BLUEs, line_name %in% tp), population = "tp")) %>%
  mutate_at(vars(location:line_name), as.factor) %>%
  group_by(trait, population) %>%
  do({
    
    df <- droplevels(.)
    
    # Table of lines by environments (i.e. plots)
    plot_table <- xtabs(formula = ~ line_name + environment, data = df)
    
    ## Harmonic means
    # Locations
    harm_env <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Reps
    harm_rep <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
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
    lrt <- fits[-1] %>% 
      map(~lr_test(model1 = fits$full, model2 = .)) %>%
      list(., names(.)) %>%
      pmap_df(~mutate(.x, term = .y)) %>%
      mutate(variation_source = str_replace_all(string = term, pattern = "\\(1\\||\\)", "")) %>%
      select(variation_source, df:p_value)
    
    # Calculate heritability
    h2 <- herit(object = fits$full, n_e = harm_env, n_r = harm_rep,
                exp = "line_name / (line_name + (environment / n_e) + (Residual / (n_e * n_r)))")
    
    # Return data_frame
    data_frame(fit = list(fits$full), lrt = list(lrt), h2 = list(h2), n_e = harm_env, n_r = harm_rep) 
    
    })


## What is the significance of each variance component?
# trait             n_e   n_r variation_source         df statistic   p_value
# 1 GrainYield       29.3     1 line_name:environment     1     4654. 0.       
# 2 GrainYield       29.3     1 environment               1    16847. 0.       
# 3 GrainYield       29.3     1 line_name                 1     1308. 2.12e-286
# 4 HeadingDate      27.9     1 line_name:environment     1     3749. 0.       
# 5 HeadingDate      27.9     1 environment               1    11437. 0.       
# 6 HeadingDate      27.9     1 line_name                 1     4798. 0.       
# 7 HeadingDateAGDD  27.9     1 line_name:environment     1     3485. 0.       
# 8 HeadingDateAGDD  27.9     1 environment               1     8294. 0.       
# 9 HeadingDateAGDD  27.9     1 line_name                 1     4841. 0.       
# 10 PlantHeight      29.0     1 line_name:environment     1     1860. 0.       
# 11 PlantHeight      29.0     1 environment               1    13568. 0.       
# 12 PlantHeight      29.0     1 line_name                 1     1991. 0. 


# Everything is significant

## Plot the proportion of variance from each source
stage_two_fits_GE_varprop <- stage_two_fits_GE %>%
  unnest(h2) %>% 
  filter(map_lgl(h2, is.data.frame)) %>% 
  unnest() %>% 
  group_by(trait, population) %>% 
  mutate(var_prop = variance / sum(variance),
         source = str_replace_all(source, c("line_name:environment" = "Genotype x Environment",
                                            "environment" = "Environment", "line_name" = "Genotype")),
         source = factor(source, levels = c("Environment", "Genotype", "Genotype x Environment",
                                            "Residual")))

## Plot
g_varprop <- stage_two_fits_GE_varprop %>% 
  filter(population == "all") %>%
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of Variance") +
  scale_fill_manual(values = umn_palette(3, 4), name = NULL) +
  theme_acs() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

# Save
ggsave(filename = "variance_components_AGDD.jpg", plot = g_varprop, path = fig_dir, width = 4, height = 4, dpi = 1000)


# Remove HeadingDate AGDD
## Plot
g_varprop <- stage_two_fits_GE_varprop %>% 
  filter(population == "all", trait != "HeadingDateAGDD") %>%
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of Variance") +
  scale_fill_manual(values = umn_palette(3, 4), name = NULL) +
  theme_acs() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

# Save
ggsave(filename = "variance_components.jpg", plot = g_varprop, path = fig_dir, width = 4, height = 4, dpi = 1000)




### Plot just the TP
g_varprop <- stage_two_fits_GE_varprop %>% 
  filter(population == "tp") %>%
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of Variance") +
  scale_fill_manual(values = umn_palette(3, 4), name = NULL) +
  theme_acs() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

# Save
ggsave(filename = "variance_components_AGDD_tp.jpg", plot = g_varprop, path = fig_dir, width = 4, height = 4, dpi = 1000)


# Remove HeadingDate AGDD
## Plot
g_varprop <- stage_two_fits_GE_varprop %>% 
  filter(population == "tp", trait != "HeadingDateAGDD") %>%
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of Variance") +
  scale_fill_manual(values = umn_palette(3, 4), name = NULL) +
  theme_acs() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

# Save
ggsave(filename = "variance_components_tp.jpg", plot = g_varprop, path = fig_dir, width = 4, height = 4, dpi = 1000)


## Look at the heritability and plot
g_herit <- stage_two_fits_GE %>% 
  mutate(h2 = map_dbl(h2, "heritability")) %>% 
  ggplot(aes(x = trait, y = h2)) + 
  geom_col(fill = "grey65") + 
  geom_text(aes(y = 0.25, label = round(h2, 2))) +
  facet_grid(~ population, labeller = labeller(population = str_to_upper)) +
  theme_acs()

ggsave(filename = "heritability.jpg", plot = g_herit, path = fig_dir, width = 4, height = 4, dpi = 1000)







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



# Group by trait and fit the multi-environment model
stage_two_fits <- bind_rows(mutate(S2_MET_BLUEs, population = "all"), mutate(filter(S2_MET_BLUEs, line_name %in% tp), population = "tp")) %>%
  mutate_at(vars(location:line_name), as.factor) %>%
  group_by(trait, population) %>%
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
## First estimate using only the TP
tp_BLUEs <- S2_MET_BLUEs %>%
  filter(line_name %in% tp_geno)

env_cor_tp <- tp_BLUEs %>% 
  split(.$trait) %>% 
  map(~select(., line_name, environment, value) %>% spread(environment, value) %>% 
        as.data.frame() %>% remove_rownames() %>% column_to_rownames("line_name") %>% 
        cor(., use = "pairwise.complete.obs"))

# Now use all data
env_cor_all <- S2_MET_BLUEs %>% 
  split(.$trait) %>% 
  map(~select(., line_name, environment, value) %>% spread(environment, value) %>% 
        as.data.frame() %>% remove_rownames() %>% column_to_rownames("line_name") %>% 
        cor(., use = "pairwise.complete.obs"))


# Convert to a data.frame for visualization 
env_cor_df <- env_cor_all %>% 
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
save("env_cor_tp", "env_cor_all", file = save_file)


















