## S2MET Phenotypic Data Summary
## 
## Author: Jeff Neyhart
## Last updated: May 2, 2018
## 

# Load libraries and directories
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

## Combine data
S2_MET_BLUEs_tomodel <- bind_rows(mutate(S2_MET_BLUEs, population = "all"), 
                                  mutate(filter(S2_MET_BLUEs, line_name %in% tp), population = "tp")) 


# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP
stage_two_fits_GE <- S2_MET_BLUEs_tomodel %>%
  mutate_at(vars(location, year, line_name), as.factor) %>%
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




### Calculate the proportion of GxE that is due to environmental genetic variance
### heterogeneity versus lack of environmental correlation

control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")

# For each environment, calculate the genetic variance via reml
env_varG <- S2_MET_BLUEs_tomodel %>% 
  group_by(population, trait, environment) %>%
  do(varG = {
    df <- .
    wts <- df$std_error^2
    
    # If more than one trial is present, add trial as a fixed effect
    if (n_distinct(df$trial) > 1) {
      formula <- value ~ 1 + trial + (1|line_name)
    } else {
      formula <- value ~ 1 + (1|line_name)
    }
    
    fit <- lmer(formula = formula, data = df, control = control, weights = wts)  
    
    as.data.frame(VarCorr(fit))[1,"vcov"]
    
  })

# Calculate the heterogeneity
env_varG_V <- env_varG %>% 
  group_by(population, trait) %>% 
  unnest() %>% 
  summarize(V = var(sqrt(varG)))


## Fit a model using all the data
prop_varcomp <- S2_MET_BLUEs_tomodel %>%
  group_by(population, trait) %>%
  do({
    df <- .
    wts <- df$std_error^2
    
    formula <- value ~ 1 + (1|line_name) + environment + (1|line_name:environment)
    fit <- lmer(formula = formula, data = df, control = control, weights = wts)  
    
    as.data.frame(VarCorr(fit))[,c("grp", "vcov")]
    
  })


# Use the estimate of varGE across all environments to calculate L
env_L <- left_join(env_varG_V, subset(prop_varcomp, grp == "line_name:environment",c(population, trait, vcov))) %>% 
  mutate(L = vcov - V)

# Use the estimate of genetic variance across all environments to calculate the 
# genetic correlation
env_r <- left_join(env_L, subset(prop_varcomp, grp == "line_name", c(population, trait, vcov)), by = c("population", "trait")) %>% 
  mutate(r_G = vcov.y / (vcov.y + L))

env_r %>% select(population, trait, r_G) %>% spread(population, r_G)
# trait             all    tp
# 1 GrainYield      0.262 0.273
# 2 HeadingDate     0.614 0.703
# 3 HeadingDateAGDD 0.614 0.697
# 4 PlantHeight     0.343 0.342


## What proportion do V and L make up of varGE?
## This is from Li et al 2018 or Cooper and DeLacey 1994
## Add to the variance component table
varGE_components <- env_r %>% 
  select(population, trait, varGE = vcov.x, V, L) %>% 
  mutate_at(vars(V, L), funs(prop = . / varGE)) 


varGE_components %>%
  mutate(heterogeneity = str_c(round(V, 3), " (", round(V_prop, 2) * 100, "%)"), 
         lackCorrelation = str_c(round(L, 3), " (", round(L_prop, 2) * 100, "%)")) %>% 
  select(population, trait, heterogeneity, lackCorrelation) %>% 
  gather(grp, value, -trait, -population) %>% 
  spread(grp, value)


# population trait           heterogeneity  lackCorrelation 
# 1 all        GrainYield      28018.464 (9%) 276021.075 (91%)
# 2 all        HeadingDate     0.98 (12%)     6.982 (88%)     
# 3 all        HeadingDateAGDD 1110.434 (11%) 8749.777 (89%)  
# 4 all        PlantHeight     1.975 (9%)     18.863 (91%)    
# 5 tp         GrainYield      26800.803 (9%) 277236.683 (91%)
# 6 tp         HeadingDate     1.214 (18%)    5.603 (82%)     
# 7 tp         HeadingDateAGDD 1342.66 (16%)  7200.932 (84%)  
# 8 tp         PlantHeight     2.018 (10%)    18.461 (90%) 


# Plot
g_varGE_comp <- varGE_components %>% 
  select(population, trait, GeneticHeterogen = V_prop, LackCorrelation = L_prop) %>% 
  gather(group, proportion, -population, -trait) %>% 
  mutate(group = factor(group, levels = rev(unique(group)))) %>%
  ggplot(data = ., aes(x = trait, y = proportion, fill = group)) + 
  geom_col() +
  geom_text(aes(y = proportion / 2, label = round(proportion, 2))) +
  scale_fill_discrete(name = NULL) +
  facet_grid(~ population) +
  theme_acs() +
  theme(legend.position = c(0.87, 0.75))

ggsave(filename = "varGE_components.jpg", plot = g_varGE_comp, width = 6, height = 3, path = fig_dir, dpi = 1000)






# Mixed model formula - everything is random!
# Decompose environments into locations and years
random_terms <- c("(1|line_name)", "(1|location)", "(1|year)", "(1|line_name:location)", "(1|line_name:year)", "(1|line_name:location:year)")


# Create formulas
reduced_forms <- combn(x = random_terms, m = length(random_terms) - 1, paste, simplify = FALSE, collapse = " + ") %>%
  map(~str_c("value ~ ", .) %>% as.formula) %>%
  set_names(rev(random_terms))

# Add the full formula
forms <- c(full = as.formula(str_c("value ~ ", str_c(random_terms, collapse = " + "))),
           reduced_forms)

## Combine data
S2_MET_BLUEs_tomodel <- bind_rows(mutate(S2_MET_BLUEs, population = "all"), 
                                  mutate(filter(S2_MET_BLUEs, line_name %in% tp), population = "tp")) 

# Lmer control
lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")

# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP
stage_two_fits_GE <- S2_MET_BLUEs_tomodel %>%
  mutate_at(vars(location, year, line_name), as.factor) %>%
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
    harm_year <- apply(X = plot_table, MARGIN = c(2,3), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Reps
    harm_rep <- apply(X = plot_table, MARGIN = c(1,2,3), sum) %>% 
      harm_mean()
    
    # Get the weights
    wts <- df$std_error^2
    
    # Fit the models and return
    fits <- fit_with(data = df, lme4::lmer, forms, weights = wts, control = lmer_control)
    
    ## Likelihood ratio tests
    lrt <- fits[-1] %>% 
      # map(~structure(., class = "lmerMod")) %>%
      map(~lr_test(model1 = fits$full, model2 = .)) %>%
      list(., names(.)) %>%
      pmap_df(~mutate(.x, term = .y)) %>%
      mutate(variation_source = str_replace_all(string = term, pattern = "\\(1\\||\\)", "")) %>%
      select(variation_source, df:p_value)
    
    # Calculate heritability
    h2 <- herit(object = fits$full, n_l = harm_loc, n_y = harm_year, n_r = harm_rep,
                exp = "line_name / (line_name + (line_name:year / n_y) + (line_name:location / n_l) + 
                (line_name:location:year / (n_l * n_y)) + (Residual / (n_l * n_y * n_r)))")
    
    # Return data_frame
    data_frame(fit = list(fits$full), lrt = list(lrt), h2 = list(h2), n_l = harm_loc, 
               n_y = harm_year, n_r = harm_rep)
    
  })

stage_two_fits_GE %>% distinct(trait, population, lrt) %>% unnest() %>% select(-df, -statistic) %>% spread(variation_source, p_value)

## Plot the LRT results
# trait           population line_name `line_name:location` `line_name:location:year` `line_name:year` location      year
# 1 GrainYield      all        2.03e- 61                1                     0.                   0.974        0 6.78e- 27
# 2 GrainYield      tp         3.44e- 52                1                     0.                   1            0 1.78e- 18
# 3 HeadingDate     all        5.70e-146                1                     0.                   1            0 2.45e-109
# 4 HeadingDate     tp         1.92e-137                1                     0.                   0.932        0 2.57e-113
# 5 HeadingDateAGDD all        6.56e-157                1                     0.                   1            0 3.36e- 40
# 6 HeadingDateAGDD tp         2.16e-149                0.929                 0.                   0.884        0 3.30e- 56
# 7 PlantHeight     all        9.06e- 80                0.950                 5.28e-235            0.950        0 1.32e-111
# 8 PlantHeight     tp         2.25e- 64                0.945                 3.08e-171            0.931        0 3.10e- 94

# Line name, GxLxY, L, and Y are significant
# 

stage_two_fits_GE <- stage_two_fits_GE %>% mutate(h2 = list(h2)) %>% distinct()


## Plot the proportion of variance from each source
stage_two_fits_GE_varprop <- stage_two_fits_GE %>%
  unnest(h2) %>% 
  filter(map_lgl(h2, is.data.frame)) %>% 
  unnest() %>% 
  group_by(trait, population) %>% 
  mutate(var_prop = variance / sum(variance), source = str_replace_all(source, ":", " x "),
         source = factor(source, levels = c("line_name", "location", "year", "line_name x location", "line_name x year",
                                            "line_name x location x year", "Residual")))

## Plot
g_varprop <- stage_two_fits_GE_varprop %>% 
  # filter(population == "all") %>%
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  facet_grid(~ population, labeller = labeller(population = str_to_upper)) +
  ylab("Proportion of Variance") +
  scale_fill_brewer(palette = "Set2") +
  theme_acs() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

# Save
ggsave(filename = "variance_components_decompose_env.jpg", plot = g_varprop, path = fig_dir, width = 6, height = 4, dpi = 1000)


## Look at the heritability and plot
g_herit <- stage_two_fits_GE %>% 
  mutate(h2 = map_dbl(h2, 1)) %>% 
  ggplot(aes(x = trait, y = h2)) + 
  geom_col(fill = "grey65") + 
  geom_text(aes(y = 0.25, label = round(h2, 2))) +
  facet_grid(~ population, labeller = labeller(population = str_to_upper)) +
  theme_acs()

ggsave(filename = "heritability_decompose_env.jpg", plot = g_herit, path = fig_dir, width = 4, height = 4, dpi = 1000)

    

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


















