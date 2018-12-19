## S2MET Phenotypic Data Summary
## 
## Author: Jeff Neyhart
## Last updated: November 15, 2018
## 



# The head directory
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load libraries and directories
library(broom)
library(modelr)
library(pbr)
library(rrBLUP)
library(ggridges)
# Load a new optimizer
library(optimx)
library(lme4qtl)
library(cowplot)
## ggrepel
library(ggrepel)


# Load the distance matrices
load(file.path(result_dir, "distance_method_results.RData"))


## significance level
alpha <- 0.05


## Filter the BLUEs for the environments of interest
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  filter(environment %in% tp_vp_env,
         trait %in% traits)




## Basic Summaries

## Number of lines per environment per trait

## Find the total number of possible line x environment combinations and find
## the proportion that are observed for each trait
## If a trait was not observed in an entire environment, that environment is not
## included in these calculations
(prob_observed <- S2_MET_BLUEs_use %>% 
    distinct(trait, environment, line_name) %>%
    mutate_at(vars(line_name, environment), as.factor) %>%
    mutate(observed = TRUE) %>%
    split(.$trait) %>%
    map_df(~{
      droplevels(.) %>%
        complete(trait, environment, line_name, fill = list(observed = FALSE)) %>%
        summarize(trait = unique(trait), prop_obs = mean(observed))
    }))

# trait           prop_obs
# 1 GrainYield         0.987
# 2 HeadingDate        0.991
# 3 HeadingDateAGDD    0.991
# 4 PlantHeight        0.990


## Find the proportion of unbalance in the dataset for the TP (tp and tp_geno)
(prob_observed <- list(tp = tp, tp_geno = tp_geno) %>% 
    map(~{
       filter(S2_MET_BLUEs_use, line_name %in% .) %>%
        distinct(trait, environment, line_name) %>%
        mutate_at(vars(line_name, environment), as.factor) %>%
        mutate(observed = TRUE) %>%
        split(.$trait) %>%
        map_df(~{
          droplevels(.) %>%
            complete(trait, environment, line_name, fill = list(observed = FALSE)) %>%
            summarize(trait = unique(trait), prop_obs = mean(observed))
        })
    }))


# $`tp`
# trait           prop_obs
# 1 GrainYield         0.985
# 2 HeadingDate        0.989
# 3 HeadingDateAGDD    0.989
# 4 PlantHeight        0.988
# 
# $tp_geno
# trait           prop_obs
# 1 GrainYield         0.985
# 2 HeadingDate        0.989
# 3 HeadingDateAGDD    0.989
# 4 PlantHeight        0.988

## Create a table with the proportions of observed combinations for the TP and VP
prop_observed_table <- data_frame(population = c("tp", "vp"), line_name = list(tp, vp)) %>%
  mutate(prop_obs = map(line_name, ~{
    filter(S2_MET_BLUEs_use, line_name %in% .) %>%
      distinct(trait, environment, line_name) %>%
      mutate_at(vars(line_name, environment), as.factor) %>%
      mutate(observed = TRUE) %>%
      split(.$trait) %>%
      map_df(~{
        droplevels(.) %>%
          complete(trait, environment, line_name, fill = list(observed = FALSE)) %>%
          summarize(trait = unique(trait), prop_obs = mean(observed))
      }) }) ) %>%
  unnest(prop_obs)

prop_observed_table %>% 
  spread(population, prop_obs) %>% 
  write_csv(x = ., path = file.path(fig_dir, "population_prop_obs.csv"))



## Range of heritabilities
stage_one_data %>% 
  filter(environment %in% tp_vp_env) %>% 
  group_by(trait) %>% 
  summarize_at(vars(heritability), funs(min, max, mean))

# trait              min   max  mean
# 1 GrainYield      0      0.856 0.441
# 2 HeadingDate     0.563  0.978 0.846
# 3 HeadingDateAGDD 0.569  0.973 0.843
# 4 PlantHeight     0.0849 0.883 0.52




## Visualization of distributions

# Sort on grain yield environmental mean
env_order <- S2_MET_BLUEs_use %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(environment, trait) %>% 
  mutate(env_mean = mean(value, na.rm = TRUE)) %>% 
  filter(trait == "HeadingDate") %>% 
  complete(environment) %>%
  arrange(env_mean) %>%
  pull(environment) %>% 
  unique()

S2_MET_BLUEs_toplot <- S2_MET_BLUEs_use %>%
  mutate(environment = parse_factor(environment, levels = env_order))
  

## Plot
g_met_dist <- S2_MET_BLUEs_toplot %>%
  ggplot(aes(x = value, y = environment, fill = environment)) +
  geom_density_ridges() +
  facet_grid(. ~ trait, scales = "free_x") +
  scale_fill_discrete(guide = FALSE) +
  ylab("Environment") +
  xlab("Phenotypic value") +
  theme_acs()

# Save it
ggsave(filename = "met_trait_dist.jpg", plot = g_met_dist, path = fig_dir, width = 4.5, height = 5, dpi = 1000)






## Stage-Two analysis


## Combine data
S2_MET_BLUEs_tomodel <- bind_rows(mutate(S2_MET_BLUEs_use, population = "all"), 
                                  mutate(filter(S2_MET_BLUEs_use, line_name %in% tp), population = "tp"),
                                  mutate(filter(S2_MET_BLUEs_use, line_name %in% vp), population = "vp"))

# Boot reps
boot_reps <- 100

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
    
    # # Lmer control
    # lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", calc.derivs = FALSE,
    #                             optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))
    
    lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
    
    # Get the weights
    wts <- pull(df, std_error)^2
    
    ## Fit the full model
    fit <- lmer(formula = value ~ 1 + (1|line_name) + (1|environment) + (1|line_name:environment), 
                data = df, control = lmer_control, weights = wts)
    
    ## Likelihood ratio tests
    lrt <- ranova(fit) %>% 
      tidy() %>% 
      filter(!str_detect(term, "none")) %>% 
      mutate(term = str_remove(term, "\\(1 \\| ") %>% 
               str_remove("\\)")) %>% 
      select(term, LRT, df, p.value)
    
    ## Calculate heritability
    exp <- "line_name / (line_name + (line_name:environment / n_e) + (Residual / (n_e * n_r)))"
    
    ## Use bootstrapping to calculate a confidence interval
    # Generate bootstrapping samples and calculate heritability using the bootMer function
    h2_boot <- bootMer(x = fit, nsim = boot_reps, FUN = function(x) 
      herit(object = x, exp = exp, n_e = harm_env, n_r = harm_rep)$heritability)
    
    h2 <- herit(object = fit, n_e = harm_env, n_r = harm_rep, exp = exp)
    
    # Add the bootstrapping results with a confidence interval
    h2$heritability <- tidy(h2_boot) %>% 
      cbind(., t(quantile(h2_boot$t, probs = c(alpha / 2, 1 - (alpha / 2))))) %>% 
      rename_at(vars(4, 5), ~c("lower", "upper"))
    
    # Return data_frame
    data_frame(fit = list(fit), lrt = list(lrt), h2 = list(h2), n_e = harm_env, n_r = harm_rep) 
    
  }) %>% ungroup()


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
  mutate(h2 = map(h2, "var_comp")) %>%
  unnest(h2) %>% 
  group_by(trait, population) %>% 
  mutate(var_prop = variance / sum(variance),
         source = str_replace_all(source, c("line_name:environment" = "Genotype x Environment",
                                            "environment" = "Environment", "line_name" = "Genotype")),
         source = factor(source, levels = c("Environment", "Genotype", "Genotype x Environment",
                                            "Residual"))) %>%
  ungroup()

## Plot both populations and all traits
g_varprop <- stage_two_fits_GE_varprop %>% 
  mutate(population = str_to_upper(population)) %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of Variance") +
  scale_fill_manual(values = umn_palette(3, 4), name = NULL) +
  facet_grid(~ population) + 
  theme_acs() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

# Save
ggsave(filename = "variance_components.jpg", plot = g_varprop, path = fig_dir, width = 8, height = 4, dpi = 1000)


## Plot just the TP and remove AGDD heading date
g_varprop_tp <- stage_two_fits_GE_varprop %>% 
  filter(population == "tp") %>%
  mutate(population = str_to_upper(population)) %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of Variance") +
  scale_fill_manual(values = umn_palette(3, 4), name = NULL) +
  theme_acs() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

# Save
ggsave(filename = "variance_components_tp.jpg", plot = g_varprop_tp, path = fig_dir, width = 4, height = 4, dpi = 1000)


## Look at the heritability and plot
g_herit <- stage_two_fits_GE %>% 
  mutate(h2 = map(h2, "heritability"), population = str_to_upper(population)) %>% 
  unnest(h2) %>% 
  mutate(statistic = statistic + bias) %>%
  ggplot(aes(x = trait, y = statistic, ymin = lower, ymax = upper, fill = population)) + 
  geom_col(position = position_dodge(0.9)) + 
  geom_errorbar(position = position_dodge(0.9), width = 0.5) +
  geom_text(aes(y = 0.60, label = round(statistic, 2)), position = position_dodge(0.9)) +
  # geom_hline(aes(yintercept = unbiased_statistic)) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  ylab("Heritability") +
  xlab("Trait") +
  labs(caption = paste0("Error bars represent a 95% confidence interval\ncalculated using ", boot_reps, " bootsrapping replicates.")) + 
  theme_acs() +
  theme(legend.position = c(0.20, 0.85), legend.direction = "horizontal")

ggsave(filename = "heritability.jpg", plot = g_herit, path = fig_dir, width = 6, height = 4, dpi = 1000)


# Plot just the TP
g_herit_tp <- stage_two_fits_GE %>% 
  filter(population == "tp") %>%
  mutate(h2 = map(h2, "heritability"), population = str_to_upper(population)) %>% 
  unnest(h2) %>% 
  mutate(statistic = statistic + bias) %>%
  ggplot(aes(x = trait, y = statistic, ymin = lower, ymax = upper, fill = population)) + 
  geom_col(position = position_dodge(0.9)) + 
  geom_errorbar(position = position_dodge(0.9), width = 0.5) +
  geom_text(aes(y = 0.60, label = round(statistic, 2)), position = position_dodge(0.9)) +
  scale_fill_brewer(palette = "Set2", name = "Population", guide = FALSE) +
  ylab("Heritability") +
  xlab("Trait") +
  labs(caption = paste0("Error bars represent a 95% confidence interval\ncalculated using ", boot_reps, " bootsrapping replicates.")) + 
  theme_acs() +
  theme(legend.position = c(0.15, 0.20))

ggsave(filename = "heritability_tp.jpg", plot = g_herit_tp, path = fig_dir, width = 4, height = 4, dpi = 1000)







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
    
  }) %>% ungroup() %>%
  mutate(varG = unlist(varG))

# Calculate the genetic heterogeneity, V
env_varG_V <- env_varG %>% 
  group_by(population, trait) %>% 
  summarize(V = var(sqrt(varG))) %>%
  ungroup()


## Use the variance components estinated in the previous random model
prop_varcomp <- stage_two_fits_GE %>%
  mutate(varcomp = map(h2, "var_comp")) %>% 
  unnest(varcomp)


# Use the estimate of varGE across all environments to calculate L
env_L <- left_join(env_varG_V, subset(prop_varcomp, source == "line_name:environment", c(population, trait, variance))) %>% 
  rename(varGE = variance) %>%
  mutate(L = varGE - V)

# Use the estimate of genetic variance across all environments to calculate the 
# genetic correlation
env_r <- left_join(env_L, subset(prop_varcomp, source == "line_name", c(population, trait, variance)), by = c("population", "trait")) %>% 
  rename(varG = variance) %>%
  mutate(r_G = varG / (varG + L))

env_r %>% 
  select(population, trait, r_G) %>% 
  spread(population, r_G)

# trait             all    tp
# 1 GrainYield      0.262 0.273
# 2 HeadingDate     0.614 0.703
# 3 HeadingDateAGDD 0.614 0.697
# 4 PlantHeight     0.343 0.342


## What proportion do V and L make up of varGE?
## This is from Li et al 2018 or Cooper and DeLacey 1994
## Add to the variance component table
varGE_components <- env_r %>% 
  mutate_at(vars(V, L), funs(prop = . / varGE)) 


varGE_components %>%
  mutate(heterogeneity = str_c(round(V, 3), " (", round(V_prop, 2) * 100, "%)"), 
         lackCorrelation = str_c(round(L, 3), " (", round(L_prop, 2) * 100, "%)")) %>% 
  select(population, trait, heterogeneity, lackCorrelation) %>% 
  gather(grp, value, -trait, -population) %>% 
  spread(grp, value)


# population trait           heterogeneity  lackCorrelation 
# 1 all        GrainYield      28018.464 (9%) 276020.293 (91%)
# 2 all        HeadingDate     0.98 (12%)     6.982 (88%)     
# 3 all        HeadingDateAGDD 1110.434 (11%) 8748.997 (89%)  
# 4 all        PlantHeight     1.975 (9%)     18.862 (91%)    
# 5 tp         GrainYield      26800.803 (9%) 277236.568 (91%)
# 6 tp         HeadingDate     1.214 (18%)    5.603 (82%)     
# 7 tp         HeadingDateAGDD 1342.66 (16%)  7200.925 (84%)  
# 8 tp         PlantHeight     2.018 (10%)    18.461 (90%) 


# Plot
g_varGE_comp <- varGE_components %>% 
  select(population, trait, GeneticHeterogen = V_prop, LackCorrelation = L_prop) %>% 
  gather(group, proportion, -population, -trait) %>% 
  mutate(group = factor(group, levels = rev(unique(group))),
         population = str_to_upper(population)) %>%
  ggplot(data = ., aes(x = trait, y = proportion, fill = group)) + 
  geom_col() +
  geom_text(aes(y = proportion / 2, label = round(proportion, 2))) +
  scale_fill_brewer(name = NULL, palette = "Set2") +
  facet_grid(~ population) +
  ylab(expression("Proportion of"~sigma[GE]^2)) +
  xlab("Trait") +
  theme_acs() +
  theme(legend.position = c(0.87, 0.75))

ggsave(filename = "varGE_components.jpg", plot = g_varGE_comp, width = 8, height = 3, path = fig_dir, dpi = 1000)



## Combine the plots for all variance components with the GxE variance components
# Plot
g_varGE_comp_tp <- varGE_components %>% 
  filter(population == "tp") %>%
  select(population, trait, GeneticHeterogen = V_prop, LackCorrelation = L_prop) %>% 
  gather(group, proportion, -population, -trait) %>% 
  mutate(group = factor(group, levels = rev(unique(group))),
         population = str_to_upper(population)) %>%
  ggplot(data = ., aes(x = trait, y = proportion, fill = group)) + 
  geom_col() +
  # geom_text(aes(y = proportion / 2, label = round(proportion, 2))) +
  scale_fill_manual(name = NULL, values = c(neyhart_palette("umn1")[3], neyhart_palette("umn2", 8)[8])) +
  ylab(expression("Proportion of"~sigma[GE]^2)) +
  xlab("Trait") +
  theme_acs() +
  theme(legend.position = "bottom")

g_varcomp_combine <- plot_grid(g_varprop_tp + theme(legend.position = "top"), 
                               g_varGE_comp_tp, ncol = 1, rel_heights = c(1, 0.4), axis = "lr", align = "v")

ggsave(filename = "variance_components_combined.jpg", plot = g_varcomp_combine, path = fig_dir,
       height = 8, width = 5, dpi = 1000)


## Output a table of all variance components for all populations
prop_varcomp1 <- prop_varcomp %>% 
  select(trait, population, source, variance) %>% 
  group_by(trait, population) %>% 
  mutate(proportion = variance / sum(variance),
         source = str_replace_all(source, ":", " x "), 
         source = str_replace_all(source, "line_name", "Genotype"), 
         source = str_to_title(source))
  

varGE_components1 <- varGE_components %>% 
  select(trait, population, V, L) %>% 
  gather(source, variance, V, L) %>%
  group_by(trait, population) %>% 
  mutate(source = ifelse(source == "V", "GeneticHeterogeneity", "LackOfCorrelation"),  
         proportion = variance / sum(variance))

# Component order
comp_order <- c("Genotype", "Environment", "Genotype X Environment", "GeneticHeterogeneity", "LackOfCorrelation", "Residual")

## Combine and output
var_comp_table <- bind_rows(prop_varcomp1, varGE_components1) %>%
  ungroup() %>%
  mutate(annotation = paste0(formatC(x = variance, digits = 3), " (", formatC(proportion * 100, digits = 2, width = 2), "%)"),
         annotation = str_trim(annotation)) %>%
  select(trait, population, source, annotation) %>%
  # spread(source, annotation) %>%
  # select(Trait = trait, Population = population, Genotype, Environment, `Genotype X Environment`, GeneticHeterogeneity, LackOfCorrelation,
  #        Residual)
  spread(trait, annotation) %>%
  mutate(source = factor(source, levels = comp_order)) %>% arrange(population, source) %>%
  rename(Population = population, Source = source)
  
  

write_csv(x = var_comp_table, path = file.path(fig_dir, "population_variance_components.csv"))

## Just the TP
var_comp_table %>% filter(Population == "tp") %>% select(-Population) %>% 
  write_csv(x = ., path = file.path(fig_dir, "tp_population_variance_components.csv"))






### Analyze the components of GxE

## Combine data
S2_MET_BLUEs_tomodel <- bind_rows(mutate(S2_MET_BLUEs_use, population = "all"), 
                                  mutate(filter(S2_MET_BLUEs_use, line_name %in% tp), population = "tp"),
                                  mutate(filter(S2_MET_BLUEs_use, line_name %in% vp), population = "vp")) 

# Lmer control
lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")

# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP
stage_two_fits_GYL <- S2_MET_BLUEs_tomodel %>%
  # filter(location %in% c("St_Paul", "Crookston", "Fargo", "Arlington", "Madison")) %>%
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
    
    ## fit the full model
    fit <- lmer(formula = value ~ 1 + (1|line_name) + (1|location) + (1|year) + (1|line_name:location) + (1|line_name:year) + (1|line_name:location:year),
                data = df, control = lmer_control, weights = wts)
    
    
    ## Likelihood ratio tests
    lrt <- ranova(fit) %>% 
      tidy() %>% 
      filter(!str_detect(term, "none")) %>% 
      mutate(term = str_remove(term, "\\(1 \\| ") %>% 
               str_remove("\\)")) %>% 
      select(term, LRT, df, p.value)
    
    ## Calculate heritability
    # Expression for heritability
    exp <- "line_name / (line_name + (line_name:location / n_l) + (line_name:year / n_y) + (line_name:location:year / (n_l * n_y)) + (Residual / (n_y * n_l * n_r)))"
    
    ## Use bootstrapping to calculate a confidence interval
    # Generate bootstrapping samples and calculate heritability using the bootMer function
    h2_boot <- bootMer(x = fit, nsim = boot_reps, FUN = function(x) 
      herit(object = x, exp = exp, n_l = harm_loc, n_y = harm_year, n_r = harm_rep)$heritability)
    
    h2 <- herit(object = fit, exp = exp, n_l = harm_loc, n_y = harm_year, n_r = harm_rep)
    
    # Add the bootstrapping results with a confidence interval
    h2$heritability <- tidy(h2_boot) %>% 
      cbind(., t(quantile(h2_boot$t, probs = c(alpha / 2, 1 - (alpha / 2))))) %>% 
      rename_at(vars(4, 5), ~c("lower", "upper"))
    
    
    # Return data_frame
    data_frame(fit = list(fit), lrt = list(lrt), h2 = list(h2), n_l = harm_loc, 
               n_y = harm_year, n_r = harm_rep)
    
  }) %>% ungroup()

stage_two_fits_GYL %>% 
  select(trait, population, lrt) %>% 
  unnest() %>% 
  select(-df, -LRT) %>% 
  spread(term, p.value)

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


## Plot the proportion of variance from each source
stage_two_fits_GYL_varprop <- stage_two_fits_GYL %>%
  mutate(h2 = map(h2, "var_comp")) %>%
  unnest(h2) %>% 
  mutate(source = map(source, ~str_split(., pattern = ":", simplify = FALSE) %>% map(str_to_title) %>% .[[1]]) %>% 
           map_chr(~paste(., collapse = " x ")),
         source = str_replace_all(source, "Line_name", "Genotype"),
         source = factor(source, levels = c("Genotype", "Location", "Year", "Genotype x Location", "Genotype x Year",
                                            "Genotype x Location x Year", "Residual"))) %>%
  group_by(trait, population) %>% 
  mutate(var_prop = variance / sum(variance)) %>%
  ungroup()

## Plot both populations and all traits
g_varprop <- stage_two_fits_GYL_varprop %>% 
  mutate(population = str_to_upper(population)) %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of Variance") +
  scale_fill_manual(values = umn_palette(3, 7), name = NULL) +
  facet_grid(~ population) + 
  theme_acs() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

# Save
ggsave(filename = "variance_components_expanded.jpg", plot = g_varprop, path = fig_dir, width = 8, height = 4, dpi = 1000)


## Plot just the TP and remove AGDD heading date
g_varprop_tp <- stage_two_fits_GYL_varprop %>% 
  filter(population == "tp") %>%
  mutate(population = str_to_upper(population)) %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of Variance") +
  scale_fill_manual(values = umn_palette(3, 7), name = NULL, guide = guide_legend(nrow = 3)) +
  theme_acs() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

# Save
ggsave(filename = "variance_components_expanded_tp.jpg", plot = g_varprop_tp, path = fig_dir, width = 4, height = 5, dpi = 1000)


## Look at the heritability and plot
g_herit <- stage_two_fits_GYL %>% 
  mutate(h2 = map(h2, "heritability"), population = str_to_upper(population)) %>% 
  unnest(h2) %>% 
  mutate(statistic = statistic + bias) %>%
  ggplot(aes(x = trait, y = statistic, ymin = lower, ymax = upper, fill = population)) + 
  geom_col(position = position_dodge(0.9)) + 
  geom_errorbar(position = position_dodge(0.9), width = 0.5) +
  geom_text(aes(y = 0.40, label = round(statistic, 2)), position = position_dodge(0.9)) +
  # geom_hline(aes(yintercept = unbiased_statistic)) +
  scale_fill_brewer(palette = "Set1", name = "Population") +
  ylab("Heritability") +
  xlab("Trait") +
  labs(caption = paste0("Error bars represent a 95% confidence interval\ncalculated using ", boot_reps, " bootsrapping replicates.")) + 
  theme_acs() +
  theme(legend.position = c(0.15, 0.85))

ggsave(filename = "heritability_expanded.jpg", plot = g_herit, path = fig_dir, width = 4, height = 4, dpi = 1000)


# Plot just the TP
g_herit <- stage_two_fits_GYL %>% 
  filter(population == "tp") %>%
  mutate(h2 = map(h2, "heritability"), population = str_to_upper(population)) %>% 
  unnest(h2) %>% 
  mutate(statistic = statistic + bias) %>%
  ggplot(aes(x = trait, y = statistic, ymin = lower, ymax = upper, fill = population)) + 
  geom_col(position = position_dodge(0.9)) + 
  geom_errorbar(position = position_dodge(0.9), width = 0.5) +
  geom_text(aes(y = 0.40, label = round(statistic, 2)), position = position_dodge(0.9)) +
  scale_fill_brewer(palette = "Set1", name = "Population") +
  ylab("Heritability") +
  xlab("Trait") +
  labs(caption = paste0("Error bars represent a 95% confidence interval\ncalculated using ", boot_reps, " bootsrapping replicates.")) + 
  theme_acs() +
  theme(legend.position = c(0.15, 0.85))

ggsave(filename = "heritability_expanded_tp.jpg", plot = g_herit, path = fig_dir, width = 4, height = 4, dpi = 1000)




## Combine and output
var_comp_table_gyl <- stage_two_fits_GYL_varprop %>%
  ungroup() %>%
  mutate(annotation = paste0(formatC(x = variance, digits = 3), " (", formatC(var_prop * 100, digits = 2, width = 2), "%)"),
         annotation = str_trim(annotation)) %>%
  select(trait, population, source, annotation) %>%
  # spread(source, annotation) %>%
  # select(Trait = trait, Population = population, Genotype, Environment, `Genotype X Environment`, GeneticHeterogeneity, LackOfCorrelation,
  #        Residual)
  spread(trait, annotation) %>%
  rename(Population = population, Source = source)



write_csv(x = var_comp_table_gyl, path = file.path(fig_dir, "population_variance_components_decomposed.csv"))

## Just the TP
var_comp_table_gyl %>% filter(Population == "tp") %>% select(-Population) %>% 
  write_csv(x = ., path = file.path(fig_dir, "tp_population_variance_components_decomposed.csv"))







    






## Environmental Correlations

## Estimate genetic correlations using the BLUEs from each environment
## First estimate using only the TP with genotypes
tp_BLUEs <- S2_MET_BLUEs %>%
  filter(line_name %in% tp,
         trait %in% traits)

tp_BLUEs_use <- tp_BLUEs %>%
  select(line_name, environment, trait, value)

## Cross all environments
env_shared_tp <- tp_BLUEs_use %>% 
  split(.$trait) %>%
  map(~{
    df <- .
    crossing(environment = unique(df$environment), environment2 = environment) %>%
      filter(environment != environment2) %>%
      mutate(shared = map2_dbl(.x = environment, .y = environment2, .f = ~{
        length(intersect(subset(df, environment == .x, line_name, drop = T), subset(df, environment == .y, line_name, drop = T))) }))
  })
    
env_cor_tp <- tp_BLUEs %>% 
  split(.$trait) %>% 
  map(~select(., line_name, environment, value) %>% spread(environment, value) %>% 
        as.data.frame() %>% remove_rownames() %>% column_to_rownames("line_name") %>% 
        cor(., use = "pairwise.complete.obs"))

## VP
env_cor_vp <- S2_MET_BLUEs %>%
  filter(line_name %in% vp, trait %in% traits) %>%
  select(line_name, environment, trait, value) %>% 
  split(.$trait) %>% 
  map(~select(., line_name, environment, value) %>% spread(environment, value) %>% 
        as.data.frame() %>% remove_rownames() %>% column_to_rownames("line_name") %>% 
        cor(., use = "pairwise.complete.obs"))





# Now use all data
env_cor_all <- S2_MET_BLUEs %>%
  filter(trait %in% traits) %>%
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

env_cor_tp_df <- env_cor_tp %>% 
  map(~as.data.frame(.) %>% rownames_to_column("environment1") %>% 
        gather(environment2, correlation, -environment1) %>% 
        filter(environment1 != environment2)) %>%
  list(., names(.)) %>%
  pmap_df(~mutate(.x, trait = .y))

env_cor_vp_df <- env_cor_vp %>% 
  map(~as.data.frame(.) %>% rownames_to_column("environment1") %>% 
        gather(environment2, correlation, -environment1) %>% 
        filter(environment1 != environment2)) %>%
  list(., names(.)) %>%
  pmap_df(~mutate(.x, trait = .y))



# Plot
# One trait at a time
g_env_cor_list <- env_cor_df %>%
  split(.$trait) %>%
  map(~{
    df <- .
    
    # Use the heatmap function to construct a dendrogram
    mat <- df %>% 
      select(-trait) %>% 
      spread(environment2, correlation) %>% 
      as.data.frame() %>% 
      column_to_rownames("environment1") %>% 
      as.matrix()
    
    hm <- heatmap(mat)
    
    # Extract the ordered factor for environments
    eLevels <- sort(unique(df$environment1))[hm$rowInd]
    
    df %>% 
      mutate_at(vars(contains("environment")), funs(factor(., levels = eLevels))) %>%
      ggplot(., aes(x = environment1, y = environment2, fill = correlation)) +
      geom_tile() +
      scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "transparent", 
                           breaks = seq(-1, 1, by = 0.5), limits = c(-1, 1), 
                           guide = guide_colorbar(title = "Correlation")) +
      ylab("Environment 2") +
      xlab("Environment 1") +
      facet_wrap(~ trait) +
      theme_acs() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
            axis.text.y = element_text(size = 5),
            # legend.position = "bottom", legend.direction = "horizontal")
            legend.position = "right", legend.direction = "vertical")
    
    
  })

## Combine plots
g_env_cor <- plot_grid(plotlist = g_env_cor_list %>% map(~. + theme(legend.position = "none")), nrow = 1)
g_env_cor1 <- plot_grid(g_env_cor, get_legend(g_env_cor_list[[1]]), nrow = 1, rel_widths = c(1, 0.1))


# Save this plot
# ggsave(filename = "environmental_correlation.jpg", path = fig_dir, plot = g_env_cor1, height = 9, width = 3.5, dpi = 1000)
ggsave(filename = "environmental_correlation.jpg", path = fig_dir, plot = g_env_cor1, height = 3, width = 10, dpi = 1000)



# Plot
# One trait at a time
g_env_cor_list_tp <- env_cor_tp_df %>%
  split(.$trait) %>%
  map(~{
    df <- .
    
    # Use the heatmap function to construct a dendrogram
    mat <- df %>% 
      select(-trait) %>% 
      spread(environment2, correlation) %>% 
      as.data.frame() %>% 
      column_to_rownames("environment1") %>% 
      as.matrix()
    
    hm <- heatmap(mat)
    
    # Extract the ordered factor for environments
    eLevels <- sort(unique(df$environment1))[hm$rowInd]
    
    df %>% 
      mutate_at(vars(contains("environment")), funs(factor(., levels = eLevels))) %>%
      ggplot(., aes(x = environment1, y = environment2, fill = correlation)) +
      geom_tile() +
      scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "transparent", 
                           breaks = seq(-1, 1, by = 0.5), limits = c(-1, 1), 
                           guide = guide_colorbar(title = "Correlation")) +
      ylab("Environment 2") +
      xlab("Environment 1") +
      facet_wrap(~ trait) +
      theme_acs() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
            axis.text.y = element_text(size = 5),
            # legend.position = "bottom", legend.direction = "horizontal")
            legend.position = "right", legend.direction = "vertical")
    
    
  })

## Combine plots
g_env_cor <- plot_grid(plotlist = g_env_cor_list_tp %>% map(~. + theme(legend.position = "none")), nrow = 1)
g_env_cor1 <- plot_grid(g_env_cor, get_legend(g_env_cor_list_tp[[1]]), nrow = 1, rel_widths = c(1, 0.1))


# Save this plot
ggsave(filename = "environmental_correlation_tp.jpg", path = fig_dir, plot = g_env_cor1, height = 3, width = 10, dpi = 1000)




## Histogram of pairwise correlations
env_cor_tp_df %>%
  ggplot(aes(x = correlation)) + 
  geom_histogram() + 
  facet_grid(~trait) +
  theme_acs()

## Min, max, mean of correlations
bind_rows(mutate(env_cor_df, population = "all"), mutate(env_cor_tp_df, population = "tp"), mutate(env_cor_vp_df, population = "vp")) %>% 
  group_by(population, trait) %>%
  summarize_at(vars(correlation), funs(min, max, mean), na.rm = T) %>%
  arrange(trait)




# population trait           min   max  mean
# 1 all        GrainYield  -0.389  0.616 0.196
# 2 tp         GrainYield  -0.193  0.616 0.256
# 3 vp         GrainYield  -0.389  0.609 0.144
# 4 all        HeadingDate -0.262  0.878 0.552
# 5 tp         HeadingDate  0.366  0.892 0.708
# 6 vp         HeadingDate -0.262  0.825 0.390
# 7 all        PlantHeight -0.132  0.716 0.342
# 8 tp         PlantHeight -0.0243 0.617 0.329
# 9 vp         PlantHeight -0.260  0.771 0.359




## Save the correlation matrices for further use
save_file <- file.path(result_dir, "environmental_genetic_correlations.RData")
save("env_cor_tp", "env_cor_all", file = save_file)


## Create a two-way table of genotypes and environments
## 
## We will fill in the missing TP geno points using the kronecker
## product of the realized relationship matrix and the environmental
## correlation matrix
## 

# split the TP blues by trait
tp_BLUEs <- S2_MET_BLUEs %>%
  filter(line_name %in% tp_geno, trait %in% traits) 

tp_BLUEs_split <- tp_BLUEs %>%
  split(.$trait) %>%
  map(~mutate_at(., vars(line_name, environment), as.factor))

# Subset the K matrix
G <- K[tp_geno, tp_geno]


# Iterate and fit models
tp_twoway_fit <- map2(.x = tp_BLUEs_split, .y = env_cor_tp[traits], ~{
  phenos <- .x
  E <- .y
  
  # Create the kronecker product
  Kge <- kronecker(E, G, make.dimnames = TRUE)
  
  
  # Model frame
  mf <- model.frame(value ~ line_name + environment, data = phenos)
  y <- model.response(mf)
  # # Try just the mean as the fixed effect - NO
  # X <- model.matrix(~ 1, data = mf)
  
  # Fixed effects of environments
  X <- model.matrix(~ 1 + environment + line_name, data = droplevels(mf), contrasts.arg = list(line_name = "contr.sum", environment = "contr.sum"))
  Z <- model.matrix(~ -1 + line_name:environment, data = mf)
  
  
  # Fit
  fit_complete <- mixed.solve(y = y, Z = Z, K = Kge, X = X)
  
  
  
  ## Remove data from 2017 and re-fit
  # Model frame
  mf <- model.frame(value ~ line_name + environment, data = droplevels(filter(phenos, year < 2017)))
  y <- model.response(mf)

  # Fixed effects of environments
  X <- model.matrix(~ 1 + environment + line_name, data = droplevels(mf), contrasts.arg = list(line_name = "contr.sum", environment = "contr.sum"))
  Z <- model.matrix(~ -1 + line_name:environment, data = mf)
  
  # Create the kronecker product
  Kge <- kronecker(E[levels(mf$environment), levels(mf$environment)], G, make.dimnames = TRUE)
  
  fit_realistic <- mixed.solve(y = y, Z = Z, K = Kge, X = X)
  
  ## Create a df
  fit_df <- data_frame(set = c("complete", "realistic"), model = list(fit_complete, fit_realistic))
  
  
  mf_list <- list(phenos, droplevels(filter(phenos, year < 2017)))
  
  ## Add fixed effects and random effects together
  fit_df1 <- fit_df %>%
    mutate(predictions = map2(.x = model, .y = mf_list, ~{
      fit <- .x
      mf <- .y
  
      mu <- fit$beta[[1]]
      env_eff <- fit$beta[2:n_distinct(mf$environment)]
      gen_eff <- fit$beta[(n_distinct(mf$environment)+1):length(fit$beta)]
      
      env_eff_df <- data_frame(environment = levels(mf$environment), effect = c(env_eff, -sum(env_eff)))
      gen_eff_df <- data_frame(line_name = levels(mf$line_name), effect = c(gen_eff, -sum(gen_eff)))
      
      ## Random effects
      ge_eff_df <- data_frame(term = names(fit$u), effect = fit$u) %>% 
        separate(term, c("environment", "line_name"), sep = ":")
      
      fitted_df <- left_join(ge_eff_df, gen_eff_df, by = "line_name") %>% 
        left_join(., env_eff_df, by = "environment") %>% 
        mutate(value = effect.x + effect.y + effect + mu) %>%
        select(line_name, environment, value)
      
      # Return the predictions and the GxE BLUPs
      list(predictions = fitted_df, varE = fit$Ve)
      
    }))
  
  fit_df1 %>% 
    mutate( varE = map_dbl(predictions, "varE"), predictions = map(predictions, "predictions")) %>% 
    select(-model)
  
})



## Create two datasets - one using all environments and one without 2017
training_sets_twoway <- tp_twoway_fit %>% 
  map2_df(.x = ., .y = names(.), ~mutate(.x, trait = .y)) %>% 
  select(trait, set, data = predictions, varE)


## Perform the AMMI analysis
# Create a two-way matrix
ammi_out <- training_sets_twoway %>%
  group_by(set, trait) %>%
  do(ammi = {
    df <- .$data[[1]]
    varE <- .$varE
    
    ## Fit a model
    fit <- aov(value ~ line_name + environment, data = df)
    # Get the effects
    effects <- model.tables(fit)
    effects_df <- effects$tables %>% 
      map(~as.matrix(.) %>% fix_data_frame(newnames = "effect"))
    
    # Add the residuals (GxE effects) to the df
    df1 <- add_residuals(df, fit, var = "ge")
    
    # Create a two-way table
    ge <- xtabs(ge ~ line_name + environment, df1)
    
    ## SVD
    ge_svd <- svd(ge)
    # Get the d values for each factor
    dm <- svd(ge)$d
    
    ## Calculate sums of squares for each factor
    SSd <- dm^2
    # Calculate the proportion of variance explained by each factor
    varprop <- SSd / sum(SSd)
    # Cumulative proportion
    cumprop <- cumsum(varprop)
    
    # Calculate the degrees of freedom for each factor
    J <- nrow(ge)
    K <- ncol(ge)
    df_d <- (J + K - 1 - (2 * seq_along(SSd)))
    
    # Calculate the mean squares
    MSd <- SSd / df_d
    
    ## Use the residual variance from the previous model as the expected mean squares of the residuals
    ## Calculate an F statistic
    stat <- MSd / varE
    
    # Calculate p-values
    p_value_Ftest <- pf(q = stat, df1 = df_d, df2 = (J * K), lower.tail = FALSE)
    

    
    # ## Use resampling to determine the null distribution of MSf / MSe
    # # First fit a main effects model, then simulate from that model to calculate ge terms
    # fit <- aov(value ~ line_name + environment, df)
    # mf <- model.frame(fit)
    # 
    # # Simulate from the model
    # sims <- simulate(fit, nsim = 1000)
    # # Update the model
    # fit_sims <- apply(X = sims, MARGIN = 2, FUN = function(data) {
    #   mf$value <- data
    #   fit1 <- lm(value ~ line_name + environment, mf)
    #   
    #   # Get the residuals / GE
    #   ge1 <- matrix(data = residuals(fit1), nrow = J, ncol = K)
    #   svd(ge1)$d
    # })
    # 
    # ## Use the matrix to determine the significance of each factor
    # SS_sims <- fit_sims^2
    # p_value_sim <- map_dbl(seq_along(SSd), ~mean(SS_sims[.,] >= SSd[.]))

    p_value_sim <- rep(NA, length(SSd))
    
    # Create a table of sums of squares
    out <- data.frame(term = c("ge", str_c("PC", seq_along(SSd)), "Residuals"),
               sum_squares = c(sum(ge^2), SSd, varE * (J * K)),
               prop_var = c(NA, varprop, NA),
               cumprop = c(NA, cumprop, NA),
               F_value = c(NA, stat, NA),
               p_value_F = c(NA, p_value_Ftest, NA),
               p_value_sim = c(NA, p_value_sim, NA))
    
    ## Select the significant environmental scores based on the F-test
    which_scores <- which(p_value_Ftest < alpha)
    
    escores <- ge_svd$v[,which_scores, drop = FALSE]
    dimnames(escores) <- list(colnames(ge), paste0("PC", seq(ncol(escores))))
    
    escores_df <- effects_df$environment %>%
      rename(environment = term) %>%
      left_join(., fix_data_frame(escores, newcol = "environment"), by = "environment") %>%
      gather(PC, eigen, -environment, -effect) %>%
      left_join(., data_frame(PC = paste0("PC", seq_along(dm)), d = dm), by = "PC") %>%
      mutate(score = eigen * sqrt(d)) %>%
      select(-d)
    
    gscores <- ge_svd$u[,which_scores, drop = FALSE]
    dimnames(gscores) <- list(rownames(ge), paste0("PC", seq(ncol(gscores))))
    
    gscores_df <- effects_df$line_name %>%
      rename(line_name = term) %>%
      left_join(., fix_data_frame(gscores, newcol = "line_name"), by = "line_name") %>%
      gather(PC, eigen, -line_name, -effect) %>%
      left_join(., data_frame(PC = paste0("PC", seq_along(dm)), d = dm), by = "PC") %>%
      mutate(score = eigen * sqrt(d)) %>%
      select(-d)
    
    # Return
    list(results = out, escores = escores_df, gscores = gscores_df) %>% map(as_data_frame)
    
  }) %>% ungroup()
 

## What total proportion of variation do the AMMI axes explain?
ammi_out %>% 
  mutate(results = map(ammi, "results")) %>% 
  unnest(results) %>%
  filter(p_value_F < alpha) %>% 
  group_by(set, trait) %>% 
  do(tail(., 1)) %>% 
  select(term, cumprop)

# set       trait       term  cumprop
# 1 complete  GrainYield  PC2     0.622
# 2 complete  HeadingDate PC6     0.776
# 3 complete  PlantHeight PC1     0.540
# 4 realistic GrainYield  PC1     0.551
# 5 realistic HeadingDate PC4     0.707
# 6 realistic PlantHeight PC1     0.601

## What proportion of variation do the first IPCA scores explain?
ammi_out %>% 
  mutate(results = map(ammi, "results")) %>% 
  unnest(results) %>%
  filter(p_value_F < alpha) %>% 
  group_by(set, trait) %>% 
  slice(1) %>% 
  select(term, cumprop)

# set       trait       term  cumprop
# 1 complete  GrainYield  PC1     0.470
# 2 complete  HeadingDate PC1     0.213
# 3 complete  PlantHeight PC1     0.540
# 4 realistic GrainYield  PC1     0.551
# 5 realistic HeadingDate PC1     0.277
# 6 realistic PlantHeight PC1     0.601



year_levels <- c(unique(tp_BLUEs$year), "Genotype")

## Plot the AMMI axes and the proportion of variance explained
g_ammi <- ammi_out %>%
  # filter(set == "complete") %>%
  group_by(set, trait) %>%
  do(plot = {
    df <- .
  
    out <- df$ammi[[1]]
    trait <- df$trait
    
    ## Combine the g and e scores into one DF
    scores_df <- map_df(out[-1], ~mutate(., group = names(.)[1]) %>% `names<-`(., c("term", names(.)[-1])))
    
    # Create a renaming vector
    propvar_rename <- setNames(paste0(str_replace_all(out$results$term, "PC", "IPCA"), " (", round(out$results$prop_var, 2) * 100, "%)"), 
                               as.character(out$results$term))
    # Renaming function
    label_rename <- function(x) str_replace_all(x, propvar_rename)
    
    scores_toplot <- scores_df %>%
      select(-eigen) %>%
      spread(PC, score) %>%
      # Add year
      mutate(year = ifelse(group == "environment", paste0("20", str_extract(term, "[0-9]{2}")), "Genotype"),
             trait = trait)
    
    ## Year color scheme
    year_color <- setNames( c(neyhart_palette("umn1", 5)[3:5], "grey"), year_levels)
    
    ## Figure out the range for x and y axis for the trait
    effect_breaks <- subset(ammi_out, trait == df$trait, ammi, drop = T) %>% 
      map(~.[-1]) %>% map(., ~map(., "effect") %>% unlist()) %>% 
      unlist() %>% range()
    
    scores_breaks <- subset(ammi_out, trait == df$trait, ammi, drop = T) %>% 
      map(~.[-1]) %>% map(., ~map(., "score") %>% unlist()) %>% 
      unlist() %>% range()
    
    
    # Plot
    scores_toplot %>% 
      ggplot(aes(x = effect, y = PC1, color = year)) + 
      geom_point() +
      geom_text_repel(data = filter(scores_toplot, group == "environment"), aes(label = term)) +
      scale_color_manual(values = year_color, name = NULL) + 
      scale_y_continuous(breaks = pretty, limits = range(scores_breaks)) +
      scale_x_continuous(breaks = pretty, limits = range(effect_breaks)) +
      ylab(propvar_rename[2]) +
      xlab(expression("Effect (deviation from"~mu*")")) +
      facet_grid(~ trait) + 
      theme_presentation2() +
      theme(legend.position = "bottom")
      # theme(legend.position = "right", legend.direction = "vertical")
    
  })


# Plot in a grid
plot_list_complete <- filter(g_ammi, set == "complete")$plot
# ammi_grid <- plot_grid(plotlist = plot_list %>% map(~. + theme(legend.position = "none")), ncol = 1)
ammi_grid_complete <- plot_grid(ncol = 1, rel_heights = c(0.9, 0.9, 1),
                                plot_list_complete[[1]] + theme(legend.position = "none", axis.title.x = element_blank()),
                                plot_list_complete[[2]] + theme(legend.position = "none", axis.title.x = element_blank()),
                                plot_list_complete[[3]] + theme(legend.position = "none"))
  
# Add legend
ammi_grid1 <- plot_grid(ammi_grid_complete, get_legend(plot_list_complete[[1]]), ncol = 1, rel_heights = c(1, 0.05))

# Save
ggsave(filename = "ammi_biplots_complete.jpg", plot = ammi_grid1, path = fig_dir, height = 11, width = 5, dpi = 1000)
  

## Realistic
# Plot in a grid
plot_list_realistic <- filter(g_ammi, set == "realistic")$plot
ammi_grid_realistic <- plot_grid(ncol = 1, rel_heights = c(0.9, 0.9, 1),
                                 plot_list_realistic[[1]] + theme(legend.position = "none", axis.title.x = element_blank()),
                                 plot_list_realistic[[2]] + theme(legend.position = "none", axis.title.x = element_blank()),
                                 plot_list_realistic[[3]] + theme(legend.position = "none"))
# Add legend
ammi_grid1 <- plot_grid(ammi_grid_realistic, get_legend(plot_list_realistic[[1]]), ncol = 1, rel_heights = c(1, 0.05))

# Save
ggsave(filename = "ammi_biplots_realistic.jpg", plot = ammi_grid1, path = fig_dir, height = 11, width = 5, dpi = 1000)

    



g_ammi <- ammi_out %>%
  group_by(trait) %>%
  do(plot = {
    df <- .
    
    trait <- df$trait
    
    ## Combine the g and e scores into one DF
    scores_df <- df %>% 
      mutate(out = map(ammi, ~map_df(.[-1], ~mutate(., group = names(.)[1]) %>% `names<-`(., c("term", names(.)[-1]))))) %>%
      unnest(out)
    
    propvar_rename_df <- df %>% 
      mutate(propvar = map(ammi, "results") %>% map(~mutate(., propvar = paste0(term, " (", formatC(prop_var * 100, 2, 2), "%)")) %>% select(term, propvar))) %>%
      unnest(propvar)
    
    scores_toplot <- scores_df %>%
      filter(PC == "PC1") %>% 
      left_join(., propvar_rename_df, by = c("PC" = "term", "set", "trait")) %>%
      # spread(PC, score) %>%
      # Add year
      mutate(year = ifelse(group == "environment", paste0("20", str_extract(term, "[0-9]{2}")), "Genotype"),
             set = str_to_title(set))
             
    
    ## Year color scheme
    year_color <- setNames( c(neyhart_palette("umn2")[c(3,4,7)], "grey80"), year_levels)
    
    # Plot
    g <- scores_toplot %>% 
      filter(group == "line_name") %>%
      ggplot(aes(x = effect, y = score, color = year)) + 
      geom_point(size = 0.5) +
      geom_point(data = filter(scores_toplot, group == "environment"), size = 0.5) +
      geom_text_repel(data = filter(scores_toplot, group == "environment"), aes(label = term), size = 2) +
      scale_color_manual(values = year_color, name = NULL) + 
      scale_y_continuous(breaks = pretty) +
      scale_x_continuous(breaks = pretty) +
      ylab("IPCA1 score") +
      xlab(expression("Effect (deviation from"~mu*")")) +
      theme_presentation2(base_size = 10) +
      theme(legend.position = "bottom")
    # theme(legend.position = "right", legend.direction = "vertical")
    
    if (trait == levels(as.factor(ammi_out$trait))[1]) {
      g + facet_grid(trait ~ set + propvar, switch = "y") 
    } else {
      g + facet_grid(trait ~ propvar, switch = "y") 
    }
    
  })


# Plot in a grid
plot_list <- g_ammi$plot

## Plot in a single column
ammi_grid <- plot_grid(plotlist = plot_list %>% map(~ . + theme(legend.position = "none", axis.title.x = element_blank())), 
                       ncol = 1, rel_heights = c(1, 0.9, 0.9)) %>%
  add_sub(label = expression("Effect (deviation from"~mu*")"), size = 10) %>% 
  ggdraw()

ammi_grid1 <- plot_grid(ammi_grid, get_legend(plot_list[[1]]), ncol = 1, rel_heights = c(1, 0.05))

## Save
ggsave(filename = "ammi_biplots_combined.jpg", plot = ammi_grid1, path = fig_dir, height = 7, width = 5, dpi = 1000)




## Save all of this
save_file <- file.path(result_dir, "genotype_environment_phenotypic_analysis.RData")
save("training_sets_twoway", "ammi_out", "stage_two_fits_GYL", "stage_two_fits_GE", file = save_file)

    














