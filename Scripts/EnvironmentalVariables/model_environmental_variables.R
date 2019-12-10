## Exhaustive search for environmental variables that are correlated with the environmental mean
## 
## This script will first fit a model to get the fixed effect of environments. Then, using information
## on planting times, it will search for environmental variables that are correlated. The output of this script
## will include the covariable matrix for environments.
## 
## Author: Jeff Neyhart
## Last modified: October 04, 2018
## 


# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load other libraries
library(lubridate)
library(lmerTest)
library(modelr)
library(broom)
library(pbr)
library(cowplot)
library(optimx)

# Load the environmental variables
load(file.path(data_dir, "environmental_data_compiled.RData"))
# Load the AMMI results
load(file.path(result_dir, "genotype_environment_phenotypic_analysis.RData"))

# Significance level
alpha <- 0.05


## Character replacement for covariates
ec_type_replace <- c("PPT" = "Prcp", "TAVG" = "AvgTemp", "TMAX" = "MaxTemp", "TMIN" = "MinTemp", "annual_TRANGE_max" = "AnnoRangeTemp",
                     "annual_TRANGE" = "AvgRangeTemp", "TRANGE" = "RangeTemp", "TSEASON" = "SeasTemp", "isothermality" = "IsoTemp", "subsoil" = "SS", "topsoil" = "TS", 
                     "claytotal_r_" = "Clay", "om_r_" = "OrgMat", "sandtotal_r_" = "Sand", "silttotal_r_" = "Silt",
                     "ph1to1h2o_r_" = "SoilpH", "latitude" = "Latitude", "longitude" = "Longitude", "elevation" = "Elevation")
## Vector to replace dates/timeframes. Use abbreviations for growth stages
ec_date_replace <- c("annual_" = "Anno", "max_" = "Max", "min_" = "Min", "interval_" = "Int",
                     "vegetative" = "VG", "flowering" = "FL", "grainfill" = "GF")














### Fit a joint regression model
##
## Models are fitted using all data or using data excluding each year
## 


# # First create a list of data split by whether all lines are used or just the TP
# data_to_model <- list(tp = tp, all = c(tp, vp)) %>% 
#   map(~filter(S2_MET_BLUEs, line_name %in% .)) %>%
#   list(., names(.)) %>%
#   pmap_df(~mutate(.x, population = .y))

data_to_model <- training_sets_twoway %>%
  filter(trait %in% traits) %>%
  select(-varE) %>% 
  unnest() %>%
  distinct(set, trait, environment, line_name) %>%
  left_join(., S2_MET_BLUEs)


# Fit the model per trait, extract the environmental means
full_model_fits <- data_to_model %>%
  mutate_at(., vars(line_name, environment), as.factor) %>%
  group_by(trait, set) %>%
  do(model = {
    
    df <- droplevels(.)
  
    # Control and weights
    control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", calc.derivs = FALSE)
    wts <- df$std_error^2
    
    # Formula
    form <- value ~ (1|line_name) + environment
    # Fit the model
    lmer(formula = form, data = df, control = control, weights = wts, contrasts = list(environment = "contr.sum"))
    
  })




    
# Extract the environmental means
env_means_all <- full_model_fits %>%
  group_by(trait, set) %>%
  do({
    
    df <- .
    mod <- df$model[[1]]
    
    mf <- model.frame(mod)
    
    # Get the fixed coefficients
    env_eff <- fixef(mod)
    # Tidy
    data_frame(environment = levels(mf$environment), 
               h = c(env_eff[-1], -sum(env_eff[-1])))  
    
    }) %>% ungroup()


## Correlate the environment means from the two sets
env_means_all %>% 
  group_by(trait, environment) %>% 
  filter(n() == n_distinct(env_means_all$set) - 1) %>% 
  spread(set, h) %>% 
  split(.$trait) %>%
  map(~cor(.[,-1:-2], use = "pairwise.complete.obs"))
  


## High correlations of environmental means between sets

# This is not surprising

    





## Calculate regression coefficients for each genotype, then test the signficance of differences among regression coefficients
data_to_model %>% 
  left_join(., env_means_all) %>%
  group_by(trait, set) %>%
  do(fwr(formula = value ~ h, data = .))

# trait           set       regression         Fstat   pvalue
# 1 GrainYield  complete      <tibble [175 x 4]> 2.51  1.49e-22
# 2 GrainYield  realistic2015 <tibble [175 x 4]> 1.71  6.13e- 8
# 3 GrainYield  realistic2016 <tibble [175 x 4]> 2.64  1.08e-23
# 4 GrainYield  realistic2017 <tibble [175 x 4]> 2.35  2.61e-18
# 5 HeadingDate complete      <tibble [175 x 4]> 1.44  1.91e- 4
# 6 HeadingDate realistic2015 <tibble [175 x 4]> 1.39  7.42e- 4
# 7 HeadingDate realistic2016 <tibble [175 x 4]> 0.793 9.76e- 1
# 8 HeadingDate realistic2017 <tibble [175 x 4]> 1.85  6.83e-10
# 9 PlantHeight complete      <tibble [175 x 4]> 1.40  5.55e- 4
# 10 PlantHeight realistic2015 <tibble [175 x 4]> 1.20  3.88e- 2
# 11 PlantHeight realistic2016 <tibble [175 x 4]> 1.25  1.59e- 2
# 12 PlantHeight realistic2017 <tibble [175 x 4]> 1.28  9.26e- 3

## As expected, all are significant
## Since all traits show genotype-specific reactions to the environment mean, genotypes should also
# demonstrate different responses to variables that are correlated with the mean












## Correlate environmental variables with the mean

## First combine one-year with multi-year ECs
environ_covariate_data1 <- environ_covariate_data %>%
  unnest() %>%
  ## Adjust names
  mutate(variable_newname = str_replace_all(variable, ec_type_replace) %>% str_replace_all(ec_date_replace)) %>%
  ## Remove pre-emergence and maturity ECs
  filter(!str_detect(variable, "maturity|pre-emergence"))


## A separate DF of variable names and new names
ec_names <- distinct(environ_covariate_data1, variable, variable_newname)

## Write a csv
write_csv(x = ec_names, path = file.path(fig_dir, "env_covariate_names.csv"))

## Remove the interval covariates
ec_names <- environ_covariate_data1 %>%
  filter(timeframe != "interval") %>%
  distinct(variable, variable_newname)

write_csv(x = ec_names, path = file.path(fig_dir, "env_covariate_names_growthstage.csv"))




## A function to scale a vector to dot product sum = 1
scale_length <- function(x) {
  # First center at 0
  x1 <- scale(x, center = TRUE, scale = FALSE)
  # Next scale to obtain 1 as the squared vector length
  scale(x1, center = FALSE, scale = sqrt(sum(x1^2)))
}
  



## Only use growth stage information moving forward
# Center each covariate and scale to achive the squared length of the vector = 1
environ_covariate_data2 <- environ_covariate_data1 %>%
  filter(timeframe == "growth_stage") %>%
  split(list(.$ec_group, .$variable, .$timeframe)) %>%
  map_df(~mutate(., scaled_value = scale_length(x = value),
                 center = attr(scaled_value, "scaled:center"),
                 scale = attr(scaled_value, "scaled:scale"),
                 scaled_value = as.numeric(scaled_value)))


## Look at histogram of all covariates
# Determine the number of rows/cols
nplots <- environ_covariate_data2 %>% group_by(ec_group, timeframe, variable_newname) %>% n_groups()
nrow <- 5
ncol <- 5
npage <- ceiling(nplots / (nrow * ncol))

# Open a PDF
pdf(file = file.path(fig_dir, "env_covariate_distribution.pdf"), onefile = TRUE)

# Iterate over pages
for (p in seq(npage)) {
  
  gg_ec_dist <- environ_covariate_data2 %>%
    ggplot(aes(x = value)) +
    geom_density() +
    # facet_grid(~ ec_group + timeframe + variable_newname, scales = "free", 
    #                              ncol = ncol, nrow = nrow) +
    ggforce::facet_wrap_paginate(~ ec_group + timeframe + variable_newname, scales = "free",
                                 ncol = ncol, nrow = nrow, page = p) +
    theme_acs()
  
  print(gg_ec_dist)
}

dev.off()



## Re-organization and removal of some covariates
environ_covariate_data3 <- environ_covariate_data2 %>%
  # filter(ec_group == "multiyear") %>% # Only look at the 10year average
  left_join(env_means_all, .) %>%
  mutate(EC_type = "summary") %>%
  filter(! variable %in% c("elevation", "latitude", "longitude")) %>%
  # # Remove some variables that are not well distributed
  filter(!variable %in% c("claytotal_r_topsoil", "sandtotal_r_subsoil", "sandtotal_r_topsoil", "silttotal_r_subsoil",
                          "om_r_subsoil")) %>%
  ## Filter out growth stage covariates that do not make sense for a trait
  filter(!(trait == "HeadingDate" & str_detect(variable, "flowering|grainfill")),
         !(trait == "PlantHeight" & str_detect(variable, "grainfill")))
  
  
  



## Correlate covariates with the environmental mean
## Here we used the actual value of the covariates (later it will be scaled)
env_mean_cor <- environ_covariate_data3 %>%
  group_by(set, trait, ec_group, variable, timeframe) %>% 
  do(test = cor.test(.$h, .$value)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         df = map_dbl(test, "parameter")) %>%
  select(-test)


## Plot the top n correlations per group
n = 5
env_mean_topn_cor <- env_mean_cor %>%
  group_by(trait, set, ec_group, timeframe) %>%
  top_n(x = ., n = n, wt = abs(cor)) %>%
  filter(set == "complete") %>%
  arrange(desc(abs(cor))) %>%
  ungroup() %>%
  mutate(test = "env_mean") %>%
  left_join(., environ_covariate_data3)

env_mean_cor_plot_list <- env_mean_topn_cor %>%
  split(.$trait) %>%
  map(~ggplot(data = ., aes(x = value, y = h)) +
        geom_point() +
        geom_label(aes(label = paste0("r = ", round(cor, 3)), x = Inf, y = -Inf, hjust = 1, vjust = -0.5), size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        facet_grid(ec_group ~ timeframe + variable, scales = "free", switch = "y") +
        theme_presentation2(base_size = 12))


# Remove some outliers and re-run
environ_covariate_data4 <- environ_covariate_data3 %>%
  filter(!(trait == "GrainYield" & h > 3000),
         !(trait == "PlantHeight" & h > 25))

env_mean_cor1 <- environ_covariate_data4 %>%
  group_by(set, trait, ec_group, variable, timeframe) %>% 
  do(test = cor.test(.$h, .$value)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         df = map_dbl(test, "parameter")) %>%
  select(-test)



env_mean_topn_cor1 <- env_mean_cor1 %>%
  group_by(trait, set, ec_group, timeframe) %>%
  top_n(x = ., n = n, wt = abs(cor)) %>%
  filter(set %in% c("complete", "realistic2017")) %>%
  arrange(desc(abs(cor))) %>%
  ungroup() %>%
  mutate(test = "env_mean") %>%
  left_join(., environ_covariate_data4)

env_mean_cor_plot_list1 <- env_mean_topn_cor1 %>%
  split(.$trait) %>%
  map(~ggplot(data = ., aes(x = value, y = h)) +
        geom_point() +
        geom_label(aes(label = paste0("r = ", round(cor, 3)), x = Inf, y = -Inf, hjust = 1, vjust = -0.5), size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        facet_grid(set + ec_group ~ timeframe + variable, scales = "free", switch = "y") +
        theme_presentation2(base_size = 12))

## Top 5 coviarates for each trait
env_mean_top5_cor <- env_mean_cor1 %>%
  group_by(set, trait, ec_group) %>% 
  top_n(n = 5, wt = abs(cor)) %>%
  arrange(set, ec_group, trait, desc(abs(cor))) %>%
  mutate(top = seq(5)) %>%
  ungroup() %>%
  mutate(test = "env_mean") %>%
  left_join(., ec_names)

## 






## Try multiple regression, where ECs are added in order of correlation
env_mean_mr <- environ_covariate_data4 %>%
  distinct(set, trait, ec_group, timeframe) %>%
  mutate(out = list(NULL))
  
for (i in seq(nrow(env_mean_mr))) {
  
  st <- env_mean_mr$set[i] # Set/scenario
  tr <- env_mean_mr$trait[i] # Trait
  grp <- env_mean_mr$ec_group[i] # 10-yr or 1-yr
  tf <- env_mean_mr$timeframe[i] # growth stage or interval

  df <- environ_covariate_data4 %>% 
    filter(trait == tr, set == st, ec_group == grp, timeframe == tf)

  # Extract EC names in arbitrary order
  ecs <- distinct(df, variable)$variable
  
  # Formula
  base_formula <- h ~ 1
  full_formula <- formula(paste("h ~ ", paste(ecs, collapse = " + ")))
  
  
  # Create a model frame
  df1 <- df %>% 
    select(trait, environment, h, variable, value) %>% 
    spread(variable, value) %>%
    # Center the ECs on 0
    mutate_at(vars(-trait:-h), ~as.numeric(scale(., scale = FALSE)))

  # Base model
  fit <- lm(formula = base_formula, data = df1)
  # Stepwise selected model
  fit_step <- step(object = fit, scope = full_formula, direction = "both", trace = 0)
  # Tidy
  fit_step_tidy <- tidy(anova(fit_step)) %>%
    filter(term != "Residuals") %>%
    select(variable = term, sumsq, p_value = p.value)
  
  ## Get the significant ECs
  sig_ecs <- attr(terms(fit_step), "term.labels")
  
  ## If the number of sig_ecs, is 0, choose the most correlated ec
  if (length(sig_ecs) == 0) {
    sig_ecs <- summarize_at(df1, vars(-trait:-h), list(~cor(., h))) %>% 
      gather(variable, correlation) %>% 
      top_n(x = ., n = 1, wt = abs(correlation)) %>%
      pull(variable)
    
    fit_step <- lm(formula(paste("h ~ ", paste(sig_ecs, collapse = " + "))), data = df1)
    
    fit_step_tidy <- tidy(anova(fit_step)) %>%
      filter(term != "Residuals") %>%
      select(variable = term, sumsq, p_value = p.value)
      
  }
  
  
  ### Fit a model with genotype, environment, and genotype-specific slopes to each ECs
  ## Get the relevant phenotypic data
  pheno_df <- data_to_model %>%
    filter(trait == tr, set == st) %>%
    # Add ECs
    left_join(., select(df1, trait, environment, sig_ecs), by = c("trait", "environment"))
  
  # Base formula
  base_formula2 <- value ~ line_name + environment
  
  
  ## If no sig_ecs, skip the next step
  if (length(sig_ecs) > 0) {
  
    # Formula
    full_formula2 <- add_predictors(base_formula2, as.formula(paste0("~", paste0("line_name:", sig_ecs, collapse = " + "))))
    
    # Fit just the base model
    base_fit2 <- lm(base_formula2, pheno_df)
    # Fit the full model with interactions
    full_fit2 <- lm(full_formula2, pheno_df)
    
  } else {
    base_fit2 <- lm(base_formula2, pheno_df)
    full_fit2 <- NULL
  }
  
  
  # ## Mixed model version
  # # Formula
  # base_formula2 <- value ~ (1|line_name) + environment
  # full_formula2 <- add_predictors(base_formula2, as.formula(paste0("~", paste0("(0 + ", sig_ecs, "|line_name)", collapse = " + "))))
  # 
  # # Fit just the base model
  # base_fit2 <- lmer(base_formula2, pheno_df)
  # # Fit the full model with interactions
  # full_fit2 <- lmer(full_formula2, pheno_df)
  

  ## Return a tidy df of the model
  env_mean_mr$out[[i]] <- data_frame(
    fit = list(fit_step),
    tidy = list(fit_step_tidy),
    fit2_base = list(base_fit2),
    fit2_full = list(full_fit2))
  
}
    

## Unnest the DF
env_mean_mr1 <- unnest(env_mean_mr)

## Mean square of covariates compared with mean square of interactions
env_mean_mr1_hyptest <- env_mean_mr1 %>%
  filter(!map_lgl(fit2_full, is.null)) %>%
  mutate(interaction_anova = map(fit2_base, ~subset(tidy(anova(.)), term == "Residuals", c(meansq, df))),
         factorial_reg_anova = map(fit2_full, ~subset(tidy(anova(.)), str_detect(term, "line_name:"), c(term, meansq, df)))) %>%
  unnest(interaction_anova) %>%
  unnest(factorial_reg_anova) %>%
  # Calculate F stat and perform test
  mutate(Fstat = meansq1 / meansq,
         pvalue = pf(q = Fstat, df1 = df1, df2 = df, lower.tail = FALSE))


## For each analysis, how many associated ECs were also significant?
env_mean_mr1_hyptest %>% 
  group_by(set, trait, ec_group, timeframe) %>% 
  summarize(nECs = n_distinct(term), nSigECs = sum(pvalue <= 0.05)) %>%
  as.data.frame()



## Unnest the DF and extract R2, the significant ECs, and the FW p-value
env_mean_mr2 <- env_mean_mr1 %>%
  mutate(fit_R2 = map_dbl(fit, ~summary(.)$adj.r.squared),
         ecs = map(fit, ~attr(terms(.), "term.labels")))



## How many ECs per trait and set
env_mean_mr2 %>% 
  filter(trait %in% traits,
         timeframe == "growth_stage",
         ec_group == "multiyear") %>%
  mutate(nECs = map_dbl(ecs, length)) %>% 
  xtabs(nECs ~ set + trait, data = .)

## How many unique ECs
env_mean_mr2 %>% 
  filter(trait %in% traits,
         timeframe == "growth_stage",
         ec_group == "multiyear") %>% 
  pull(ecs) %>% 
  unlist() %>%
  unique()

##  27


## ECs common to traits within sets
env_mean_mr2 %>% 
  filter(trait %in% traits,
         timeframe == "growth_stage",
         ec_group == "multiyear") %>%
  select(set, trait, ecs) %>%
  crossing(., .) %>%
  filter(set != set1, trait == trait1) %>%
  mutate(common_ec = map2(ecs, ecs1, intersect)) %>%
  mutate_at(vars(ecs, ecs1, common_ec), ~map_dbl(., length)) %>%
  as.data.frame()


##

# set       trait ecs          set1      trait1 ecs1 common_ec
# 1       complete  GrainYield   4 realistic2015  GrainYield    6         2
# 2       complete  GrainYield   4 realistic2016  GrainYield   11         2
# 3       complete  GrainYield   4 realistic2017  GrainYield    3         1
# 4  realistic2015  GrainYield   6      complete  GrainYield    4         2
# 5  realistic2015  GrainYield   6 realistic2016  GrainYield   11         4
# 6  realistic2015  GrainYield   6 realistic2017  GrainYield    3         1
# 7  realistic2016  GrainYield  11      complete  GrainYield    4         2
# 8  realistic2016  GrainYield  11 realistic2015  GrainYield    6         4
# 9  realistic2016  GrainYield  11 realistic2017  GrainYield    3         2
# 10 realistic2017  GrainYield   3      complete  GrainYield    4         1
# 11 realistic2017  GrainYield   3 realistic2015  GrainYield    6         1
# 12 realistic2017  GrainYield   3 realistic2016  GrainYield   11         2
# 13      complete HeadingDate   2 realistic2015 HeadingDate    2         2
# 14      complete HeadingDate   2 realistic2016 HeadingDate    1         0
# 15      complete HeadingDate   2 realistic2017 HeadingDate    4         2
# 16 realistic2015 HeadingDate   2      complete HeadingDate    2         2
# 17 realistic2015 HeadingDate   2 realistic2016 HeadingDate    1         0
# 18 realistic2015 HeadingDate   2 realistic2017 HeadingDate    4         2
# 19 realistic2016 HeadingDate   1      complete HeadingDate    2         0
# 20 realistic2016 HeadingDate   1 realistic2015 HeadingDate    2         0
# 21 realistic2016 HeadingDate   1 realistic2017 HeadingDate    4         0
# 22 realistic2017 HeadingDate   4      complete HeadingDate    2         2
# 23 realistic2017 HeadingDate   4 realistic2015 HeadingDate    2         2
# 24 realistic2017 HeadingDate   4 realistic2016 HeadingDate    1         0
# 25      complete PlantHeight   7 realistic2015 PlantHeight    3         2
# 26      complete PlantHeight   7 realistic2016 PlantHeight   14         5
# 27      complete PlantHeight   7 realistic2017 PlantHeight   13         2
# 28 realistic2015 PlantHeight   3      complete PlantHeight    7         2
# 29 realistic2015 PlantHeight   3 realistic2016 PlantHeight   14         3
# 30 realistic2015 PlantHeight   3 realistic2017 PlantHeight   13         2
# 31 realistic2016 PlantHeight  14      complete PlantHeight    7         5
# 32 realistic2016 PlantHeight  14 realistic2015 PlantHeight    3         3
# 33 realistic2016 PlantHeight  14 realistic2017 PlantHeight   13         6
# 34 realistic2017 PlantHeight  13      complete PlantHeight    7         2
# 35 realistic2017 PlantHeight  13 realistic2015 PlantHeight    3         2
# 36 realistic2017 PlantHeight  13 realistic2016 PlantHeight   14         6

   

## # 


## For each set, trait, and ec_group, determine the overlap in ECs that were determined
## to be significant

# (overlap <- env_mean_mr2 %>% 
#   select(set:timeframe, ecs) %>% 
#   spread(timeframe, ecs) %>%
#   mutate(ec_overlap = map2(.x = growth_stage, .y = interval, intersect)))


## Determine EC overlap between LOYO years
(overlap <- env_mean_mr2 %>%
    select(set:timeframe, ecs) %>%
    filter(set != "complete") %>%
    spread(set, ecs) %>%
    mutate(ec_overlap = select(., contains("realistic")) %>% pmap(., ~reduce(., intersect))) %>%
    mutate_at(vars(contains("realistic")), list(~map2(., ec_overlap, setdiff))) )


## Plot the significant ECs
env_mean_mr_plot_df <- env_mean_mr2 %>% 
  mutate(ecs = map(ecs, ~tibble(variable = .))) %>% 
  unnest(ecs) %>% 
  left_join(., environ_covariate_data4) %>%
  group_by(trait, set, ec_group, timeframe) %>%
  nest(variable_newname, environment, h, value)
  
env_mean_mr_plot_list <- vector("list", nrow(env_mean_mr_plot_df))
ncol <- 5
nrow <- 4


# Open a PDF
pdf(file = file.path(fig_dir, "env_mean_significant_covariates.pdf"), onefile = TRUE)

## Iterate over each row and make a plot
for (i in seq_along(env_mean_mr_plot_list)) {
  
  row <- env_mean_mr_plot_df[i,]
  df <- unnest(row)
  
  ## Determine pages
  nplots <- n_distinct(df$variable_newname)
  npage <- ceiling(nplots / (ncol * nrow))
  
  ## Create a plot title
  title <- paste0(row[,1:4], collapse = " : ")
  
  for (p in seq(npage)) {
    
    gg <- df %>%
      ggplot(aes(x = value, y = h)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      scale_x_continuous(name = "EC value", breaks = pretty) +
      scale_y_continuous(name = "Environmental effect", breaks = pretty) +
      ggforce::facet_wrap_paginate(~ variable_newname, ncol = ncol, nrow = nrow, scales = "free_x") +
      labs(title = title, caption = paste0("Page ", p, "/", npage)) +
      theme_acs()
  
    print(gg)
    
  }
}

# Close the PDF
dev.off()



#### Repeat, but include only pages that are relevant to paper
env_mean_mr_plot_df_paper <- env_mean_mr_plot_df %>%
  filter(trait %in% traits,
         # set %in% c("complete", "realistic2017"),
         ec_group == "multiyear",
         timeframe == "growth_stage")

env_mean_mr_plot_list <- vector("list", nrow(env_mean_mr_plot_df_paper))

# Open a PDF
pdf(file = file.path(fig_dir, "env_mean_significant_covariates_paper.pdf"), onefile = TRUE)

## Iterate over each row and make a plot
for (i in seq_along(env_mean_mr_plot_list)) {
  
  row <- env_mean_mr_plot_df_paper[i,]
  df <- unnest(row)
  
  ## Determine pages
  nplots <- n_distinct(df$variable_newname)
  npage <- ceiling(nplots / (ncol * nrow))
  
  ## Create a plot title
  title <- paste0(row[,1:4], collapse = " : ")
  
  for (p in seq(npage)) {
    
    gg <- df %>%
      ggplot(aes(x = value, y = h)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      scale_x_continuous(name = "EC value", breaks = pretty) +
      scale_y_continuous(name = "Environmental effect", breaks = pretty) +
      ggforce::facet_wrap_paginate(~ variable_newname, ncol = ncol, nrow = nrow, scales = "free_x") +
      labs(title = title, caption = paste0("Page ", p, "/", npage)) +
      theme_acs()
    
    print(gg)
    
  }
}

# Close the PDF
dev.off()




## Select those variables that are significant
ec_env_mean_sig <- env_mean_mr_plot_df_paper %>% 
  unnest(data)

## Create a similar data.frame for all variables - this will be used later.
## This will include 
env_combine <- environ_covariate_data3 %>%
  filter(timeframe == "growth_stage", ec_group == "multiyear") %>%
  distinct(trait, set, ec_group, timeframe, variable_newname, environment, value)












##################
## Now correlated ECs with the AMMI environmental IPCA scores
##################


## This one does not remove the outliers - 
ec_score_df_use <- ammi_out %>%
  mutate(escore = map(ammi, "escores")) %>% 
  unnest(escore) %>% 
  filter(PC == "PC1") %>%
  left_join(environ_covariate_data3, ., by = c("trait", "set", "environment"))
    



## Estimate the correlation between the IPCA score and the scaled EC value
env_ipca_cor1 <- ec_score_df_use %>% 
  group_by(set, ec_group, trait, timeframe, variable, PC) %>% 
  do(test = cor.test(.$score, .$scaled_value)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         df = map_dbl(test, "parameter")) %>%
  select(-test)



## Plot the top n correlations per group
n = 5
env_ipca_topn_cor <- env_ipca_cor1 %>%
  group_by(trait, set, ec_group, timeframe) %>%
  top_n(x = ., n = n, wt = abs(cor)) %>%
  filter(set %in% c("complete", "realistic2017")) %>%
  arrange(desc(abs(cor))) %>%
  ungroup() %>%
  mutate(test = "env_ipca") %>%
  left_join(., ec_score_df_use)

env_ipca_cor_plot_list <- env_ipca_topn_cor %>%
  split(.$trait) %>%
  map(~ggplot(data = ., aes(x = score, y = scaled_value)) +
        geom_point() +
        geom_label(aes(label = paste0("r = ", round(cor, 3)), x = Inf, y = -Inf, hjust = 1, vjust = -0.5), size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        facet_grid(set + ec_group ~ timeframe + variable, scales = "free", switch = "y") +
        theme_presentation2(base_size = 12))

## 






## Try multiple regression, where ECs are added in order of correlation
env_ipca_mr <- ec_score_df_use %>%
  group_by(set, trait, ec_group, timeframe) %>%
  nest() %>%
  mutate(out = list(NULL))

for (i in seq(nrow(env_ipca_mr))) {
  
  row <- env_ipca_mr[i,]
  
  st <- row$set[1] # Set/scenario
  tr <- row$trait[1] # Trait
  grp <- row$ec_group[1] # 10-yr or 1-yr EC
  tf <- row$timeframe[1] # Interval or growth stage
  
  df <- unnest(row, data)
  
  # Extract EC names in arbitrary order
  ecs <- distinct(df, variable)$variable
  
  # Formula
  base_formula <- score ~ 1
  full_formula <- formula(paste("score ~ ", paste(ecs, collapse = " + ")))
  
  
  # Create a model frame
  df1 <- df %>% 
    select(trait, environment, score, variable, scaled_value) %>% 
    spread(variable, scaled_value)
  
  # Base model
  fit <- lm(formula = base_formula, data = df1)
  # Stepwise selected model
  fit_step <- step(object = fit, scope = full_formula, direction = "both", trace = 0)
  # Tidy
  fit_step_tidy <- tidy(anova(fit_step)) %>%
    # filter(term != "Residuals") %>%
    select(variable = term, sumsq, p_value = p.value) %>%
    ## Proportion of sums sqs
    mutate(sumsq_prop = sumsq / sum(sumsq))
  
  ## Get the significant ECs
  sig_ecs <- attr(terms(fit_step), "term.labels")
  
  ## If the number of sig_ecs, is 0, choose the most correlated ec
  if (length(sig_ecs) == 0) {
    sig_ecs <- summarize_at(df1, vars(-trait:-score), list(~cor(., score))) %>% 
      gather(variable, correlation) %>% 
      top_n(x = ., n = 1, wt = abs(correlation)) %>%
      pull(variable)
    
    fit_step <- lm(formula(paste("score ~ ", paste(sig_ecs, collapse = " + "))), data = df1)
    
    fit_step_tidy <- tidy(anova(fit_step)) %>%
      filter(term != "Residuals") %>%
      select(variable = term, sumsq, p_value = p.value)
    
  }
  
  
  ## New DF that includes the centered values of the ECs, instead of the scaled values
  df1 <- df %>% 
    select(trait, environment, h, variable, value) %>% 
    spread(variable, value) %>%
    # Center the ECs on 0
    mutate_at(vars(-trait:-h), ~as.numeric(scale(., scale = FALSE)))
  
  ### Fit a model with genotype, environment, and genotype-specific slopes to each ECs
  ## Get the relevant phenotypic data
  pheno_df <- data_to_model %>%
    filter(trait == tr, set == st) %>%
    # Add ECs
    left_join(., select(df1, trait, environment, sig_ecs), by = c("trait", "environment"))
  
  ## If no sig_ecs, skip the next step
  if (length(sig_ecs) > 0) {
    
    # Formula
    base_formula2 <- value ~ line_name + environment
    full_formula2 <- add_predictors(base_formula2, as.formula(paste0("~", paste0("line_name:", sig_ecs, collapse = " + "))))
    
    # Fit just the base model
    base_fit2 <- lm(base_formula2, pheno_df)
    # Fit the full model with interactions
    full_fit2 <- lm(full_formula2, pheno_df)
    
  } else {
    base_formula2 <- value ~ line_name + environment
    base_fit2 <- lm(base_formula2, pheno_df)
    full_fit2 <- NULL
  }
  
  
  # ## Mixed model version
  # # Formula
  # base_formula2 <- value ~ (1|line_name) + environment
  # full_formula2 <- add_predictors(base_formula2, as.formula(paste0("~", paste0("(0 + ", sig_ecs, "|line_name)", collapse = " + "))))
  # 
  # # Fit just the base model
  # base_fit2 <- lmer(base_formula2, pheno_df)
  # # Fit the full model with interactions
  # full_fit2 <- lmer(full_formula2, pheno_df)
  
  
  ## Return a tidy df of the model
  env_ipca_mr$out[[i]] <- data_frame(
    fit = list(fit_step),
    tidy = list(fit_step_tidy),
    fit2_base = list(base_fit2),
    fit2_full = list(full_fit2))
  
}


## Unnest the DF
env_ipca_mr1 <- unnest(env_ipca_mr, out)

## Mean square of covariates compared with mean square of interactions
env_ipca_mr1_hyptest <- env_ipca_mr1 %>%
  filter(!map_lgl(fit2_full, is.null)) %>%
  mutate(interaction_anova = map(fit2_base, ~subset(tidy(anova(.)), term == "Residuals", c(meansq, df))),
         factorial_reg_anova = map(fit2_full, ~subset(tidy(anova(.)), str_detect(term, "line_name:"), c(term, meansq, df)))) %>%
  unnest(interaction_anova) %>%
  unnest(factorial_reg_anova) %>%
  # Calculate F stat and perform test
  mutate(Fstat = meansq1 / meansq,
         pvalue = pf(q = Fstat, df1 = df1, df2 = df, lower.tail = FALSE))


## For each analysis, how many associated ECs were also significant?
env_ipca_mr1_hyptest %>% 
  group_by(set, trait, ec_group, timeframe) %>% 
  summarize(nECs = n_distinct(term), nSigECs = sum(pvalue <= 0.05)) %>%
  as.data.frame()



## Unnest the DF and extract R2, the significant ECs, and the FW p-value
env_ipca_mr2 <- env_ipca_mr1 %>%
  mutate(fit_R2 = map_dbl(fit, ~summary(.)$adj.r.squared),
         ecs = map(fit, ~attr(terms(.), "term.labels")))



## How many ECs per trait and set
env_ipca_mr2 %>% 
  filter(trait %in% traits,
         timeframe == "growth_stage",
         ec_group == "multiyear") %>%
  mutate(nECs = map_dbl(ecs, length)) %>% 
  xtabs(nECs ~ set + trait, data = .)

## How many unique ECs
(unique_ec_ipca <- env_ipca_mr2 %>% 
  filter(trait %in% traits,
         timeframe == "growth_stage",
         ec_group == "multiyear") %>% 
  pull(ecs) %>% 
  unlist() %>%
  unique())

length(unique_ec_ipca)
##  20


## df to use for discovering unique ECs
env_ipca_mr2_use <- env_ipca_mr2 %>% 
  filter(trait %in% traits,
         timeframe == "growth_stage",
         ec_group == "multiyear") %>%
  select(set, trait, ecs)


## First find ECs that are shared between at least two traits, in any scenario
env_ipca_across_trait_shared <- env_ipca_mr2_use %>%
  crossing(., .) %>% 
  filter(trait != trait1) %>%
  mutate(common_ec = map2(ecs, ecs1, intersect))
  
# Print
env_ipca_across_trait_shared_list <- env_ipca_across_trait_shared$common_ec %>% 
  unlist() %>% 
  unique()
length(env_ipca_across_trait_shared_list)

# 7


## Next find ECs that are shared between within a trait
env_ipca_within_trait_shared <- env_ipca_mr2_use %>%
  crossing(., .) %>% 
  filter(trait == trait1, set != set1) %>%
  mutate(common_ec = map2(ecs, ecs1, intersect))

# Print
env_ipca_within_trait_shared_list <- env_ipca_within_trait_shared$common_ec %>% 
  unlist() %>% 
  unique()
length(env_ipca_within_trait_shared_list)

# 5


## Intersect these lists
(env_ipca_within_across_shared_list <- intersect(env_ipca_within_trait_shared_list, env_ipca_across_trait_shared_list))

# 3

## Find the ECs unique to a trait-set combo
env_ipca_mr2_use %>% 
  mutate(ecs = map(ecs, ~setdiff(., union(env_ipca_within_trait_shared_list, env_ipca_across_trait_shared_list)))) %>% 
  pull(ecs) %>% 
  unlist() %>% 
  unique()

# 11



## Plot the significant ECs
env_ipca_mr_plot_df <- env_ipca_mr2 %>% 
  mutate(ecs = map(ecs, ~tibble(variable = .))) %>% 
  unnest(ecs) %>% 
  left_join(., environ_covariate_data4) %>%
  group_by(trait, set, ec_group, timeframe) %>%
  nest(variable_newname, environment, h, value)

env_ipca_mr_plot_list <- vector("list", nrow(env_ipca_mr_plot_df))
ncol <- 5
nrow <- 4


# Open a PDF
pdf(file = file.path(fig_dir, "env_ipca_significant_covariates.pdf"), onefile = TRUE)

## Iterate over each row and make a plot
for (i in seq_along(env_ipca_mr_plot_list)) {
  
  row <- env_ipca_mr_plot_df[i,]
  df <- unnest(row)
  
  ## Determine pages
  nplots <- n_distinct(df$variable_newname)
  npage <- ceiling(nplots / (ncol * nrow))
  
  ## Create a plot title
  title <- paste0(row[,1:4], collapse = " : ")
  
  for (p in seq(npage)) {
    
    gg <- df %>%
      ggplot(aes(x = value, y = h)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      scale_x_continuous(name = "EC value", breaks = pretty) +
      scale_y_continuous(name = "Environmental effect", breaks = pretty) +
      ggforce::facet_wrap_paginate(~ variable_newname, ncol = ncol, nrow = nrow, scales = "free_x") +
      labs(title = title, caption = paste0("Page ", p, "/", npage)) +
      theme_acs()
    
    print(gg)
    
  }
}

# Close the PDF
dev.off()



#### Repeat, but include only pages that are relevant to paper
env_ipca_mr_plot_df_paper <- env_ipca_mr_plot_df %>%
  filter(trait %in% traits,
         # set %in% c("complete", "realistic2017"),
         ec_group == "multiyear",
         timeframe == "growth_stage")

env_ipca_mr_plot_list <- vector("list", nrow(env_ipca_mr_plot_df_paper))

# Open a PDF
pdf(file = file.path(fig_dir, "env_ipca_significant_covariates_paper.pdf"), onefile = TRUE)

## Iterate over each row and make a plot
for (i in seq_along(env_ipca_mr_plot_list)) {
  
  row <- env_ipca_mr_plot_df_paper[i,]
  df <- unnest(row)
  
  ## Determine pages
  nplots <- n_distinct(df$variable_newname)
  npage <- ceiling(nplots / (ncol * nrow))
  
  ## Create a plot title
  title <- paste0(row[,1:4], collapse = " : ")
  
  for (p in seq(npage)) {
    
    gg <- df %>%
      ggplot(aes(x = value, y = h)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      scale_x_continuous(name = "EC value", breaks = pretty) +
      scale_y_continuous(name = "Environmental effect", breaks = pretty) +
      ggforce::facet_wrap_paginate(~ variable_newname, ncol = ncol, nrow = nrow, scales = "free_x") +
      labs(title = title, caption = paste0("Page ", p, "/", npage)) +
      theme_acs()
    
    print(gg)
    
  }
}

# Close the PDF
dev.off()




## Select those variables that are significant
ec_env_ipca_sig <- env_ipca_mr_plot_df_paper %>% 
  unnest(data)









## Create a table with significant ECs
ec_env_sig_combine <- bind_rows(mutate(ec_env_mean_sig, group = "EC_Mean"), mutate(ec_env_ipca_sig, group = "EC_IPCA")) %>%
  filter(ec_group == "multiyear", trait %in% traits, timeframe == "growth_stage") %>%
  distinct(set, trait, variable_newname, group) %>% 
  # filter(set %in% c("complete", "realistic2017")) %>%
  group_by(set, group, trait) %>% 
  mutate(n = seq(n())) %>%
  ungroup() %>%
  mutate(set = str_replace_all(set, set_replace),
         set = str_replace_all(set, "([0-9]{4})", " \\(\\1\\)"),
         trait = str_add_space(trait))

ec_env_sig_combine_towrite <- ec_env_sig_combine %>% 
  spread(trait, variable_newname) %>%
  filter(group != "EC_Mean") %>%
  select(-n)

write_csv(x = ec_env_sig_combine_towrite, path = file.path(fig_dir, "significant_ecs.csv"), na = " ")





## Kable extra
library(kableExtra)
## Officer and flextable
library(officer)
library(flextable)

## Microsoft word table
ec_env_sig_combine_towrite1 <- ec_env_sig_combine_towrite %>%
  select(-group)

typology <- data.frame(
  col_keys = names(ec_env_sig_combine_towrite1),
  trait = c("Scenario", "Grain Yield", "Heading Date", "Plant Height"),
  stringsAsFactors = FALSE )

## Find rows/column to italicize (shared ecs within trait)
italicize <- ec_env_sig_combine_towrite1 %>%
  mutate_all(~. %in% names(which(table(.) > 1))) %>%
  mutate(set = FALSE) %>%
  as.matrix() %>%
  which()
italicize_coord <- tibble(j = ceiling(italicize / nrow(ec_env_sig_combine_towrite1))) %>%
  mutate(i = italicize - ((j - 1) * nrow(ec_env_sig_combine_towrite1)), italic = TRUE)

## Find rows/columns to underline (shared ECs across traits)
underline <- ec_env_sig_combine_towrite1 %>%
  mutate_all(~. %in% subset(ec_names, variable %in% env_ipca_across_trait_shared_list, variable_newname, drop = TRUE)) %>%
  as.matrix() %>%
  which()
underline_coord <- tibble(j = ceiling(underline / nrow(ec_env_sig_combine_towrite1))) %>%
  mutate(i = underline - ((j - 1) * nrow(ec_env_sig_combine_towrite1)), underlined = TRUE)

## Combine into a style df
style_df <- full_join(italicize_coord, underline_coord) %>% 
  mutate_at(vars(italic, underlined), ~ifelse(is.na(.), FALSE, .))



## Create the flextable
ft <- flextable(ec_env_sig_combine_towrite1) %>%
  set_header_df(mapping = typology, key = "col_keys") %>%
  # Merge cells with identical scenarios
  theme_booktabs() %>%
  border_remove() %>%
  ## Add borders
  border(i = 1, border.top = fp_border(width = 1.5), border.bottom = fp_border(width = 1.5), part = "header") %>%
  border(i = nrow(ec_env_sig_combine_towrite1), border.bottom = fp_border(width = 1.5), part = "body") %>%
  # Add borders to separate scenarios
  border(i = map_dbl(map(unique(ec_env_sig_combine_towrite1$set), ~which(ec_env_sig_combine_towrite1$set == .)), first)[-1],
         border.top = fp_border(width = 1)) %>% 
  autofit() %>%
  ## Align first column to left, subsequent column centered
  align(j = 1, align = "left", part = "all") %>%
  align(j = seq(2, ncol(ec_env_sig_combine_towrite1)), align = "center", part = "all") %>%
  # Vertical align the first column to top
  valign(x = ., j = seq(ncol(ec_env_sig_combine_towrite1)), valign = "top") %>%
  merge_v(x = ., j = 1)

## Loop to italicize and underline
fti <- ft
for (i in seq(nrow(style_df))) {
  fti <- style(x = fti, i = style_df$i[i], j = style_df$j[i], pr_t = fp_text(italic = style_df$italic[i], underlined = style_df$underlined[i]))
}

ft1 <- fti

## Output to MS word
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft1)
print(doc, target = file.path(fig_dir, "significant_ipca_ecs_paper.docx"))

    














## Group by trait and find number of common ECs between env_mean and IPCA
bind_rows(env_mean_topn_cor1, env_ipca_topn_cor) %>%
  filter(ec_group == "multiyear", trait %in% traits, timeframe == "growth_stage") %>%
  distinct(set, trait, test, variable) %>%
  group_by(set, trait, variable) %>%
  filter(n() > 1)

# set           trait       test     variable        
# 1 realistic2017 HeadingDate env_mean flowering_TRANGE
# 2 realistic2017 HeadingDate env_ipca flowering_TRANGE

## No overlap for GY or PH

## What about overall?



## Common significant ECs
bind_rows(mutate(ec_env_mean_sig, test = "EC_Mean"), mutate(ec_env_ipca_sig, test = "EC_IPCA")) %>%
  filter(ec_group == "multiyear", trait %in% traits, timeframe == "growth_stage") %>% 
  distinct(trait, set, variable_newname, test) %>%
  group_by(variable_newname) %>% 
  summarize(n = n(), n_trait = n_distinct(trait), n_test = n_distinct(test)) %>%
  arrange(desc(n_trait), desc(n_test))

bind_rows(mutate(ec_env_mean_sig, test = "EC_Mean"), mutate(ec_env_ipca_sig, test = "EC_IPCA")) %>%
  filter(ec_group == "multiyear", trait %in% traits, timeframe == "growth_stage") %>% 
  distinct(trait, set, variable_newname, test) %>%
  group_by(variable_newname, trait) %>% 
  summarize(n = n(), n_trait = n_distinct(trait), n_test = n_distinct(test)) %>%
  arrange(desc(n_test))


## What is the most frequently seen variable?
bind_rows(env_mean_topn_cor1, env_ipca_topn_cor) %>%
  filter(ec_group == "multiyear", trait %in% traits, timeframe == "growth_stage") %>%
  distinct(set, trait, test, variable)  %>%
  group_by(variable) %>%
  summarize(n = n(), nTrait = n_distinct(trait)) %>%
  arrange(desc(nTrait))

# variable                n nTrait
# 1 om_r_topsoil            4      3
# 2 annual_TAVG             3      2
# 3 annual_TMAX             2      2
# 4 annual_TSEASON          3      2
# 5 claytotal_r_subsoil     2      2
# 6 flowering_TMAX          2      2
# 7 flowering_TMIN          2      2
# 8 flowering_TRANGE        3      2
# 9 flowering_TSEASON       2      2
# 10 grainfill_TAVG          2      2



# ### Poster figures
# 
# 
# ## For each trait, plot the ECs most correlated with the mean or the IPCA score
# env_mean_cor_top <- env_mean_cor_sig %>%
#   filter(set == "complete", ec_group == "multi_year", variable != "ph1to1h2o_r_subsoil") %>%
#   left_join(., ec_env_df3) %>%
#   mutate(variable1 = case_when(
#     variable == "interval_2_PPT" ~ "Month 2 Precipitation",
#     variable == "interval_2_TMIN" ~ "Month 2 Average Min. Temp.",
#     variable == "annual_TSEASON" ~ "Annual Temperature Seasonality (Max. - Min.)"
#   )) %>%
#   group_by(trait) %>% 
#   top_n(n = 1, wt = -p_value) %>%
#   ggplot(aes(x = value, y = h)) + 
#   geom_smooth(method = "lm", se = FALSE) + 
#   geom_point() + 
#   facet_wrap( ~ trait + variable, scale = "free", labeller = labeller(trait = str_add_space), ncol = 1) + 
#   ylab(expression("Environmental effect (deviation from"~mu*")")) +
#   xlab("Covariate value") +
#   scale_x_continuous(breaks = pretty) +
#   scale_y_continuous(breaks = pretty) +
#   labs(title = expression(Example~italic(Mean-EC)~covariables)) +
#   theme_presentation2()
# 
# ggsave(filename = "env_mean_cor_top_poster.jpg", plot = env_mean_cor_top, path = fig_dir, width = 5, height = 10, dpi = 1000)
# 
# 
# ## Subset the most important covariates
# env_ipca_cor_sig_plot <-  env_ipca_cor_sig %>%
#   filter(set == "complete", ec_group == "multi_year") %>%
#   group_by(trait) %>% 
#   top_n(n = 1, wt = -p_value) %>%
#   ungroup() %>%
#   left_join(., filter(ec_score_df, PC == "PC1")) %>%
#   left_join(., unnest(env_ipca_mr1, cors)) %>%
#   mutate(variable1 = case_when(
#     # variable == "max_PPT" ~ "Max. Monthly Precipitation",
#     variable == "interval_2_PPT" ~ "Month 2 Precip.",
#     variable == "interval_1_TMAX" ~ "Month 1 Average Max. Temp.",
#     variable == "isothermality" ~ "Annual Isothermality"
#   )) 
# 
# 
# 
# 
# env_ipca_cor_top <-env_ipca_cor_sig_plot %>%
#   ggplot(aes(x = value, y = score)) + 
#   geom_smooth(method = "lm", se = FALSE) + 
#   geom_point() + 
#   # geom_text(aes(x = Inf, y = -Inf, label = paste0("r = ", formatC(correlation, 2))), hjust = 1, vjust = -1) + 
#   facet_wrap( ~ trait + variable, scale = "free", labeller = labeller(trait = str_add_space), ncol = 1) + 
#   ylab("Environmental IPCA score") +
#   xlab("Covariate value") +
#   scale_x_continuous(breaks = pretty) +
#   scale_y_continuous(breaks = pretty) +
#   labs(title = expression(Example~italic(IPCA-EC)~covariables)) +
#   theme_presentation2()
# 
# ggsave(filename = "env_ipca_cor_top_poster.jpg", plot = env_ipca_cor_top, path = fig_dir, width = 5, height = 10, dpi = 1000)















## Data.frame with all ECs for all sets
ec_all <- env_combine %>%
  filter(set == "complete", trait %in% traits) %>%
  select(-set) %>%
  crossing(., set = unique(ec_env_ipca_sig$set))

## Data.frame of EC data for all environments
ec_combine_use <- unnest(tibble(trait = traits, environment = tp_vp_env_trait[traits])) %>%
  left_join(., environ_covariate_data2) %>%
  select(trait, ec_group, timeframe, variable_newname, environment, value)

## Create a data.frame with each variable chosen for each set, trait, ec_group, timeframe, and test
env_variable_combine <- bind_rows(
  distinct(ec_env_mean_sig, trait, set, ec_group,  timeframe, variable_newname) %>% mutate(group = "EC_Mean"),
  distinct(ec_env_ipca_sig, trait, set, ec_group,  timeframe, variable_newname) %>% mutate(group = "EC_IPCA"),
  distinct(ec_all, trait, ec_group, timeframe, variable_newname, set) %>% mutate(group = "EC_All")
  # Combine with the large EC dataset
) %>% left_join(., ec_combine_use)


## Compare the overlap of variables that are significantly correlated with the mean
## and correlated with the IPCA score
env_variable_combine %>%
  filter(group != "EC_All") %>%
  distinct(set, ec_group, trait, group, timeframe, variable_newname) %>%
  split(.$group) %>%
  map(~select(., -group)) %>%
  reduce(dplyr::intersect)

# set           ec_group  trait       timeframe    variable_newname
# 1 realistic2017 multiyear GrainYield  growth_stage FL_Prcp         
# 2 realistic2017 multiyear HeadingDate growth_stage VG_MinTemp      
# 3 realistic2017 multiyear PlantHeight growth_stage OrgMatTS        
# 4 realistic2017 multiyear PlantHeight growth_stage GF_Prcp



ec_mats <- env_variable_combine %>%
  group_by(set, ec_group, timeframe, trait, group) %>% 
  do(mat = {
    df <- .
    df %>% 
      select(environment, variable_newname, value) %>% 
      spread(variable_newname, value) %>% 
      as.data.frame() %>% 
      column_to_rownames("environment") %>% 
      as.matrix()
  }) %>% ungroup()


## Create similarity matrices according to Malosetti 2016
rel_mat_df <- ec_mats %>%
  mutate(sim_mat = map(mat, ~{
    
    mat <- .x
    
    ## Apply over columns
    rel_mat <- apply(X = mat, MARGIN = 2, FUN = function(z) {
      # calculate relationship matrix
      list(as.matrix(dist(z)) / diff(range(z)))
    }) %>% map(1) %>%
      # Sum
      reduce(.x = ., .f = `+`)
      
    # subtract the values from 1
    1 - rel_mat
    
  }))
  


## Save the distance matrices
save_file <- file.path(result_dir, "environmental_covariable_distance_mat.RData")
save("ec_mats", "rel_mat_df", file = save_file)




# ################################
# ## Phenotypic analysis
# ################################
# 
# # Packages
# library(sommer)
# 
# # Phenotype data to model
# pheno_tomodel <- S2_MET_BLUEs %>%
#   filter(trait %in% traits) %>%
#   filter(line_name %in% tp) %>%
#   group_by(trait) %>%
#   nest() %>%
#   left_join(., filter(rel_mat_df, set == "complete"))
#   
# 
# # Use the covariates to model GxE
# ec_gxe_modeling <- pheno_tomodel %>%
#   group_by(trait, set, ec_group, timeframe, group) %>%
#   do({
#     row <- .
#     df <- row$data[[1]] %>%
#       mutate_at(vars(line_name, environment), as.factor) %>%
#       mutate(ge = interaction(line_name, environment, drop = TRUE, sep = ":"))
#     
#     # Environment relationship
#     Emat <- row$sim_mat[[1]]
#     # Genomic relationship
#     K <- diag(nlevels(df$line_name)) %>%
#       `dimnames<-`(., list(levels(df$line_name), levels(df$line_name)))
#     KE <- kronecker(X = K, Y = Emat, make.dimnames = TRUE)
#     KE <- KE[levels(df$ge), levels(df$ge)]
#     
#     ## Fit a model - use sommer
#     fit <- mmer(fixed = value ~ 1, random = ~ line_name + environment +
#                   vs(ge, Gu = KE), data = df)
#     
#     
#     
#     fit <- relmatLmer(value ~ (1|line_name) + (1|environment) + (1|ge) +,
#                       data = df, relmat = list(ge = KE))
#     
#     # Calculate variance explained
#     
#     
#     
#   })












