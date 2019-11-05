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
  

# # Remove intervals with insufficient observations
# environ_covariate_data2 <- environ_covariate_data1 %>% 
#   filter(str_detect(variable, "interval")) %>%
#   mutate(interval = str_extract(variable, "interval_[0-9]"), 
#          variable = str_remove(variable, "interval_[0-9]_")) %>% 
#   group_by(variable, interval, ec_group) %>% 
#   filter(n() == n_distinct(.$environment)) %>%
#   ungroup() %>%
#   unite(variable, interval, variable, sep = "_") %>%
#   bind_rows(., filter(environ_covariate_data1, !str_detect(variable, "interval"))) %>%
#   # Center each covariate and scale to achive the squared length of the vector = 1
#   split(list(.$ec_group, .$variable, .$timeframe)) %>%
#   map_df(~mutate(., scaled_value = scale_length(x = value),
#                  center = attr(scaled_value, "scaled:center"),
#                  scale = attr(scaled_value, "scaled:scale"),
#                  scaled_value = as.numeric(scaled_value)))


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
                          "om_r_subsoil"))

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
  
  st <- env_mean_mr$set[i]
  tr <- env_mean_mr$trait[i]
  grp <- env_mean_mr$ec_group[i]
  tf <- env_mean_mr$timeframe[i]

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



## Number of unique ECs
common_ecs <- env_mean_mr2 %>% 
  filter(set %in% c("complete", "realistic2017")) %>%
  filter(trait %in% traits) %>%
  filter(timeframe == "growth_stage") %>%
  pull(ecs) %>%
  unlist() %>%
  unique()
  
length(common_ecs)         

## 25 with AGDD
## 21 without


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




# ##### Use partial least squares
# library(pls)
# 
# ## This helps avoid the issue of collinearity between ECs
# 
# 
# 
# ## Try multiple regression, where ECs are added in order of correlation
# env_mean_pls <- ec_env_df3 %>%
#   distinct(set, trait, ec_group) %>%
#   mutate(out = list(NULL))
# 
# # Number of permutations
# perm_rep <- 5000
# 
# for (i in seq(nrow(env_mean_pls))) {
#   
#   st <- env_mean_pls$set[i]
#   tr <- env_mean_pls$trait[i]
#   grp <- env_mean_pls$ec_group[i]
#   
#   df <- ec_env_df3 %>% 
#     filter(trait == tr, set == st, ec_group == grp)
#   
#   ## Run correlations and rank the ECs by correlation
#   cors <- df %>% 
#     group_by(variable) %>% 
#     do(test = cor.test(.$h, .$value)) %>% 
#     ungroup() %>%
#     mutate(correlation = map_dbl(test, "estimate"),
#            pvalue = map_dbl(test, "p.value")) %>%
#     arrange(desc(abs(correlation)))
#   
#   # Extract EC names in order of decreasing corelation
#   ecs <- ordered(cors$variable)
#   
#   # Formula
#   base_formula <- formula(paste("h ~ ", paste(ecs[1], collapse = " + ")))
#   full_formula <- formula(paste("h ~ ", paste(ecs, collapse = " + ")))
#   
# 
#   # Create a model frame
#   df1 <- df %>% 
#     select(trait, environment, h, variable, value) %>% 
#     spread(variable, value)
#   
#   ## Create matrices of response and predictors
#   Y <- df1 %>% 
#     select(environment, h) %>% 
#     as.data.frame() %>% 
#     column_to_rownames("environment") %>% 
#     as.matrix()
#   
#   X <- df1 %>% 
#     select(-h, -trait) %>% 
#     as.data.frame() %>% 
#     column_to_rownames("environment") %>% 
#     as.matrix()
#   
#   ## fit the pls model
#   fit <- plsr(Y ~ X, validation = "LOO")
#   # Determine the number of components by minimizing the RMSE
#   # This is minimum RMSE when the number of components is >= 1
#   ncomp_rmse <- which.min(RMSEP(fit)$val["CV",,-1])
#   
#   # Refit the model
#   fit1 <- plsr(Y ~ X, validation = "LOO", ncomp = ncomp_rmse)
#   # Coefficients
#   fit1_coef <- coef(fit1)
#   
#   ## Permutation
#   # Randomize the response and predictors, get coefficients, establish a null distribution for coefficients
#   perm_out <- matrix(data = NA, nrow = ncol(X), ncol = perm_rep, dimnames = list(colnames(X), paste0("rep", seq(perm_rep))))
#   
#   for (r in seq(perm_rep)) {
#     
#     Y_perm <- apply(X = Y, MARGIN = 2, FUN = sample)
#     X_perm <- apply(X = X, MARGIN = 2, FUN = sample)
#     
#     # Fit
#     fit_perm <- plsr(Y_perm ~ X_perm,  validation = "LOO", ncomp = ncomp_rmse)
#     # Coefficients
#     fit_perm_coef <- coef(fit_perm)
#     
#     perm_out[,r] <- fit_perm_coef
#     
#   }
#   
#   ## Determine what ECs are significant
#   perm_out1 <- cbind(fit = abs(fit1_coef), perm_out)
#   plsr_pvalue <- apply(X = perm_out1, MARGIN = 1, FUN = function(p) mean(p[-1] >= p[1]))
#   
#   ## Return the ECs that are significant
#   env_mean_pls$out[[i]] <- names(which(plsr_pvalue <= alpha))
#   
#   
# }








## Now correlated ECs with the AMMI environmental IPCA scores

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
  
  st <- row$set[1]
  tr <- row$trait[1]
  grp <- row$ec_group[1]
  tf <- row$timeframe[1]
  
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



## Number of unique ECs
common_ecs <- env_ipca_mr2 %>% 
  filter(set %in% c("complete", "realistic2017")) %>%
  filter(trait %in% traits) %>%
  filter(timeframe == "growth_stage") %>%
  pull(ecs) %>%
  unlist() %>%
  unique()

length(common_ecs)         

## 25 with AGDD
## 21 without


# ## For each set, trait, and ec_group, determine the overlap in ECs that were determined
# ## to be significant
# (overlap <- env_ipca_mr2 %>% 
#     select(set:timeframe, ecs) %>% 
#     spread(timeframe, ecs) %>%
#     mutate(ec_overlap = map2(.x = growth_stage, .y = interval, intersect)))

## Determine EC overlap between LOYO years
(overlap <- env_ipca_mr2 %>%
    select(set:timeframe, ecs) %>%
    filter(set != "complete") %>%
    spread(set, ecs) %>%
    mutate(ec_overlap = select(., contains("realistic")) %>% pmap(., ~reduce(., intersect))) %>%
    mutate_at(vars(contains("realistic")), list(~map2(., ec_overlap, setdiff))) )



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









## Combine mean-correlated and ipca-correlated variables and output a table
## Also note which ones were used in MR
ec_top5_cor_towrite <- bind_rows(env_mean_topn_cor1, env_ipca_topn_cor) %>%
  filter(ec_group == "multiyear", trait %in% traits, timeframe == "growth_stage") %>%
  distinct(set, trait, variable_newname, test, cor) %>%
  group_by(set, trait, test) %>%
  mutate(top = order(abs(cor), decreasing = T)) %>%
  ungroup() %>%
  mutate(annotation = paste0(variable_newname, " (", round(cor, 2), ")")) %>%
  select(set, trait, top, test, annotation) %>%
  spread(test, annotation) %>%
  mutate(set = str_replace_all(set, set_replace),
         trait = str_add_space(trait)) %>%
  rename(Set = set, Trait = trait, Rank =top, IPCA = env_ipca, Mean = env_mean)

write_csv(x = ec_top5_cor_towrite, path = file.path(fig_dir, "top5_correlated_ecs.csv"))


## Create a table with significant ECs
ec_env_sig_combine <- bind_rows(mutate(ec_env_mean_sig, group = "EC_Mean"), mutate(ec_env_ipca_sig, group = "EC_IPCA")) %>%
  filter(ec_group == "multiyear", trait %in% traits, timeframe == "growth_stage") %>%
  distinct(set, trait, variable_newname, group) %>% 
  # filter(set %in% c("complete", "realistic2017")) %>%
  group_by(set, group, trait) %>% 
  mutate(n = seq(n())) %>%
  ungroup() %>%
  mutate(set = str_replace_all(set, set_replace),
         trait = str_add_space(trait))

ec_env_sig_combine_towrite <- ec_env_sig_combine %>% 
  spread(trait, variable_newname) %>%
  select(-n)

write_csv(x = ec_env_sig_combine_towrite, path = file.path(fig_dir, "significant_ecs.csv"), na = " ")



## Kable extra
library(kableExtra)

## Officer and flextable
library(officer)
library(flextable)

ec_top5_cor_towrite1 <- ec_top5_cor_towrite %>% 
  gather(test, annotation, -Set:-Rank) %>%   
  unite("test", test, Trait, sep = ".") %>%
  spread(test, annotation)

## Microsoft word table

typology <- data.frame(
  col_keys = names(ec_top5_cor_towrite1),
  what = c("", " ", "IPCA-EC", "IPCA-EC", "IPCA-EC", "Mean-EC", "Mean-EC", "Mean-EC"),
  trait = c("Scenario", "Top", "Grain Yield", "Heading Date", "Plant Height", "Grain Yield", "Heading Date", "Plant Height"),
  stringsAsFactors = FALSE )

ft <- flextable(ec_top5_cor_towrite1) %>%
  set_header_df(mapping = typology, key = "col_keys") %>%
  merge_h(part = "header") %>%
  theme_booktabs() %>%
  border_remove() %>%
  # border(border.top = fp_border(), part = "header") %>%
  border(i = 1, border.top = fp_border(width = 1.5), part = "header") %>%
  border(i = 2, j = 3:8, border.top = fp_border(width = 1.5), part = "header") %>%
  border(i = 2, border.bottom = fp_border(width = 1.5), part = "header") %>%
  border(i = nrow(ec_top5_cor_towrite1), border.bottom = fp_border(width = 1.5), part = "body") %>%
  autofit() %>%
  align(j = 1, align = "left", part = "all") %>%
  align(j = 2:8, align = "center", part = "all")

## Output to MS word
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = file.path(fig_dir, "top5_correlated_ecs_paper.docx"))

    





# ## Alternative
# 
# ec_top5_cor_towrite <- bind_rows(env_mean_top5_cor, env_ipca_top5_cor) %>%
#   filter(ec_group == "multi_year") %>%
#   select(set, trait, top, variable_newname, test, cor) %>%
#   # mutate(annotation = test) %>%
#   mutate(annotation = paste0(variable_newname, " (", round(cor, 2), ")"),
#          set = str_replace_all(set, set_replace),
#          test = ifelse(test == "env_mean", "Mean-EC", "IPCA-EC")) %>%
#   select(Scenario = set, trait, Distance.method = test, Top = top, annotation) %>%
#   # unite("test", test, trait, sep = ".") %>%
#   spread(trait, annotation) %>%
#   arrange(Scenario, Distance.method)
# 
# typology <- data.frame(
#   col_keys = names(ec_top5_cor_towrite),
#   trait = c("Scenario", "Distance method", "Top", "Grain Yield", "Heading Date", "Plant Height"),
#   stringsAsFactors = FALSE )
# 
# ft <- flextable(ec_top5_cor_towrite) %>%
#   set_header_df(mapping = typology, key = "col_keys") %>%
#   merge_h(part = "header") %>%
#   theme_booktabs() %>%
#   border_remove() %>%
#   # border(border.top = fp_border(), part = "header") %>%
#   border(i = 1, border.top = fp_border(width = 1.5), border.bottom = fp_border(width = 1.5), part = "header") %>%
#   # border(i = 2, j = 3:8, border.top = fp_border(width = 1.5), part = "header") %>%
#   # border(i = 1, border.bottom = fp_border(width = 1.5), part = "header") %>%
#   border(i = nrow(ec_top5_cor_towrite), border.bottom = fp_border(width = 1.5), part = "body") %>%
#   autofit() %>%
#   align(j = 1, align = "left", part = "all") %>%
#   align(j = seq(2, ncol(ec_top5_cor_towrite)), align = "center", part = "all")
# 
# ## Output to MS word
# doc <- read_docx()
# doc <- body_add_flextable(doc, value = ft)
# print(doc, target = file.path(fig_dir, "top5_correlated_ecs_paper_alt.docx"))



















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
# 1 interval_1_PPT          8      3
# 2 interval_3_TMAX         5      3
# 3 latitude                5      3
# 4 min_PPT                 5      3
# 5 om_r_topsoil            6      3
# 6 silttotal_r_topsoil     3      3
# 7 annual_TRANGE_max       2      2
# 8 annual_TSEASON          3      2
# 9 interval_1_TSEASON      4      2
# 10 interval_2_PPT          7      2



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

















# ## Generate random samples of environments. For each sample, calculate the same
# ## correlation as above. How often is the same EC strongly correlated?
# n_sample <- 1000
# f_sample <- 0.6
# 
# set.seed(500)
# env_ipca_samples <- env_score_combine %>%
#   group_by(ec_group, trait, IPCA) %>%
#   do({
#     df <- .
#     # Generate samples of environments
#     cor_sample <- replicate(n_sample, sample(unique(df$environment), size = ceiling(f_sample * n_distinct(df$environment))), simplify = FALSE) %>% 
#       map(~filter(df, environment %in% .)) %>% 
#       map_df(~group_by(., variable) %>% summarize(cor = cor(score, value)) %>% mutate(rank = rank(-abs(cor))))
#     
#     # Return the average rank and average correlation of the variables
#     cor_sample %>% 
#       group_by(variable) %>% 
#       summarise_at(vars(cor, rank), mean) %>% 
#       arrange(rank)
#     
#   })
# 
# 
# 
# ## Check the selected covariables from above against the rank generated
# ## from the samples
# env_ipca_cor_sig1 %>% 
#   left_join(., rename(env_ipca_samples, sample_cor = cor)) %>%
#   as.data.frame()



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
sim_mat_df <- ec_mats %>%
  group_by(set, ec_group, trait, group) %>%
  do(sim_mat = {
    df <- .
    mat <- df$mat[[1]]

    # Calculate a similiarity matrix for each covariate
    sim_mats <- seq(ncol(mat)) %>%
      map(~mat[,.]) %>%
    # sim_mats <- seq(nrow(mat)) %>%
    #   map(~mat[.,]) %>%
      map(~{
        # First calculate the difference in the range
        range_diff <- diff(range(.))
        # Then calculate the distance between environments standardized to that range
        sim <- (outer(., ., `-`)) / range_diff
      })
    
    ## Sum over similarity matrices
    1 - (apply(X = simplify2array(x = sim_mats), MARGIN = c(1,2), FUN = sum) / length(sim_mats))
  
  }) %>% ungroup()


## Create similarity matrices according to Jarquin 2014
sim_mat_df1 <- ec_mats %>% 
  mutate(sim_mat = map(mat, ~tcrossprod(.) / ncol(.)))

## Create a different similarity matrix using the standarized squared difference between
## environments
sim_mat_df2 <- ec_mats %>%
  group_by(set, ec_group, trait, group) %>%
  do(sim_mat = {
    df <- .
    mat <- df$mat[[1]]
    
    # Calculate a similiarity matrix for each covariate
    sim_mats <- seq(ncol(mat)) %>%
      map(~mat[,.]) %>%
      # sim_mats <- seq(nrow(mat)) %>%
      #   map(~mat[.,]) %>%
      map(~{
        # First calculate the difference in the range
        range_diff <- diff(range(.))
        # Then calculate the distance between environments standardized to that range
        sim <- (outer(., ., `-`)^2) / range_diff
      })
    
    ## Sum over similarity matrices
    1 - (apply(X = simplify2array(x = sim_mats), MARGIN = c(1,2), FUN = sum) / length(sim_mats))
    
  }) %>% ungroup()



## Save the distance matrices
save_file <- file.path(result_dir, "environmental_covariable_distance_mat.RData")
save("ec_mats", "sim_mat_df", "sim_mat_df1", "sim_mat_df2", file = save_file)






















##### Appendix


## Comment-out the interval variables. We might revisit this.

# ## What interval variables are highly correlated
# one_year_env_interval_df1 <- one_year_daily_summary_interval %>% 
#   map(unnest) %>%
#   map(select, -variable) %>% 
#   reduce(full_join) %>% 
#   gather(variable, value, -begin:-environment) %>% 
#   left_join(env_means_all, .) %>%
#   # Remove NA
#   filter(!is.na(value))
# 
# env_interval_cor <- one_year_env_interval_df1 %>% 
#   group_by(population, trait, variable, begin, end) %>% 
#   do(test = cor.test(.$h, .$value)) %>%
#   ungroup() %>%
#   mutate(cor = map_dbl(test, "estimate"),
#          pvalue = map_dbl(test, "p.value"),
#          df = map_dbl(test, "parameter")) %>%
#   select(-test)
# 
# ## Plot
# env_interval_cor_plot_list <- env_interval_cor %>% 
#   filter(population == "all") %>%
#   split(list(.$trait, .$variable)) %>% 
#   map(~ggplot(data = ., aes(x = begin, y = end, fill = cor)) +
#         geom_tile() +
#         scale_fill_gradient2() + 
#         facet_wrap(~ trait + variable) +
#         theme_acs() + theme(legend.position = c(0.85, 0.35)))
# 
# # Cowplot
# env_interval_cor_plot <- plot_grid(plotlist = env_interval_cor_plot_list, ncol = n_distinct(env_interval_cor$trait))
# ggsave(filename = "one_year_interval_variable_correlation_space.jpg", plot = env_interval_cor_plot,
#        path = fig_dir, width = 12, height = 12, dpi = 1000)
# 
# 
# ## Plot
# env_interval_cor_plot_list <- env_interval_cor %>% 
#   filter(population == "tp") %>%
#   split(list(.$trait, .$variable)) %>% 
#   map(~qplot(x = begin, y = end, fill = cor, geom = "tile", data = .) +
#         scale_fill_gradient2() + facet_wrap(~ trait + variable) +
#         theme_acs() + theme(legend.position = c(0.85, 0.35)))
# 
# # Cowplot
# env_interval_cor_plot <- plot_grid(plotlist = env_interval_cor_plot_list, ncol = n_distinct(env_interval_cor$trait))
# ggsave(filename = "one_year_interval_variable_correlation_space_tp.jpg", plot = env_interval_cor_plot,
#        path = fig_dir, width = 12, height = 12, dpi = 1000)
# 
# 
# 
# 
# 
# 
# 
# 
# ## Find the intervals most correlated with the environmental mean
# ## Weigh the expected direction of the correlation appropriately
# env_interval_cor_best <- env_interval_cor %>% 
#   mutate(wt_cor = ifelse(str_detect(trait, "HeadingDate"), cor * -1, cor)) %>%
#   # Filter for appropriate correlations
#   filter(sign(wt_cor) == 1) %>%
#   group_by(population, trait, variable) %>%
#   # filter(variable %in% c("AGDD", "GDD")) %>%
#   # filter(pvalue <= 0.05) %>%
#   top_n(n = 10, wt = wt_cor) # Take the top 25
# 
# 
# 
# 
# ## Create extra variables to use
# env_interval_use <- env_interval_cor_best %>% 
#   left_join(one_year_env_interval_df1) %>%
#   unite(variable, variable, begin, end, sep = "_")

##