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
library(lme4)
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
ec_date_replace <- c("annual_" = "Anno", "max_" = "Max", "min_" = "Min", "interval_" = "Int")



### Fit a joint regression model
##
## Models are fitted using all data or using just 2015 - 2016 data
## 


# # First create a list of data split by whether all lines are used or just the TP
# data_to_model <- list(tp = tp, all = c(tp, vp)) %>% 
#   map(~filter(S2_MET_BLUEs, line_name %in% .)) %>%
#   list(., names(.)) %>%
#   pmap_df(~mutate(.x, population = .y))

data_to_model <- training_sets_twoway %>%
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
  filter(n() == 2) %>% 
  spread(set, h) %>% 
  group_by(trait) %>% 
  summarize(set_cor = cor(complete, realistic))

# trait           set_cor
# 1 GrainYield        1.000
# 2 HeadingDate       1.000
# 3 HeadingDateAGDD   1.000
# 4 PlantHeight       1.000

# This is not surprising

    





## Calculate regression coefficients for each genotype, then test the signficance of differences among regression coefficients
data_to_model %>% 
  left_join(., env_means_all) %>%
  group_by(trait, set) %>%
  do(fwr(formula = value ~ h, data = .))

# trait           set       regression         Fstat   pvalue
# 1 GrainYield      complete  <tibble [175 x 4]>  3.33 1.75e-41
# 2 GrainYield      realistic <tibble [175 x 4]>  2.62 6.21e-24
# 3 HeadingDate     complete  <tibble [175 x 4]>  1.57 4.09e- 6
# 4 HeadingDate     realistic <tibble [175 x 4]>  2.01 1.15e-12
# 5 HeadingDateAGDD complete  <tibble [175 x 4]>  1.76 7.15e- 9
# 6 HeadingDateAGDD realistic <tibble [175 x 4]>  1.59 2.82e- 6
# 7 PlantHeight     complete  <tibble [175 x 4]>  1.52 1.88e- 5
# 8 PlantHeight     realistic <tibble [175 x 4]>  1.50 4.88e- 5

## As expected, all are significant, but GY is most significant and HD is least
## Since all traits show genotype-specific reactions to the environment mean, genotypes should also
# demonstrate different responses to variables that are correlated with the mean



## Correlate environmental variables with the mean

## First combine one-year with multi-year ECs
ec_env_df <- bind_rows(mutate(one_year_env_df, ec_group = "one_year"),
                       mutate(multi_year_env_df, ec_group = "multi_year")) %>%
  ## Adjust names
  mutate(variable_newname = str_replace_all(variable, ec_type_replace) %>% str_replace_all(ec_date_replace))


## A separate DF of variable names and new names
ec_names <- distinct(ec_env_df, variable, variable_newname)

## Write a csv
write_csv(x = ec_names, path = file.path(fig_dir, "env_covariate_names.csv"))


## A functio to scale a vector to dot product sum = 1
scale_length <- function(x) {
  # First center at 0
  x1 <- scale(x, center = TRUE, scale = FALSE)
  # Next scale to obtain 1 as the squared vector length
  scale(x1, center = FALSE, scale = sqrt(sum(x1^2)))
}
  


# Remove intervals with insufficient observations
ec_env_df1 <- ec_env_df %>% 
  filter(str_detect(variable, "interval")) %>%
  mutate(interval = str_extract(variable, "interval_[0-9]"), 
         variable = str_remove(variable, "interval_[0-9]_")) %>% 
  group_by(variable, interval, ec_group) %>% 
  filter(n() == n_distinct(.$environment)) %>%
  ungroup() %>%
  unite(variable, interval, variable, sep = "_") %>%
  bind_rows(., filter(ec_env_df, !str_detect(variable, "interval"))) %>%
  # Center each covariate and scale to achive the squared length of the vector = 1
  split(list(.$ec_group, .$variable)) %>%
  map_df(~mutate(., scaled_value = scale_length(x = value),
              center = attr(scaled_value, "scaled:center"),
              scale = attr(scaled_value, "scaled:scale"),
              scaled_value = as.numeric(scaled_value)))



## Re-organization and removal of some covariates
ec_env_df2 <- ec_env_df1 %>%
  left_join(env_means_all, .) %>%
  mutate(EC_type = "summary") %>%
  # filter(! variable %in% c("elevation", "latitude", "longitude")) %>%
  filter(! variable %in% c("elevation")) %>%
  # Remove some variables that are not well distributed
  filter(!variable %in% c("claytotal_r_topsoil", "sandtotal_r_subsoil", "sandtotal_r_topsoil", "silttotal_r_subsoil",
                          "om_r_subsoil"))

## Correlate covariates with the environmental mean
## Here we used the actual value of the covariates
env_mean_cor <- ec_env_df2 %>%
  group_by(set, trait, ec_group, variable) %>% 
  do(test = cor.test(.$h, .$value)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         df = map_dbl(test, "parameter")) %>%
  select(-test)





# Remove some outliers and re-run
ec_env_df3 <- ec_env_df2 %>%
  filter(!(trait == "GrainYield" & h > 3000),
         !(trait == "PlantHeight" & h > 25))

env_mean_cor1 <- ec_env_df3 %>%
  group_by(set, trait, ec_group, variable) %>% 
  do(test = cor.test(.$h, .$value)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         df = map_dbl(test, "parameter")) %>%
  select(-test)



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
env_mean_mr <- ec_env_df3 %>%
  distinct(set, trait, ec_group) %>%
  mutate(out = list(NULL))
  
for (i in seq(nrow(env_mean_mr))) {
  
  st <- env_mean_mr$set[i]
  tr <- env_mean_mr$trait[i]
  grp <- env_mean_mr$ec_group[i]

  df <- ec_env_df3 %>% 
    filter(trait == tr, set == st, ec_group == grp)

  ## Run correlations and rank the ECs by correlation
  cors <- df %>% 
    group_by(variable) %>% 
    summarize(correlation = cor(h, value)) %>% 
    arrange(desc(abs(correlation)))
  
  # Extract EC names in order of decreasing corelation
  ecs <- ordered(cors$variable)
  
  # Formula
  base_formula <- formula(paste("h ~ ", paste(ecs[1], collapse = " + ")))
  full_formula <- formula(paste("h ~ ", paste(ecs, collapse = " + ")))
  
  
  # Create a model frame
  df1 <- df %>% 
    select(trait, environment, h, variable, value) %>% 
    spread(variable, value)

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
  
  ## Get predictions of the environmental mean
  df2 <- df1 %>% 
    add_predictions(fit_step) %>% 
    select(trait, environment, h, pred)
  
  ## Perform FW regression and test for heterogenity of slopes
  fit_step_fw <- fwr(formula = value ~ pred, data = left_join(df2, data_to_model, by = c("trait", "environment")))
  
  ## Return a tidy df of the model
  env_mean_mr$out[[i]] <- data_frame(
    fit = list(fit_step),
    tidy = list(fit_step_tidy),
    fw = list(fit_step_fw))
  
}
    

env_mean_mr1 <- unnest(env_mean_mr) %>%
  mutate(fit_R2 = map_dbl(fit, ~summary(.)$adj.r.squared))


# ## Run multiple testing correction, then select the variables that are still significantly correlated
# env_mean_cor_sig <- env_mean_cor1 %>% 
#   group_by(set, ec_group, trait) %>% 
#   # mutate(ntest = n(), padj = p.adjust(p = pvalue, method = "bonf")) %>% 
#   mutate(padj = pvalue) %>% # Skip multiple test correction
#   filter(padj <= alpha) %>%
#   ungroup()
# 
# 
# ## Summarize
# env_mean_cor_sig %>% 
#   group_by(set, ec_group, trait) %>% 
#   summarize(n_var = n_distinct(variable), min_cor = min(cor), max_cor = max(cor))

env_mean_cor_sig <- env_mean_mr1 %>%
  unnest(tidy) %>%
  left_join(ec_names) %>%
  left_join(., env_mean_cor1)

## Summarize number of ECs per trait
env_mean_cor_sig %>% 
  filter(ec_group == "multi_year") %>% 
  mutate(nEC = n_distinct(variable)) %>% 
  group_by(trait, set, nEC) %>% 
  summarize(nEC_group = n_distinct(variable))


## Note if the top 5 




## Plot the significant results
g_mean_cor <- env_mean_cor_sig %>% 
  split(., list(.$set, .$ec_group)) %>%
  map(~{
    
    df <- .
    
    print(unique(df$ec_group))
    
    set <- unique(df$set)
    group <- unique(df$ec_group)
    
    ## Determine the image size
    n_col <- 4
    n_row <- ceiling(nrow(df) / n_col)

    ## image height / width
    width <- 8
    height <- 6 * (n_row / (6/2))

    ## Create a facet annotation
    df1 <- mutate(df, annotation = paste0(variable_newname, " (r = ", round(cor, 2), ")"))

    g <- df1 %>%
      select(trait, set, ec_group, variable, annotation) %>%
      left_join(., ec_env_df3, by = c("set", "trait", "ec_group", "variable")) %>%
      ggplot(aes(x = value, y = h)) + 
      geom_smooth(method = "lm", se = FALSE) + 
      geom_point() + 
      facet_wrap( ~ trait + annotation, scale = "free", nrow = n_row, ncol = n_col) + 
      scale_y_continuous(breaks = pretty, name = "Environmental mean") +
      scale_x_continuous(breaks = pretty, name = "Covariate value") +
      theme_acs()

      ggsave(filename = paste0("env_mean_cor_sig_", group, "_", set, ".jpg"), path = fig_dir, plot = g, 
             width = width, height = height, dpi = 1000)
    
  })
    

    



## Find overlapping variables between sets and plot
env_mean_cor_set_overlap <- env_mean_cor_sig %>% 
  split(.$set) %>% 
  map(~select(., trait, ec_group, variable, variable_names)) %>% 
  reduce(., dplyr::intersect)

# trait       ec_group   variable          
# 1 GrainYield  multi_year interval_1_TSEASON
# 2 GrainYield  one_year   interval_2_PPT    
# 3 HeadingDate multi_year interval_1_TSEASON
# 4 HeadingDate multi_year om_r_topsoil      
# 5 PlantHeight multi_year interval_2_TMAX   
# 6 PlantHeight multi_year annual_TSEASON    
# 7 PlantHeight one_year   interval_1_TRANGE 
# 8 PlantHeight one_year   interval_1_TSEASON






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


## Combine the summary variables with the interval variables
## Then select those variables that are significant
env_combine <- ec_env_df2 %>%
  # bind_rows(env_interval_use) %>%
  select(set, trait, ec_group, environment, h, variable_newname, value)

## Merge not based on set - we want the ECs that were determined significant using a subset of data
ec_env_mean_sig <- select(env_mean_cor_sig, set, trait, ec_group, variable_newname) %>%
  left_join(., distinct(env_combine, trait, ec_group, environment, variable_newname, value), by = c("trait", "ec_group", "variable_newname"))



# ## Generate random samples of environments. For each sample, calculate the same
# ## correlation as above. How often is the same EC strongly correlated?
# n_sample <- 1000
# f_sample <- 0.6
# 
# set.seed(500)
# env_samples <- env_combine %>%
#   group_by(population, ec_group, trait) %>%
#   do({
#     df <- .
#     # Generate samples of environments
#     cor_sample <- replicate(n_sample, sample(unique(df$environment), size = ceiling(f_sample * n_distinct(df$environment))), simplify = FALSE) %>% 
#       map(~filter(df, environment %in% .)) %>% 
#       map_df(~group_by(., variable) %>% summarize(cor = cor(h, value)) %>% mutate(rank = rank(-abs(cor))))
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
# 
# 
# ## Check the selected covariables from above against the rank generated
# ## from the samples
# env_mean_cor_sig %>% 
#   left_join(., rename(env_samples, sample_cor = cor)) %>%
#   as.data.frame()








## Now correlated ECs with the AMMI environmental IPCA scores

ec_score_df <- ammi_out %>% # Remember to use the data with outliers removed
  mutate(escore = map(ammi, "escores")) %>% 
  unnest(escore) %>% 
  left_join(ec_env_df3, ., by = c("trait", "set", "environment"))

## This one does not remove the outliers - 
ec_score_df_use <- ammi_out %>%
  mutate(escore = map(ammi, "escores")) %>% 
  unnest(escore) %>% 
  left_join(ec_env_df2, ., by = c("trait", "set", "environment"))
    
env_ipca_cor1 <- ec_score_df %>% 
  group_by(set, ec_group, trait, variable, PC) %>% 
  do(test = cor.test(.$score, .$scaled_value)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         df = map_dbl(test, "parameter")) %>%
  select(-test)



env_ipca_top5_cor <- env_ipca_cor1 %>%
  filter(PC == "PC1") %>%
  group_by(set, trait, ec_group) %>%
  top_n(n = 5, wt = abs(cor)) %>%
  arrange(set, ec_group, trait, desc(abs(cor))) %>%
  mutate(top = seq(5)) %>%
  ungroup() %>%
  mutate(test = "env_ipca") %>%
  left_join(., ec_names)
  






## Try multiple regression, where ECs are added in order of correlation
env_ipca_mr <- ec_score_df %>%
  distinct(set, trait, ec_group) %>%
  mutate(out = list(NULL))


for (i in seq(nrow(env_ipca_mr))) {
  
  st <- env_ipca_mr$set[i]
  tr <- env_ipca_mr$trait[i]
  grp <- env_ipca_mr$ec_group[i]
  
  df <- ec_score_df %>% 
    filter(trait == tr, set == st, ec_group == grp, PC == "PC1")
  
  ## Run correlations and rank the ECs by correlation
  cors <- df %>% 
    group_by(variable) %>% 
    summarize(correlation = cor(score, scaled_value)) %>% 
    arrange(desc(abs(correlation)))
  
  # Extract EC names in order of decreasing corelation
  ecs <- ordered(cors$variable)
  
  # Formula
  base_formula <- formula(paste("score ~ ", paste(ecs[1], collapse = " + ")))
  full_formula <- formula(paste("score ~ ", paste(ecs, collapse = " + ")))
  
  
  # Create a model frame
  df1 <- df %>% 
    filter(PC == "PC1") %>%
    select(trait, environment, score, variable, scaled_value) %>% 
    spread(variable, scaled_value)
  
  # Base model
  fit <- lm(formula = base_formula, data = df1)
  # Stepwise selected model
  fit_step <- step(object = fit, scope = full_formula, direction = "both", trace = 0)
  # Tidy
  fit_step_tidy <- tidy(anova(fit_step)) %>%
    filter(term != "Residuals") %>%
    select(variable = term, sumsq, p_value = p.value)
  
  ## Get predictions of the environmental score
  df2 <- df1 %>% 
    add_predictions(fit_step) %>% 
    select(trait, environment, score, pred)
  
  ## Return a tidy df of the model
  env_ipca_mr$out[[i]] <- data_frame(
    fit = list(fit_step),
    tidy = list(fit_step_tidy),
    cors = list(cors) )
  
}


env_ipca_mr1 <- unnest(env_ipca_mr) %>%
  mutate(fit_R2 = map_dbl(fit, ~summary(.)$adj.r.squared))


env_ipca_cor_sig <- env_ipca_mr1 %>%
  unnest(tidy) %>% 
  left_join(ec_names) %>%
  left_join(., filter(env_ipca_cor1, PC == "PC1"))


## Summarize number of ECs per trait
env_ipca_cor_sig %>% 
  filter(ec_group == "multi_year") %>% 
  mutate(nEC = n_distinct(variable)) %>% 
  group_by(trait, set, nEC) %>% 
  summarize(nEC_group = n_distinct(variable))




## Plot the significant results
g_ipca_cor <- env_ipca_cor_sig %>% 
  split(., list(.$set, .$ec_group)) %>%
  map(~{
    
    df <- .
    
    print(unique(df$ec_group))
    
    set <- unique(df$set)
    group <- unique(df$ec_group)
    
    ## Determine the image size
    n_col <- 4
    n_row <- ceiling(nrow(df) / n_col)
    
    ## image height / width
    width <- 8
    height <- 6 * (n_row / (6/2))
    
    ## Create a facet annotation
    df1 <- mutate(df, annotation = paste0(variable_newname, " (r = ", round(cor, 2), ")"))
    
    g <- df1 %>%
      select(trait, set, ec_group, variable, annotation) %>%
      left_join(., filter(ec_score_df, PC == "PC1"), by = c("set", "trait", "ec_group", "variable")) %>%
      ggplot(aes(x = value, y = score)) + 
      geom_smooth(method = "lm", se = FALSE) + 
      geom_point() + 
      facet_wrap( ~ trait + annotation, scale = "free", nrow = n_row, ncol = n_col) + 
      scale_y_continuous(breaks = pretty, name = "Environmental mean") +
      scale_x_continuous(breaks = pretty, name = "Covariate value") +
      theme_acs()
    
    ggsave(filename = paste0("env_ipca_cor_sig_", group, "_", set, ".jpg"), path = fig_dir, plot = g, 
           width = width, height = height, dpi = 1000)
    
  })







## Find overlapping variables between sets and plot
env_ipca_cor_set_overlap <- env_ipca_cor_sig %>% 
  split(.$set) %>% 
  map(~select(., trait, ec_group, variable)) %>% 
  reduce(., dplyr::intersect)




## Combine the summary variables with the interval variables
## Then select those variables that are significant
env_ipca_combine <- ec_score_df_use %>%
  filter(PC == "PC1") %>%
  # bind_rows(env_interval_use) %>%
  select(set, trait, ec_group, environment, score, variable, variable_newname, value)

## Merge not based on set - we want the ECs that were determined significant using a subset of data
ec_env_ipca_sig <- select(env_ipca_cor_sig, set, trait, ec_group, variable) %>%
  left_join(., distinct(env_ipca_combine, trait, ec_group, environment, variable, value), by = c("trait", "ec_group", "variable")) %>%
  left_join(ec_names)




## Combine mean-correlated and ipca-correlated variables and output a table
## Also note which ones were used in MR
ec_top5_cor_towrite <- bind_rows(env_mean_top5_cor, env_ipca_top5_cor) %>%
  filter(ec_group == "multi_year") %>%
  select(set, trait, top, variable_newname, test, cor) %>%
  # mutate(annotation = test) %>%
  mutate(annotation = paste0(variable_newname, " (", round(cor, 2), ")")) %>%
  select(set, trait, top, test, annotation) %>%
  spread(test, annotation) %>%
  rename(Set = set, Trait = trait, Rank = top, IPCA = env_ipca, Mean = env_mean)

write_csv(x = ec_top5_cor, path = file.path(fig_dir, "top5_correlated_ecs.csv"))


## Kable extra
library(kableExtra)

ec_top5_cor_towrite <- bind_rows(env_mean_top5_cor, env_ipca_top5_cor) %>%
  filter(ec_group == "multi_year") %>%
  select(set, trait, top, variable_newname, test, cor) %>%
  # mutate(annotation = test) %>%
  mutate(annotation = paste0(variable_newname, " (", round(cor, 2), ")")) %>%
  select(set, trait, top, test, annotation) %>%
  unite("test", trait, test, sep = ":") %>%
  spread(test, annotation) %>%
  rename_all(~str_replace_all(string = ., pattern = "([A-Za-z]*:)([a-z_]*)", "\\2"))

kable(ec_top5_cor_towrite, escape = TRUE, format = "html") %>%
  kable_styling("hover", full_width = F) %>%
  add_header_above(c(" " = 2, "Grain Yield" = 2, "Heading Date" = 2, "Plant Height" = 2)) %>%
  column_spec(column = 3:8, width = "10em") %>%
  # Save
  save_kable(file = file.path(fig_dir, "top5_correlated_ecs.html"), self_contained = T)




## Group by trait and find number of common ECs between env_mean and IPCA
bind_rows(env_mean_top5_cor, env_ipca_top5_cor) %>%
  filter(ec_group == "multi_year") %>%
  select(set, trait, test, variable) %>%
  group_by(set, trait, variable) %>%
  filter(n() > 1)

# set       trait      test     variable      
# 1 complete  GrainYield env_mean max_PPT       
# 2 realistic GrainYield env_mean max_PPT       
# 3 realistic GrainYield env_mean interval_2_PPT
# 4 complete  GrainYield env_ipca max_PPT       
# 5 realistic GrainYield env_ipca interval_2_PPT
# 6 realistic GrainYield env_ipca max_PPT

## No overlap for HD or PH

## What about overall?

## Group by trait and find number of common ECs between env_mean and IPCA
bind_rows(mutate(env_mean_cor_sig, test = "mean"), mutate(env_ipca_cor_sig, test = "ipca")) %>%
  filter(ec_group == "multi_year") %>%
  select(trait, test, variable) %>%
  group_by(trait, variable) %>%
  filter(n() > 1)


## Common significant ECs
bind_rows(mutate(env_mean_cor_sig, test = "mean"), mutate(env_ipca_cor_sig, test = "ipca"))  %>% 
  filter(ec_group == "multi_year") %>% 
  group_by(variable) %>% 
  summarize(n = n(), n_trait = n_distinct(trait), n_test = n_distinct(test)) %>%
  arrange(desc(n_trait), desc(n_test))

bind_rows(mutate(env_mean_cor_sig, test = "mean"), mutate(env_ipca_cor_sig, test = "ipca"))  %>% 
  filter(ec_group == "multi_year") %>% 
  group_by(variable, trait) %>% 
  summarize(n = n(), n_trait = n_distinct(trait), n_test = n_distinct(test)) %>%
  arrange(desc(n_test))


## What is the most frequently seen variable?
bind_rows(env_mean_top5_cor, env_ipca_top5_cor) %>%
  filter(ec_group == "multi_year") %>%
  select(set, trait, test, variable) %>%
  group_by(variable) %>%
  summarize(n = n(), nTrait = n_distinct(trait)) %>%
  arrange(desc(nTrait))

bind_rows(env_mean_top5_cor, env_ipca_top5_cor) %>%
  filter(ec_group == "multi_year") %>%
  filter(variable %in% c("min_PPT", "max_PPT")) %>%
  arrange(variable)



### Poster figures


## For each trait, plot the ECs most correlated with the mean or the IPCA score
env_mean_cor_top <- env_mean_cor_sig %>%
  filter(set == "complete", ec_group == "multi_year", variable != "ph1to1h2o_r_subsoil") %>%
  left_join(., ec_env_df3) %>%
  mutate(variable1 = case_when(
    variable == "interval_2_PPT" ~ "Month 2 Precipitation",
    variable == "interval_2_TMIN" ~ "Month 2 Average Min. Temp.",
    variable == "annual_TSEASON" ~ "Annual Temperature Seasonality (Max. - Min.)"
  )) %>%
  group_by(trait) %>% 
  top_n(n = 1, wt = -p_value) %>%
  ggplot(aes(x = value, y = h)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  geom_point() + 
  facet_wrap( ~ trait + variable, scale = "free", labeller = labeller(trait = str_add_space), ncol = 1) + 
  ylab(expression("Environmental effect (deviation from"~mu*")")) +
  xlab("Covariate value") +
  scale_x_continuous(breaks = pretty) +
  scale_y_continuous(breaks = pretty) +
  labs(title = expression(Example~italic(Mean-EC)~covariables)) +
  theme_presentation2()

ggsave(filename = "env_mean_cor_top_poster.jpg", plot = env_mean_cor_top, path = fig_dir, width = 5, height = 10, dpi = 1000)


## Subset the most important covariates
env_ipca_cor_sig_plot <-  env_ipca_cor_sig %>%
  filter(set == "complete", ec_group == "multi_year") %>%
  group_by(trait) %>% 
  top_n(n = 1, wt = -p_value) %>%
  ungroup() %>%
  left_join(., filter(ec_score_df, PC == "PC1")) %>%
  left_join(., unnest(env_ipca_mr1, cors)) %>%
  mutate(variable1 = case_when(
    # variable == "max_PPT" ~ "Max. Monthly Precipitation",
    variable == "interval_2_PPT" ~ "Month 2 Precip.",
    variable == "interval_1_TMAX" ~ "Month 1 Average Max. Temp.",
    variable == "isothermality" ~ "Annual Isothermality"
  )) 




env_ipca_cor_top <-env_ipca_cor_sig_plot %>%
  ggplot(aes(x = value, y = score)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  geom_point() + 
  # geom_text(aes(x = Inf, y = -Inf, label = paste0("r = ", formatC(correlation, 2))), hjust = 1, vjust = -1) + 
  facet_wrap( ~ trait + variable, scale = "free", labeller = labeller(trait = str_add_space), ncol = 1) + 
  ylab("Environmental IPCA score") +
  xlab("Covariate value") +
  scale_x_continuous(breaks = pretty) +
  scale_y_continuous(breaks = pretty) +
  labs(title = expression(Example~italic(IPCA-EC)~covariables)) +
  theme_presentation2()

ggsave(filename = "env_ipca_cor_top_poster.jpg", plot = env_ipca_cor_top, path = fig_dir, width = 5, height = 10, dpi = 1000)

















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
  filter(set == "complete") %>%
  bind_rows(., mutate(., set = "realistic"))



env_variable_combine <- bind_rows(
  mutate(ec_env_mean_sig, group = "EC_Mean"),
  mutate(ec_env_ipca_sig, group = "EC_IPCA"),
  ec_all %>% select(-h) %>% mutate(group = "EC_All")
  ) %>%
  ## Remove one-year covariates with the realistic scenario
  filter(!(set == "realistic" & ec_group == "one_year"))


## Compare the overlap of variables that are significantly correlated with the mean
## and correlated with the IPCA score
env_variable_combine %>% 
  filter(group != "EC_All") %>% 
  distinct(set, ec_group, trait, group, variable_newname) %>% 
  split(.$group) %>%
  map(~select(., -group)) %>%
  reduce(dplyr::intersect)





ec_mats <- env_variable_combine %>% 
  group_by(set, ec_group, trait, group) %>% 
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
save("ec_mats", "sim_mat_df", "sim_mat_df1", "sim_mat_df2", "env_ipca_mr1", "env_mean_mr1", file = save_file)


## Use PCA to visualize the environments
ec_mats_pca <- ec_mats %>%
  mutate(pca = map(mat, ~prcomp(.) %>% tidy()))



## Plot
g_ec_pca <- ec_mats_pca %>%
  split(.$ec_group) %>%
  map(~{
    unnest(., pca) %>%
      filter(PC %in% c(1,2)) %>%
      mutate(year = paste0("20", str_extract(row, "[0-9]{2}")),
             PC = paste0("PC", PC)) %>%
      spread(PC, value) %>%
      ## Plot
      ggplot(aes(x = PC1, y = PC2, color = year)) +
      geom_point() +
      facet_wrap(trait ~ group, scales = "free", ncol = 3) +
      theme_acs()
  })

for(i in seq_along(g_ec_pca)) {
  filename <- str_c(names(g_ec_pca)[i], "_ec_mat_pca.jpg")
  ggsave(filename = filename, plot = g_ec_pca[[i]], path = fig_dir, width = 5, height = 5, dpi = 1000)
}



## PCA of the similarity matrix
## Use PCA to visualize the environments
sim_mat_pca <- sim_mat_df %>%
  mutate(pca = map(sim_mat, ~prcomp(.) %>% tidy()))
  # mutate(pca = map(sim_mat, ~cmdscale(.) %>% tidy() %>% rename(row = .rownames, x = X1, y = X2)))



## PCA of the similarity matrix
## Use PCA to visualize the environments
sim_mat_pca <- mutate(sim_mat_df2, pca = map(sim_mat, ~prcomp(.) %>% tidy()))



## Plot
g_sim_pca <- sim_mat_pca %>%
  filter(population == "tp") %>%
  split(.$ec_group) %>%
  map(~{
    unnest(., pca) %>%
      filter(PC %in% c(1,2)) %>%
      mutate(year = paste0("20", str_extract(row, "[0-9]{2}")),
             PC = paste0("PC", PC)) %>%
      spread(PC, value) %>%
      ## Plot
      ggplot(aes(x = PC1, y = PC2, color = year)) +
      geom_point() +
      facet_wrap(trait ~ group, scales = "free", ncol = 3) +
      theme_acs()
  })




## Using the different similarity matrices, rank the distance from each environment to every other environment. 
## Are the rankings consistent?
sim_mat_df_rank <- sim_mat_df$sim_mat %>% map(~apply(X = ., MARGIN = 1, FUN = order, decreasing = TRUE)[-1,])
sim_mat_df1_rank <- sim_mat_df1$sim_mat %>% map(~apply(X = ., MARGIN = 1, FUN = order, decreasing = TRUE)[-1,])
sim_mat_df2_rank <- sim_mat_df2$sim_mat %>% map(~apply(X = ., MARGIN = 1, FUN = order, decreasing = TRUE)[-1,])

