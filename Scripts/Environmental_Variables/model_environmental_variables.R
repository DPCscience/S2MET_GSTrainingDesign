## Exhaustive search for environmental variables that are correlated with the environmental mean
## 
## This script will first fit a model to get the fixed effect of environments. Then, using information
## on planting times, it will search for environmental variables that are correlated. This script will
## also generate distance matrices to be used in further analysis.
## 
## Author: Jeff Neyhart
## Last modified: July 31, 2018
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


# Significance level
alpha <- 0.05


### Fit a joint regression model


# Fit the model per trait, extract the environmental means
full_model_fits <- S2_MET_BLUEs %>%
  mutate_at(vars(line_name, environment), as.factor) %>%
  group_by(trait) %>%
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
  group_by(trait) %>%
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



    
# 
# 
# ## Calculating using only the tp
# env_means_tp <- S2_MET_BLUEs %>%
#   filter(line_name %in% tp) %>%
#   mutate_at(vars(line_name, environment), as.factor) %>%
#   group_by(trait) %>%
#   do({
#     
#     df <- droplevels(.)
#     
#     print(unique(df$trait))
#     
#     # Control and weights
#     control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
#     wts <- df$std_error^2
#     
#     # Formula
#     form <- value ~ (1|line_name) + environment
#     # Fit the model
#     fit <- lmer(formula = form, data = df, control = control, weights = wts, contrasts = list(environment = "contr.sum"))
#     
#     # Get the fixed coefficients
#     env_eff <- fixef(fit)
#     # Tidy
#     data_frame(environment = levels(df$environment), 
#                h = c(env_eff[-1], -sum(env_eff[-1])))  })
# 
# ## Calculate the remaining environments using just the vp
# ## Calculating using only the tp
# env_means_vp <- S2_MET_BLUEs %>% 
#   group_by(environment) %>%
#   filter(!any(line_name %in% tp)) %>%
#   ungroup() %>%
#   mutate_at(vars(line_name, environment), as.factor) %>%
#   group_by(trait) %>%
#   do({
#     
#     df <- droplevels(.)
#     
#     print(unique(df$trait))
#     
#     # Control and weights
#     control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
#     wts <- df$std_error^2
#     
#     # Formula
#     form <- value ~ (1|line_name) + environment
#     # Fit the model
#     fit <- lmer(formula = form, data = df, control = control, weights = wts, contrasts = list(environment = "contr.sum"))
#     
#     # Get the fixed coefficients
#     env_eff <- fixef(fit)
#     # Tidy
#     data_frame(environment = levels(df$environment), 
#                h = c(env_eff[-1], -sum(env_eff[-1])))  })
#     
# 
# ## Calculate correlations
# env_means_tp %>% 
#   left_join(env_means_all, by = c("trait", "environment")) %>% 
#   group_by(trait) %>% 
#   summarize(env_mean_cor = cor(h.x, h.y))
# 
# env_means_tp %>% 
#   left_join(env_means_all, by = c("trait", "environment")) %>%
#   ggplot(aes(x = h.x, y = h.y)) +
#   geom_point() +
#   facet_wrap(~ trait, scales = "free", ncol = 1) +
#   ylab("Environmental mean (n = 233)") +
#   xlab("Environmental mean (n = 183)")
# 
# 
# env_means_vp %>% 
#   left_join(env_means_all, by = c("trait", "environment")) %>% 
#   group_by(trait) %>% 
#   summarize(env_mean_cor = cor(h.x, h.y))
# 
# env_means_vp %>% 
#   left_join(env_means_all, by = c("trait", "environment")) %>%
#   ggplot(aes(x = h.x, y = h.y)) +
#   geom_point() +
#   facet_wrap(~ trait, scales = "free", ncol = 1) +
#   ylab("Environmental mean (n = 233)") +
#   xlab("Environmental mean (n = 183)")
# 
# 




## Calculate regression coefficients for each genotype, then test the signficance of differences among regression coefficients
gen_b_fit <- S2_MET_BLUEs %>% 
  filter(line_name %in% tp) %>% # Just the TP
  left_join(., env_means_all) %>% 
  group_by(trait) %>%
  do(fwr(formula = value ~ h, data = .))

## As expected, all are significant, but GY is most significant and HD is least
## Since all traits show genotype-specific reactions to the environment mean, genotypes should also
# demonstrate different responses to variables that are correlated with the mean




## Correlate environmental variables with the mean

# For each summarized environmental covariable, find the correlation with the environmental mean

## One-year
one_year_env_df1 <- one_year_env_df %>%
  filter(!str_detect(variable, "day")) %>%
  left_join(env_means_all, .)

env_mean_cor <- one_year_env_df1 %>% 
  group_by(trait, variable) %>% 
  do(test = cor.test(.$h, .$value)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         df = map_dbl(test, "parameter")) %>%
  select(-test)


## What interval variables are highly correlated
one_year_env_interval_df1 <- one_year_daily_summary_interval %>% 
  map(unnest) %>%
  map(select, -variable) %>% 
  reduce(full_join) %>% 
  gather(variable, value, -begin:-environment) %>% 
  left_join(env_means_all, .) %>%
  # Remove NA
  filter(!is.na(value))

env_interval_cor <- one_year_env_interval_df1 %>% 
  group_by(trait, variable, begin, end) %>% 
  do(test = cor.test(.$h, .$value)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         df = map_dbl(test, "parameter")) %>%
  select(-test)

## Plot
env_interval_cor_plot_list <- env_interval_cor %>% 
  split(list(.$trait, .$variable)) %>% 
  map(~qplot(x = begin, y = end, fill = cor, geom = "tile", data = .) +
        scale_fill_gradient2() + facet_wrap(~ trait + variable) +
        theme_acs() + theme(legend.position = c(0.85, 0.35)))

# Cowplot
env_interval_cor_plot <- plot_grid(plotlist = env_interval_cor_plot_list, ncol = 3)
ggsave(filename = "one_year_interval_variable_correlation_space.jpg", plot = env_interval_cor_plot,
       path = fig_dir, width = 10, height = 12, dpi = 1000)


## Find the intervals most correlated with the environmental mean
## Weigh the expected direction of the correlation appropriately
env_interval_cor_best <- env_interval_cor %>% 
  mutate(wt_cor = ifelse(trait == "HeadingDate", cor * -1, cor)) %>%
  # Filter for appropriate correlations
  filter(sign(wt_cor) == 1) %>%
  group_by(trait, variable) %>%
  # Or filter by the top correlations
  # top_n(n = 3, wt = wt_cor)
  # Only use AGDD and GDD
  filter(variable %in% c("AGDD", "GDD"), pvalue <= 0.05) %>%
  top_n(n = 10, wt = wt_cor) # Take only the top 10

  

## Visualize the correlations
env_interval_cor_best_plot <- env_interval_cor_best %>% 
  group_by(trait, variable) %>% 
  top_n(n = 1, wt = wt_cor) %>% 
  left_join(one_year_env_interval_df1) %>%
  # Plot
  qplot(x = value, y = h, data = .) + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ variable + trait + begin + end, scales = "free", ncol = 2,
             labeller = labeller(.multi_line = FALSE)) +
  theme_acs()
  
## Save
ggsave(filename = "one_year_interval_variable_top_correlation.jpg", plot = env_interval_cor_best_plot,
       path = fig_dir, width = 5, height = 5, dpi = 1000)


## Create extra variables to use
env_interval_use <- env_interval_cor_best %>% 
  left_join(one_year_env_interval_df1) %>%
  unite(variable, variable, begin, end, sep = "_")




## What summary variables are highly correlated?
env_mean_summary_cor_top <- env_mean_cor %>% 
  group_by(trait) %>% 
  # do(top_n(x = ., n = n, wt = abs(.$cor))) %>%
  # filter(abs(cor) >= 0.2) %>%
  filter(pvalue <= 0.10) %>%
  arrange(trait, desc(abs(cor))) %>%
  # What variables have positive versus negative correlations?
  mutate(sign = sign(cor))


## Combine the summary variables with the interval variables
env_mean_cor_top <- bind_rows(env_interval_use, left_join(env_mean_summary_cor_top, one_year_env_df1)) %>% 
  select(-df, -sign, -wt_cor) %>%
  ungroup()


# Create a model matrix
env_mean_tomodel <- env_mean_cor_top %>%
  select(-cor, -pvalue) %>%
  split(.$trait) %>%
  map(~mutate(., variable = factor(variable, levels = unique(.$variable)))) %>%
  map(~spread(., variable, value))




# Fit a linear model
env_mean_fit <- env_mean_tomodel %>%
  map_df(~{
    df <- .
    df1 <- select(df, -trait, -environment)
    
    fit <- lm(h ~ ., df1)
    # Reduce by backwards elimination
    fit_step <- step(object = fit, direction = "backward")
    
    distinct(df, trait) %>%
      mutate(fit_orig = list(fit), fit_red = list(fit_step))
    
  })


# Extract the r-squared and pvalues
env_mean_fit_summ <- env_mean_fit %>% 
  gather(type, mod, -trait) %>%
  mutate(summary = map(mod, ~summary(.)[c("r.squared", "adj.r.squared")] %>% as.data.frame()),
         terms = map_dbl(mod, ~terms(.) %>% attr("term.labels") %>% length),
         df_res = map_dbl(mod, df.residual)) %>%
  unnest(summary)


# Pull out the datasets of the final models
env_mean_fit_mf <- env_mean_fit %>%
  mutate(red_data = map(fit_red, model.frame)) %>%
  select(-contains("fit")) %>%
  split(.$trait) %>%
  map(unnest) %>%
  map2(.y = ., .x = map(env_mean_tomodel, ~select(., trait:h)), .f = left_join)
  

## Combine this with model matrices that have been reduced by MR
pc_data_tomodel <- map(env_mean_fit_mf, ~select(., -trait, -h) %>% as.data.frame() %>% column_to_rownames("environment")) %>% 
  data_frame(trait = names(.), model = "mr_red", data = .)

# Add to this the data from non- mr reduced variables
pc_data_tomodel1 <- map(env_mean_tomodel, ~select(., -trait, -h) %>% as.data.frame() %>% column_to_rownames("environment")) %>% 
  data_frame(trait = names(.), model = "top_cor", data = .) %>%
  bind_rows(., pc_data_tomodel)

# Add to this data from all environmental covariables
pc_data_tomodel2 <- one_year_env_interval_df1 %>%
  select(-h) %>% 
  filter(variable %in% c("AGDD", "GDD")) %>% 
  unite(variable, variable, begin, end, sep = "_") %>% 
  spread(variable, value) %>% 
  # Combine with the summary covariables
  full_join(., spread(one_year_env_df, variable, value)) %>% 
  split(.$trait) %>% 
  map(~select(., -trait) %>% as.data.frame() %>% column_to_rownames("environment")) %>% 
  data_frame(trait = names(.), model = "all_ec", data = .) %>%
  bind_rows(., pc_data_tomodel1)



## Try principal component regression of the final variables
# Use permutation to determine what PCs to keep
env_mean_pc_fit <- pc_data_tomodel2  %>%
  mutate(pca_summ = map(data, ~prcomp_perm(x = ., permutations = 1000)))
      

## Determine the PCs to keep
env_mean_pc_fit_tokeep <- env_mean_pc_fit %>%
  mutate(n_final_PC = map(pca_summ, "summary") %>% map_dbl(~subset(., p_value <= alpha, PC, drop = TRUE) %>% parse_number() %>% max),
         tidy = map(pca_summ, "tidy")) %>%
  unnest(tidy) %>%
  group_by(trait, model) %>%
  filter(PC <= n_final_PC) %>%
  rename(environment = row)
  
env_mean_pc_fit_tokeep %>% summarize(n_PC = n_distinct(PC))


# Correlation with mean
# Add the significant PCs together to form the index
env_mean_pc_fit_index <- env_mean_pc_fit_tokeep %>% 
  nest(environment:value) %>% 
  mutate(data = map(data, ~spread(., PC, value) %>% mutate(index = rowSums(select(., -environment))) %>% 
                      select(environment, index)),
         variable = "PC") %>%
  unnest()


## Select the one EC that is most correlated with the mean. This could be an interval EC or a summary EC
env_mean_cor_top_one <- env_mean_cor_top %>% 
  group_by(trait) %>% 
  top_n(n = 1, wt = abs(cor)) %>%
  mutate(model = "top_one_cor") %>%
  select(trait, model, variable, environment, index = value)


## Combine with the PCs
env_mean_ec_index <- bind_rows(env_mean_pc_fit_index, env_mean_cor_top_one)

# Df to plot
env_mean_ec_index_toplot <- env_mean_ec_index %>%
  left_join(., env_means_all)

env_mean_ec_index_toplot %>% 
  group_by(trait, model) %>% 
  do(test = cor.test(.$index, .$h)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"))

# trait       model       test            cor   pvalue
# 1 GrainYield  all_ec      <S3: htest>  0.0450 0.801    
# 2 GrainYield  mr_red      <S3: htest>  0.649  0.0000331
# 3 GrainYield  top_cor     <S3: htest>  0.295  0.0900   
# 4 GrainYield  top_one_cor <S3: htest>  0.594  0.000216 
# 5 HeadingDate all_ec      <S3: htest> -0.0773 0.679    
# 6 HeadingDate mr_red      <S3: htest> -0.597  0.000393 
# 7 HeadingDate top_cor     <S3: htest> -0.500  0.00419  
# 8 HeadingDate top_one_cor <S3: htest> -0.687  0.0000195
# 9 PlantHeight all_ec      <S3: htest> -0.155  0.396    
# 10 PlantHeight mr_red      <S3: htest> -0.397  0.0243   
# 11 PlantHeight top_cor     <S3: htest>  0.232  0.201    
# 12 PlantHeight top_one_cor <S3: htest>  0.397  0.0243  



## Test genotypes for their reaction to these indices
gen_b_ec_fit <- S2_MET_BLUEs %>% 
  filter(line_name %in% tp) %>% # Just the TP
  left_join(., env_mean_ec_index_toplot) %>% 
  group_by(trait, model) %>%
  do(fwr(formula = value ~ index, data = .))

## Not surprisingly, the top correlated variable is the most significant, but not the PC indices
# trait       model       regression          Fstat   pvalue
# 1 GrainYield  all_ec      <tibble [183 x 4]>  0.106 1     
# 2 GrainYield  mr_red      <tibble [183 x 4]>  0.255 1     
# 3 GrainYield  top_cor     <tibble [183 x 4]>  0.158 1     
# 4 GrainYield  top_one_cor <tibble [183 x 4]> 14.4   0     
# 5 HeadingDate all_ec      <tibble [183 x 4]>  0.140 1     
# 6 HeadingDate mr_red      <tibble [183 x 4]>  0.325 1     
# 7 HeadingDate top_cor     <tibble [183 x 4]>  0.307 1     
# 8 HeadingDate top_one_cor <tibble [183 x 4]> 48.6   0     
# 9 PlantHeight all_ec      <tibble [183 x 4]>  0.150 1     
# 10 PlantHeight mr_red      <tibble [183 x 4]>  0.236 1     
# 11 PlantHeight top_cor     <tibble [183 x 4]>  0.167 1     
# 12 PlantHeight top_one_cor <tibble [183 x 4]>  1.22  0.0247



# Plot
env_mean_index_plot <- env_mean_ec_index_toplot %>%
  ggplot(aes(x = index, y = h, label = environment)) +
  geom_point() + 
  facet_wrap(~ trait + model + variable, nrow = 3, scales = "free") +
  geom_smooth(method = "lm", se = FALSE) + 
  theme_acs()

## Filter expression
filter_exp <- expression((trait == "GrainYield" & h > 2000) | (trait == "PlantHeight" & h > 20))

env_mean_index_plot_label <- env_mean_index_plot + 
  geom_text(data = subset(env_mean_ec_index_toplot, eval(filter_exp)), 
            size = 1.5, check_overlap = T, hjust = 1)

ggsave(filename = "one_year_env_variable_index.jpg", plot = env_mean_index_plot_label, path = fig_dir,
       height = 6, width = 8, dpi = 1000)



## Create distance matrices
env_mean_pc_fit_dist <- env_mean_pc_fit_tokeep %>% 
  nest(environment:value) %>%
  mutate(dist = map(data, ~spread(., PC, value) %>% as.data.frame() %>% column_to_rownames("environment") %>% 
                      as.matrix() %>% dist() ))

env_mean_cor_top_one_dist <- env_mean_cor_top_one %>%
  group_by(trait, model, variable) %>%
  nest(environment, index) %>%
  mutate(dist = map(data, ~as.data.frame(.) %>% column_to_rownames("environment") %>% 
                      as.matrix() %>% dist() ))

## Combine
one_year_ec_dist <- bind_rows(env_mean_pc_fit_dist, env_mean_cor_top_one_dist) %>% 
  select(trait, model, variable, dist)


# ### Test the reaction of genotypes to all summary variables and interval variables
# # First break-up the interval data to prevent a gigantic df from forming
# one_year_env_interval_df2 <- one_year_env_interval_df1 %>%
#   rename(index = value) %>%
#   group_by(trait, variable) %>%
#   nest(environment, begin:index)
# 
# pheno_use <- S2_MET_BLUEs %>% 
#   filter(line_name %in% tp)
# 
# for (i in seq(nrow(one_year_env_interval_df2))) {
#   df <- one_year_env_interval_df2$data[[i]] %>%
#     left_join(pheno_use, ., by = "environment")
#   
#   gen_b_fit_out <- df %>%
#     group_by(trait, begin, end) %>%
#     do(fwr(formula = value ~ index, data = .))
#   
# }












#### Multi-year ####

multi_year_env_df1 <- multi_year_env_df %>%
  filter(!str_detect(variable, "day")) %>%
  left_join(env_means_all, .)


env_mean_cor <- multi_year_env_df1 %>% 
  group_by(trait, variable) %>% 
  do(test = cor.test(.$h, .$value)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         df = map_dbl(test, "parameter")) %>%
  select(-test)


## What interval variables are highly correlated
multi_year_env_interval_df1 <- multi_year_daily_summary_interval %>% 
  map(unnest) %>%
  map(select, -variable) %>% 
  reduce(full_join) %>% 
  gather(variable, value, -begin:-environment) %>% 
  left_join(env_means_all, .) %>%
  # Remove NA
  filter(!is.na(value))

env_interval_cor <- multi_year_env_interval_df1 %>% 
  group_by(trait, variable, begin, end) %>% 
  do(test = cor.test(.$h, .$value)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         df = map_dbl(test, "parameter")) %>%
  select(-test)

## Plot
env_interval_cor_plot_list <- env_interval_cor %>% 
  split(list(.$trait, .$variable)) %>% 
  map(~qplot(x = begin, y = end, fill = cor, geom = "tile", data = .) +
        scale_fill_gradient2() + facet_wrap(~ trait + variable) +
        theme_acs() + theme(legend.position = c(0.85, 0.35)))

# Cowplot
env_interval_cor_plot <- plot_grid(plotlist = env_interval_cor_plot_list, ncol = 3)
ggsave(filename = "multi_year_interval_variable_correlation_space.jpg", plot = env_interval_cor_plot,
       path = fig_dir, width = 10, height = 12, dpi = 1000)


## Find the intervals most correlated with the environmental mean
## Weigh the expected direction of the correlation appropriately
env_interval_cor_best <- env_interval_cor %>% 
  mutate(wt_cor = ifelse(trait == "HeadingDate", cor * -1, cor)) %>%
  # Filter for appropriate correlations
  filter(sign(wt_cor) == 1) %>%
  group_by(trait, variable) %>%
  # filter(variable %in% c("AGDD", "GDD"), wt_cor >= quantile(x = wt_cor, probs = 0.995))
  # Only use AGDD and GDD
  filter(variable %in% c("AGDD", "GDD"), pvalue <= 0.05) %>%
  top_n(n = 10, wt = wt_cor) # Take only the top 10


## Visualize the correlations
env_interval_cor_best_plot <- env_interval_cor_best %>% 
  group_by(trait, variable) %>% 
  top_n(n = 1, wt = wt_cor) %>% 
  left_join(multi_year_env_interval_df1) %>%
  # Plot
  qplot(x = value, y = h, data = .) + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ variable + trait + begin + end, scales = "free", ncol = 3,
             labeller = labeller(.multi_line = FALSE)) +
  theme_acs()

## Save
ggsave(filename = "multi_year_interval_variable_top_correlation.jpg", plot = env_interval_cor_best_plot,
       path = fig_dir, width = 8, height = 3, dpi = 1000)


## Create extra variables to use
env_interval_use <- env_interval_cor_best %>% 
  left_join(multi_year_env_interval_df1) %>%
  unite(variable, variable, begin, end, sep = "_")




## What summary variables are highly correlated?
env_mean_summary_cor_top <- env_mean_cor %>% 
  group_by(trait) %>% 
  # do(top_n(x = ., n = n, wt = abs(.$cor))) %>%
  # filter(abs(cor) >= 0.2) %>%
  filter(pvalue <= 0.10) %>%
  arrange(trait, desc(abs(cor))) %>%
  # What variables have positive versus negative correlations?
  mutate(sign = sign(cor))


## Combine the summary variables with the interval variables
env_mean_cor_top <- bind_rows(env_interval_use, left_join(env_mean_summary_cor_top, multi_year_env_df1)) %>% 
  select(-df, -sign, -wt_cor) %>%
  ungroup()


# Create a model matrix
env_mean_tomodel <- env_mean_cor_top %>%
  select(-cor, -pvalue) %>%
  split(.$trait) %>%
  map(~mutate(., variable = factor(variable, levels = unique(.$variable)))) %>%
  map(~spread(., variable, value))




# Fit a linear model
env_mean_fit <- env_mean_tomodel %>%
  map_df(~{
    df <- .
    df1 <- select(df, -trait, -environment)
    
    fit <- lm(h ~ ., df1)
    # Reduce by backwards elimination
    fit_step <- step(object = fit, direction = "backward")
    
    distinct(df, trait) %>%
      mutate(fit_orig = list(fit), fit_red = list(fit_step))
    
  })


# Extract the r-squared and pvalues
env_mean_fit_summ <- env_mean_fit %>% 
  gather(type, mod, -trait) %>%
  mutate(summary = map(mod, ~summary(.)[c("r.squared", "adj.r.squared")] %>% as.data.frame()),
         terms = map_dbl(mod, ~terms(.) %>% attr("term.labels") %>% length),
         df_res = map_dbl(mod, df.residual)) %>%
  unnest(summary)

## Note for heading date only one variable remains. PC will be complicated.

# Pull out the datasets of the final models
env_mean_fit_mf <- env_mean_fit %>%
  mutate(red_data = map(fit_red, model.frame)) %>%
  select(-contains("fit")) %>%
  split(.$trait) %>%
  map(unnest) %>%
  map2(.y = ., .x = map(env_mean_tomodel, ~select(., trait:h)), .f = left_join)

## 


## Combine this with model matrices that have been reduced by MR
pc_data_tomodel <- map(env_mean_fit_mf, ~select(., -trait, -h) %>% as.data.frame() %>% column_to_rownames("environment")) %>% 
  data_frame(trait = names(.), model = "mr_red", data = .)

# Add to this the data from non- mr reduced variables
pc_data_tomodel1 <- map(env_mean_tomodel, ~select(., -trait, -h) %>% as.data.frame() %>% column_to_rownames("environment")) %>% 
  data_frame(trait = names(.), model = "top_cor", data = .) %>%
  bind_rows(., pc_data_tomodel)

# Add to this data from all environmental covariables
pc_data_tomodel2 <- multi_year_env_interval_df1 %>%
  select(-h) %>% 
  filter(variable %in% c("AGDD", "GDD")) %>% 
  unite(variable, variable, begin, end, sep = "_") %>% 
  spread(variable, value) %>% 
  # Combine with the summary covariables
  full_join(., spread(one_year_env_df, variable, value)) %>% 
  split(.$trait) %>% 
  map(~select(., -trait) %>% as.data.frame() %>% column_to_rownames("environment")) %>% 
  data_frame(trait = names(.), model = "all_ec", data = .) %>%
  bind_rows(., pc_data_tomodel1)



## Try principal component regression of the final variables
# Use permutation to determine what PCs to keep
env_mean_pc_fit <- pc_data_tomodel2 %>%
  mutate(pca_summ = map(data, ~{
    if (ncol(.) < 2) {
      list(tidy = NA, summary = NA)
    } else {
      prcomp_perm(x = ., permutations = 1000)
    } }))


## Determine the PCs to keep

env_mean_pc_fit_tokeep <- env_mean_pc_fit %>%
  mutate(keep = map_lgl(pca_summ, ~all(!is.na(.)))) %>%
  filter(keep) %>%
  mutate(n_final_PC = map(pca_summ, "summary") %>% map_dbl(~subset(., p_value <= alpha, PC, drop = TRUE) %>% parse_number() %>% max),
         n_final_PC = ifelse(is.infinite(n_final_PC), 1, n_final_PC),
         tidy = map(pca_summ, "tidy")) %>%
  unnest(tidy) %>%
  group_by(trait, model) %>%
  filter(PC <= n_final_PC) %>%
  rename(environment = row)



# Correlation with mean
# Add the significant PCs together to form the index
env_mean_pc_fit_index <- env_mean_pc_fit_tokeep %>% 
  nest(environment:value) %>% 
  mutate(data = map(data, ~spread(., PC, value) %>% mutate(index = rowSums(select(., -environment))) %>% 
                      select(environment, index)),
         variable = "PC") %>%
  unnest()


## Select the one EC that is most correlated with the mean. This could be an interval EC or a summary EC
env_mean_cor_top_one <- env_mean_cor_top %>% 
  group_by(trait) %>% 
  top_n(n = 1, wt = abs(cor)) %>%
  mutate(model = "top_one_cor") %>%
  select(trait, model, variable, environment, index = value)


## Combine with the PCs
env_mean_ec_index <- bind_rows(env_mean_pc_fit_index, env_mean_cor_top_one)

# Df to plot
env_mean_ec_index_toplot <- env_mean_ec_index %>%
  left_join(., env_means_all)

env_mean_ec_index_toplot %>% 
  group_by(trait, model) %>% 
  do(test = cor.test(.$index, .$h)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"))

# trait       model       test            cor    pvalue
# 1 GrainYield  all_ec      <S3: htest>  0.101 0.569    
# 2 GrainYield  mr_red      <S3: htest> -0.587 0.000265 
# 3 GrainYield  top_cor     <S3: htest> -0.302 0.0828   
# 4 GrainYield  top_one_cor <S3: htest> -0.663 0.0000191
# 5 HeadingDate all_ec      <S3: htest>  0.633 0.000133 
# 6 HeadingDate top_cor     <S3: htest>  0.602 0.000342 
# 7 HeadingDate top_one_cor <S3: htest> -0.679 0.0000272
# 8 PlantHeight all_ec      <S3: htest>  0.188 0.302    
# 9 PlantHeight mr_red      <S3: htest> -0.434 0.0131   
# 10 PlantHeight top_cor     <S3: htest>  0.480 0.00546  
# 11 PlantHeight top_one_cor <S3: htest>  0.398 0.0242  


## Test genotypes for their reaction to these indices
gen_b_ec_fit <- S2_MET_BLUEs %>% 
  filter(line_name %in% tp) %>% # Just the TP
  left_join(., env_mean_ec_index_toplot) %>% 
  group_by(trait, model) %>%
  do(fwr(formula = value ~ index, data = .))

## Not surprisingly, the top correlated variable is the most significant, but not the PC indices
gen_b_ec_fit

# trait       model       regression          Fstat   pvalue
# 1 GrainYield  all_ec      <tibble [183 x 4]>  0.106   1.00e+ 0
# 2 GrainYield  mr_red      <tibble [183 x 4]>  0.228   1.00e+ 0
# 3 GrainYield  top_cor     <tibble [183 x 4]>  0.0981  1.00e+ 0
# 4 GrainYield  top_one_cor <tibble [183 x 4]>  1.91    6.02e-12
# 5 HeadingDate all_ec      <tibble [183 x 4]>  0.376  10.00e- 1
# 6 HeadingDate top_cor     <tibble [183 x 4]>  0.292   1.00e+ 0
# 7 HeadingDate top_one_cor <tibble [183 x 4]> 37.4     0.      
# 8 PlantHeight all_ec      <tibble [183 x 4]>  0.165   1.00e+ 0
# 9 PlantHeight mr_red      <tibble [183 x 4]>  0.255   1.00e+ 0
# 10 PlantHeight top_cor     <tibble [183 x 4]>  0.179   1.00e+ 0
# 11 PlantHeight top_one_cor <tibble [183 x 4]> 16.2     0.



env_mean_ec_index_toplot %>% 
  group_by(trait, model) %>% 
  do(test = cor.test(.$index, .$h)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"))

# Plot
env_mean_index_plot <- env_mean_ec_index_toplot %>%
  ggplot(aes(x = index, y = h, label = environment)) +
  geom_point() + 
  facet_wrap(~ trait + model + variable, nrow = 3, scales = "free") +
  geom_smooth(method = "lm", se = FALSE) + 
  theme_acs()

## Filter expression
filter_exp <- expression((trait == "GrainYield" & h > 2000) | (trait == "PlantHeight" & h > 20))

env_mean_index_plot_label <- env_mean_index_plot + 
  geom_text(data = subset(env_mean_ec_index_toplot, eval(filter_exp)), 
            size = 1.5, check_overlap = T, hjust = 1)

ggsave(filename = "multi_year_env_variable_index.jpg", plot = env_mean_index_plot_label, path = fig_dir,
       height = 6, width = 8, dpi = 1000)



## Create distance matrices
env_mean_pc_fit_dist <- env_mean_pc_fit_tokeep %>% 
  nest(environment:value) %>%
  mutate(dist = map(data, ~spread(., PC, value) %>% as.data.frame() %>% column_to_rownames("environment") %>% 
                      as.matrix() %>% dist() ))

env_mean_cor_top_one_dist <- env_mean_cor_top_one %>%
  group_by(trait, model, variable) %>%
  nest(environment, index) %>%
  mutate(dist = map(data, ~as.data.frame(.) %>% column_to_rownames("environment") %>% 
                      as.matrix() %>% dist() ))

## Combine
multi_year_ec_dist <- bind_rows(env_mean_pc_fit_dist, env_mean_cor_top_one_dist) %>% 
  select(trait, model, variable, dist)


# 
# 
# ### Test the reaction of genotypes to all summary variables and interval variables
# # First break-up the interval data to prevent a gigantic df from forming
# multi_year_env_interval_df2 <- multi_year_env_interval_df1 %>%
#   rename(index = value) %>%
#   group_by(trait, variable) %>%
#   nest(environment, begin:index)
# 
# pheno_use <- S2_MET_BLUEs %>% 
#   filter(line_name %in% tp)
# 
# for (i in seq(nrow(multi_year_env_interval_df2))) {
#   df <- multi_year_env_interval_df2$data[[i]] %>%
#     left_join(pheno_use, ., by = "environment")
#   
#   gen_b_fit_out <- df %>%
#     group_by(trait, begin, end) %>%
#     do(fwr(formula = value ~ index, data = .))
#   
# }




## Save the distance matrices
save_file <- file.path(result_dir, "environmental_covariable_distance_mat.RData")
save("one_year_ec_dist", "multi_year_ec_dist", file = save_file)




## Use multi-dimensional scaling to visualize the environments
one_year_ec_mds <- one_year_ec_dist %>% 
  mutate(mds = map(dist, ~cmdscale(.) %>% as.data.frame() %>% `names<-`(., c("x", "y")) %>% rownames_to_column("environment"))) %>% 
  unnest(mds)

# Plot
one_year_ec_mds_plot <- one_year_ec_mds %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_text(aes(label = environment), size = 2) + 
  # geom_point() + 
  # geom_text_repel(aes(label = environment), size = 2, ) + 
  facet_wrap(trait ~ model, scales = "free") + 
  theme_acs()

ggsave(filename = "one_year_env_ec_mds.jpg", plot = one_year_ec_mds_plot, path = fig_dir,
       height = 6, width = 8, dpi = 1000)




multi_year_ec_mds <- multi_year_ec_dist %>% 
  mutate(mds = map(dist, ~cmdscale(.) %>% as.data.frame() %>% `names<-`(., c("x", "y")) %>% rownames_to_column("environment"))) %>% 
  unnest(mds)

# Plot
multi_year_ec_mds_plot <- multi_year_ec_mds %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_text(aes(label = environment), size = 2) + 
  # geom_point() + 
  # geom_text_repel(aes(label = environment), size = 2, ) + 
  facet_wrap(trait ~ model, scales = "free") + 
  theme_acs()

ggsave(filename = "multi_year_env_ec_mds.jpg", plot = multi_year_ec_mds_plot, path = fig_dir,
       height = 6, width = 8, dpi = 1000)







