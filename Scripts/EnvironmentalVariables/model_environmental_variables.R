## Exhaustive search for environmental variables that are correlated with the environmental mean
## 
## This script will first fit a model to get the fixed effect of environments. Then, using information
## on planting times, it will search for environmental variables that are correlated. This script will
## also generate distance matrices to be used in further analysis.
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


# Significance level
alpha <- 0.05


### Fit a joint regression model
# First create a list of data split by whether all lines are used or just the TP
data_to_model <- list(tp = tp, all = c(tp, vp)) %>% 
  map(~filter(S2_MET_BLUEs, line_name %in% .)) %>%
  list(., names(.)) %>%
  pmap_df(~mutate(.x, population = .y))


# Fit the model per trait, extract the environmental means
full_model_fits <- data_to_model %>%
  mutate_at(., vars(line_name, environment), as.factor) %>%
  group_by(trait, population) %>%
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
  group_by(trait, population) %>%
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
S2_MET_BLUEs %>% 
  left_join(., env_means_all) %>% 
  group_by(trait, population) %>%
  do(fwr(formula = value ~ h, data = .))

# trait           regression         Fstat   pvalue
# 1 GrainYield      all        <tibble [233 x 4]>  3.95 7.41e-76
# 2 GrainYield      tp         <tibble [233 x 4]>  2.91 1.43e-41
# 3 HeadingDate     all        <tibble [233 x 4]>  1.61 2.02e- 8
# 4 HeadingDate     tp         <tibble [233 x 4]>  1.28 3.22e- 3
# 5 HeadingDateAGDD all        <tibble [233 x 4]>  2.44 2.98e-28
# 6 HeadingDateAGDD tp         <tibble [233 x 4]>  1.73 1.25e-10
# 7 PlantHeight     all        <tibble [233 x 4]>  1.65 3.37e- 9
# 8 PlantHeight     tp         <tibble [233 x 4]>  1.50 2.35e- 6

## As expected, all are significant, but GY is most significant and HD is least
## Since all traits show genotype-specific reactions to the environment mean, genotypes should also
# demonstrate different responses to variables that are correlated with the mean




## Correlate environmental variables with the mean

# For each summarized environmental covariable, find the correlation with the environmental mean

## One-year
one_year_env_df1 <- one_year_env_df %>%
  filter(!str_detect(variable, "day")) %>%
  left_join(env_means_all, .) %>%
  mutate(EC_type = "summary")

env_mean_cor <- one_year_env_df1 %>% 
  group_by(population, trait, variable) %>% 
  do(test = cor.test(.$h, .$value)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         df = map_dbl(test, "parameter")) %>%
  select(-test)

env_mean_cor %>% group_by(population, trait) %>% top_n(n = 2, wt = abs(cor))

# trait           variable               cor   pvalue    df
# 1 all        GrainYield      month_4_PPT         -0.546 0.000843    32
# 2 all        GrainYield      ph1to1h2o_r_topsoil  0.594 0.000216    32
# 3 all        HeadingDate     month_4_TMAX         0.344 0.0580      29
# 4 all        HeadingDate     temp_seasonality    -0.356 0.0491      29
# 5 all        HeadingDateAGDD om_r_topsoil         0.206 0.266       29
# 6 all        HeadingDateAGDD ph1to1h2o_r_topsoil -0.229 0.214       29
# 7 all        PlantHeight     month_6_PPT          0.397 0.0246      30
# 8 all        PlantHeight     ph1to1h2o_r_topsoil  0.391 0.0269      30
# 9 tp         GrainYield      month_4_PPT         -0.560 0.00158     27
# 10 tp         GrainYield      ph1to1h2o_r_topsoil  0.594 0.000680    27
# 11 tp         HeadingDate     month_4_TMAX         0.429 0.0228      26
# 12 tp         HeadingDate     temp_seasonality    -0.375 0.0492      26
# 13 tp         HeadingDateAGDD month_4_PPT          0.416 0.0275      26
# 14 tp         HeadingDateAGDD ph1to1h2o_r_topsoil -0.303 0.116       26
# 15 tp         PlantHeight     isothermality       -0.477 0.00883     27
# 16 tp         PlantHeight     month_6_PPT          0.432 0.0191      27

# Plot
env_mean_cor %>%
  filter(population == "all") %>%
  group_by(trait) %>%
  top_n(n = 5, wt = abs(cor)) %>%
  select(population:variable) %>%
  left_join(., one_year_env_df1) %>%
  ggplot(aes(x = value, y = h)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  geom_point() + 
  facet_wrap( ~ trait + variable, scale = "free", ncol = 5) + 
  ylab("Environmental mean") +
  xlab("Covariate") +
  theme_acs()



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



## Combine the summary variables with the interval variables
env_combine <- one_year_env_df1 %>%
  # bind_rows(env_interval_use) %>%
  select(population, trait, environment, h, variable, value)



## Generate random samples of environments. For each sample, calculate the same
## correlation as above. How often is the same EC strongly correlated?
n_sample <- 1000
f_sample <- 0.6

set.seed(500)
one_year_env_samples <- env_combine %>%
  group_by(population, trait) %>%
  do({
    df <- .
    # Generate samples of environments
    cor_sample <- replicate(n_sample, sample_frac(tbl = distinct(df, environment), size = f_sample), simplify = FALSE) %>% 
      map(~left_join(., df, by = "environment")) %>% 
      map_df(~group_by(., variable) %>% summarize(cor = cor(h, value)) %>% mutate(rank = rank(-abs(cor))))
    
    # Return the average rank and average correlation of the variables
    cor_sample %>% 
      group_by(variable) %>% 
      summarise_at(vars(cor, rank), mean) %>% 
      arrange(rank)
    
  })

# Number of covariables to select
n_ec <- 5

top_n(one_year_env_samples, n_ec, -rank) %>% as.data.frame()


## Just summary variables
# population           trait             variable        cor   rank
# 1         all      GrainYield  ph1to1h2o_r_topsoil  0.5718817  2.423
# 2         all      GrainYield          month_4_PPT -0.5422179  2.689
# 3         all      GrainYield              min_PPT -0.5135731  3.115
# 4         all      GrainYield  ph1to1h2o_r_subsoil  0.3844677  6.994
# 5         all      GrainYield annual_precipitation -0.2726273 12.972
# 6         all     HeadingDate     temp_seasonality -0.3531687  7.442
# 7         all     HeadingDate         month_4_TMAX  0.3430303  9.213
# 8         all     HeadingDate        isothermality  0.2905725 12.315
# 9         all     HeadingDate         month_6_TMAX  0.2705173 13.744
# 10        all     HeadingDate          month_8_PPT -0.2707158 14.382
# 11        all HeadingDateAGDD  ph1to1h2o_r_topsoil -0.2310684 12.543
# 12        all HeadingDateAGDD         month_7_TAVG -0.1583984 16.126
# 13        all HeadingDateAGDD         om_r_topsoil  0.2212285 16.492
# 14        all HeadingDateAGDD          month_4_PPT  0.1963070 16.548
# 15        all HeadingDateAGDD         month_7_TMIN -0.1654573 17.342
# 16        all     PlantHeight          month_6_PPT  0.4017050  5.223
# 17        all     PlantHeight  ph1to1h2o_r_topsoil  0.3734711  8.667
# 18        all     PlantHeight              max_PPT  0.3490671  9.962
# 19        all     PlantHeight         month_6_TMIN  0.2993061 10.882
# 20        all     PlantHeight        isothermality -0.2960932 11.391
# 21         tp      GrainYield          month_4_PPT -0.5553062  2.536
# 22         tp      GrainYield  ph1to1h2o_r_topsoil  0.5669364  3.070
# 23         tp      GrainYield  ph1to1h2o_r_subsoil  0.4597698  5.389
# 24         tp      GrainYield              min_PPT -0.4194168  6.388
# 25         tp      GrainYield          month_3_PPT -0.3332770 10.931
# 26         tp     HeadingDate         month_4_TMAX  0.4236365  6.560
# 27         tp     HeadingDate     temp_seasonality -0.3681706  8.604
# 28         tp     HeadingDate          month_8_PPT -0.3385865 11.122
# 29         tp     HeadingDate         month_4_TAVG  0.3440702 11.236
# 30         tp     HeadingDate         month_5_TMAX -0.3057176 13.161
# 31         tp HeadingDateAGDD          month_4_PPT  0.4139577  4.931
# 32         tp HeadingDateAGDD             max_TMAX -0.2729342 11.789
# 33         tp HeadingDateAGDD  ph1to1h2o_r_topsoil -0.3018567 12.439
# 34         tp HeadingDateAGDD         month_7_TAVG -0.2708759 12.866
# 35         tp HeadingDateAGDD         month_7_TMAX -0.2366124 13.464
# 36         tp     PlantHeight        isothermality -0.4756837  4.085
# 37         tp     PlantHeight          month_6_PPT  0.4314312  7.567
# 38         tp     PlantHeight annual_diurnal_range -0.3852985  9.617
# 39         tp     PlantHeight         month_3_TMAX -0.3640840 10.528
# 40         tp     PlantHeight         month_8_TMAX -0.3678424 10.935


## Instead of correlation, calculate an Fstat based on the linear reaction of genotypes to the covariable
ec_fwr <- S2_MET_BLUEs %>% 
    left_join(., rename(env_combine, ec_value = value)) %>% 
    group_by(population, trait, variable) %>%
    do(fwr(formula = value ~ ec_value, data = .)) %>%
    select(-regression) %>%
    as.data.frame()

# Adjust the p values and select the significant covariables
ec_fwr_sig <- ec_fwr %>% 
  group_by(population, trait) %>% 
  mutate(p_adj = p.adjust(pvalue, "bonf")) %>% 
  filter(p_adj <= alpha)







## Plot the correlations for the n highest correlated variables
# g_env_mean_variable <- top_n(one_year_env_samples, n_ec, -rank) %>% 
g_env_mean_variable <- top_n(ec_fwr_sig, n_ec, Fstat) %>% 
  filter(population == "all") %>%
  left_join(., env_combine) %>% 
  ggplot(aes(x = value, y = h)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  geom_point() + 
  facet_wrap(trait ~ variable, scale = "free", ncol = n_ec) + 
  ylab("Environmental mean") +
  xlab("Covariate") +
  theme_acs()

ggsave(filename = "one_year_env_variable_mean.jpg", plot = g_env_mean_variable, path = fig_dir, width = 8, height = 6, dpi = 1000)


# g_env_mean_variable <- top_n(one_year_env_samples, n_ec, -rank) %>% 
g_env_mean_variable <- top_n(ec_fwr_sig, n_ec, Fstat) %>% 
  filter(population == "tp") %>%
  left_join(., env_combine) %>% 
  ggplot(aes(x = value, y = h)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  geom_point() + 
  facet_wrap(trait ~ variable, scale = "free", ncol = n_ec) + 
  ylab("Environmental mean") +
  xlab("Covariate") +
  theme_acs()

ggsave(filename = "one_year_env_variable_mean_tp.jpg", plot = g_env_mean_variable, path = fig_dir, width = 8, height = 6, dpi = 1000)







## Examine some outliers
env_combine %>%
  filter((trait == "PlantHeight" & h > 25) | (trait == "GrainYield" & h > 3000) | (trait == "HeadingDateAGDD" & h > 500)) %>%
  distinct(population,trait, environment, h)

# Remove outliers
env_combine_outliers <- env_combine %>%
  filter(!(trait == "PlantHeight" & environment == "CRM17")) %>%
  filter(!(trait == "GrainYield" & h > 3000)) %>%
  filter(!(trait == "HeadingDateAGDD" & h > 500))


set.seed(500)
# Remove outliers and re-adjust
one_year_env_samples_outliers <- env_combine_outliers %>%
  group_by(population, trait) %>%
  do({
    df <- .
    # Generate samples of environments
    cor_sample <- replicate(n_sample, sample_frac(tbl = distinct(df, environment), size = f_sample), simplify = FALSE) %>% 
      map(~left_join(., df, by = "environment")) %>% 
      map_df(~group_by(., variable) %>% summarize(cor = cor(h, value)) %>% mutate(rank = rank(-abs(cor))))
    
    # Return the average rank and average correlation of the variables
    cor_sample %>% 
      group_by(variable) %>% 
      summarise_at(vars(cor, rank), mean) %>% 
      arrange(rank)
    
  })


## Instead of correlation, calculate an Fstat based on the linear reaction of genotypes to the covariable
ec_fwr_outlier <- S2_MET_BLUEs %>% 
  inner_join(., rename(env_combine_outliers, ec_value = value)) %>% 
  group_by(population, trait, variable) %>%
  do(fwr(formula = value ~ ec_value, data = .)) %>%
  select(-regression) %>%
  as.data.frame()

# Adjust the p values and select the significant covariables
ec_fwr_outlier_sig <- ec_fwr_outlier %>% 
  group_by(population, trait) %>% 
  mutate(p_adj = p.adjust(pvalue, "bonf")) %>% 
  filter(p_adj <= alpha)




# Re-plot
## Plot the correlations for the 3 highest correlated variables
# g_env_mean_variable <- top_n(one_year_env_samples_outliers, n_ec, -rank) %>% 
g_env_mean_variable <- top_n(ec_fwr_outlier_sig, n_ec, Fstat) %>% 
  slice(1:n_ec) %>%
  filter(population == "all") %>%
  left_join(., env_combine_outliers) %>% 
  ggplot(aes(x = value, y = h)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  geom_point() + 
  facet_wrap(trait ~ variable, scale = "free", ncol = n_ec) + 
  ylab("Environmental mean") +
  xlab("Covariate") +
  theme_acs()

ggsave(filename = "one_year_env_variable_mean_outliers.jpg", plot = g_env_mean_variable, path = fig_dir, width = 8, height = 6, dpi = 1000)


# g_env_mean_variable <- top_n(one_year_env_samples_outliers, n_ec, -rank) %>% 
g_env_mean_variable <- top_n(ec_fwr_outlier_sig, n_ec, Fstat) %>%
  slice(1:n_ec) %>%
  filter(population == "tp") %>%
  left_join(., env_combine_outliers) %>% 
  ggplot(aes(x = value, y = h)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  geom_point() + 
  facet_wrap(trait ~ variable, scale = "free", ncol = n_ec) + 
  ylab("Environmental mean") +
  xlab("Covariate") +
  theme_acs()

ggsave(filename = "one_year_env_variable_mean_outliers_tp.jpg", plot = g_env_mean_variable, path = fig_dir, width = 8, height = 6, dpi = 1000)







## Create df for the top 1 and top 5 correlated ECs, or all ECs
env_mean_cor_top1 <- one_year_env_samples_outliers %>% 
  group_by(population, trait) %>% 
  top_n(1, -rank) %>% 
  slice(1) %>%
  ungroup() %>%
  left_join(env_combine) %>% 
  select(population, trait, environment, variable, value) %>%
  mutate(group = "Top1EC_Cor")

env_mean_cor_top5 <- one_year_env_samples_outliers %>% 
  group_by(population, trait) %>% 
  top_n(5, -rank) %>% 
  slice(1:5) %>%
  ungroup() %>%
  left_join(env_combine) %>% 
  select(population, trait, environment, variable, value) %>%
  mutate(group = "Top5EC_Cor")

## select on Fstat
env_mean_Fstat_top1 <- ec_fwr_outlier_sig %>% 
  group_by(population, trait) %>% 
  top_n(1, Fstat) %>% 
  slice(1) %>%
  ungroup() %>%
  left_join(env_combine) %>% 
  select(population, trait, environment, variable, value) %>%
  mutate(group = "Top1EC_Fstat")

env_mean_Fstat_top5 <- ec_fwr_outlier_sig %>% 
  group_by(population, trait) %>% 
  top_n(5, Fstat) %>% 
  slice(1:5) %>%
  ungroup() %>%
  left_join(env_combine) %>% 
  select(population, trait, environment, variable, value) %>%
  mutate(group = "Top5EC_Fstat")

env_mean_all <- one_year_env_samples_outliers %>% 
  left_join(env_combine) %>% 
  select(population, trait, environment, variable, value) %>%
  mutate(group = "AllEC")

env_variable_combine <- bind_rows(env_mean_cor_top1, env_mean_cor_top5, env_mean_Fstat_top1, env_mean_Fstat_top5, env_mean_all)



## Combine and create distance matrices using all data or the top correlated EC
env_variable_cov_mat <- env_variable_combine %>%
  group_by(population, trait, group) %>% 
  do(cov = {
    df <- .
    select(df, environment:value) %>%
      spread(variable, value) %>% 
      as.data.frame() %>%
      mutate_at(vars(-environment), scale) %>% 
      column_to_rownames("environment") %>% 
      as.matrix() %>% 
      {tcrossprod(.) / ncol(.)} 
    }) %>%
  # Calculate distance
  ungroup() %>%
  mutate(dist = map(cov, dist))


## Rename
one_year_ec_dist <- env_variable_cov_mat


















#### Multi-year ####


## Correlate environmental variables with the mean

# For each summarized environmental covariable, find the correlation with the environmental mean

multi_year_env_df1 <- multi_year_env_df %>%
  filter(!str_detect(variable, "day")) %>%
  left_join(env_means_all, .)

env_mean_cor <- multi_year_env_df1 %>% 
  group_by(population, trait, variable) %>% 
  do(test = cor.test(.$h, .$value)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         df = map_dbl(test, "parameter")) %>%
  select(-test)

env_mean_cor %>% group_by(population, trait) %>% top_n(n = 2, wt = abs(cor))

# population trait           variable                    cor     pvalue    df
# 1 all        GrainYield      min_PPT                  -0.663 0.0000192     32
# 2 all        GrainYield      ph1to1h2o_r_topsoil       0.594 0.000216      32
# 3 all        HeadingDate     isothermality             0.392 0.0294        29
# 4 all        HeadingDate     temp_seasonality         -0.404 0.0243        29
# 5 all        HeadingDateAGDD month_3_PPT               0.215 0.246         29
# 6 all        HeadingDateAGDD ph1to1h2o_r_topsoil      -0.229 0.214         29
# 7 all        PlantHeight     annual_temperature_range  0.338 0.0588        30
# 8 all        PlantHeight     ph1to1h2o_r_topsoil       0.391 0.0269        30
# 9 tp         GrainYield      annual_temperature_range  0.727 0.00000781    27
# 10 tp         GrainYield      temp_seasonality          0.673 0.0000630     27
# 11 tp         HeadingDate     isothermality             0.438 0.0197        26
# 12 tp         HeadingDate     temp_seasonality         -0.421 0.0258        26
# 13 tp         HeadingDateAGDD month_3_PPT               0.332 0.0848        26
# 14 tp         HeadingDateAGDD month_7_TMAX             -0.317 0.100         26
# 15 tp         PlantHeight     annual_temperature_range  0.439 0.0171        27
# 16 tp         PlantHeight     temp_seasonality          0.532 0.00297       27


## Plot these
env_mean_cor %>% 
  group_by(population, trait) %>%
  top_n(n = 2, wt = abs(cor)) %>%
  left_join(., multi_year_env_df1) %>%
  ggplot(aes(x = value, y = h)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~  trait + population + variable, scales = "free", ncol = 4) +
  theme_acs()
  





# ## What interval variables are highly correlated
# multi_year_env_interval_df1 <- multi_year_daily_summary_interval %>% 
#   map(unnest) %>%
#   map(select, -variable) %>% 
#   reduce(full_join) %>% 
#   gather(variable, value, -begin:-environment) %>% 
#   left_join(env_means_all, .) %>%
#   # Remove NA
#   filter(!is.na(value))
# 
# env_interval_cor <- multi_year_env_interval_df1 %>% 
#   group_by(population, trait, variable, begin, end) %>% 
#   do(test = cor.test(.$h, .$value)) %>%
#   ungroup() %>%
#   mutate(cor = map_dbl(test, "estimate"),
#          pvalue = map_dbl(test, "p.value"),
#          df = map_dbl(test, "parameter")) %>%
#   select(-test)
# 
# 
# 
# ## Plot
# env_interval_cor_plot_list <- env_interval_cor %>% 
#   filter(population == "all") %>%
#   split(list(.$trait, .$variable)) %>% 
#   map(~qplot(x = begin, y = end, fill = cor, geom = "tile", data = .) +
#         scale_fill_gradient2() + facet_wrap(~ trait + variable) +
#         theme_acs() + theme(legend.position = c(0.85, 0.35)))
# 
# # Cowplot
# env_interval_cor_plot <- plot_grid(plotlist = env_interval_cor_plot_list, ncol = n_distinct(env_interval_cor$trait))
# ggsave(filename = "multi_year_interval_variable_correlation_space.jpg", plot = env_interval_cor_plot,
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
# ggsave(filename = "multi_year_interval_variable_correlation_space_tp.jpg", plot = env_interval_cor_plot,
#        path = fig_dir, width = 12, height = 12, dpi = 1000)
# 
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
#   left_join(multi_year_env_interval_df1) %>%
#   unite(variable, variable, begin, end, sep = "_")



## Combine the summary variables with the interval variables
env_combine <- multi_year_env_df1 %>%
  # bind_rows(env_interval_use) %>%
  select(population, trait, environment, h, variable, value)

# Number of covariables to select
n_ec <- 5

## Generate random samples of environments. For each sample, calculate the same
## correlation as above. How often is the same EC strongly correlated?
n_sample <- 1000
f_sample <- 0.6

set.seed(500)
multi_year_env_samples <- env_combine %>%
  group_by(population, trait) %>%
  do({
    df <- .
    # Generate samples of environments
    cor_sample <- replicate(n_sample, sample_frac(tbl = distinct(df, environment), size = f_sample), simplify = FALSE) %>% 
      map(~left_join(., df, by = "environment")) %>% 
      map_df(~group_by(., variable) %>% summarize(cor = cor(h, value)) %>% mutate(rank = rank(-abs(cor))))
    
    # Return the average rank and average correlation of the variables
    cor_sample %>% 
      group_by(variable) %>% 
      summarise_at(vars(cor, rank), mean) %>% 
      arrange(rank)
    
  })


top_n(multi_year_env_samples, n_ec, -rank)


# population           trait                 variable        cor   rank
# 1         all      GrainYield                  min_PPT -0.6575531  1.928
# 2         all      GrainYield      ph1to1h2o_r_topsoil  0.5718817  4.258
# 3         all      GrainYield              month_3_PPT -0.5564224  4.373
# 4         all      GrainYield              month_7_PPT -0.4847908  7.049
# 5         all      GrainYield              month_4_PPT -0.4423563  8.966
# 6         all     HeadingDate         temp_seasonality -0.3982210  5.743
# 7         all     HeadingDate             month_3_TMAX  0.3739004  6.090
# 8         all     HeadingDate            isothermality  0.3876106  6.102
# 9         all     HeadingDate             month_3_TAVG  0.3378879  8.569
# 10        all     HeadingDate             month_3_TMIN  0.2803799 13.326
# 11        all HeadingDateAGDD      ph1to1h2o_r_topsoil -0.2310684 12.450
# 12        all HeadingDateAGDD              month_3_PPT  0.2176072 13.207
# 13        all HeadingDateAGDD         temp_seasonality -0.2106442 14.003
# 14        all HeadingDateAGDD             month_7_TAVG -0.1739534 15.604
# 15        all HeadingDateAGDD annual_temperature_range -0.1844010 15.997
# 16        all     PlantHeight      ph1to1h2o_r_topsoil  0.3734711  7.274
# 17        all     PlantHeight annual_temperature_range  0.3266239  8.820
# 18        all     PlantHeight         temp_seasonality  0.3370670  9.720
# 19        all     PlantHeight                 min_TMIN -0.2576557 13.116
# 20        all     PlantHeight             month_7_TAVG  0.2401431 13.661
# 21         tp      GrainYield annual_temperature_range  0.7123123  1.620
# 22         tp      GrainYield         temp_seasonality  0.6650965  2.694
# 23         tp      GrainYield                  min_PPT -0.5857913  4.818
# 24         tp      GrainYield              month_3_PPT -0.5792967  5.121
# 25         tp      GrainYield      ph1to1h2o_r_topsoil  0.5669364  5.862
# 26         tp     HeadingDate            isothermality  0.4373452  4.858
# 27         tp     HeadingDate             month_3_TMAX  0.4046446  5.774
# 28         tp     HeadingDate         temp_seasonality -0.4166541  5.910
# 29         tp     HeadingDate             month_3_TAVG  0.3528061  9.147
# 30         tp     HeadingDate     annual_diurnal_range  0.2981868 12.832
# 31         tp HeadingDateAGDD              month_3_PPT  0.3309195  8.028
# 32         tp HeadingDateAGDD             month_7_TMAX -0.3107832 10.833
# 33         tp HeadingDateAGDD      ph1to1h2o_r_topsoil -0.3018567 11.213
# 34         tp HeadingDateAGDD                 max_TMAX -0.2765676 12.517
# 35         tp HeadingDateAGDD             month_8_TMAX -0.2632006 13.392
# 36         tp     PlantHeight         temp_seasonality  0.5276867  2.851
# 37         tp     PlantHeight            isothermality -0.4206109  6.228
# 38         tp     PlantHeight annual_temperature_range  0.4328569  6.738
# 39         tp     PlantHeight      ph1to1h2o_r_topsoil  0.3697951 10.552
# 40         tp     PlantHeight                 min_TMIN -0.3443040 10.902


## Print the top 3 per trait
# population           trait                 variable        cor   rank
# 1         all      GrainYield                  min_PPT -0.6624807  2.672
# 2         all      GrainYield phototherm_daily_139_139  0.5922994  4.088
# 3         all      GrainYield phototherm_daily_138_139  0.5547831  6.928
# 4         all     HeadingDate         phototherm_69_69 -0.7067192  4.860
# 5         all     HeadingDate         phototherm_69_70 -0.7065902  5.228
# 6         all     HeadingDate         phototherm_68_69 -0.7063274  5.836
# 7         all HeadingDateAGDD phototherm_daily_139_139 -0.3631174  7.004
# 8         all HeadingDateAGDD              GDD_139_139 -0.3531900  8.484
# 9         all HeadingDateAGDD phototherm_daily_138_139 -0.3283787 12.580
# 10        all     PlantHeight phototherm_daily_139_140  0.4286537 12.700
# 11        all     PlantHeight phototherm_daily_140_140  0.4074949 15.108
# 12        all     PlantHeight phototherm_daily_144_144  0.3900052 17.272
# 13         tp      GrainYield annual_temperature_range  0.7074483  3.360
# 14         tp      GrainYield phototherm_daily_139_139  0.6547078  3.736
# 15         tp      GrainYield phototherm_daily_138_139  0.6258182  5.480
# 16         tp     HeadingDate         phototherm_70_70 -0.7799304  4.216
# 17         tp     HeadingDate         phototherm_69_70 -0.7798161  4.600
# 18         tp     HeadingDate         phototherm_69_71 -0.7797201  4.852
# 19         tp HeadingDateAGDD phototherm_daily_139_139 -0.4529366  9.944
# 20         tp HeadingDateAGDD phototherm_daily_136_139 -0.4476957 10.176
# 21         tp HeadingDateAGDD phototherm_daily_137_139 -0.4402978 11.128
# 22         tp     PlantHeight         temp_seasonality  0.5319591  7.084
# 23         tp     PlantHeight phototherm_daily_144_144  0.4420357 13.864
# 24         tp     PlantHeight phototherm_daily_139_140  0.4397592 15.148


## Instead of correlation, calculate an Fstat based on the linear reaction of genotypes to the covariable
ec_fwr <- S2_MET_BLUEs %>% 
  left_join(., rename(env_combine, ec_value = value)) %>% 
  group_by(population, trait, variable) %>%
  do(fwr(formula = value ~ ec_value, data = .)) %>%
  select(-regression) %>%
  as.data.frame()

# Adjust the p values and select the significant covariables
ec_fwr_sig <- ec_fwr %>% 
  group_by(population, trait) %>% 
  mutate(p_adj = p.adjust(pvalue, "bonf")) %>% 
  filter(p_adj <= alpha)








## Plot the correlations for the 3 highest correlated variables
# g_env_mean_variable <- top_n(multi_year_env_samples, n_ec, -rank) %>% 
g_env_mean_variable <- top_n(ec_fwr_sig, n_ec, Fstat) %>% 
  split(.$population) %>%
  map(~left_join(., env_combine) %>% 
        ggplot(aes(x = value, y = h)) + 
        geom_smooth(method = "lm", se = FALSE) + 
        geom_point() + 
        facet_wrap(trait ~ variable, scale = "free", ncol = n_ec) + 
        ylab("Environmental mean") +
        xlab("Covariate") +
        theme_acs() )

ggsave(filename = "multi_year_env_variable_mean.jpg", plot = g_env_mean_variable$all, path = fig_dir, width = 8, height = 6, dpi = 1000)

ggsave(filename = "multi_year_env_variable_mean_tp.jpg", plot = g_env_mean_variable$tp, path = fig_dir, width = 8, height = 6, dpi = 1000)








## Examine some outliers
env_combine %>%
  filter((trait == "PlantHeight" & h > 25) | (trait == "GrainYield" & h > 3000) | (trait == "HeadingDateAGDD" & h > 500)) %>%
  distinct(population, trait, environment, h)

# Remove outliers
env_combine_outliers <- env_combine %>%
  filter(!(trait == "PlantHeight" & environment == "CRM17")) %>%
  filter(!(trait == "GrainYield" & h > 3000)) %>%
  filter(!(trait == "HeadingDateAGDD" & h > 500))


set.seed(500)
# Remove outliers and re-adjust
multi_year_env_samples_outliers <- env_combine_outliers %>%
  group_by(population, trait) %>%
  do({
    df <- .
    # Generate samples of environments
    cor_sample <- replicate(n_sample, sample_frac(tbl = distinct(df, environment), size = f_sample), simplify = FALSE) %>% 
      map(~left_join(., df, by = "environment")) %>% 
      map_df(~group_by(., variable) %>% summarize(cor = cor(h, value)) %>% mutate(rank = rank(-abs(cor))))
    
    # Return the average rank and average correlation of the variables
    cor_sample %>% 
      group_by(variable) %>% 
      summarise_at(vars(cor, rank), mean) %>% 
      arrange(rank)
    
  })

## Instead of correlation, calculate an Fstat based on the linear reaction of genotypes to the covariable
ec_fwr_outlier <- S2_MET_BLUEs %>% 
  inner_join(., rename(env_combine_outliers, ec_value = value)) %>% 
  group_by(population, trait, variable) %>%
  do(fwr(formula = value ~ ec_value, data = .)) %>%
  select(-regression) %>%
  as.data.frame()

# Adjust the p values and select the significant covariables
ec_fwr_outlier_sig <- ec_fwr_outlier %>% 
  group_by(population, trait) %>% 
  mutate(p_adj = p.adjust(pvalue, "bonf")) %>% 
  filter(p_adj <= alpha)



# Re-plot
## Plot the correlations for the 3 highest correlated variables
# g_env_mean_variable <- top_n(multi_year_env_samples_outliers, n_ec, -rank) %>% 
g_env_mean_variable <- top_n(ec_fwr_outlier_sig, n_ec, Fstat) %>% 
  slice(1:5) %>%
  split(.$population) %>%
  map(~left_join(., env_combine_outliers) %>% 
        ggplot(aes(x = value, y = h)) + 
        geom_smooth(method = "lm", se = FALSE) + 
        geom_point() + 
        facet_wrap(trait ~ variable, scale = "free", ncol = n_ec) + 
        ylab("Environmental mean") +
        xlab("Covariate") +
        theme_acs() )

ggsave(filename = "multi_year_env_variable_mean_outliers.jpg", plot = g_env_mean_variable$all, path = fig_dir, width = 8, height = 6, dpi = 1000)

ggsave(filename = "multi_year_env_variable_mean_outliers_tp.jpg", plot = g_env_mean_variable$tpex, path = fig_dir, width = 8, height = 6, dpi = 1000)





## Create df for the top 1 and top 5 correlated ECs, or all ECs
env_mean_cor_top1 <- multi_year_env_samples_outliers %>% 
  group_by(population, trait) %>% 
  top_n(1, -rank) %>% 
  slice(1) %>%
  ungroup() %>%
  left_join(env_combine) %>% 
  select(population, trait, environment, variable, value) %>%
  mutate(group = "Top1EC_Cor")

env_mean_cor_top5 <- multi_year_env_samples_outliers %>% 
  group_by(population, trait) %>% 
  top_n(5, -rank) %>% 
  slice(1:5) %>%
  ungroup() %>%
  left_join(env_combine) %>% 
  select(population, trait, environment, variable, value) %>%
  mutate(group = "Top5EC_Cor")

## select on Fstat
env_mean_Fstat_top1 <- ec_fwr_outlier_sig %>% 
  group_by(population, trait) %>% 
  top_n(1, Fstat) %>% 
  slice(1) %>%
  ungroup() %>%
  left_join(env_combine) %>% 
  select(population, trait, environment, variable, value) %>%
  mutate(group = "Top1EC_Fstat")

env_mean_Fstat_top5 <- ec_fwr_outlier_sig %>% 
  group_by(population, trait) %>% 
  top_n(5, Fstat) %>% 
  slice(1:5) %>%
  ungroup() %>%
  left_join(env_combine) %>% 
  select(population, trait, environment, variable, value) %>%
  mutate(group = "Top5EC_Fstat")

env_mean_all <- multi_year_env_samples_outliers %>% 
  left_join(env_combine) %>% 
  select(population, trait, environment, variable, value) %>%
  mutate(group = "AllEC")

env_variable_combine <- bind_rows(env_mean_cor_top1, env_mean_cor_top5, env_mean_Fstat_top1, env_mean_Fstat_top5, env_mean_all)




## Combine and create distance matrices using all data or the top correlated EC
env_variable_cov_mat <- env_variable_combine %>%
  group_by(population, trait, group) %>% 
  do(cov = {
    df <- .
    select(df, environment:value) %>%
      spread(variable, value) %>% 
      as.data.frame() %>%
      mutate_at(vars(-environment), scale) %>% 
      column_to_rownames("environment") %>% 
      as.matrix() %>% 
      {tcrossprod(.) / ncol(.)} 
  }) %>%
  # Calculate distance
  ungroup() %>%
  mutate(dist = map(cov, dist))


## Rename
multi_year_ec_dist <- env_variable_cov_mat




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
  facet_wrap(trait ~ group, scales = "free", ncol = 3) + 
  theme_acs()

ggsave(filename = "one_year_env_ec_mds.jpg", plot = one_year_ec_mds_plot, path = fig_dir,
       height = 8, width = 6, dpi = 1000)




multi_year_ec_mds <- multi_year_ec_dist %>% 
  mutate(mds = map(dist, ~cmdscale(.) %>% as.data.frame() %>% `names<-`(., c("x", "y")) %>% rownames_to_column("environment"))) %>% 
  unnest(mds)

# Plot
multi_year_ec_mds_plot <- multi_year_ec_mds %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_text(aes(label = environment), size = 2) + 
  # geom_point() + 
  # geom_text_repel(aes(label = environment), size = 2, ) + 
  facet_wrap(trait ~ group, scales = "free", ncol = 3) + 
  theme_acs()

ggsave(filename = "multi_year_env_ec_mds.jpg", plot = multi_year_ec_mds_plot, path = fig_dir,
       height = 8, width = 6, dpi = 1000)







