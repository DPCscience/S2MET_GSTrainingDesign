## S2MET Predictions
## 
## Growth stage prediction testing
## 
## Author: Jeff Neyhart
## Last modified: 29 July 2019
## 
## This R script will predict the environmental average growth stage
## based on GDD.
## 

# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load other libraries
library(lubridate)
library(rvest)
library(modelr)
library(broom)

# Load phenotypic data - trial BLUEs
load("C:/GoogleDrive/BarleyLab/Breeding/PhenotypicData/Final/MasterPhenotypes/S2_tidy_BLUE.RData")

## Filter BLUEs for relevant traits and for appropriate trials
S2_MET_BLUEs <- s2_tidy_BLUE %>% 
  filter(trait %in% c(traits, "MaturityDate"), 
         trial %in% trial_info$trial,
         environment %in% tp_vp_env,
         line_name %in% c(tp, vp)) %>%
  rename(std_error = std.error) %>%
  droplevels()


## Load the weather dataset
load(file.path(data_dir, "Climate_Data/NOAA_Data/noaa_stations_trial_data.RData"))

# Load the GDD information and growth stages
load(file.path(data_dir, "agdd_growth_staging.RData"))
## Load phenotypic data
load(file.path(result_dir, "genotype_environment_phenotypic_analysis.RData"))


### Test Haun versus Zadoks growth stage

env_test <- one_year_growth_staging %>%
  filter(environment == environment[1])

# Predict growth stage using polynomial model (Bauer)
agdd <- env_test$agdd


## Fit a different model to haun growth stage
ndawn_table <- html_session(url = "https://ndawn.ndsu.nodak.edu/help-barley-growing-degree-days.html") %>%
  html_table()

ndawn_haun <- ndawn_table[[1]] %>%
  tbl_df() %>%
  filter(!is.na(GDDRequired)) %>%
  mutate(HaunStage = parse_number(HaunStage)) %>%
  rename_all(make.names) %>%
  rename(gdd_far = GDDRequired, agdd_far = AccumulatedGDD.Req) %>%
  mutate(gdd_cent = gdd_far / 1.8, agdd_cent = cumsum(gdd_cent))

# Model
haun_fit <- lm(HaunStage ~ agdd_far, data = ndawn_haun)
haun_fit2 <- lm(HaunStage ~ agdd_cent, data = ndawn_haun)


## Get the haun/zadoks/feekes growth stage table
growth_stage_table <- html_session(url = "https://www.ag.ndsu.edu/crops/spring-wheat-articles/small-grains-development-using-zadoks-feekes-haun") %>%
  html_table() %>% 
  map_df(~.) %>%
  tbl_df() %>%
  mutate_at(vars(-Description), parse_number) %>%
  filter(!is.na(Zadoks))


## Add predictions
# env_test1 <- env_test %>%
growth_staging <- bind_rows(
  mutate(one_year_growth_staging, timescale = "oneyear"),
  select(mutate(multiyear_growth_staging, timescale = "multiyear"), environment, dap, timescale, gdd = avg_gdd, agdd = avg_agdd)
) %>%
  select(-contains("growth_stage")) %>%
  group_by(environment, timescale) %>%
  mutate(bauer_gs = stage_growth(agdd, method = "haun.bauer"),
         enz_gs = stage_growth(agdd, method = "haun.enz"),
         montana_gs = stage_growth(agdd, method = "montana")) %>%
  ungroup()
    






### Validate by comparing the predicted flowering time via AGDD with the average
### heading date in an environment

data_to_model <- training_sets_twoway %>%
  filter(set %in% c("complete", "realistic2017")) %>%
  select(-varE) %>% 
  unnest() %>%
  distinct(set, environment) %>%
  inner_join(., S2_MET_BLUEs)


# Fit the model per trait, extract the environmental means
full_model_fits <- data_to_model %>%
  mutate_at(., vars(line_name, environment), as.factor) %>%
  group_by(trait, set) %>%
  do({
    
    df <- droplevels(.)
    # print(unique(df$trait))
    df1 <- df
    ## Adjust contrasts of environments
    contrasts(df1$environment) <- `colnames<-`(contr.sum(levels(df1$environment)), head(levels(df1$environment), -1))
    
    # Control and weights
    control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", calc.derivs = FALSE)
    wts <- df1$std_error^2
    
    # Formula
    form <- value ~ (1|line_name) + environment
    # Fit the model
    fit <- lmer(formula = form, data = df1, control = control, weights = wts)
    
    ## Get coefficients for environments
    env_eff <- fixef(fit)
    # Tidy
    env_means <- tibble(environment = levels(df1$environment), h = c(env_eff[-1], -sum(env_eff[-1])), mean = h + env_eff[1]) 
    
    tibble(env_means = list(env_means), model = list(fit))
    
  }) %>% ungroup()

# Unnest the environmental mean DFs
env_means_all <- full_model_fits %>%
  unnest(env_means)


###########################
## Heading Date
###########################


## Use the growth staging to determine the dap when flowering begins
flowering_agdd <- growth_staging %>%
  select(environment, dap, timescale, contains("gs")) %>%
  gather(model, stage, contains("gs")) %>%
  filter(stage == "flowering") %>%
  group_by(environment, model, timescale) %>%
  summarize(dap = mean(dap)) %>%
  ungroup()



## Compare average flowering with predicted
compare_flowering_oneyear <- env_means_all %>% 
  filter(trait == "HeadingDate") %>%
  inner_join(., flowering_agdd) %>%
  # Calculate bias
  mutate(bias = dap - mean) %>%
  # Prepare for plotting
  mutate(set = str_remove_all(set, "2017"),
         set = str_replace_all(set, set_replace))

aggregate(bias ~ set + model + timescale, compare_flowering_oneyear, mean)


## The enz method has the lowest bias

## Summarize
compare_flowering_oneyear %>%
  group_by(set, model, timescale) %>%
  do(test = cor.test(.$mean, .$dap)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"))

# Plot
limits <- select(compare_flowering_oneyear, mean, dap) %>% 
  {c(min(.), max(.))}

g_compare_flowering_oneyear <- compare_flowering_oneyear %>%
  split(.$timescale) %>%
  map(~{
    df <- .
    
    ggplot(data = ., aes(x = dap, y = mean)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_point() +
      ggrepel::geom_text_repel(data = subset(df, abs(bias) - 5 > 5), aes(label = environment), size = 3) +
      facet_grid(model ~ set) +
      scale_y_continuous(name = "Average observed heading date", breaks = pretty, limits = limits) +
      scale_x_continuous(name = "GDD predicted flowering date", breaks = pretty, limits = limits) +
      theme_presentation2() })

## Save
ggsave(filename = "one_year_gdd_predicted_FT.jpg", plot = g_compare_flowering_oneyear$oneyear,
       path = fig_dir, width = 6, height = 4, dpi = 1000)
ggsave(filename = "multiyear_gdd_predicted_FT.jpg", plot = g_compare_flowering_oneyear$multiyear,
       path = fig_dir, width = 6, height = 4, dpi = 1000)





###########################
## Maturity Date
###########################



## Use the growth staging to determine the dap when flowering begins
maturity_agdd <- growth_staging %>%
  select(environment, dap, timescale, contains("gs")) %>%
  gather(model, stage, contains("gs")) %>%
  filter(stage == "maturity") %>%
  group_by(environment, model, timescale) %>%
  summarize(dap = min(dap)) %>%
  ungroup()



## Compare average flowering with predicted
compare_maturity <- env_means_all %>% 
  filter(trait == "MaturityDate") %>%
  inner_join(., flowering_agdd) %>%
  # Calculate bias
  mutate(bias = dap - mean) %>%
  # Prepare for plotting
  mutate(set = str_remove_all(set, "2017"),
         set = str_replace_all(set, set_replace))

aggregate(bias ~ set + model + timescale, compare_maturity, mean)


## The enz method has the lowest bias

## Summarize
compare_maturity %>%
  group_by(set, model, timescale) %>%
  do(test = cor.test(.$mean, .$dap)) %>%
  ungroup() %>%
  mutate(cor = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"))

# Plot
limits <- select(compare_maturity, mean, dap) %>% 
  {c(min(.), max(.))}

g_compare_maturity <- compare_maturity %>%
  split(.$timescale) %>%
  map(~{
    df <- .
    
    ggplot(data = ., aes(x = dap, y = mean)) +
      geom_abline(slope = 1, intercept = 0) +
      geom_point() +
      ggrepel::geom_text_repel(data = subset(df, abs(bias) - 5 > 5), aes(label = environment), size = 3) +
      facet_grid(model ~ set) +
      scale_y_continuous(name = "Average observed maturity date", breaks = pretty, limits = limits) +
      scale_x_continuous(name = "GDD predicted maturity date", breaks = pretty, limits = limits) +
      theme_presentation2() })




## Try to use heading and maturity date to calibrate genotype-specific models


# Subset heading and maturity info
phenology_phenos <- S2_MET_BLUEs %>%
  filter(trait %in% c("HeadingDate", "MaturityDate")) %>%
  # Express these traits in AGDD
  left_join(., subset(growth_staging, timescale == "oneyear", select = c(environment, dap, agdd))) %>%
  # Round hd and md to whole days
  mutate(value = round(value),
         stage = ifelse(trait == "HeadingDate", 10.5, 16)) %>%
  filter(value == dap) %>%
  group_by(environment) %>%
  filter(n_distinct(trait) == 2) %>%
  ungroup()

## Filter for environments that have both HD and MD info
phenology_phenos1 <- phenology_phenos %>%
  # distinct(environment, line_name) %>%
  # crossing(., data.frame(trait = "PlantingDate", stage = 0, agdd = 0, stringsAsFactors = FALSE)) %>%
  # full_join(phenology_phenos, ., by = c("environment", "line_name", "trait", "agdd", "stage")) %>%
  mutate(line_name = as.factor(line_name))
  
## Fit a model per genotype
phenology_models <- phenology_phenos1 %>%
  group_by(line_name) %>%
  do(fit = lm(stage ~ 1 + agdd, data = .) ) %>%
  ungroup()

## Extract coefficients
phenology_coef <- phenology_models %>%
  mutate(coefs = map(fit, ~tibble(intercept = coef(.)[1], agdd = coef(.)[2]))) %>%
  unnest(coefs) %>%
  select(-fit)


## Plot regression lines
phenology_phenos1 %>%
  ggplot(aes(x = agdd, y = stage, shape = trait, color = line_name, group = line_name)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_discrete(guide = FALSE)
    

## Fit a model using all data (test for variation among slopes)
phenology_full_model <- lm(stage ~ -1 + agdd + line_name:agdd, data = phenology_phenos1)
anova(phenology_full_model)

## Anova indicates no difference among slopes

# Should we allow diffence in intercept?
# No, because all plants should start at the same stage

# Get effects
phenology_full_coef <- tidy(phenology_full_model) %>% 
  mutate(line_name = levels(phenology_phenos1$line_name)) %>% 
  select(line_name, slope = estimate) %>% 
  mutate(slope = ifelse(line_name == line_name[1], slope, slope + slope[1]),
         intercept = 0)


## Use the models to predict growth stage for each genotype in all environments
predicted_growth_stage <- subset(growth_staging, timescale == "oneyear", select = c(environment, dap, agdd)) %>%
  crossing(., phenology_full_coef) %>%
  mutate(stage = intercept + (slope * agdd))

## Determine predicted heading date
predicted_hd <- predicted_growth_stage %>%
  group_by(environment, line_name) %>%
  top_n(x = ., n = 1, wt = -abs(stage - 10.5)) %>%
  ungroup() %>%
  select(environment, line_name, predicted_hd = dap, agdd, stage)
  


## Compare with observed
predicted_obs_hd <- S2_MET_BLUEs %>%
  filter(trait %in% c("HeadingDate")) %>%
  left_join(., predicted_hd)

## Correlation by environment
predicted_obs_hd %>%
  group_by(environment) %>%
  summarize(acc = cor(value, predicted_hd))

# Overall correlation
cor(predicted_obs_hd$value, predicted_obs_hd$predicted_hd)

# Plot
limits <- range(unlist(predicted_obs_hd[,c("predicted_hd", "value")]))
predicted_obs_hd %>%
  ggplot(aes(x = value, y = predicted_hd, color = environment)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_discrete(guide = FALSE) +
  scale_x_continuous(breaks = pretty, limits = limits) +
  scale_y_continuous(breaks = pretty, limits = limits)
  


  























