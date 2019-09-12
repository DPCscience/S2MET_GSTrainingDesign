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
library(lmerTest)
library(cowplot)
## ggrepel
library(ggrepel)


# Load the distance matrices
load(file.path(result_dir, "distance_method_results.RData"))
# Load the tidy S2 data
load("C:/GoogleDrive/BarleyLab/Breeding/PhenotypicData/Final/MasterPhenotypes/S2_tidy_pheno.RData")



## significance level
alpha <- 0.05


## Filter the BLUEs for the environments of interest
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  filter(environment %in% tp_vp_env)


## Determine the number of environments per trait
S2_MET_BLUEs_use %>% 
  group_by(trait) %>% 
  summarize_at(vars(location, year, environment), n_distinct)

# trait       location  year environment
# 1 GrainYield        13     3          23
# 2 HeadingDate       12     3          26
# 3 PlantHeight       13     3          27



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

# trait       prop_obs
# 1 GrainYield     0.986
# 2 HeadingDate    0.991
# 3 PlantHeight    0.991


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
            summarize(trait = unique(trait), prop_obs = mean(observed), n_obs = sum(observed), n_total = n(), n_pred = n_total - n_obs)
        })
    }))

# $`tp`
# # A tibble: 3 x 4
# trait       prop_obs n_obs n_total
# 1 GrainYield     0.984  4140    4209
# 2 HeadingDate    0.989  4704    4758
# 3 PlantHeight    0.989  4888    4941
# 
# $tp_geno
# # A tibble: 3 x 4
# trait       prop_obs n_obs n_total
# 1 GrainYield     0.984  3960    4025
# 2 HeadingDate    0.989  4500    4550
# 3 PlantHeight    0.989  4675    4725



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
left_join(distinct(S2_MET_BLUEs_use, trait, trial, environment), stage_one_data) %>% 
  group_by(trait) %>% 
  summarize_at(vars(heritability), funs(min, max, mean))

# trait         min   max  mean
# 1 GrainYield  0.180 0.856 0.542
# 2 HeadingDate 0.563 0.978 0.846
# 3 PlantHeight 0.106 0.883 0.532


## Output a table
intra_env_herit <- left_join(distinct(S2_MET_BLUEs_use, trait, trial, environment), stage_one_data)  %>%
  filter(!str_detect(trial, "S2C1")) %>%
  distinct(environment, trait, heritability) %>% 
  mutate(heritability = formatC(heritability, digits = 2, width = 2, format = "g", flag = "#")) %>%
  spread(trait, heritability)

write_csv(x = intra_env_herit, path = file.path(fig_dir, "intra_environment_heritability.csv"))







## Plot a map

# First create a df to plot
trial_info_toplot <- trial_info %>% 
  filter(is.na(notes)) %>%
  group_by(location, latitude, longitude) %>%
  summarize(n_years = n_distinct(year)) %>%
  ungroup()

## Order environments on latitude
# Sort on latitute
loc_order <- trial_info_toplot %>% 
  arrange(latitude, location) %>% 
  pull(location)


# Create a gradient of colors
f_colors <- colorRampPalette(colors = umn_palette(2)[3:5])
colors_use <- set_names(f_colors(length(loc_order)), loc_order)


# Get the map data for canada
canada <- map_data("world", "Canada")

# Download map data for US by county
usa_county <- map_data(map = "county")
# Download state data
usa_state <- map_data(map = "state")

# Adjust the groups in the states
usa_state <- usa_state %>%
  mutate(group = group + max(canada$group))

# Adjust the groups in the counties
usa_county <- usa_county %>%
  mutate(group = group + max(usa_state$group))

# Tidy and combine
north_america <- bind_rows(usa_state, usa_county, canada)


## Now create labels for each location (as a combination of all environments)
loc_info_toplot <- trial_info_toplot %>%
  mutate(n_years = factor(n_years, levels = sort(unique(n_years))))

## Now create labels for each location (as a combination of all environments)
use_loc_info_toplot <-  trial_info %>% 
  filter(environment %in% tp_vp_env) %>%
  filter(is.na(notes)) %>%
  group_by(location, latitude, longitude) %>%
  summarize(n_years = n_distinct(year)) %>%
  ungroup() %>%
  mutate(n_years = factor(n_years, levels = sort(unique(n_years))),
         location = str_to_title(str_replace_all(location, "_", " ")))



## Plot only the locations used in this study

# Map
g_map <- ggplot(data = north_america, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = canada, fill = NA, color = "grey50", lwd = 0.3) + # Add canada
  geom_polygon(data = usa_state, aes(x = long, y = lat, group = group), fill = NA, color = "grey50", lwd = 0.3) +
  # geom_point(data = use_loc_info_toplot, aes(x = longitude, y = latitude, group = location, color = location, size = n_years), alpha = 0.70) +
  geom_point(data = use_loc_info_toplot, aes(x = longitude, y = latitude, size = n_years, group = location), alpha = 0.70, color = neyhart_palette("umn2")[3]) +
  coord_fixed(ratio = 1.5, xlim = c(-125, -65), ylim = c(35, 50), ) +
  scale_color_manual(guide = FALSE, values = colors_use) + 
  scale_size_discrete(name = "Number of years observed\nat location", guide = guide_legend(title.position = "left", order = 1, title.hjust = 1)) +
  scale_x_continuous(breaks = pretty, name = "Longitude") + 
  scale_y_continuous(breaks = pretty, name = "Latitude") +
  theme_presentation2(base_size = 8) + 
  theme(panel.background = element_blank(), panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center", plot.title = element_text(hjust = 0.5)) +
  ggsn::north(location = "bottomright", symbol = 12, x.min = -125, x.max = -60, 
              y.min = 35, y.max = 50) + 
  ggsn::scalebar(x.min = -125, x.max = -60, y.min = 35, y.max = 50, dist = 500, 
                 dd2km = TRUE, model = "WGS84", location = "bottomleft", st.size = 1.5, 
                 st.dist = 0.05, st.bottom = FALSE)


# Save the figure
ggsave(filename = "site_map_paper.jpg", plot = g_map, path = fig_dir, width = 3.5, height = 3, dpi = 1000)


## Collapse Ithaca
use_loc_info_toplot1 <- use_loc_info_toplot %>% 
  mutate(location = str_remove(location, "[0-9]{1}"), 
         n_years = parse_number(as.character(n_years))) %>% 
  group_by(location) %>% 
  mutate(n_trials = sum(n_years)) %>% 
  slice(1) %>%
  ungroup()
  

# Coordinate limits
long_limit <- c(-115, -70)
lat_limit <- c(37, 50)

## Different version of the map - grey
g_map_alt <- ggplot(data = north_america, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "grey85") +
  geom_polygon(data = canada, fill = NA, color = "white", lwd = 0.3) + # Add canada
  geom_polygon(data = usa_state, aes(x = long, y = lat, group = group), fill = NA, color = "white", lwd = 0.3) +
  # geom_point(data = use_loc_info_toplot, aes(x = longitude, y = latitude, group = location, color = location, size = n_years), alpha = 0.70) +
  geom_point(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location), size = 3.5) +
  ggrepel::geom_text_repel(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, label = location), 
                           size = 2, hjust = 0.5, nudge_x = -1, segment.size = 0.2, point.padding = unit(2, "pt"),
                           min.segment.length = 1) +
  geom_text(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, label = n_trials), size = 2, 
            color = ifelse(use_loc_info_toplot1$location == "Arlington", "black", "white")) +
  coord_fixed(ratio = 1.5, xlim = long_limit, ylim = lat_limit) +
  scale_color_manual(guide = FALSE, values = colors_use) + 
  scale_x_continuous(breaks = NULL, name = NULL, labels = NULL) + 
  scale_y_continuous(breaks = NULL, name = NULL, labels = NULL) +
  theme_classic() +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill = alpha("white", 0)), axis.line = element_blank()) +
  ggsn::north(location = "bottomleft", symbol = 12, x.min = min(long_limit) - 2, x.max = max(long_limit), 
              y.min = min(lat_limit) + 2, y.max = max(lat_limit)) + 
  ggsn::scalebar(x.min = min(long_limit), x.max = max(long_limit), y.min = min(lat_limit), y.max = max(lat_limit), dist = 500, 
                 dd2km = TRUE, model = "WGS84", location = "bottomleft", st.size = 1.5,
                 st.dist = 0.05, st.bottom = FALSE)


# Save the figure
ggsave(filename = "site_map_paper_alt.jpg", plot = g_map_alt, path = fig_dir, width = 3.5, height = 1.75, dpi = 1000)





## Now create labels for each location (as a combination of all environments)
loc_info_toplot <-  trial_info %>% 
  filter(is.na(notes)) %>%
  group_by(location, latitude, longitude) %>%
  summarize(n_years = n_distinct(year)) %>%
  ungroup() %>%
  mutate(n_years = factor(n_years, levels = sort(unique(n_years))))
  


# Map
g_map1 <- ggplot(data = north_america, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = canada, fill = NA, color = "grey50", lwd = 0.5) + # Add canada
  geom_polygon(data = usa_state, aes(x = long, y = lat, group = group), fill = NA, color = "grey50", lwd = 0.5) +
  geom_point(data = loc_info_toplot, aes(x = longitude, y = latitude, group = location, color = location, size = n_years), alpha = 0.70) +
  coord_map(projection = "albers", lat0 = 35, lat1 = 50,  xlim = c(-125, -65), ylim = c(35, 50)) +
  scale_color_manual(guide = FALSE, values = colors_use) + 
  scale_size_discrete(name = "Number of years observed\nat location", guide = guide_legend(title.position = "left", order = 1, title.hjust = 1)) +
  scale_x_continuous(breaks = pretty) + 
  scale_y_continuous(breaks = pretty) +
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_presentation2(base_size = 10) + 
  theme(panel.background = element_blank(), panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center", plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())


# Save the figure
ggsave(filename = "site_map1.jpg", plot = g_map1, path = fig_dir,
       width = 4, height = 3, dpi = 1000)


## Project map - no color
g_map_nocolor_base <- ggplot(data = north_america, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = canada, fill = NA, color = "grey50", lwd = 0.3) + # Add canada
  geom_polygon(data = usa_state, aes(x = long, y = lat, group = group), fill = NA, color = "grey50", lwd = 0.3) +
  coord_map(projection = "albers", lat0 = 35, lat1 = 50,  xlim = c(-125, -65), ylim = c(35, 48.1)) +
  # scale_color_manual(guide = FALSE, values = colors_use) + 
  scale_size_discrete(name = "Number of trial years at a location", guide = guide_legend(title.position = "left", order = 1, title.hjust = 1)) +
  scale_x_continuous(breaks = pretty) + 
  scale_y_continuous(breaks = pretty) +
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_presentation2(base_size = 10) + 
  theme(panel.background = element_blank(), panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center", plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

# Save the figure
ggsave(filename = "site_map_nocolor_blank.jpg", plot = g_map_nocolor_base, path = fig_dir,
       width = 5, height = 2.5, dpi = 1000)

g_map_nocolor <- g_map_nocolor_base +
  geom_point(data = loc_info_toplot, aes(x = longitude, y = latitude, group = location, size = n_years), 
             color = neyhart_palette("barley")[1], alpha = 0.70)

# Save the figure
ggsave(filename = "site_map_nocolor.jpg", plot = g_map_nocolor, path = fig_dir,
       width = 5, height = 4, dpi = 1000)


g_map_nocolor1 <- g_map_nocolor_base + 
  geom_point(data = use_loc_info_toplot, aes(x = longitude, y = latitude, group = location, size = n_years), 
             color = neyhart_palette("barley")[1], alpha = 0.70)

# Save the figure
ggsave(filename = "site_map_nocolor_filter.jpg", plot = g_map_nocolor1, path = fig_dir,
       width = 5, height = 4, dpi = 1000)











## Visualization of distributions
env_order <- S2_MET_BLUEs_use %>% 
  distinct(environment, location, year) %>%
  left_join(., data_frame(location = names(colors_use), color = colors_use)) %>%
  mutate(location = factor(location, levels = loc_order)) %>% 
  arrange(location, year) %>% 
  {set_names(x = .$color, .$environment)}

## Sort on Lat/Long and year
S2_MET_BLUEs_toplot <- S2_MET_BLUEs_use %>%
  mutate(environment = parse_factor(environment, levels = env_order))


## Plot
g_met_dist <- S2_MET_BLUEs_toplot %>%
  ggplot(aes(x = value, y = environment, fill = environment)) +
  geom_density_ridges() +
  facet_grid(. ~ trait, scales = "free_x") +
  scale_fill_manual(values = env_order, guide = FALSE) +
  ylab("Environment") +
  xlab("Phenotypic value") +
  theme_presentation2(base_size = 10)

# Save it
ggsave(filename = "met_trait_dist.jpg", plot = g_met_dist, path = fig_dir, width = 4.5, height = 5, dpi = 1000)


## Combine map with distributions
g_map_dist_combine <- plot_grid(g_map1, g_met_dist, ncol = 1, rel_heights = c(0.65, 1))
ggsave(filename = "map_and_trait_dist.jpg", plot = g_map_dist_combine, path = fig_dir, width = 5, height = 6, dpi = 1000)












## Remove the maps package so it doesn't clobber purrr
detach("package:maps")




## Stage-Two analysis


## Combine data
S2_MET_BLUEs_tomodel <- bind_rows(mutate(S2_MET_BLUEs_use, population = "all"), 
                                  mutate(filter(S2_MET_BLUEs_use, line_name %in% tp), population = "tp"),
                                  mutate(filter(S2_MET_BLUEs_use, line_name %in% vp), population = "vp"))

S2_MET_BLUEs_tomodel1 <- S2_MET_BLUEs %>%
  filter(environment %in% tp_vp_env)

S2_MET_BLUEs_tomodel1 <- bind_rows(mutate(S2_MET_BLUEs_tomodel1, population = "all"), 
                                  mutate(filter(S2_MET_BLUEs_tomodel1, line_name %in% tp), population = "tp"),
                                  mutate(filter(S2_MET_BLUEs_tomodel1, line_name %in% vp), population = "vp"))

S2_MET_BLUEs_tomodel <- S2_MET_BLUEs_tomodel1


# Boot reps
boot_reps <- 10

# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP

#### Note
#### 
#### I am not sure whether to weight the residuals. It seems to result in a more realistic estimate of residual
#### variance, but it may not be correct.
#### 


# stage_two_fits_GE <- S2_MET_BLUEs_tomodel %>%
#   mutate_at(vars(location, year, line_name), as.factor) %>%
#   group_by(trait, population) %>%
#   do({
#     
#     df <- droplevels(.)
#     
#     # Table of lines by environments (i.e. plots)
#     plot_table <- xtabs(formula = ~ line_name + environment, data = df)
#     
#     ## Harmonic means
#     # Locations
#     harm_env <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
#       ifelse(. > 1, 1, .) %>%
#       rowSums() %>% 
#       harm_mean()
#     
#     # Reps
#     harm_rep <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
#       harm_mean()
#     
#     # # Lmer control
#     # lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", calc.derivs = FALSE,
#     #                             optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))
#     
#     lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
#     
#     # Get the weights
#     wts <- pull(df, std_error)^2
# 
#     ## Fit the full model
#     ## Environment fixed
#     fit <- lmer(formula = value ~ 1 + environment + (1|line_name) + (1|line_name:environment),
#                 data = df, control = lmer_control, weights = wts)
#     
#     # ## Fit the full model
#     # fit <- lmer(formula = value ~ 1 + (1|line_name) + (1|environment) + (1|line_name:environment), 
#     #             data = df, control = lmer_control)
#     
#     ## Likelihood ratio tests
#     lrt <- ranova(fit) %>% 
#       tidy() %>% 
#       filter(!str_detect(term, "none")) %>% 
#       mutate(term = str_remove(term, "\\(1 \\| ") %>% 
#                str_remove("\\)")) %>% 
#       select(term, LRT, df, p.value)
#     
#     ## Calculate heritability
#     exp <- "line_name / (line_name + (line_name:environment / n_e) + (Residual / (n_e * n_r)))"
#     
#     ## Use bootstrapping to calculate a confidence interval
#     # Generate bootstrapping samples and calculate heritability using the bootMer function
#     h2_boot <- bootMer(x = fit, nsim = boot_reps, FUN = function(x) 
#       herit(object = x, exp = exp, n_e = harm_env, n_r = harm_rep)$heritability)
#     
#     h2 <- herit(object = fit, n_e = harm_env, n_r = harm_rep, exp = exp)
#     
#     # Add the bootstrapping results with a confidence interval
#     h2$heritability <- tidy(h2_boot) %>% 
#       cbind(., t(quantile(h2_boot$t, probs = c(alpha / 2, 1 - (alpha / 2))))) %>% 
#       rename_at(vars(4, 5), ~c("lower", "upper"))
#     
#     # Return data_frame
#     data_frame(fit = list(fit), lrt = list(lrt), h2 = list(h2), n_e = harm_env, n_r = harm_rep) 
#     
#   }) %>% ungroup()
# 
# 
# 
# 
# # Everything is significant
# 
# ## Plot the proportion of variance from each source
# stage_two_fits_GE_varprop <- stage_two_fits_GE %>%
#   mutate(h2 = map(h2, "var_comp")) %>%
#   unnest(h2) %>% 
#   group_by(trait, population) %>% 
#   mutate(var_prop = variance / sum(variance),
#          source = str_replace_all(source, c("line_name:environment" = "Genotype x Environment",
#                                             "environment" = "Environment", "line_name" = "Genotype")),
#          source = factor(source, levels = c("Environment", "Genotype", "Genotype x Environment",
#                                             "Residual"))) %>%
#   ungroup()
# 
# ## Plot both populations and all traits
# g_varprop <- stage_two_fits_GE_varprop %>% 
#   mutate(population = str_to_upper(population)) %>% 
#   ggplot(aes(x = trait, y = var_prop, fill = source)) + 
#   geom_col() +
#   ylab("Proportion of Variance") +
#   scale_fill_manual(values = umn_palette(3, 4), name = NULL) +
#   facet_grid(~ population) + 
#   theme_acs() +
#   theme(axis.title.x = element_blank(),
#         legend.position = "bottom")
# 
# # Save
# # ggsave(filename = "variance_components.jpg", plot = g_varprop, path = fig_dir, width = 8, height = 4, dpi = 1000)
# ggsave(filename = "variance_components_withAGDD.jpg", plot = g_varprop, path = fig_dir, width = 8, height = 4, dpi = 1000)
# 
# 
# 
# ## Plot just the TP and remove AGDD heading date
# g_varprop_tp <- stage_two_fits_GE_varprop %>% 
#   filter(population == "tp") %>%
#   mutate(population = str_to_upper(population)) %>% 
#   ggplot(aes(x = trait, y = var_prop, fill = source)) + 
#   geom_col() +
#   ylab("Proportion of Variance") +
#   scale_fill_manual(values = umn_palette(3, 4), name = NULL) +
#   theme_acs() +
#   theme(axis.title.x = element_blank(),
#         legend.position = "bottom")
# 
# # Save
# # ggsave(filename = "variance_components_tp.jpg", plot = g_varprop_tp, path = fig_dir, width = 4, height = 4, dpi = 1000)
# ggsave(filename = "variance_components_tp_withAGDD.jpg", plot = g_varprop_tp, path = fig_dir, width = 4, height = 4, dpi = 1000)
# 
# 
# 
# ## Look at the heritability and plot
# g_herit <- stage_two_fits_GE %>% 
#   filter(trait %in% traits) %>%
#   mutate(h2 = map(h2, "heritability"), population = str_to_upper(population)) %>% 
#   unnest(h2) %>% 
#   mutate(statistic = statistic + bias) %>%
#   ggplot(aes(x = trait, y = statistic, ymin = lower, ymax = upper, fill = population)) + 
#   geom_col(position = position_dodge(0.9)) + 
#   geom_errorbar(position = position_dodge(0.9), width = 0.5) +
#   geom_text(aes(y = 0.60, label = round(statistic, 2)), position = position_dodge(0.9)) +
#   # geom_hline(aes(yintercept = unbiased_statistic)) +
#   scale_fill_brewer(palette = "Set2", name = NULL) +
#   ylab("Heritability") +
#   xlab("Trait") +
#   labs(caption = paste0("Error bars represent a 95% confidence interval\ncalculated using ", boot_reps, " bootsrapping replicates.")) + 
#   theme_acs() +
#   theme(legend.position = c(0.20, 0.85), legend.direction = "horizontal")
# 
# ggsave(filename = "heritability1.jpg", plot = g_herit, path = fig_dir, width = 6, height = 4, dpi = 1000)
# # ggsave(filename = "heritability_withADGG1.jpg", plot = g_herit, path = fig_dir, width = 6, height = 4, dpi = 1000)
# 
# 
# 
# # Plot just the TP
# g_herit_tp <- stage_two_fits_GE %>% 
#   filter(population == "tp") %>%
#   mutate(h2 = map(h2, "heritability"), population = str_to_upper(population)) %>% 
#   unnest(h2) %>% 
#   mutate(statistic = statistic + bias) %>%
#   ggplot(aes(x = trait, y = statistic, ymin = lower, ymax = upper, fill = population)) + 
#   geom_col(position = position_dodge(0.9)) + 
#   geom_errorbar(position = position_dodge(0.9), width = 0.5) +
#   geom_text(aes(y = 0.60, label = round(statistic, 2)), position = position_dodge(0.9)) +
#   scale_fill_brewer(palette = "Set2", name = "Population", guide = FALSE) +
#   ylab("Heritability") +
#   xlab("Trait") +
#   labs(caption = paste0("Error bars represent a 95% confidence interval\ncalculated using ", boot_reps, " bootsrapping replicates.")) + 
#   theme_acs() +
#   theme(legend.position = c(0.15, 0.20))
# 
# # ggsave(filename = "heritability_tp.jpg", plot = g_herit_tp, path = fig_dir, width = 4, height = 4, dpi = 1000)
# ggsave(filename = "heritability_tp_withAGDD.jpg", plot = g_herit_tp, path = fig_dir, width = 4, height = 4, dpi = 1000)





### Analyze the components of GxE

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
    fit <- lmer(formula = value ~ 1 + (1|line_name) + (1|location) + (1|year) + (1|location:year) + (1|line_name:location) + 
                  (1|line_name:year) + (1|line_name:location:year),
                data = df, control = lmer_control, weights = wts)
                
    # fit <- lmer(formula = value ~ 1 + (1|line_name) + (1|location) + (1|year) + (1|line_name:location) + (1|line_name:year) + (1|line_name:location:year),
    #             data = df, control = lmer_control)
    
    
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



## Create a table for output
stage_two_fits_GYL_varcomp <- stage_two_fits_GYL %>%
  mutate(var_comp = map(h2, "var_comp"),
         var_comp = map2(var_comp, lrt, ~full_join(.x, .y, by = c("source" = "term")))) %>% 
  unnest(var_comp)


## Plot the proportion of variance from each source
stage_two_fits_GYL_varprop <- stage_two_fits_GYL_varcomp %>% 
  mutate(source = map(source, ~str_split(., pattern = ":", simplify = FALSE) %>% map(str_to_title) %>% .[[1]]) %>% 
           map_chr(~paste(., collapse = " x ")),
         source = str_replace_all(source, "Line_name", "Genotype"),
         source = factor(source, levels = c("Genotype", "Location", "Year", "Location x Year", "Genotype x Location", "Genotype x Year",
                                            "Genotype x Location x Year", "Residual")))

stage_two_fits_GYL_varprop1 <- stage_two_fits_GYL_varprop %>%
  group_by(trait, population) %>% 
  mutate(var_prop = variance / sum(variance)) %>%
  ungroup()

## Colors
var_comp_colors <- set_names(umn_palette(3, 8), levels(stage_two_fits_GYL_varprop$source))

## Plot both populations and all traits
g_varprop <- stage_two_fits_GYL_varprop1 %>% 
  mutate(population = str_to_upper(population)) %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of variance") +
  scale_fill_manual(values = var_comp_colors, name = NULL) +
  facet_grid(~ population) + 
  theme_presentation2(base_size = 10) +
  theme(axis.title.x = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1))

# Save
ggsave(filename = "variance_components_expanded_withAGDD.jpg", plot = g_varprop, path = fig_dir, width = 9, height = 6, dpi = 1000)


## Remove AGDD
g_varprop <- stage_two_fits_GYL_varprop1 %>% 
  filter(trait %in% traits) %>%
  mutate(population = str_to_upper(population)) %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of variance") +
  scale_fill_manual(values = var_comp_colors, name = NULL) +
  facet_grid(~ population) + 
  theme_presentation2(base_size = 10) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

# Save
ggsave(filename = "variance_components_expanded.jpg", plot = g_varprop, path = fig_dir, width = 7, height = 5, dpi = 1000)



#### Heading Date Analysis
#### 
#### Remove potential outlier environments and re-analyze
#### 

# Group by trait and fit the multi-environment model
# Fit models in the TP and the TP + VP
stage_two_fits_HD <- S2_MET_BLUEs_tomodel %>%
  filter(trait == "HeadingDate", ! environment %in% c("EON17", "CRM15", "HNY15")) %>%
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
    fit <- lmer(formula = value ~ 1 + (1|line_name) + (1|location) + (1|year) + (1|location:year) + (1|line_name:location) + 
                  (1|line_name:year) + (1|line_name:location:year),
                data = df, control = lmer_control, weights = wts)
    
    # fit <- lmer(formula = value ~ 1 + (1|line_name) + (1|location) + (1|year) + (1|line_name:location) + (1|line_name:year) + (1|line_name:location:year),
    #             data = df, control = lmer_control)
    
    
    ## Likelihood ratio tests
    lrt <- ranova(fit) %>% 
      tidy() %>% 
      filter(!str_detect(term, "none")) %>% 
      mutate(term = str_remove(term, "\\(1 \\| ") %>% 
               str_remove("\\)")) %>% 
      select(term, LRT, df, p.value)

    # Return data_frame
    data_frame(fit = list(fit), lrt = list(lrt), n_l = harm_loc, 
               n_y = harm_year, n_r = harm_rep)
    
  }) %>% ungroup()

## Create a table for output
stage_two_fits_GYL_varcomp <- stage_two_fits_HD %>%
  mutate(var_comp = map(fit, ~as.data.frame(VarCorr(.)) %>% select(source = grp, variance = vcov)),
         var_comp = map2(var_comp, lrt, ~full_join(.x, .y, by = c("source" = "term")))) %>% 
  unnest(var_comp) %>%
  mutate(source = map(source, ~str_split(., pattern = ":", simplify = FALSE) %>% map(str_to_title) %>% .[[1]]) %>% 
           map_chr(~paste(., collapse = " x ")),
         source = str_replace_all(source, "Line_name", "Genotype"),
         source = factor(source, levels = c("Genotype", "Location", "Year", "Location x Year", "Genotype x Location", "Genotype x Year",
                                            "Genotype x Location x Year", "Residual"))) %>%
  group_by(trait, population) %>% 
  mutate(var_prop = variance / sum(variance)) %>%
  ungroup()



#######
#######
#######










## Look at genetic components
stage_two_fits_GYL_varprop2 <- stage_two_fits_GYL_varprop %>%
  filter(str_detect(source, "Genotype")) %>%
  group_by(trait, population) %>% 
  mutate(var_prop = variance / sum(variance)) %>%
  ungroup()



## Plot both populations and all traits
g_varprop_gen <- stage_two_fits_GYL_varprop2 %>% 
  mutate(population = str_to_upper(population)) %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col() +
  ylab("Proportion of genetic-related variance") +
  scale_fill_manual(values = var_comp_colors, name = NULL) +
  facet_grid(~ population) + 
  theme_acs() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom")

# Save
ggsave(filename = "variance_components_genetic_expanded.jpg", plot = g_varprop_gen, path = fig_dir, width = 8, height = 4, dpi = 1000)





## Look at the heritability and plot
g_herit <- stage_two_fits_GYL %>% 
  mutate(h2 = map(h2, "heritability"), population = str_to_upper(population)) %>% 
  unnest(h2) %>% 
  mutate(statistic = statistic + bias) %>%
  ggplot(aes(x = trait, y = statistic, ymin = lower, ymax = upper, fill = population)) + 
  geom_col(position = position_dodge(0.9)) + 
  geom_errorbar(position = position_dodge(0.9), width = 0.5) +
  geom_text(aes(y = 0.40, label = round(statistic, 2)), position = position_dodge(0.9), size = 3) +
  # geom_hline(aes(yintercept = unbiased_statistic)) +
  scale_fill_brewer(palette = "Set1", name = "Population") +
  ylab("Heritability") +
  xlab("Trait") +
  labs(caption = paste0("Error bars represent a 95% confidence interval\ncalculated using ", boot_reps, " bootsrapping replicates.")) + 
  theme_acs() +
  theme(legend.position = "bottom")

ggsave(filename = "heritability_expanded_withAGDD.jpg", plot = g_herit, path = fig_dir, width = 5, height = 4, dpi = 1000)


# Remove AGDD
g_herit <- stage_two_fits_GYL %>% 
  filter(trait %in% traits) %>%
  mutate(h2 = map(h2, "heritability"), population = str_to_upper(population)) %>% 
  unnest(h2) %>% 
  mutate(statistic = statistic + bias) %>%
  ggplot(aes(x = trait, y = statistic, ymin = lower, ymax = upper, fill = population)) + 
  geom_col(position = position_dodge(0.9)) + 
  geom_errorbar(position = position_dodge(0.9), width = 0.5) +
  geom_text(aes(y = 0.40, label = round(statistic, 2)), position = position_dodge(0.9), size = 3) +
  scale_fill_brewer(palette = "Set1", name = "Population") +
  ylab("Heritability") +
  xlab("Trait") +
  labs(caption = paste0("Error bars represent a 95% confidence interval\ncalculated using ", boot_reps, " bootsrapping replicates.")) + 
  theme_acs() +
  theme(legend.position = "bottom")

ggsave(filename = "heritability_expanded.jpg", plot = g_herit, path = fig_dir, width = 4, height = 4, dpi = 1000)











### Calculate the proportion of GxE that is due to environmental genetic variance
### heterogeneity versus lack of environmental correlation

# For each environment, calculate the genetic variance via reml
env_varG <- S2_MET_BLUEs_tomodel %>% 
  group_by(population, trait, environment) %>%
  do(varG = {
    df <- .
    wts <- df$std_error^2
    
    # If more than one trial is present, average over trials
    if (n_distinct(df$trial) > 1) {
      formula <- value ~ 1 + trial + (1|line_name)
    } else {
      formula <- value ~ 1 + (1|line_name)
    }
    
    fit <- lmer(formula = formula, data = df, control = lmer_control, weights = wts)  
    
    as.data.frame(VarCorr(fit))[1,"vcov"]
    
  }) %>% ungroup() %>%
  mutate(varG = unlist(varG))

# Calculate the genetic heterogeneity, V
env_varG_V <- env_varG %>% 
  group_by(population, trait) %>% 
  summarize(V = var(sqrt(varG))) %>%
  ungroup()


## Use the variance components estinated in the previous random model
# prop_varcomp <- stage_two_fits_GE %>%
prop_varcomp <- stage_two_fits_GYL %>%
  mutate(varcomp = map(h2, "var_comp")) %>% 
  unnest(varcomp)


# Use the estimate of varGE across all environments to calculate L
env_L <- env_varG_V %>%
  # left_join(., subset(prop_varcomp, source == "line_name:environment", c(population, trait, variance))) %>% 
  left_join(., subset(prop_varcomp, source == "line_name:location:year", c(population, trait, variance))) %>% 
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

# # A tibble: 4 x 4
# trait             all    tp    vp
# 1 GrainYield      0.303 0.308 0.222
# 2 HeadingDate     0.642 0.733 0.490
# 3 HeadingDateAGDD 0.647 0.725 0.485
# 4 PlantHeight     0.360 0.363 0.387


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
# 1 all        GrainYield      22853.898 (9%) 233134.396 (91%)
# 2 all        HeadingDate     1.057 (15%)    6.127 (85%)     
# 3 all        HeadingDateAGDD 1301.759 (15%) 7646.485 (85%)  
# 4 all        PlantHeight     2.034 (11%)    16.958 (89%)    
# 5 tp         GrainYield      19983.738 (8%) 239570.652 (92%)
# 6 tp         HeadingDate     1.271 (21%)    4.78 (79%)      
# 7 tp         HeadingDateAGDD 1444.957 (19%) 6300.711 (81%)  
# 8 tp         PlantHeight     2 (11%)        16.474 (89%)    
# 9 vp         GrainYield      19347.096 (9%) 190971.903 (91%)
# 10 vp         HeadingDate     1.181 (22%)    4.135 (78%)     
# 11 vp         HeadingDateAGDD 1205.73 (18%)  5435.906 (82%)  
# 12 vp         PlantHeight     2.41 (13%)     15.976 (87%)


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

# ggsave(filename = "varGE_components.jpg", plot = g_varGE_comp, width = 8, height = 3, path = fig_dir, dpi = 1000)
ggsave(filename = "varGE_components_withAGDD.jpg", plot = g_varGE_comp, width = 8, height = 3, path = fig_dir, dpi = 1000)



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
comp_order <- unique(stage_two_fits_GYL_varprop1$source) %>% 
  {.[order(str_count(., "x"))]} %>%
  {c(str_subset(., "[^Residual]"), "Genetic Heterogeneity", "Lack Of Correlation", "Residual")}

## Combine and output
var_comp_table <- select(stage_two_fits_GYL_varprop1, trait, population, source, variance, proportion = var_prop, p.value) %>%
  bind_rows(., varGE_components1) %>%
  ungroup() %>%
  mutate(significance = case_when(p.value < 0.001 ~ "***",
                                  p.value < 0.01 ~ "**",
                                  p.value < 0.05 ~ "*",
                                  TRUE ~ ""),
         variance = signif(variance, 3) %>% formatC(x = ., digits = 3, format = "f") %>% str_remove(., "[0]*$") %>% str_remove(., "\\.$"),
         proportion = signif(proportion * 100, 3) %>% formatC(x = ., digits = 3, format = "f") %>% str_remove(., "[0]*$") %>% str_remove(., "\\.$"),
         annotation = paste0(variance, significance, " (", proportion, "%)"),
         annotation = str_trim(annotation)) %>%
  select(trait, population, source, annotation) %>%
  mutate(source = str_add_space(source),
         trait = str_add_space(trait),
         source = factor(source, levels = comp_order),
         population = toupper(population)) %>% 
  arrange(population, source) %>%
  rename_all(str_to_title) %>%
  spread(Population, Annotation)


write_csv(x = var_comp_table, path = file.path(fig_dir, "population_variance_components_decomposed_withAGDD.csv"))

# Remove AGDD
write_csv(x = filter(var_comp_table, Trait %in% str_add_space(traits)), path = file.path(fig_dir, "population_variance_components_decomposed.csv"))


## Just the TP
var_comp_table %>% 
  select(Trait, Source, TP) %>% 
  spread(Trait, TP) %>% 
  write_csv(x = ., path = file.path(fig_dir, "tp_population_variance_components_decomposed_withAGDD.csv"))

var_comp_table %>% 
  filter(Trait %in% str_add_space(traits)) %>%
  select(Trait, Source, TP) %>% 
  spread(Trait, TP) %>% 
  write_csv(x = ., path = file.path(fig_dir, "tp_population_variance_components_decomposed.csv"))




##






## Environmental Correlations

## Estimate genetic correlations using the BLUEs from each environment
## First estimate using only the TP with genotypes
tp_BLUEs <- S2_MET_BLUEs_use %>%
  filter(line_name %in% tp) # %>% filter(trait %in% traits)

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
env_cor_vp <- S2_MET_BLUEs_use %>%
  filter(line_name %in% vp) %>%
  # filter(trait %in% traits) %>%
  select(line_name, environment, trait, value) %>% 
  split(.$trait) %>% 
  map(~select(., line_name, environment, value) %>% spread(environment, value) %>% 
        as.data.frame() %>% remove_rownames() %>% column_to_rownames("line_name") %>% 
        cor(., use = "pairwise.complete.obs"))


# Now use all data
env_cor_all <- S2_MET_BLUEs_use %>%
  # filter(trait %in% traits) %>%
  split(.$trait) %>% 
  map(~select(., line_name, environment, value) %>% spread(environment, value) %>% 
        as.data.frame() %>% remove_rownames() %>% column_to_rownames("line_name") %>% 
        cor(., use = "pairwise.complete.obs"))


# Convert to a data.frame for visualization 
env_cor_df <- list(all = env_cor_all, train = env_cor_tp, test = env_cor_vp) %>% 
  map(~{
    map(., ~as.data.frame(.) %>% rownames_to_column("environment1") %>% 
          gather(environment2, correlation, -environment1) %>% 
          filter(environment1 != environment2)) %>%
      list(., names(.)) %>%
      pmap_df(~mutate(.x, trait = .y))
  }) %>%
  map2_df(.x = ., names(.), ~mutate(.x, population = .y)) %>%
  mutate(population = factor(population, levels = c("all", "train", "test")))




# Plot
## Colors
heat_colors <- wesanderson::wes_palette("Zissou1")

# One trait at a time
g_env_cor_df <- env_cor_df %>%
  group_by(trait) %>%
  do(plot = {
  
    df <- .
    
    # Use the heatmap function to construct a dendrogram
    # Use only data from the whole population
    mat <- df %>% 
      filter(population == "all") %>%
      select(-trait, -population) %>% 
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
      scale_fill_gradientn(colors = heat_colors[c(1,3,5)], na.value = "transparent", 
                           breaks = seq(-1, 1, by = 0.5), limits = c(-1, 1), 
                           guide = guide_colorbar(title = "Phenotypic\ncorrelation")) +
      ylab("Environment 2") +
      xlab("Environment 1") +
      facet_wrap(~ trait + population, labeller = labeller(trait = str_add_space, population = str_to_title)) +
      theme_presentation2(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
            axis.text.y = element_text(size = 5),
            # legend.position = "bottom", legend.direction = "horizontal")
            legend.position = "right", legend.direction = "vertical")
    
    
  }) %>% ungroup()




## Combine plots
g_env_cor <- plot_grid(plotlist = subset(g_env_cor_df, trait %in% traits, plot, drop = T) %>% 
                         map(~. + theme(legend.position = "none", axis.title = element_blank())), 
                       nrow = length(subset(g_env_cor_df, trait != "HeadingDateAGDD", trait, drop = T)), scale = 0.9)

g_env_cor1 <- g_env_cor + 
  draw_label("Environment 1", x = 0.5, y = 0, vjust=-0.5, angle= 0) +
  draw_label("Environment 2", x = 0, y = 0.5, vjust= 1.5, angle=90)

g_env_cor2 <- plot_grid(g_env_cor1, get_legend(g_env_cor_df$plot[[1]]), nrow = 1, rel_widths = c(1, 0.1))


# Save this plot
ggsave(filename = "environmental_correlation.jpg", path = fig_dir, plot = g_env_cor2, height = 8, width = 8, dpi = 1000)







## Included heading date AGDD
g_env_cor <- plot_grid(plotlist = g_env_cor_df$plot %>% map(~. + theme(legend.position = "none", axis.title = element_blank())), 
                       nrow = n_distinct(g_env_cor_df$trait), scale = 0.9)

g_env_cor1 <- g_env_cor + 
  draw_label("Environment 1", x = 0.5, y = 0, vjust=-0.5, angle= 0) +
  draw_label("Environment 2", x = 0, y = 0.5, vjust= 1.5, angle=90)

g_env_cor2 <- plot_grid(g_env_cor1, get_legend(g_env_cor_list[[1]]), nrow = 1, rel_widths = c(1, 0.1))

ggsave(filename = "environmental_correlation_withAGDD.jpg", path = fig_dir, plot = g_env_cor2, height = 12, width = 9, dpi = 1000)



#######


## Are the correlation matrices significantly different between populations?
env_cor_mantel <- map2(.x = env_cor_tp, .y = env_cor_vp, ~ade4::mantel.rtest(m1 = dist(.x), m2 = dist(.y), nrepet = 5000))
# P-value is for null hypothesis that the two matrices are unrelated
map_dbl(env_cor_mantel, "pvalue")

## Are the phenotypic distance matrices significantly different between populations?
# First calculate the distance based on each environment
ge_mean_D <- S2_MET_BLUEs %>%
  filter(trait %in% traits) %>%
  mutate(population = ifelse(line_name %in% tp, "TP", "VP")) %>%
  # Split by trait
split(.$trait) %>%
  map2(.x = ., .y = tp_vp_env_trait[names(.)], ~filter(.x, environment %in% .y)) %>%
  map(~split(., .$population)) %>%
  map(~map(., function(df){
    dist1 <- dist_env(x = df, gen.col = "line_name", env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 }))


env_PD_mantel <- map(ge_mean_D, ~ade4::mantel.rtest(m1 = .[[1]], m2 = .[[2]], nrepet = 5000))
# P-value is for null hypothesis that the two matrices are unrelated
map_dbl(env_PD_mantel, "pvalue")



## Simulate two populations with the same correlation between environments
nEnv <- 10
nPop <- c(200, 200)
covDist <- as.dist(matrix(0, nrow = nEnv, ncol = nEnv))

pvalueNull <- replicate(100, {
  covDist <- as.dist(matrix(0, nrow = nEnv, ncol = nEnv))
  covDist[seq_along(covDist)] <- runif(n = length(covDist), min = 0.3, max = 1)
  covMat <- as.matrix(covDist); diag(covMat) <- 1
  
  # Simulate genotypic values from a multi-variate distribution
  GList <- map(nPop, ~mvtnorm::rmvnorm(n = ., mean = rep(0, nrow(covMat)), sigma = covMat)) %>%
    map(cor) %>%
    map(dist)
  
  # Test using mantel
  mantel_out <- ade4::mantel.rtest(m1 = GList[[1]], m2 = GList[[2]], nrepet = 500)
  mantel_out$pvalue
})

## True discovery rate
mean(pvalueNull <= 0.05)

## Simulate two populations with different correlations between environments
pvalueAlt <- replicate(100, {
  covDist <- replicate(2, as.dist(matrix(0, nrow = nEnv, ncol = nEnv)), simplify = FALSE) %>%
    map(function(d) {d[seq_along(d)] <- runif(n = length(d), min = 0.3, max = 1); d})
  covMatList <- map(covDist, ~as.matrix(.) %>% `diag<-`(., 1))
  
  # Simulate genotypic values from a multi-variate distribution
  GList <- map2(nPop, covMatList, ~mvtnorm::rmvnorm(n = .x, mean = rep(0, nrow(.y)), sigma = .y)) %>%
    map(cor) %>%
    map(dist)
  
  # Test using mantel
  mantel_out <- ade4::mantel.rtest(m1 = GList[[1]], m2 = GList[[2]], nrepet = 500)
  mantel_out$pvalue
})

## FDR?
mean(pvalueAlt <= 0.05)


## What about smaller populations? Different sizes? More environments?
nEnv <- 25
nPop <- c(175, 50)
covDist <- as.dist(matrix(0, nrow = nEnv, ncol = nEnv))

pvalueNull <- replicate(100, {
  covDist <- as.dist(matrix(0, nrow = nEnv, ncol = nEnv))
  covDist[seq_along(covDist)] <- runif(n = length(covDist), min = 0.3, max = 1)
  covMat <- as.matrix(covDist); diag(covMat) <- 1
  
  # Simulate genotypic values from a multi-variate distribution
  GList <- map(nPop, ~mvtnorm::rmvnorm(n = ., mean = rep(0, nrow(covMat)), sigma = covMat)) %>%
    map(cor) %>%
    map(dist)
  
  # Test using mantel
  mantel_out <- ade4::mantel.rtest(m1 = GList[[1]], m2 = GList[[2]], nrepet = 500)
  mantel_out$pvalue
})

## True discovery rate
mean(pvalueNull <= 0.05)

## Simulate two populations with different correlations between environments
pvalueAlt <- replicate(100, {
  covDist <- replicate(2, as.dist(matrix(0, nrow = nEnv, ncol = nEnv)), simplify = FALSE) %>%
    map(function(d) {d[seq_along(d)] <- runif(n = length(d), min = 0.3, max = 1); d})
  covMatList <- map(covDist, ~as.matrix(.) %>% `diag<-`(., 1))
  
  # Simulate genotypic values from a multi-variate distribution
  GList <- map2(nPop, covMatList, ~mvtnorm::rmvnorm(n = .x, mean = rep(0, nrow(.y)), sigma = .y)) %>%
    map(cor) %>%
    map(dist)
  
  # Test using mantel
  mantel_out <- ade4::mantel.rtest(m1 = GList[[1]], m2 = GList[[2]], nrepet = 500)
  mantel_out$pvalue
})

## FDR?
mean(pvalueAlt <= 0.05)





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
      scale_fill_gradientn(colors = heat_colors[c(1,3,5)], na.value = "transparent", 
                           breaks = seq(-1, 1, by = 0.5), limits = c(-1, 1), 
                           guide = guide_colorbar(title = "Correlation")) +
      ylab("Environment 2") +
      xlab("Environment 1") +
      facet_wrap(~ trait) +
      theme_presentation2(base_size = 10) +
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
  arrange(population, trait)




# population trait           min   max  mean
# 1 all        GrainYield  -0.148  0.572 0.253
# 2 all        HeadingDate  0.134  0.878 0.636
# 3 all        PlantHeight -0.0473 0.659 0.330
# 4 tp         GrainYield  -0.171  0.615 0.259
# 5 tp         HeadingDate  0.366  0.892 0.707
# 6 tp         PlantHeight -0.0243 0.617 0.333
# 7 vp         GrainYield  -0.264  0.595 0.180
# 8 vp         HeadingDate  0.0683 0.825 0.459
# 9 vp         PlantHeight -0.260  0.771 0.351


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
tp_BLUEs <- S2_MET_BLUEs_use %>%
  filter(line_name %in% tp_geno) 

tp_BLUEs_split <- tp_BLUEs %>%
  split(.$trait) %>%
  map(~mutate_at(., vars(line_name, environment), as.factor))

# Subset the K matrix
G <- K[tp_geno, tp_geno]


data_to_model <- full_join(
  data_frame(trait = names(tp_BLUEs_split), phenos = tp_BLUEs_split),
  data_frame(trait = names(env_cor_tp), env_cor = env_cor_tp)
)
  
  
# Iterate and fit models
tp_twoway_fit <- data_to_model %>%
  group_by(trait) %>%
  do({
    
    phenos <- .$phenos[[1]]
    E <- .$env_cor[[1]]
    
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
    
    
    ## Empty list of model fits
    fit_list <- set_names(vector("list", length = n_distinct(phenos$year)), unique(phenos$year)) 
    
    ## Drop each year and re-fit
    for (yr in unique(phenos$year)) {
      
      ## Remove data from that year and re-fit
      # Model frame
      mf <- model.frame(value ~ line_name + environment, data = droplevels(filter(phenos, year != yr)))
      y <- model.response(mf)
      
      # Fixed effects of environments
      X <- model.matrix(~ 1 + environment + line_name, data = droplevels(mf), contrasts.arg = list(line_name = "contr.sum", environment = "contr.sum"))
      Z <- model.matrix(~ -1 + line_name:environment, data = mf)
      
      # Create the kronecker product
      Kge <- kronecker(E[levels(mf$environment), levels(mf$environment)], G, make.dimnames = TRUE)
      
      fit_list[[as.character(yr)]] <- mixed.solve(y = y, Z = Z, K = Kge, X = X)
      
    }
    
    ## Create a df
    fit_df <- data_frame(set = c("complete", paste0("realistic", names(fit_list))), model = c(list(fit_complete), fit_list))
    
    mf_list <- c(list(phenos), map(names(fit_list), ~droplevels(filter(phenos, year != .))))
    
    
    
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
    
  }) %>% ungroup()



## Create two datasets - one using all environments and one without 2017
training_sets_twoway <- tp_twoway_fit %>% 
  ungroup() %>%
  rename(data = predictions)


## Perform the AMMI analysis
# Create a two-way matrix
ammi_out <- training_sets_twoway %>%
  group_by(set, trait) %>%
  do(ammi = {
    df <- .$data[[1]]
    # varE <- .$varE
    
    ## Use the bi-linear package
    fit <- Bilinear::bilinear(x = df, G = "line_name", E = "environment", y = "value", model = "AMMI", alpha = 0.05, 
                              test = "bootstrap", B = 1000)
    
    ## Create escores and gscores
    gscores <- fit$scores$Gscores
    escores <- fit$scores$Escores
    
    ## Create df of eigenvalues
    geigen <- fit$svdE$u[,seq(1, fit$sigPC), drop = FALSE]
    dimnames(geigen) <- dimnames(gscores)
    geigen_df <- rownames_to_column(as.data.frame(geigen), "line_name") %>% gather(PC, eigen, -line_name)
    
    eeigen <- fit$svdE$v[,seq(1, fit$sigPC), drop = FALSE]
    dimnames(eeigen) <- dimnames(escores)
    eeigen_df <- rownames_to_column(as.data.frame(eeigen), "environment") %>% gather(PC, eigen, -environment)
    
    
    gscores_df <- left_join(rownames_to_column(data.frame(effect = fit$Geffect), "line_name"), rownames_to_column(as.data.frame(gscores), "line_name")) %>%
      gather(PC, score, -line_name, -effect) %>%
      left_join(., geigen_df)
    escores_df <- left_join(rownames_to_column(data.frame(effect = fit$Eeffect), "environment"), rownames_to_column(as.data.frame(escores), "environment")) %>%
      gather(PC, score, -environment, -effect) %>%
      left_join(., eeigen_df)
    
    # Return
    map(list(results = rownames_to_column(fit$ANOVA, "term"), escores = escores_df, gscores = gscores_df), as_data_frame)
    
    
    
    
    
    # ## Fit a model
    # fit <- aov(value ~ line_name + environment, data = df)
    # # Get the effects
    # effects <- model.tables(fit)
    # effects_df <- effects$tables %>% 
    #   map(~as.matrix(.) %>% fix_data_frame(newnames = "effect"))
    # 
    # # Add the residuals (GxE effects) to the df
    # df1 <- add_residuals(df, fit, var = "ge")
    # 
    # # Create a two-way table of GxE + residuals
    # ge <- xtabs(ge ~ line_name + environment, df1)
    # 
    # ## SVD
    # ge_svd <- svd(ge)
    # # Get the d values for each factor
    # dm <- svd(ge)$d
    # 
    # ## Calculate sums of squares for each factor
    # SSd <- dm^2
    # # Calculate the proportion of variance explained by each factor
    # varprop <- SSd / sum(SSd)
    # # Cumulative proportion
    # cumprop <- cumsum(varprop)
    # 
    # # Calculate the degrees of freedom for each factor
    # J <- nrow(ge)
    # K <- ncol(ge)
    # df_d <- (J + K - 1 - (2 * seq_along(SSd)))
    # 
    # # Calculate the mean squares
    # MSd <- SSd / df_d
    # 
    # ## Use the residual variance from the previous model as the expected mean squares of the residuals
    # ## Calculate an F statistic
    # stat <- MSd / varE
    # 
    # # Calculate p-values
    # p_value_Ftest <- pf(q = stat, df1 = df_d, df2 = (J * K), lower.tail = FALSE)
    # 
    # 
    # 
    # # ## Use resampling to determine the null distribution of MSf / MSe
    # # # First fit a main effects model, then simulate from that model to calculate ge terms
    # # fit <- aov(value ~ line_name + environment, df)
    # # mf <- model.frame(fit)
    # # 
    # # # Simulate from the model
    # # sims <- simulate(fit, nsim = 1000)
    # # # Update the model
    # # fit_sims <- apply(X = sims, MARGIN = 2, FUN = function(data) {
    # #   mf$value <- data
    # #   fit1 <- lm(value ~ line_name + environment, mf)
    # #   
    # #   # Get the residuals / GE
    # #   ge1 <- matrix(data = residuals(fit1), nrow = J, ncol = K)
    # #   svd(ge1)$d
    # # })
    # # 
    # # ## Use the matrix to determine the significance of each factor
    # # SS_sims <- fit_sims^2
    # # p_value_sim <- map_dbl(seq_along(SSd), ~mean(SS_sims[.,] >= SSd[.]))
    # 
    # p_value_sim <- rep(NA, length(SSd))
    # 
    # # Create a table of sums of squares
    # out <- data.frame(term = c("ge", str_c("PC", seq_along(SSd)), "Residuals"),
    #            sum_squares = c(sum(ge^2), SSd, varE * (J * K)),
    #            prop_var = c(NA, varprop, NA),
    #            cumprop = c(NA, cumprop, NA),
    #            F_value = c(NA, stat, NA),
    #            p_value_F = c(NA, p_value_Ftest, NA),
    #            p_value_sim = c(NA, p_value_sim, NA))
    # 
    # ## Select the significant environmental scores based on the F-test
    # which_scores <- which(p_value_Ftest < alpha)
    # 
    # escores <- ge_svd$v[,which_scores, drop = FALSE]
    # dimnames(escores) <- list(colnames(ge), paste0("PC", seq(ncol(escores))))
    # 
    # escores_df <- effects_df$environment %>%
    #   rename(environment = term) %>%
    #   left_join(., fix_data_frame(escores, newcol = "environment"), by = "environment") %>%
    #   gather(PC, eigen, -environment, -effect) %>%
    #   left_join(., data_frame(PC = paste0("PC", seq_along(dm)), d = dm), by = "PC") %>%
    #   mutate(score = eigen * sqrt(d)) %>%
    #   select(-d)
    # 
    # gscores <- ge_svd$u[,which_scores, drop = FALSE]
    # dimnames(gscores) <- list(rownames(ge), paste0("PC", seq(ncol(gscores))))
    # 
    # gscores_df <- effects_df$line_name %>%
    #   rename(line_name = term) %>%
    #   left_join(., fix_data_frame(gscores, newcol = "line_name"), by = "line_name") %>%
    #   gather(PC, eigen, -line_name, -effect) %>%
    #   left_join(., data_frame(PC = paste0("PC", seq_along(dm)), d = dm), by = "PC") %>%
    #   mutate(score = eigen * sqrt(d)) %>%
    #   select(-d)
    # 
    # # Return
    # list(results = out, escores = escores_df, gscores = gscores_df) %>% map(as_data_frame)
    
  }) %>% ungroup()
 

## What total proportion of variation do the significant AMMI axes explain?
ammi_out %>% 
  mutate(results = map(ammi, ~mutate(.$results, maxPC = tail(.$escores$PC, 1))),
         results = map(results, ~filter(., str_detect(term, "PC"))),
         results = map(results, ~mutate(., cumprop = SS / sum(SS))),
         results = map(results, ~filter(., parse_number(term) <= parse_number(maxPC)))) %>%
  unnest(results) %>%
  group_by(set, trait) %>%
  summarize_at(vars(cumprop), sum)

## Way too much

## What proportion of variation do the first IPCA scores explain?
ammi_out1 <- ammi_out %>% 
  mutate(results = map(ammi, ~mutate(.$results, maxPC = tail(.$escores$PC, 1))),
         results = map(results, ~filter(., str_detect(term, "PC"))),
         results = map(results, ~mutate(., cumprop = SS / sum(SS))),
         results = map(results, ~filter(., parse_number(term) <= parse_number(maxPC))))

# set       trait       term     Df         SS        MS Pvalue  ` `   maxPC cumprop
# 1 complete      GrainYield      PC1     195 146319165. 750355.   < 0.001 ***   PC20    0.508
# 2 complete      HeadingDate     PC1     198      2727.     13.8  < 0.001 ***   PC22    0.222
# 3 complete      HeadingDateAGDD PC1     198   3752348.  18951.   < 0.001 ***   PC21    0.243
# 4 complete      PlantHeight     PC1     199      7943.     39.9  < 0.001 ***   PC21    0.542
# 5 realistic2015 GrainYield      PC1     192 144361163. 751881.   < 0.001 ***   PC17    0.565
# 6 realistic2015 HeadingDate     PC1     195      2440.     12.5  < 0.001 ***   PC20    0.258
# 7 realistic2015 HeadingDateAGDD PC1     195   2874885.  14743.   < 0.001 ***   PC19    0.245
# 8 realistic2015 PlantHeight     PC1     196      6713.     34.3  < 0.001 ***   PC19    0.545
# 9 realistic2016 GrainYield      PC1     185  62145501. 335922.   < 0.001 ***   PC10    0.516
# 10 realistic2016 HeadingDate     PC1     186      2067.     11.1  < 0.001 ***   PC11    0.307
# 11 realistic2016 HeadingDateAGDD PC1     186   2951084.  15866.   < 0.001 ***   PC11    0.343
# 12 realistic2016 PlantHeight     PC1     188      4266.     22.7  < 0.001 ***   PC13    0.577
# 13 realistic2017 GrainYield      PC1     185  79435444. 429381.   < 0.001 ***   PC8     0.546
# 14 realistic2017 HeadingDate     PC1     187      1689.      9.03 < 0.001 ***   PC11    0.275
# 15 realistic2017 HeadingDateAGDD PC1     187   2295333.  12275.   < 0.001 ***   PC13    0.292
# 16 realistic2017 PlantHeight     PC1     186      3200.     17.2  < 0.001 ***   PC12    0.593


## Use fitted values of AMMI model to cluster environments based on the best lines in each
ammi_fitted <- ammi_out %>%
  group_by(set, trait) %>% 
  do({
    df <- .
    
    ## Filter for PC1
    scores <- df$ammi[[1]][c("escores", "gscores")] %>%
      map(filter, PC == "PC1")
    
    # Create matrices of genotype and environment effects
    g_effect <- scores$gscores %>% 
      select(line_name, effect) %>% 
      as.data.frame() %>% 
      column_to_rownames("line_name") %>% 
      as.matrix()
    
    e_effect <- scores$escores %>% 
      select(environment, effect) %>% 
      as.data.frame() %>% 
      column_to_rownames("environment") %>% 
      as.matrix()
    
    
    g_effect1 <- replicate(nrow(e_effect), g_effect[,1])
    e_effect1 <- replicate(nrow(g_effect), e_effect[,1])
    
    main_effect <- g_effect1 + t(e_effect1)
    colnames(main_effect) <- row.names(e_effect1)
    
    # Create interaction scores from the PC
    g_scores <- scores$gscores %>% 
      select(line_name, score) %>%
      as.data.frame() %>% 
      column_to_rownames("line_name") %>% 
      as.matrix()
    
    e_scores <- scores$escores %>% 
      select(environment, score) %>% 
      as.data.frame() %>% 
      column_to_rownames("environment") %>% 
      as.matrix()
    
    int_effect <- tcrossprod(g_scores, e_scores)
    
    ## Calculate fitted values
    fitted_effect <- main_effect + int_effect
    
    ## Determine the best (and worst) in each environment
    top <- apply(X = fitted_effect, MARGIN = 2, FUN = function(x) row.names(fitted_effect)[order(x, decreasing = TRUE)[1:5]])
    bottom <- apply(X = fitted_effect, MARGIN = 2, FUN = function(x) row.names(fitted_effect)[order(x, decreasing = FALSE)[1:5]])
    
    ## Return the fitted values, the top, and bottom
    fitted_effect1 <- as.data.frame(fitted_effect) %>% rownames_to_column("line_name") %>% gather(environment, fitted, -line_name)
    tibble(fitted = list(fitted_effect1), top = list(top), bottom = list(bottom))
    
  })


## Assign clusters
nTop <- 1
ammi_fitted_clusters <- ammi_fitted %>% 
  mutate(favorable = ifelse(trait == "GrainYield", top, bottom),
         clusters = map(favorable, ~.[seq_len(nTop),,drop = F] %>% apply(X = ., MARGIN = 2, FUN = paste, collapse = "/") %>% 
                          as.data.frame() %>% rownames_to_column("environment") %>% rename_at(vars(-environment), ~"line") %>% 
                          mutate(cluster = as.numeric(as.factor(line))) %>% arrange(cluster) ))





year_levels <- c(sort(unique(S2_MET_BLUEs$year)), "Genotype")

## Year color scheme
year_color <- setNames( c(neyhart_palette("umn1", 5)[3:5], "grey"), year_levels)



# Rename an object
ammi_out <- ammi_out1

## Plot the AMMI axes and the proportion of variance explained
g_ammi <- ammi_out %>%
  left_join(., ammi_fitted_clusters) %>%
  # filter(set == "complete") %>%
  group_by(set, trait) %>%
  do(plot = {
    df <- .
  
    out <- df$ammi[[1]]
    res <- df$results[[1]]
    trait <- df$trait
    
    ## Combine the g and e scores into one DF
    scores_df <- map_df(out[c("escores", "gscores")], ~mutate(., group = names(.)[1]) %>% `names<-`(., c("term", names(.)[-1])))
    
    # Create a renaming vector
    res1 <- mutate(res, annotation = paste0(str_replace_all(term, "PC", "IPCA"), " (", round(cumprop, 2) * 100, "%)"))
    propvar_rename <- setNames(object = res1$term, nm = res1$annotation)
      
    # Renaming function
    label_rename <- function(x) str_replace_all(x, propvar_rename)
    
    scores_toplot <- scores_df %>%
      select(-eigen) %>%
      spread(PC, score) %>%
      # Add year
      mutate(year = ifelse(group == "environment", paste0("20", str_extract(term, "[0-9]{2}")), "Genotype"),
             trait = trait)
  
    
    ## Figure out the range for x and y axis for the trait
    effect_breaks <- subset(ammi_out, trait == df$trait, ammi, drop = T) %>% 
      map(~.[-1]) %>% map(., ~map(., "effect") %>% unlist()) %>% 
      unlist() %>% range()
    
    scores_breaks <- subset(ammi_out, trait == df$trait, ammi, drop = T) %>% 
      map(~.[-1]) %>% map(., ~map(., "score") %>% unlist()) %>% 
      unlist() %>% range()
    
    ## Units
    unit_use <- traits_unit1[df$trait]
    
    # Plot
    scores_toplot %>% 
      mutate(trait_use = paste0(str_add_space(trait), " (", unit_use, ")"),
             trait_use = str_replace_all(trait_use, " ", "~")) %>%
      ggplot(aes(x = effect, y = PC1, color = year)) + 
      geom_point(size = 0.75) +
      geom_text_repel(data = filter(scores_toplot, group == "environment"), aes(label = term), size = 2) +
      scale_color_manual(values = year_color, name = NULL, guide = guide_legend(override.aes = list(size = 3))) + 
      scale_y_continuous(breaks = pretty, limits = range(scores_breaks)) +
      scale_x_continuous(breaks = pretty, limits = range(effect_breaks)) +
      ylab(names(propvar_rename)[1]) +
      xlab(expression("Effect (deviation from"~mu*")")) +
      facet_grid(~ trait_use, labeller = labeller(trait_use = label_parsed)) + 
      theme_presentation2(base_size = 10) +
      theme(legend.position = "bottom")
      # theme(legend.position = "right", legend.direction = "vertical")
    
  })


# Plot in a grid
plot_list_complete <- filter(g_ammi, set == "complete")$plot
# ammi_grid <- plot_grid(plotlist = plot_list %>% map(~. + theme(legend.position = "none")), ncol = 1)
ammi_grid_complete <- plot_grid(ncol = 1, rel_heights = c(rep(0.9, length(plot_list_complete) - 1), 1), 
                                plotlist = c(head(plot_list_complete, -1) %>% map(~. + theme(legend.position = "none", axis.title.x = element_blank())),
                                             list(tail(plot_list_complete, 1)[[1]] + theme(legend.position = "none"))))

# Add legend
ammi_grid1 <- plot_grid(ammi_grid_complete, get_legend(plot_list_complete[[1]]), ncol = 1, rel_heights = c(1, 0.05))

# Save
ggsave(filename = "ammi_biplots_complete_withAGDD.jpg", plot = ammi_grid1, path = fig_dir, height = 8, width = 3, dpi = 1000)
  

## Remove AGDD
plot_list_complete <- subset(g_ammi, set == "complete" & trait %in% traits, plot, drop = T)
# ammi_grid <- plot_grid(plotlist = plot_list %>% map(~. + theme(legend.position = "none")), ncol = 1)
ammi_grid_complete <- plot_grid(ncol = 1, rel_heights = c(rep(0.9, length(plot_list_complete) - 1), 1), 
                                plotlist = c(head(plot_list_complete, -1) %>% map(~. + theme(legend.position = "none", axis.title.x = element_blank())),
                                             list(tail(plot_list_complete, 1)[[1]] + theme(legend.position = "none"))))

# Add legend
ammi_grid1 <- plot_grid(ammi_grid_complete, get_legend(plot_list_complete[[1]]), ncol = 1, rel_heights = c(1, 0.05))

# Save
ggsave(filename = "ammi_biplots_complete.jpg", plot = ammi_grid1, path = fig_dir, height = 6.5, width = 3, dpi = 1000)



### Plot clusters as colors, years as different shapes; keep genotype constant


## Year shape scheme
year_shape <- setNames(c(15, 17, 18, 16), names(year_color))


## Plot the AMMI axes and the proportion of variance explained
g_ammi <- ammi_out %>%
  left_join(., ammi_fitted_clusters) %>%
  # filter(set == "complete") %>%
  group_by(set, trait) %>%
  do(plot = {
    df <- .
    
    out <- df$ammi[[1]]
    res <- df$results[[1]]
    trait <- df$trait
    clus <- df$clusters[[1]] %>%
      ## Assign cluster number in descending order of size
      group_by(cluster) %>% 
      mutate(nEnv = n()) %>% 
      ungroup() %>% 
      mutate(cluster2 = factor(nEnv, levels = sort(unique(nEnv), decreasing = TRUE)), 
             cluster2 = as.numeric(cluster2))
    
    ## Combine the g and e scores into one DF
    scores_df <- map_df(out[c("escores", "gscores")], ~mutate(., group = names(.)[1]) %>% `names<-`(., c("term", names(.)[-1])))
    
    # Create a renaming vector
    res1 <- mutate(res, annotation = paste0(str_replace_all(term, "PC", "IPCA"), " (", round(cumprop, 2) * 100, "%)"))
    propvar_rename <- setNames(object = res1$term, nm = res1$annotation)
    
    # Renaming function
    label_rename <- function(x) str_replace_all(x, propvar_rename)
    
    scores_toplot <- scores_df %>%
      select(-eigen) %>%
      spread(PC, score) %>%
      # Add year
      mutate(year = ifelse(group == "environment", paste0("20", str_extract(term, "[0-9]{2}")), "Genotype"),
             trait = trait) %>%
      # Add clusters
      left_join(., select(clus, term = environment, cluster = cluster2), by = "term") %>%
      mutate(cluster = ifelse(is.na(cluster), 0, cluster),
             cluster = as.factor(cluster))
    

    ## Vector for cluster colors
    cluster_colors <- setNames(c("grey", neyhart_palette("umn1")[-1:-2][seq(n_distinct(scores_toplot$cluster) - 1)]), 
                                sort(unique(scores_toplot$cluster)))
    
    
    ## Figure out the range for x and y axis for the trait
    effect_breaks <- subset(ammi_out, trait == df$trait, ammi, drop = T) %>% 
      map(~.[-1]) %>% map(., ~map(., "effect") %>% unlist()) %>% 
      unlist() %>% range()
    
    scores_breaks <- subset(ammi_out, trait == df$trait, ammi, drop = T) %>% 
      map(~.[-1]) %>% map(., ~map(., "score") %>% unlist()) %>% 
      unlist() %>% range()
    
    ## Units
    unit_use <- traits_unit1[df$trait]
    
    ## Annotation specs for labelling cluster
    cluster_annotation_df <- scores_toplot %>% 
      distinct(cluster) %>%
      filter(cluster != "0") %>%
      mutate(clusterN = as.numeric(cluster) - 1,
             label = paste0("Cluster ", cluster),
             x = Inf, y = -Inf, hjust = 1.2,
             vjust = rev(clusterN) * -1.1,
             color = cluster_colors[-1][clusterN])
    
    # Plot
    scores_toplot %>% 
      mutate(trait_use = paste0(str_add_space(trait), " (", unit_use, ")"),
             trait_use = str_replace_all(trait_use, " ", "~")) %>%
      # ggplot(aes(x = effect, y = PC1, color = year)) + 
      ggplot(aes(x = effect, y = PC1, color = cluster, shape = year)) + 
      geom_point(size = ifelse(scores_toplot$group == "environment", 1.5, 0.5)) +
      geom_text_repel(data = filter(scores_toplot, group == "environment"), aes(label = term), size = 2) +
      scale_color_manual(values = cluster_colors, name = NULL, guide = FALSE) +
      scale_shape_manual(values = year_shape, name = NULL, guide = guide_legend(override.aes = list(size = 3))) + 
      scale_y_continuous(breaks = pretty, limits = range(scores_breaks)) +
      scale_x_continuous(breaks = pretty, limits = range(effect_breaks)) +
      ylab(names(propvar_rename)[1]) +
      xlab(expression("Effect (deviation from"~italic(mu)*")")) +
      facet_grid(~ trait_use, labeller = labeller(trait_use = label_parsed)) + 
      theme_presentation2(base_size = 10) +
      theme(legend.position = "bottom") +
      ## Annotate with cluster color
      annotate(geom = "text", size = 2, x = cluster_annotation_df$x, y = cluster_annotation_df$y, 
               label = cluster_annotation_df$label, hjust = cluster_annotation_df$hjust,
               vjust = cluster_annotation_df$vjust, color = cluster_annotation_df$color)
    
  })

## Combine plots
plot_list_complete <- subset(g_ammi, set == "complete" & trait %in% traits, plot, drop = T)
# ammi_grid <- plot_grid(plotlist = plot_list %>% map(~. + theme(legend.position = "none")), ncol = 1)
ammi_grid_complete <- plot_grid(ncol = 1, rel_heights = c(rep(0.9, length(plot_list_complete) - 1), 1), 
                                plotlist = c(head(plot_list_complete, -1) %>% map(~. + theme(legend.position = "none", axis.title.x = element_blank())),
                                             list(tail(plot_list_complete, 1)[[1]] + theme(legend.position = "none"))))

# Add legend
ammi_grid1 <- plot_grid(ammi_grid_complete, get_legend(plot_list_complete[[1]]), ncol = 1, rel_heights = c(1, 0.05))

# Save
ggsave(filename = "ammi_biplots_complete_alt1.jpg", plot = ammi_grid1, path = fig_dir, 
       height = 6.5, width = 3, dpi = 1000)









## Save all of this
save_file <- file.path(result_dir, "genotype_environment_phenotypic_analysis.RDAta")
save("training_sets_twoway", "ammi_out", "stage_two_fits_GYL", "stage_two_fits_GE", file = save_file)



## Year color scheme
year_shape <- setNames(c(15, 17, 18, 16), names(year_shape))


## Re-plot AMMI using different shaped points
## Plot the AMMI axes and the proportion of variance explained
g_ammi_alt <- ammi_out1 %>%
  left_join(., ammi_fitted_clusters) %>%
  # filter(set == "complete") %>%
  group_by(set, trait) %>%
  do(plot = {
    df <- .
    
    out <- df$ammi[[1]]
    res <- df$results[[1]]
    trait <- df$trait
    
    ## Combine the g and e scores into one DF
    scores_df <- map_df(out[c("escores", "gscores")], ~mutate(., group = names(.)[1]) %>% `names<-`(., c("term", names(.)[-1])))
    
    # Create a renaming vector
    res1 <- mutate(res, annotation = paste0(str_replace_all(term, "PC", "IPCA"), " (", round(cumprop, 2) * 100, "%)"))
    propvar_rename <- setNames(object = res1$term, nm = res1$annotation)
    
    # Renaming function
    label_rename <- function(x) str_replace_all(x, propvar_rename)
    
    scores_toplot <- scores_df %>%
      select(-eigen) %>%
      spread(PC, score) %>%
      # Add year
      mutate(year = ifelse(group == "environment", paste0("20", str_extract(term, "[0-9]{2}")), "Genotype"),
             trait = trait)

    
    ## Figure out the range for x and y axis for the trait
    effect_breaks <- subset(ammi_out, trait == df$trait, ammi, drop = T) %>% 
      map(~.[-1]) %>% map(., ~map(., "effect") %>% unlist()) %>% 
      unlist() %>% range()
    
    scores_breaks <- subset(ammi_out, trait == df$trait, ammi, drop = T) %>% 
      map(~.[-1]) %>% map(., ~map(., "score") %>% unlist()) %>% 
      unlist() %>% range()
    
    ## Units
    unit_use <- traits_unit1[df$trait]
    
    # Plot
    scores_toplot %>% 
      mutate(trait_use = paste0(str_add_space(trait), " (", unit_use, ")"),
             trait_use = str_replace_all(trait_use, " ", "~")) %>%
      ggplot(aes(x = effect, y = PC1)) + 
      geom_point(aes(color = year, shape = year), size = 1.5) +
      geom_text_repel(data = filter(scores_toplot, group == "environment"), aes(label = term), size = 2) +
      scale_color_manual(values = year_color, name = NULL, guide = guide_legend(override.aes = list(size = 3, label.size = 0))) + 
      scale_shape_manual(values = year_shape, name = NULL) + 
      scale_y_continuous(breaks = pretty, limits = range(scores_breaks)) +
      scale_x_continuous(breaks = pretty, limits = range(effect_breaks)) +
      ylab(names(propvar_rename)[1]) +
      xlab(expression("Effect (deviation from"~mu*")")) +
      facet_grid(~ trait_use, labeller = labeller(trait_use = label_parsed)) + 
      theme_presentation2(base_size = 10) +
      theme(legend.position = "bottom")
    # theme(legend.position = "right", legend.direction = "vertical")
    
  })



## Remove AGDD
plot_list_complete <- subset(g_ammi_alt, set == "complete" & trait %in% traits, plot, drop = T)
# ammi_grid <- plot_grid(plotlist = plot_list %>% map(~. + theme(legend.position = "none")), ncol = 1)
ammi_grid_complete <- plot_grid(ncol = 1, rel_heights = c(rep(0.9, length(plot_list_complete) - 1), 1), 
                                plotlist = c(head(plot_list_complete, -1) %>% map(~. + theme(legend.position = "none", axis.title.x = element_blank())),
                                             list(tail(plot_list_complete, 1)[[1]] + theme(legend.position = "none"))))

# Add legend
ammi_grid1 <- plot_grid(ammi_grid_complete, get_legend(plot_list_complete[[1]]), ncol = 1, rel_heights = c(1, 0.05))

# Save
ggsave(filename = "ammi_biplots_complete_greyscale.jpg", plot = ammi_grid1, path = fig_dir, height = 6.5, width = 3, dpi = 1000)



## Phenotypic Correlation Heatmap in greyscale
# One trait at a time
g_env_cor_df <- env_cor_df %>%
  group_by(trait) %>%
  do(plot = {
    
    df <- .
    
    # Use the heatmap function to construct a dendrogram
    # Use only data from the whole population
    mat <- df %>% 
      filter(population == "all") %>%
      select(-trait, -population) %>% 
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
      # scale_fill_gradientn(colors = heat_colors[c(1,3,5)], na.value = "transparent", 
      #                      breaks = seq(-1, 1, by = 0.5), limits = c(-1, 1), 
      #                      guide = guide_colorbar(title = "Phenotypic\ncorrelation")) +
      scale_fill_gradientn(colours = rev(grey.colors(n = 10)), na.value = "transparent", breaks = seq(-1, 1, by = 0.5), limits = c(-1, 1), 
                      guide = guide_colorbar(title = "Phenotypic\ncorrelation")) +
      ylab("Environment 2") +
      xlab("Environment 1") +
      facet_wrap(~ trait + population, labeller = labeller(trait = str_add_space, population = str_to_title)) +
      theme_presentation2(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
            axis.text.y = element_text(size = 5),
            # legend.position = "bottom", legend.direction = "horizontal")
            legend.position = "right", legend.direction = "vertical")
    
    
  }) %>% ungroup()




## Combine plots
g_env_cor <- plot_grid(plotlist = subset(g_env_cor_df, trait %in% traits, plot, drop = T) %>% 
                         map(~. + theme(legend.position = "none", axis.title = element_blank())), 
                       nrow = length(subset(g_env_cor_df, trait != "HeadingDateAGDD", trait, drop = T)), scale = 0.9)

g_env_cor1 <- g_env_cor + 
  draw_label("Environment 1", x = 0.5, y = 0, vjust=-0.5, angle= 0) +
  draw_label("Environment 2", x = 0, y = 0.5, vjust= 1.5, angle=90)

g_env_cor2 <- plot_grid(g_env_cor1, get_legend(g_env_cor_df$plot[[1]]), nrow = 1, rel_widths = c(1, 0.1))


# Save this plot
ggsave(filename = "environmental_correlation_greyscale.jpg", path = fig_dir, plot = g_env_cor2, height = 8, width = 8, dpi = 1000)












