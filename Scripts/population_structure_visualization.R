## S2MET
## 
## Visualizations
## 
## 1. Population structure
## 2. Map of trial locations
## 
## 

library(broom)
library(cowplot)

# The head directory
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))


## Plot relatedness of TP versus TP-VP relatedness
K_df <- K[tp_geno, c(tp_geno, vp_geno)] %>% 
  as.data.frame() %>% 
  rownames_to_column("line_name") %>% 
  gather(line_name2, G, -line_name) %>% 
  filter(line_name != line_name2) %>%
  mutate(pop = ifelse(line_name2 %in% vp_geno, "TPVP", "TPTP"),
         pop_x = ifelse(pop == "TPVP", "VP", "TP"),
         pop_y = "TP",
         pop = as.factor(pop))

## Model
fit <- lm(G ~ pop, data = K_df, subset = line_name != line_name2)
effects::allEffects(fit) %>% plot


## Visualize with heatmap
heat_colors <- wesanderson::wes_palette(name = "Zissou1")

g_relatedness_heat <- ggplot(data = K_df, aes(x = line_name2, y = line_name, fill = G)) +
  geom_tile() +
  scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[3], high = heat_colors[5], name = "Genomic relatedness") +
  facet_grid(pop_y ~ pop_x, scales = "free_x", space = "free_x", switch = "both") +
  theme_presentation(base_size = 10) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        panel.spacing.x = unit(0.1, "line"), legend.position = "top")

## Plot average relatedness
g_relatedness_bar <- group_by(K_df, pop_x, pop_y) %>% 
  mutate_at(vars(G), funs(mean, sd, n())) %>% 
  split(.$pop) %>%
  map_df(~mutate(., row = seq(nrow(.)))) %>%
  mutate(row_mean = round(mean(row))) %>%
  ungroup() %>% 
  mutate(se = sd / sqrt(n)) %>% 
  mutate_at(vars(mean, se), funs(ifelse(row != row_mean, NA, .))) %>%
  ggplot(aes(x = line_name2, y = mean, ymin = mean - se, ymax = mean + se, fill = pop_x)) +
  geom_col(width = 35) +
  geom_errorbar(width = 10) +
  scale_fill_manual(values = heat_colors[c(5,1)], guide = FALSE) +
  scale_y_continuous(name = "Genomic relatedness\nwith TP", breaks = pretty) +
  scale_x_discrete(expand = c(0, 0)) +
  facet_grid(~ pop_x, scales = "free_x", space = "free_x") +
  theme_presentation(base_size = 10) +
  theme(panel.grid = element_blank(), axis.ticks.y = element_line(), axis.text.x = element_blank(), axis.title.x = element_blank(), panel.spacing.x = unit(0.1, "line"))

g_relatedness <- plot_grid(g_relatedness_heat, g_relatedness_bar, ncol = 1, align = "hv", axis = "lr", rel_heights = c(1, 0.5))
ggsave(filename = "tp_vp_relatedness.jpg", plot = g_relatedness, path = fig_dir, width = 4, height = 6, dpi = 1000)








Ksvd <- prcomp(K)

# Color vector for programs
program_colors <- subset(entry_list, Line %in% tp_geno, Program, drop = T) %>%
  unique() %>%
  setNames(object = umn_palette(2, n = 7)[-1:-2], .)

# Plot
g_pop_struc <- tidy(Ksvd) %>% 
  rename(line_name = row) %>% 
  left_join(., select(entry_list, line_name = Line, program = Program)) %>% 
  filter(PC %in% 1:2) %>% 
  mutate(PC = paste0("PC", PC)) %>% 
  spread(PC, value)  %>% 
  ggplot(aes(x = PC1, y = PC2, color = program)) + 
  geom_point(size = 3) + 
  scale_color_manual(values = program_colors, name = NULL) +
  theme_poster() +
  theme(legend.position = c(0.10, 0.75), legend.key.height = unit(1.5, "lines"))

# Save
ggsave(filename = "tp_population_structure.jpg", plot = g_pop_struc, path = fig_dir, width = 6, height = 5, dpi = 1000)



## Plot the trial locations

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

# Edit the trial information
trial_info_toplot <- trial_info %>%
  filter(is.na(notes)) %>%
  group_by(location) %>%
  mutate(n_year = as.factor(n_distinct(year))) %>%
  select(location, lat = latitude, long = longitude, n_year) %>%
  distinct()

# Map
g_map <- ggplot(data = north_america, aes(x = long, y = lat)) +
  geom_polygon(fill = "white") +
  geom_polygon(data = canada, aes(group = group), fill = "grey85", color = "grey50", lwd = 0.75) + # Add canada
  geom_polygon(data = usa_state, aes(group = group), fill = "grey85", color = "grey50", lwd = 0.75) +
  geom_point(data = trial_info_toplot, aes(size = n_year)) +
  coord_map(projection = "mercator", xlim = c(-125, -60), ylim = c(25, 50)) +
  scale_size_manual(values = seq(3,7, length.out = 3), name = "Years in\nExperiment") +
  theme_void(base_size = 14) +
  theme(legend.position = c(0.85, 0.15))


# Save the figure
ggsave(filename = "trial_location_map.jpg", plot = g_map, path = fig_dir,
       width = 8, height = 5, dpi = 1000)
  












