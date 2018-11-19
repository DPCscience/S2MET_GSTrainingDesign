## S2MET
## 
## Visualizations
## 
## 1. Population structure
## 2. Map of trial locations
## 
## 

library(broom)

# The head directory
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# PCA of marker genotypes
K <- A.mat(s2_imputed_mat_use[tp_geno,])

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
  












