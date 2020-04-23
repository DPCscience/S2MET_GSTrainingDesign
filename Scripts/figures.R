## S2MET Predictions
## 
## Script for recreating all figures
## 
## 1. Population structure
## 2. Map of trial locations
## 
## 


# Packages to load
library(broom)
library(cowplot)
library(patchwork)
library(ggrepel)
library(ggsn)

# The head directory
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

## Calculate genomic relationship
genom_rel <- matrix(data = 0, nrow = nrow(K), ncol = ncol(K), dimnames = dimnames(K))



# Iterate over rows and columns
for (j in seq(ncol(K))) {
  for (i in seq(nrow(K))) {

    # Skip if i <= j
    if (i <= j) {
      next
      
    } else {
      genom_rel[i,j] <- K[i,j, drop = F] / prod(sqrt(diag(K)[c(i, j)]))
      
    }
  }
}

genom_rel1 <- as.matrix(as.dist(genom_rel))

## Convert to data.frame
genom_rel_df <- genom_rel1 %>% 
  as.data.frame() %>% 
  rownames_to_column("line_name") %>% 
  gather(line_name2, G, -line_name) %>% 
  filter(line_name != line_name2)


## Create DF of TP-TP relatedness
genome_rel_tptp_df <- genom_rel1[tp_geno, tp_geno] %>% 
  as.dist() %>%
  tidy(diagonal = F) %>%
  rename(line_name2 = item1, line_name = item2, G = distance)

## Create DF of TP-VP relatedness
genome_rel_tpvp_df <- genom_rel1[tp_geno, vp_geno] %>% 
  as.data.frame() %>% 
  rownames_to_column("line_name") %>% 
  gather(line_name2, G, -line_name) %>% 
  filter(line_name != line_name2)

## Combine
genom_rel_df1 <- bind_rows(genome_rel_tptp_df, genome_rel_tpvp_df) %>%
  mutate(pop = ifelse(line_name2 %in% vp_geno, "TPVP", "TPTP"),
         pop_x = ifelse(pop == "TPVP", "Test", "Train"),
         pop_x = factor(pop_x, levels = c("Train", "Test")),
         pop_y = "Train",
         pop = as.factor(pop))


## Resample the genomic relatedness of size n_test,
## calculate genomic relationship. Repeat
null_genome_rel_mean <- replicate(100000, sample(c(tp_geno, vp_geno), size = length(vp_geno)), simplify = FALSE) %>%
  map(~subset(genom_rel_df1, line_name2 %in% ., G, drop = T)) %>%
  map_dbl(mean)

## Calculate realized mean relatedness between tp and vp
real_genome_rel_mean <- group_by(genom_rel_df1, pop_x, pop_y) %>% 
  summarize_at(vars(G), list(~mean, ~sd, ~n())) %>% 
  mutate(se = sd / sqrt(n),
         comparison = paste0(pop_y, "-", pop_x))

## Plot with realized means
tibble(x = null_genome_rel_mean) %>% 
  ggplot(aes(x = x)) + 
  geom_histogram(fill = "grey85", color = "black") +
  geom_vline(data = real_genome_rel_mean, aes(xintercept = mean, color = comparison), lwd = 2 )

## Calculate p-value for trait-test relatedness
real_genome_rel_mean %>%
  mutate(p.value = mean(null_genome_rel_mean < mean))

# pop_x pop_y      mean    sd     n      se comparison 
# 1 Train Train  0.000780 0.188 15225 0.00152 Train-Train
# 2 Test  Train -0.0213   0.129  8400 0.00141 Train-Test




## Model

## Visualize with heatmap
heat_colors <- wesanderson::wes_palette(name = "Zissou1")

g_relatedness_heat <- ggplot(data = genom_rel_df1, aes(x = line_name2, y = line_name, fill = G)) +
  geom_tile() +
  scale_fill_gradient2(low = heat_colors[1], mid = heat_colors[3], high = heat_colors[5], name = "Genomic relationship") +
  facet_grid(pop_y ~ pop_x, scales = "free_x", space = "free_x", switch = "both") +
  theme_presentation2(base_size = 10) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), panel.border = element_blank(),
        panel.spacing.x = unit(0.1, "line"), legend.position = "top")

## Plot average relatedness
g_relatedness_bar <- genom_rel_df1 %>% 
  ## Create a variable to determine the width of the facet strip
  distinct(line_name2, pop_x) %>% 
  split(.$pop_x) %>%
  map_df(~mutate(., row = seq(nrow(.)))) %>%
  group_by(pop_x) %>%
  mutate(row_mean = round(mean(row))) %>%
  left_join(., real_genome_rel_mean) %>%
  mutate_at(vars(mean, se), funs(ifelse(row != row_mean, NA, .))) %>%
  ## Plot
  ggplot(aes(x = line_name2, y = mean, ymin = mean - se, ymax = mean + se, fill = pop_x)) +
  geom_hline(yintercept = 0, lwd = 0.5) +
  geom_col(width = 35) +
  geom_errorbar(width = 10) +
  scale_fill_manual(values = heat_colors[c(5,1)], guide = FALSE) +
  scale_y_continuous(name = "Genomic relationship\nwith Train", breaks = pretty) +
  scale_x_discrete(expand = c(0, 0)) +
  facet_grid(~ pop_x, scales = "free_x", space = "free_x", switch = "x") +
  theme_presentation(base_size = 10) +
  theme(panel.grid = element_blank(), axis.ticks.y = element_line(), axis.text.x = element_blank(), axis.title.x = element_blank(), panel.spacing.x = unit(0.1, "line"))

## Use patchwork to combine the plots
g_relatedness <- g_relatedness_heat + g_relatedness_bar + 
  plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = "a")

ggsave(filename = "tp_vp_relatedness.jpg", plot = g_relatedness, path = fig_dir, width = 4, height = 6, dpi = 1000)






Ksvd <- prcomp(K)

# Color vector for programs
program_colors <- entry_list %>%
  mutate(Program = parse_character(Program)) %>%
  filter(!is.na(Program)) %>%
  distinct(Program) %>%
  pull() %>%
  setNames(object = c(neyhart_palette("umn1", 7)[-1:-2], neyhart_palette("umn2")[4]), .)

# Plot
g_pop_struc <- tidy(Ksvd) %>% 
  rename(line_name = row) %>% 
  left_join(., select(entry_list, line_name = Line, program = Program)) %>% 
  filter(program != "M2") %>%
  filter(PC %in% 1:2) %>% 
  mutate(PC = paste0("PC", PC),
         program = factor(program, levels = names(program_colors))) %>% 
  spread(PC, value)  %>% 
  ggplot(aes(x = PC1, y = PC2, color = program)) + 
  geom_point(size = 3) + 
  scale_color_manual(values = program_colors, name = NULL, drop = FALSE) +
  theme_poster() +
  theme(legend.position = c(0.10, 0.75), legend.key.height = unit(1.5, "lines"), legend.box = "horizontal")

# Save
ggsave(filename = "population_structure.jpg", plot = g_pop_struc, path = fig_dir, width = 6, height = 5, dpi = 1000)
ggsave(filename = "population_structure_tp.jpg", plot = g_pop_struc, path = fig_dir, width = 6, height = 5, dpi = 1000)








## Greyscale heatmap

g_relatedness_heat <- ggplot(data = genom_rel_df1, aes(x = line_name2, y = line_name, fill = G)) +
  geom_tile() +
  scale_fill_gradient2(low = grey.colors(11)[11], mid = grey.colors(11)[5], high = grey.colors(11)[1], name = "Genomic relationship") +
  facet_grid(pop_y ~ pop_x, scales = "free_x", space = "free_x", switch = "both") +
  theme_presentation2(base_size = 10) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), panel.border = element_blank(),
        panel.spacing.x = unit(0.1, "line"), legend.position = "top")

## Plot average relatedness
g_relatedness_bar <- genom_rel_df1 %>% 
  distinct(line_name2, pop_x) %>% 
  split(.$pop_x) %>%
  map_df(~mutate(., row = seq(nrow(.)))) %>%
  group_by(pop_x) %>%
  mutate(row_mean = round(mean(row))) %>%
  left_join(., real_genome_rel_mean) %>%
  mutate_at(vars(mean, se), funs(ifelse(row != row_mean, NA, .))) %>%
  ggplot(aes(x = line_name2, y = mean, ymin = mean - se, ymax = mean + se, fill = pop_x)) +
  geom_col(width = 35) +
  geom_errorbar(width = 10) +
  scale_fill_manual(values = heat_colors[c(5,1)], guide = FALSE) +
  scale_y_continuous(name = "Genomic relationship\nwith Train", breaks = pretty) +
  scale_x_discrete(expand = c(0, 0)) +
  facet_grid(~ pop_x, scales = "free_x", space = "free_x", switch = "x") +
  theme_presentation(base_size = 10) +
  theme(panel.grid = element_blank(), axis.ticks.y = element_line(), axis.text.x = element_blank(), axis.title.x = element_blank(), panel.spacing.x = unit(0.1, "line"))


g_relatedness <- plot_grid(g_relatedness_heat, g_relatedness_bar, ncol = 1, align = "hv", axis = "lr", rel_heights = c(1, 0.5),
                           labels = LETTERS[1:2])
ggsave(filename = "tp_vp_relatedness_greyscale.jpg", plot = g_relatedness, path = fig_dir, width = 4, height = 6, dpi = 1000)







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


## Now create labels for each location (as a combination of all environments)
use_loc_info_toplot <-  trial_info %>% 
  filter(environment %in% tp_vp_env) %>%
  filter(is.na(notes)) %>%
  group_by(location, latitude, longitude) %>%
  summarize(n_years = n_distinct(year)) %>%
  ungroup() %>%
  mutate(n_years = factor(n_years, levels = sort(unique(n_years))),
         location = str_to_title(str_replace_all(location, "_", " ")))


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
  geom_text_repel(data = use_loc_info_toplot1, aes(x = longitude, y = latitude, group = location, label = location), 
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
  north(location = "bottomleft", symbol = 12, x.min = min(long_limit) - 2, x.max = max(long_limit), 
              y.min = min(lat_limit) + 2, y.max = max(lat_limit)) +
  # Add an example point
  geom_point(aes(x = min(long_limit) + 0.3, y = min(lat_limit) + 1), size = 3, inherit.aes = FALSE) +
  geom_text(aes(x = min(long_limit) + 0.3, y = min(lat_limit) + 1, label = "2"), color = "white", 
            size = 2, inherit.aes = FALSE) +
  geom_text(aes(x = min(long_limit) + 1.2, y = min(lat_limit) + 1, label = "No. years at location"), 
            size = 2, inherit.aes = FALSE, hjust = 0, nudge_x = 0.1)


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

















