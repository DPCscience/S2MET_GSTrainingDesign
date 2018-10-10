## S2MET
## 
## Visualize population structure
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


