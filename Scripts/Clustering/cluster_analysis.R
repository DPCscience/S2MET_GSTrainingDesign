## Analyze the clustering and heritability results
## 
## Author: Jeff Neyhart
## Last modified: March 29, 2018
## 
## This script will perform different clustering procedures on the S2MET data
## and calculate the within and across-cluster heritabilities.
## 


# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load packages and the source script
library(lme4)
library(broom)

# Load the clustering results
load(file.path(result_dir, "distance_methods_results.RData"))
# Load the heritability results
load(file.path(result_dir, "cluster_heritability_results.RData"))



## Assign traits to whether a higher value is favorable
higher_favorable_traits <- "GrainYield"

# Edit the traits to make higher value more favorable
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  mutate(value = if_else(trait %in% higher_favorable_traits, value, value * -1))


## Manipulate the cluster results and cut the trees
# What should be the minimum and maximum number of clusters?
min_k <- 2
max_k <- 10
seq_k <- seq(min_k, max_k)

# Create a name replacement vector for the distance methods
dist_method_replace <- c(
  "great_circle_dist" = "Great Circle\nDistance", "env_cor_dist" = "Environment Genetic\nCorrelation",
  "ge_mean_D" = "Phenotypic\nDistance", "ge_PCA_dist" = "GxE BLUP PCA", 
  "ec_one_PCA_dist" = "1 yr Environmental\nCovariates", "ec_multi_PCA_dist" = "10 yr Environmental\nCovariates")


## Combine the data.frames
## Replicate the clusters and add each cluster k from min_k to max_k
clust_method_df_use <- clust_method_df %>%
  rerun(.n = length(seq_k), .) %>%
  list(., seq_k) %>% 
  pmap_df(~mutate(.x, k = .y))


## Use multi-dimensional scaling to convert the distance matrices to a 2D matrix
clust_method_mds <- clust_method_df %>%
  mutate(mds = map(dist, ~cmdscale(d = ., k = 2) %>%
                     fix_data_frame(newnames = c("x", "y"), newcol = "environment"))) %>%
  unnest(mds)


# Plot function
plot_clust <- function(tr, label = TRUE) {
  g <- filter(clust_method_mds, trait == tr, population == "all") %>% 
    mutate(dist_method = str_replace_all(dist_method, dist_method_replace)) %>% 
    ggplot(aes(x = x, y = y)) + 
    geom_point() + 
    facet_wrap(~ dist_method, scales = "free") +
    labs(title = tr) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  if (label) {
    g1 <- g + geom_text(aes(label = environment), size = 2, check_overlap = TRUE, vjust = 2)
  } else {
    g1 <- g
  }
  return(g1)
}

## Plot the the results of multi-dimensional scaling
g_dist_mds_all <- set_names(traits, traits) %>%
  map(plot_clust)



## Using the value of k, cut the cluster tree and create data.frames for the environment
## and the assigned cluster
clust_method_df_tomodel <- clust_method_df_use %>% 
  mutate(env_cluster = list(cluster, k) %>% 
           pmap(~cutree(tree = .x, k = .y) %>% 
                  data.frame(environment = names(.), cluster = ., stringsAsFactors = FALSE)) )








### Cluster heritability results
# Tidy up
cluster_method_herit_df <- cluster_method_herit_out %>% 
  bind_rows()

# Extract the heritability estimates
cluster_method_herit <- cluster_method_herit_df %>%
  mutate(herit = map(out, "herit")) %>%
  unnest(herit) %>% 
  select(-out) %>%
  gather(herit_type, estimate, across, within)

# Create plots per population
## Define a common plot modifier function
plot_fun <- function(pop) {
  cluster_method_herit %>% 
    mutate(dist_method = str_replace_all(dist_method, dist_method_replace)) %>%
    filter(population == pop) %>% 
    ggplot(aes(x = k, y = estimate, shape = herit_type)) +
    geom_point() + 
    geom_line() +
    facet_grid(trait ~ dist_method) + 
    labs(title = str_c("Population: ", pop)) +
    theme_bw()
}

# Map over the population type
population <- unique(cluster_method_herit$population) %>% set_names(., .)
g_herit_list <- map(population, plot_fun)

ggsave(filename = file.path(fig_dir, "cluster_herit_pop_all.jpg"), 
       plot = g_herit_list$all, width = 10, height = 6)

ggsave(filename = file.path(fig_dir, "cluster_herit_pop_pop.jpg"), 
       plot = g_herit_list$tp, width = 10, height = 6)







