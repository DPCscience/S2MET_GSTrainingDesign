## Analyze the clustering and heritability results
## 
## Author: Jeff Neyhart
## Last modified: May 2, 2018
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
library(lubridate)

# Load the clustering results
load(file.path(result_dir, "distance_methods_results.RData"))
# Load the heritability results
load(file.path(result_dir, "cluster_heritability_results.RData"))



## Assign traits to whether a higher value is favorable
higher_favorable_traits <- "GrainYield"

# Edit the traits to make higher value more favorable
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  mutate(value = if_else(trait %in% higher_favorable_traits, value, value * -1))


## Find the correlation between the different distance methods for pairs of environments
dist_method_df <- clust_method_df %>% 
  mutate(dist_df = map(dist, ~as.matrix(.) %>% as.data.frame() %>% rownames_to_column("environment1") %>% 
                         gather(environment2, distance, -environment1) %>% filter(environment1 != environment2))) %>%
  unnest(dist_df)

# Calculate the correlation
dist_method_cor <- dist_method_df %>% 
  spread(dist_method, distance) %>% 
  group_by(trait, population) %>% 
  do(cor(select(., -trait:-environment2)) %>%
       as.data.frame() %>%
       rownames_to_column("dist_method1") %>%
       gather(dist_method2, correlation, -dist_method1))

## Plot
dist_method_cor %>% 
  ggplot(aes(x = dist_method1, y = dist_method2, fill = correlation)) + 
  geom_tile() + 
  facet_grid(trait ~ population)








## Manipulate the cluster results and cut the trees
# What should be the minimum and maximum number of clusters?
min_k <- 2
max_k <- 20
seq_k <- seq(min_k, max_k)


## Combine the data.frames
## Replicate the clusters and add each cluster k from min_k to max_k
clust_method_df_use <- clust_method_df %>%
  rerun(.n = length(seq_k), .) %>%
  list(., seq_k) %>% 
  pmap_df(~mutate(.x, k = .y))





##### Visualize Distance in 2 Dimensions

## Use multi-dimensional scaling to convert the distance matrices to a 2D matrix
clust_method_mds <- clust_method_df %>%
  filter(dist_method != "ec_one_PCA_dist") %>%
  mutate(mds = map(dist, ~cmdscale(d = ., k = 2) %>%
                     fix_data_frame(newnames = c("x", "y"), newcol = "environment"))) %>%
  unnest(mds) %>%
  mutate(dist_method = as_replaced_factor(dist_method, dist_method_replace),
         ## Extract year from environment
         year = as.factor(year(parse_date_time(x = parse_number(environment), orders = "y"))))
         




### First plot everything for each trait
g_mds_traits <- clust_method_mds %>% 
  split(.$trait) %>%
  map(~{
    tr <- unique(.$trait)
    ggplot(., aes(x = x, y = y)) +
      geom_point() +
      geom_text(aes(label = environment), size = 2, check_overlap = TRUE, vjust = 2) +
      facet_wrap(population ~ dist_method, scales = "free", ncol = 3) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      labs(title = str_c("Trait: ", tr))
  })
    

## Use this space to visualize individual traits
g_mds_traits$GrainYield
g_mds_traits$HeadingDate
g_mds_traits$PlantHeight



### Second plot everything for each distance measure
g_mds_dist <- clust_method_mds %>% 
  mutate(dist_method = str_replace_all(dist_method, pattern = "\n", replacement = " ")) %>%
  split(.$dist_method) %>%
  map(~{
    dm <- unique(.$dist_method)
    ggplot(., aes(x = x, y = y)) +
      geom_point() +
      geom_text(aes(label = environment), size = 2, check_overlap = TRUE, vjust = 2) +
      facet_wrap(population ~ trait, scales = "free", ncol = 3) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      labs(title = str_c("Distance Method: ", dm))
  })

## Use this space to visualize individual traits
g_mds_dist$`Phenotypic Distance`
g_mds_dist$`10 yr Environmental Covariates`
g_mds_dist$`Great Circle Distance`




## Plot all distance methods for each trait and each distance methods
g_clust_trait <- clust_method_mds %>% 
  filter(population == "all") %>%
  ggplot(aes(x = x, y = y, color = year)) +
  geom_point() +
  # geom_text(aes(label = environment), size = 2, check_overlap = TRUE, vjust = 2) +
  facet_wrap(trait ~ dist_method, scales = "free") +
  scale_color_discrete(name = "Year") + 
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())

ggsave(filename = "cluster_trait_mds.jpg", plot = g_clust_trait, path = fig_dir,
       height = 6, width = 6, dpi = 1000)





# 
# ### Attempt to create plots with different number of clusters
# test_dist <- subset(clust_method_df, trait == "GrainYield" & dist_method == "ge_mean_D" & population == "all", 
#                     dist, drop = T)[[1]]
# 
# test_clust <- hclust(d = test_dist)
# 
# ## Calculate the silhouette across k = 1:20 clusters
# ks <- seq(2, length(test_clust$labels) - 1)
# 
# ks_sil <- ks %>%
#   map(~cutree(test_clust, k = .)) %>%
#   map(~silhouette(x = ., dist = test_dist))
# 
# # Calculate the mean silhouette at each splitting
# ks_ss <- ks_sil %>%
#   map_dbl(~mean(.[,3]))
# 
# # Plot
# plot(ks_ss)
# 
# 
# 
# ## A function that calculates the mean distances within a set of k clusters
# within_dist <- function(d, method = "complete", k) {
#   
#   # Convert the distance object to a matrix
#   d_mat <- as.matrix(d)
#   
#   # Cluster
#   d_clust <- hclust(d = d, method = method)
#   # Cut the tree
#   tree_cut <- cutree(tree = d_clust, k = k)
#   
#   # Split the members and then iterate over the members to calculate the average distance within the group
#   split(tree_cut, tree_cut) %>%
#     map(~d_mat[names(.),names(.)]) %>%
#     map_dbl(~ifelse(is.null(dim(.)), 0, mean(.[lower.tri(.)])))
#   
# }
#   
  











## Using the value of k, cut the cluster tree and create data.frames for the environment
## and the assigned cluster
clust_method_df_tomodel <- clust_method_df_use %>% 
  mutate(env_cluster = list(cluster, k) %>% 
           pmap(~cutree(tree = .x, k = .y) %>% 
                  data.frame(environment = names(.), cluster = ., stringsAsFactors = FALSE)) ) %>%
  unnest(env_cluster)

# How many environments per cluster per k?
clust_method_df_summ <- clust_method_df_tomodel %>% 
  mutate_at(vars(k, cluster), as.factor) %>%
  group_by(trait, population, dist_method, k, cluster) %>% 
  summarize(n_env = n_distinct(environment)) %>%
  mutate(prop_env = n_env / sum(n_env)) %>%
  ungroup()

plot_fun <- function(pop, label = TRUE) {
  clust_method_df_summ %>% 
    mutate(dist_method = str_replace_all(dist_method, dist_method_replace),
           dist_method = factor(dist_method, levels = dist_method_replace)) %>% 
    filter(population == pop) %>% 
    ggplot(aes(x = k, y = prop_env, fill = cluster)) + 
    geom_col() + 
    facet_grid(trait ~ dist_method) +
    scale_x_discrete(breaks = function(x) {x1 = as.numeric(x); seq(min(x1), max(x1), 3)}) +
    labs(title = str_c("Population: ", pop)) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

g_clust_nenv_all <- plot_fun("all")

g_clust_nenv_tp <- plot_fun("tp")











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
    mutate(dist_method = str_replace_all(dist_method, dist_method_replace),
           dist_method = factor(dist_method, levels = dist_method_replace)) %>% 
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



# Examine the log likelihood of the model fits over number of clusters
cluster_method_loglik <- cluster_method_herit_df %>%
  mutate(loglik = map(out, "loglik")) %>%
  unnest(loglik) %>% 
  select(-out)


## Define a common plot modifier function
plot_fun <- function(pop) {
  cluster_method_loglik %>% 
    mutate(dist_method = str_replace_all(dist_method, dist_method_replace),
           dist_method = factor(dist_method, levels = dist_method_replace)) %>% 
    filter(population == pop) %>% 
    ggplot(aes(x = k, y = loglik)) +
    geom_point() + 
    geom_line() +
    facet_grid(trait ~ dist_method, scales = "free") + 
    labs(title = str_c("Population: ", pop)) +
    theme_bw()
}

# Map over the population type
population <- unique(cluster_method_herit$population) %>% set_names(., .)
g_loglik_list <- map(population, plot_fun)

ggsave(filename = file.path(fig_dir, "cluster_loglik_pop_all.jpg"), 
       plot = g_loglik_list$all, width = 10, height = 6)

ggsave(filename = file.path(fig_dir, "cluster_loglik_pop_pop.jpg"), 
       plot = g_loglik_list$tp, width = 10, height = 6)



## Combine the heritability and log-likelihood trends



