## Use environmental clusters to calculate heritability
## 
## Author: Jeff Neyhart
## Last modified: March 29, 2018
## 
## This script will perform different clustering procedures on the S2MET data
## and calculate the within and across-cluster heritabilities.
## 


### Run for MSI

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"
source(file.path(repo_dir, "source_MSI.R"))


## Load packages
packages <- c("lme4")
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))

# Source some scripts from pbr
source("/panfs/roc/groups/6/smithkp/neyha001/R/my_packages/pbr/R/convenience_functions.R")
source("/panfs/roc/groups/6/smithkp/neyha001/R/my_packages/pbr/R/herit.R")

### Run for local machine

# # Run the source script
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# 
# # Load packages and the source script
# library(lme4)


# Load the clustering results
load(file.path(result_dir, "distance_methods_results.RData"))


## Number of cores
n_core <- ifelse(Sys.info()["sysname"] == "Windows", 1, detectCores())


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




## Using the value of k, cut the cluster tree and create data.frames for the environment
## and the assigned cluster
clust_method_df_tomodel <- clust_method_df_use %>% 
  mutate(env_cluster = list(cluster, k) %>% 
           pmap(~cutree(tree = .x, k = .y) %>% 
                  data.frame(environment = names(.), cluster = ., stringsAsFactors = FALSE)) )

# List of models
# forms <- formulas(~ value,
#                   no_gc = ~ (1 | line_name) + (1|environment) + (1|line_name:environment) + 
#                     (1|line_name:environment:cluster),
#                   full = ~ (1 | line_name) + environment + (1|line_name:environment) + 
#                     (1|line_name:environment:cluster) + (1|line_name:cluster) )


# Breakup by core
clust_method_df_tomodel_core <- clust_method_df_tomodel %>% 
  assign_cores(n_core) %>%
  split(.$core)

## Parallelize over the core df
cluster_method_herit_out <- mclapply(X = clust_method_df_tomodel_core, FUN = function(core_df) {
 
  # Create a list to store the results
  results_list <- vector("list", nrow(core_df))
  
  # Iterate over the list
  for (i in seq_along(results_list)) {
    
    df <- unnest(core_df[i,], env_cluster)

    #Print a message
    print(paste(unique(df$trait), unique(df$dist_method), sep = ", "))
    
    # Combine the phenotypic data
    pheno_cluster <- left_join(df, S2_MET_BLUEs, by = c("trait", "environment"))
    
    # Filter all data or the tp, depending on the population
    if (unique(df$population) == "tp") {
      pheno_tomodel <- pheno_cluster %>%
        filter(line_name %in% tp)
      
    } else {
      pheno_tomodel <- pheno_cluster
      
    }
    
    # Extract the weights
    wts <- pheno_tomodel$std_error^2
    
    # lmer control
    lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore",
                                calc.derivs = FALSE)
    
    # Main formula
    form <- value ~ (1 | line_name) + (1|location) + (1|year) + (1|line_name:location) + 
      (1|line_name:year) + (1|line_name:location:cluster) + (1|line_name:year:cluster) + 
      (1|line_name:cluster) + (1|line_name:location:year:cluster)    
    
    
    # Fit the model
    fit <- lmer(formula = form, data = pheno_tomodel, control = lmer_control, weights = wts)
    
    # Calculate heritability
    herit_df <- cluster_heritability(object = fit, breakup_env = TRUE)
    
    # Return the variance components
    var_comp <- as.data.frame(VarCorr(fit)) %>%
      select(source = grp, variance = vcov)
    
    # Return a list
    results_list[[i]] <- list(herit = herit_df, var_comp = var_comp, loglik = as.numeric(logLik(fit)))
    
  }
  
  # Add the list to the core_df and return
  core_df %>%
    select(trait:dist_method, k) %>%
    mutate(out = results_list)
  
}, mc.cores = n_core)



save_file <- file.path(result_dir, "cluster_heritability_results.RData")
save("cluster_method_herit_out", file = save_file)





