## S2MET Clustering and Prediction Procedures
## 
## This script is designed to be run on a local machine or on MSI

# Is this MSI?
MSI <- TRUE

# List of packages
packages <- c("tidyverse", "stringr", "readxl", "rrBLUP", "gws", "pbr", "purrrlyr",
              "psych", "parallel")

if (MSI) {
  ## MSI directories
  # Directory to save the files
  proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"  

  # Set the directory of the R packages
  package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"
  
  # Load all packages
  invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))
  
  # Detect cores
  n_cores <- detectCores()
  
} else {
  ## Local directories
  # The head directory
  proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET/"
  
  # Load all packages
  invisible(lapply(packages, library, character.only = TRUE))
  
  # Local = 1 core
  n_cores <- 1
  
}

# Prediction directory
pred_dir <- file.path(proj_dir, "Predictions")


## Load data
if (MSI) {
  load(file.path(pred_dir, "MSI/S2MET_MSI_prediction_material.RData"))
  
} else {

  # Other directories
  geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/Genotypic Data/GBS Genotype Data/"
  env_var_dir <- file.path(proj_dir, "Environmental_Variables/")
  pheno_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET/Phenotype_Data/"
  # Directory where entry data is located
  entry_dir <- file.path(proj_dir, "Plant_Materials")

  # Load the genotypic data
  load(file.path(geno_dir, "S2_genos_mat.RData"))
  
  # Load the phenotypic data
  load(file.path(pheno_dir, "S2_MET_BLUE.RData"))
  load(file.path(pheno_dir, "S2_MET_tidy.RData"))
  
  # Load environmental data
  load(file.path(env_var_dir, "environmental_data_compiled.RData"))
  
  ## Load the entry data
  entry_list <- file.path(entry_dir, "S2MET_project_entries.xlsx") %>%
    read_excel()
  
}

# # Data to save to MSI
# to_save <- c("entry_list", "s2_imputed_mat", "S2_MET_tidy", "S2_MET_BLUE", 
#              "one_year_mat", "multi_year_mat")
# save(list = to_save, file = file.path(pred_dir, "MSI/S2MET_MSI_prediction_material.RData"))

# Separate into entries and checks
entries <- entry_list %>% 
  filter(Class %in% c("S2TP", "S2C1R"))

## Separate the TP and VP
tp <- entries %>%
  filter(Class == "S2TP", Line %in% row.names(s2_imputed_mat)) %>%
  pull(Line)

vp <- entries %>%
  filter(Class == "S2C1R", Line %in% row.names(s2_imputed_mat)) %>%
  pull(Line)

checks <- entry_list %>% 
  filter(Class == "Check") %>%
  pull(Line)

# Filter the BLUEs for the entries
S2_MET_BLUE_entries <- S2_MET_BLUE %>%
  filter(line_name %in% c(tp, vp)) %>%
  droplevels()

S2_MET_tidy_use <- S2_MET_tidy %>%
  mutate_at(vars(line_name), as.character) %>%
  mutate(line = ifelse(!line_name %in% checks, line_name, "00check"),
         check = ifelse(line_name %in% checks, line_name, "00entry")) %>%
  droplevels()



# Calculate the genomic relationship matrix
## Matrix construction
# Genomic relationship matrix
# Convert the matrix name
M <- s2_imputed_mat[c(tp, vp), ]

# Calculate the genomic relationship matrix
G_1 <- A.mat(X = M, min.MAF = 0, max.missing = 1)







# Define some functions for calculating heritability

# Functions
# Define the clustering function
hclusCut <- function(x, k, ...) {
  list(cluster = cutree(hclust(as.dist(x), method = "average", ...), k=k))
}

# Function to calculate base heritability
base_herit <- function(mod, g, e) {
  # Model frame
  mf <- model.frame(mod)
  # Unique environments
  n_e <- n_distinct(mf$environment)
  # Reps in environments
  n_r <- pbr::harm_mean(table(subset(mf, , c(g, e))))
  
  # Variance components
  varcomp <- as.data.frame(VarCorr(mod))
  varG <- subset(varcomp, grp == g, vcov, drop = TRUE)
  varGE <- subset(varcomp, grp == paste(g, e, sep = ":"), vcov, drop = TRUE)
  varR <- subset(varcomp, grp == "Residual", vcov, drop = TRUE)
  
  # Heritability 
  varG / (varG + (varGE / n_e) + (varR / (n_e * n_r)))
} # Close the function

# Function to calculate across and within cluster heritabilities
clust_herit <- function(mod, g, e, cl) {
  # Model frame
  mf <- model.frame(mod)
  # Unique environments
  n_e <- n_distinct(mf[[e]])
  # Reps in environments
  n_r <- pbr::harm_mean(table(subset(mf, , c(g, e))))
  # Clusters
  n_s <- n_distinct(mf[[cl]])
  
  varcomp <- as.data.frame(VarCorr(mod))
  varG <- subset(varcomp, grp == g, vcov, drop = TRUE)
  varGE <- subset(varcomp, grp == paste(g, e, sep = ":"), vcov, drop = TRUE)
  varR <- subset(varcomp, grp == "Residual", vcov, drop = TRUE)
  varGS <- subset(varcomp, grp == paste(g, cl, sep = ":"), vcov, drop = TRUE)
  
  # Heritability - across
  H_a <- varG / (varG + (varGS / n_s) + (varGE / n_e) + (varR / (n_e * n_r)))
  # Heritability - within
  H_w <- (varG + varGS) / ((varG + varGS) + (n_s * ( (varGE / n_e) + (varR / (n_e * n_r)) ) ))
  data.frame(selection = c("across", "within"), heritability = c(H_a, H_w))
}


# A function that takes training and test data and spits back a prediction accuracy
pred_acc_out <- function(train_df, test_df) {
  # First, if either df has no rows, return NA
  if (nrow(train_df) == 0 | nrow(test_df) == 0) {
    return(as.data.frame(NA))
  } else {
    # Create model matrices
    mf <- model.frame(value ~ line_name + environment, data = train_df)
    y <- model.response(mf)
    if (n_distinct(mf$environment) == 1) {
      X <- model.matrix(~ 1, data = mf)
    } else {
      X <- model.matrix(~ 1 + environment, data = mf)
    }
    Z <- model.matrix(~ line_name, data = mf)
    
    # Solve the mixed model
    rr_out <- mixed.solve(y = y, Z = Z, K = G_1, X = X, method = "REML")

    # Extract the predictions for the test set
    test_pred <- rr_out$u %>% 
      data_frame(line_name = row.names(.), pred = .) %>% 
      filter(line_name %in% test_df$line_name)
    
    # Measure accuracy
    full_join(test_pred, select(test_df, line_name, value)) %>% 
      do({pred_acc = boot_cor(x = .$pred, y = .$value, boot.reps = 1000)})
  }
}





### Clustering

#### Distance Measure Based on Genotypic Means

# Cluster first on all available data - including the validation set,
# then cluster only on the TP (n = 175)

# First fit a base model using the tidy data
base_mod <- S2_MET_tidy_use %>%
  group_by(trait) %>%
  do(mod = lmer(value ~ (1|line) + check + environment + (1|line:environment), data = .))

# Calculate base heritability
herit_base <- base_mod %>%
  mutate(H_base = base_herit(mod, g = "line", e = "environment"))

# Fit another model, but removing the validation population individuals
tp_mod <- S2_MET_tidy_use %>%
  filter(line_name %in% c(tp, checks)) %>%
  group_by(trait) %>%
  do(mod = lmer(value ~ (1|line) + check + environment + (1|line:environment), data = .))

# First create a data_frame of two populations: tp + vp and only tp
S2_MET_BLUE_entries_pop <- data_frame(
  population = c("all", "tp"),
  data = list(filter(S2_MET_BLUE_entries, line_name %in% c(tp, vp)), 
              filter(S2_MET_BLUE_entries, line_name %in% c(tp))))

  
# Perform the distance calculation on the line means in each environment
ge_mean_D <- S2_MET_BLUE_entries_pop %>% 
  unnest() %>%
  group_by(population, trait) %>%
  do(D = {
    dist1 <- dist_env(x = ., gen.col = "line_name", 
                      env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })



# Extract the GxE BLUPs and run PCA
# Extract the blups from the base_mod and then from the tp_mod
ge_blup <- data_frame(
  population = c("all", "tp"),
  mod_df = list(base_mod, tp_mod)) %>%
  unnest() %>%
  # Group by the population and then the trait
  group_by(population, trait) %>%
  do(mat = {
    ranef(.$mod[[1]])$`line:environment` %>%
      rownames_to_column("grp") %>% 
      separate("grp", into = c("line", "environment"), sep = ":", fill = "right") %>%
      rename(effect = `(Intercept)`) %>%
      filter(line != "00check") %>%
      spread(environment, effect) %>%
      as.data.frame() %>%
      remove_rownames() %>%
      column_to_rownames("line") %>%
      as.matrix() })

# Impute using mean over each environment then run PCA
# Number of PCs to use
n_PC <- 2

ge_blup_PCA <- ge_blup %>%
  mutate(prcomp = list({
    mat %>% 
      apply(MARGIN = 2, FUN = function(env) ifelse(is.na(env), mean(env, na.rm = T), env)) %>% 
      prcomp() })) %>%
  mutate(PCs = list(prcomp$rotation)) %>%
  mutate(PCA = list(dist(PCs[,1:n_PC])))

## Use the ge_means to perform factor analysis, extract the loadings
## matrix, then calculate the distance.
ge_mean_FA <- S2_MET_BLUE_entries_pop %>%
  unnest() %>%
  group_by(population, trait) %>%
  do(mat = {
    select(., population, trait, line_name, environment, value) %>%
      spread(., environment, value) %>%
      as.data.frame() %>%
      remove_rownames() %>%
      column_to_rownames("line_name") %>%
      select(-trait, -population) %>%
      as.matrix() }) %>%
  # Impute using environmental mean
  mutate(mat = list(apply(mat, MARGIN = 2, FUN = function(env) 
    ifelse(is.na(env), mean(env, na.rm = T), env)))) %>%
  mutate(fa_out = list(fa(r = mat, nfactors = 2))) %>%
  # Grab the loadings
  mutate(delta = list(fa_out$loadings)) %>%
  # Distance matrix
  mutate(FA = list(dist(delta)))

## Use the PCs of the environmental covariates to form clusters
EC_PCA <- as_data_frame(expand.grid(
  population = c("all", "tp"),
  trait = unique(S2_MET_tidy_use$trait),
  EC_one = list(prcomp(t(one_year_mat))$rotation[,1:n_PC]),
  EC_multi = list(prcomp(t(multi_year_mat))$rotation[,1:n_PC])
)) %>%
  rowwise() %>%
  mutate(EC_one = list(dist(EC_one)),
         EC_multi = list(dist(EC_multi)))



# Combine the relevant distance matrices into a single DF
dist_df <- list(
  ge_mean_D,
  select(ge_blup_PCA, population, trait, PCA),
  select(ge_mean_FA, population, trait, FA),
  EC_PCA
) %>%
  reduce(full_join)



### Clustering
# For all distance matrices, create cluster objects
clust_df <- dist_df %>% 
  group_by(population, trait) %>% 
  do(D_clust = hclust(.$D[[1]], method = "average"),
     PCA_clust = hclust(.$PCA[[1]], method = "average"),
     FA_clust = hclust(.$FA[[1]], method = "average"),
     EC_one_clust = hclust(.$EC_one[[1]], method = "average"),
     EC_multi_clust = hclust(.$EC_multi[[1]], method = "average"))


# For each trait and cluster strategy, determine the within and across cluster heritability
# Build a data.frame and add the appropriate number of clusters (n_e - 1)
clust_df_tomodel <- clust_df %>% 
  gather(method, clust, -trait, -population) %>% 
  # Add k's for clusters
  by_row(function(i) seq(2, length(i$clust[[1]]$order) - 1), .to = "k") %>% 
  unnest(k, .drop = FALSE)

# Iterate over clusters and calculate the heritabilities
# Do this first for all entries in the population
cluster_herit_all <- clust_df_tomodel %>%
  filter(population == "all") %>%
  by_row(function(i) {
    
    # Cluster and cut the tree
    clusters <- cutree(i$clust[[1]], k = i$k) %>%
      data_frame(environment = names(.), cluster = .)
    
    # Assign clusters to the dataset
    S2_MET_tidy_use_clust <- left_join(S2_MET_tidy_use, clusters, "environment") %>%
      mutate(cluster = as.factor(cluster)) %>%
      filter(trait == i$trait) %>%
      droplevels()
    
    # Fit a new model with clusters
    clust_mod <- lmer(value ~ (1|line) + environment + (1|line:cluster) + 
                        (1|line:environment) + (1|line:environment:cluster), 
                      data = S2_MET_tidy_use_clust)
    
    clust_herit(clust_mod, g = "line", e = "environment", cl = "cluster") }, .to = "out_df")
    

# Repeat for just the TP
cluster_herit_tp <- clust_df_tomodel %>%
  filter(population == "tp", !str_detect(method, "EC")) %>%
  by_row(function(i) {
    
    # Cluster and cut the tree
    clusters <- cutree(i$clust[[1]], k = i$k) %>%
      data_frame(environment = names(.), cluster = .)
    
    # Assign clusters to the dataset
    S2_MET_tidy_use_clust <- left_join(S2_MET_tidy_use, clusters, "environment") %>%
      mutate(cluster = as.factor(cluster)) %>%
      filter(trait == i$trait, line_name %in% c(tp, checks)) %>%
      droplevels()
    
    # Fit a new model with clusters
    clust_mod <- lmer(value ~ (1|line) + environment + (1|line:cluster) + 
                        (1|line:environment) + (1|line:environment:cluster), 
                      data = S2_MET_tidy_use_clust)
    
    clust_herit(clust_mod, g = "line", e = "environment", cl = "cluster") }, .to = "out_df")

# Save the data
save_file <- file.path(pred_dir, "Results/cluster_heritability.RData")
save("cluster_herit_all", "cluster_herit_tp", file = save_file)

# TESTING

## Prediction
# For all clusters, drop one environment and use the remaining environments to predict
# it.
# Iterate over clusters
# First do this using clusters created by all entries
cluster_pred_acc_all <- clust_df_tomodel %>%
  filter(population == "all") %>%
  # Add core
  mutate(core = sort(rep(seq(1, n_cores), length.out = nrow(.)))) %>%
  # Split by core
  split(.$core) #%>%
  mclapply(X = ., FUN = function(core_df) {

#core_df <- cluster_pred_acc_all[[1]]

    core_df %>%
      by_row(function(i) {
        
        # Cluster and cut the tree
        clusters <- cutree(i$clust[[1]], k = i$k) %>%
          data_frame(cluster_env = names(.), cluster = .) %>%
          mutate(cluster = paste("cluster", cluster, sep = ""))
        
        # Assign clusters to the dataset
        S2_MET_BLUE_entries_clust <- 
          left_join(S2_MET_BLUE_entries, clusters, c("environment" = "cluster_env")) %>%
          mutate(cluster = as.factor(cluster)) %>%
          filter(trait == as.character(i$trait)) %>%
          droplevels()
        
        # Create a df of the the training and testing data for each cluster/environmnent
        pred_df <- clusters %>% 
          group_by(cluster, cluster_env) %>% 
          do(data = filter(S2_MET_BLUE_entries_clust, cluster == .$cluster)) %>%
          mutate(
            train_df = list(filter(data, environment != cluster_env, line_name %in% tp)), 
            test_df = list(filter(data, environment == cluster_env, line_name %in% vp))) %>%
          select(-data)
        
        # Iterate and detemine prediction accuracy
        pred_df %>% 
          by_row(function(ii) pred_acc_out(train_df = ii$train_df[[1]], 
                                           test_df = ii$test_df[[1]]),
                 .to = "acc") %>%
          select(cluster, pred_env = cluster_env, acc)
        
      }, .to = "acc_df") }, mc.cores = n_cores)
        

# Next do this using clusters created by only the tp entries
# Remove the EC clustering strategies because these are identical to above
cluster_pred_acc_tp <- clust_df_tomodel %>%
  filter(population == "tp", !str_detect(method, "EC")) %>%
  # Add core
  mutate(core = sort(rep(seq(1, n_cores), length.out = nrow(.)))) %>%
  # Split by core
  split(.$core) %>%
  mclapply(X = ., FUN = function(core_df) {
    core_df %>%
      by_row(function(i) {
    
        # Cluster and cut the tree
        clusters <- cutree(i$clust[[1]], k = i$k) %>%
          data_frame(cluster_env = names(.), cluster = .) %>%
          mutate(cluster = paste("cluster", cluster, sep = ""))
        
        # Assign clusters to the dataset
        S2_MET_BLUE_entries_clust <- 
          left_join(S2_MET_BLUE_entries, clusters, c("environment" = "cluster_env")) %>%
          mutate(cluster = as.factor(cluster)) %>%
          filter(trait == as.character(i$trait)) %>%
          droplevels()
        
        # Create a df of the the training and testing data for each cluster/environmnent
        pred_df <- clusters %>% 
          group_by(cluster, cluster_env) %>% 
          do(data = filter(S2_MET_BLUE_entries_clust, cluster == .$cluster)) %>%
          mutate(
            train_df = list(filter(data, environment != cluster_env, line_name %in% tp)), 
            test_df = list(filter(data, environment == cluster_env, line_name %in% vp))) %>%
          select(-data)
        
        # Iterate and detemine prediction accuracy
        pred_df %>% 
          by_row(function(ii) pred_acc_out(train_df = ii$train_df[[1]], 
                                           test_df = ii$test_df[[1]]),
                 .to = "acc") %>%
          select(cluster, pred_env = cluster_env, acc)
        
      }, .to = "acc_df") }, mc.cores = n_cores)

# Save the data
save_file <- file.path(pred_dir, "Results/cluster_predictions.RData")
save("cluster_pred_acc_all", "cluster_pred_acc_tp", file = save_file)
 


