## S2MET Cluster Heritability
## 
## This script will perform different clustering procedures on the S2MET data
## and calculate the within and across-cluster heritabilities.

# List of packages
packages <- c("dplyr", "tidyr", "tibble", "purrr", "readr", "stringr", "readxl", "modelr", 
              "psych", "parallel", "pbr", "purrrlyr")
 
# Set the directory of the R packages
package_dir <- NULL
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET/"
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/" 

# Geno, pheno, and enviro data
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/"
pheno_dir <- file.path(proj_dir, "Phenotype_Data/")
env_var_dir <- file.path(proj_dir, "Environmental_Variables")

geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/GBS_Genos"
pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"
env_var_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Environmental_Data"


# Other directories
fig_dir <- file.path(proj_dir, "Figures/")
pred_dir <- file.path(proj_dir, "Predictions")
entry_dir <- file.path(proj_dir, "Plant_Materials")
analysis_dir <- file.path(proj_dir, "Analysis")
result_dir <- file.path(proj_dir, "Results")




# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))
# Load environmental data
load(file.path(env_var_dir, "environmental_data_compiled.RData"))

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))


# Grab the entry names that are not checks
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)

vp <- entry_list %>% 
  filter(Class == "S2C1R") %>% 
  pull(Line)

# Find the tp and vp that are genotypes
tp_geno <- intersect(tp, row.names(s2_imputed_mat))
vp_geno <- intersect(vp, row.names(s2_imputed_mat))

# Define the checks
checks <- entry_list %>% 
  filter(Class == "Check") %>% 
  pull(Line)

entries <- entry_list %>% 
  pull(Line)

# Extract the tp and vp from the G matrix
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno, vp_geno),]

# Detect cores
n_cores <- detectCores()


# Create a new data.frame to hold the different datasets
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  group_by(trait) %>% 
  nest() %>%
  mutate(data = map(data, droplevels))


# Fit a new model with compound symmetry and get the GxE BLUPs
S2_MET_BLUPs_ge <- S2_MET_BLUEs_use %>%
  group_by(trait) %>%
  do({
    
    # Create an object for the data.frame
    ## Center and scale the data
    df <- unnest(.) %>%
      mutate(value = scale(value))
    
    # Extract the standard errors
    wts <- df$std_error^2
    
    # lmer control
    lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore",
                                calc.derivs = FALSE)
    
    # Fit the model
    form <- value ~ (1|line_name) + environment + (1|line_name:environment)
    fit <- lmer(formula = form, data = df, control = lmer_control, weights = wts)
    
    # Extract the model.frame
    mf <- model.frame(fit)
    
    # Extract the BLUPs of each genotype in each environment
    blups <- ranef(fit)$`line_name:environment` %>%
      rownames_to_column("term") %>%
      separate(term, c("line_name", "environment"), sep = ":") %>%
      rename(value = "(Intercept)")
    
    # Extract the variance components
    var_comp <- as.data.frame(VarCorr(fit))
    
    # Return a data.frame
    data_frame(fit = list(fit), var_comp = list(var_comp), BLUP = list(blups)) })



# Perform the distance calculation on the line means in each environment
ge_mean_D <- S2_MET_BLUEs_use %>% 
  unnest() %>%
  group_by(trait) %>%
  do(D = {
    dist1 <- dist_env(x = ., gen.col = "line_name", 
                      env.col = "environment", pheno.col = "value")
    # Replace missing with the mean
    dist1[is.na(dist1)] <- mean(dist1, na.rm = T)
    dist1 })



# Factor analysis using BLUEs
ge_mean_FA <- S2_MET_BLUEs_use %>% 
  group_by(trait) %>%
  mutate(BLUEs = list({
    data[[1]] %>% 
      select(line_name, environment, value) %>% 
      complete(line_name, environment) %>% 
      group_by(environment) %>% 
      mutate(value = ifelse(is.na(value), mean(value, na.rm = TRUE), value)) %>% 
      spread(environment, value) %>% 
      as.data.frame() %>% 
      remove_rownames() %>% 
      column_to_rownames("line_name") %>% 
      as.matrix() })) %>%
  do(fa_out = fa(r = .$BLUEs[[1]], nfactors = 2, rotate = "varimax")) %>%
  # Grab the loadings
  mutate(delta = list(structure(fa_out$loadings, class = "matrix"))) %>%
  # Distance matrix
  mutate(FA = list(dist(delta)))


# Number of PCs to choose for GxE BLUPs and ECs
n_PC <- 2

# Extract the GxE BLUPs and run PCA
ge_blup_PCA <- S2_MET_BLUPs_ge %>% 
  do(mat = {
    .$BLUP[[1]] %>%
      spread(environment, value) %>% 
      as.data.frame() %>%
      remove_rownames() %>% 
      column_to_rownames("line_name") %>%
      as.matrix() }) %>%
  mutate(prcomp = list({
    mat %>% 
      apply(MARGIN = 2, FUN = function(env) ifelse(is.na(env), mean(env, na.rm = T), env)) %>% 
      prcomp() })) %>%
  mutate(PCs = list(prcomp$rotation)) %>%
  mutate(PCs_use = list(PCs[,1:n_PC])) %>%
  mutate(PCA = list(dist(PCs_use)))


## Use the PCs of the environmental covariates to form clusters
EC_PCA <- as_data_frame(expand.grid(
  trait = unique(S2_MET_BLUEs$trait),
  EC_one = list(prcomp(t(one_year_mat))$rotation[,1:n_PC]),
  EC_multi = list(prcomp(t(multi_year_mat))$rotation[,1:n_PC])
)) %>%
  rowwise() %>%
  mutate(EC_one_dist = list(dist(EC_one)),
         EC_multi_dist = list(dist(EC_multi)))


# Combine all matrices, including PCAs and distance matrices
clust_method_df <- list(
  ge_mean_D,
  ge_blup_PCA,
  ge_mean_FA,
  EC_PCA) %>%
  reduce(full_join)

# # Save this
# save_file <- file.path(result_dir, "S2_MET_cluster_df.RData")
# save("clust_method_df", file = save_file)


# Combine the relevant distance matrices into a single DF
dist_df <- clust_method_df %>%
  select(trait, D, PCA, FA, EC_one_dist, EC_multi_dist)



# Mapping function
map_hclust <- function(x) map(x, hclust, method = "average")

# For all distance matrices, create cluster objects
clust_df <- dist_df %>% 
  group_by(trait) %>% 
  summarize_all(funs(map_hclust))

# Set a maximum for the number of clusters
max_K <- 20




# For each trait and cluster strategy, determine the within and across cluster heritability
# Build a data.frame and add the appropriate number of clusters (n_e - 1)
clust_df_tomodel <- clust_df %>% 
  gather(method, clust, -trait) %>% 
  # Add k's for clusters
  by_row(function(i) seq(2, max_K), .to = "k") %>% 
  unnest(k, .drop = FALSE)

cluster_herit_all <- clust_df_tomodel %>%
  mutate(core = sort(rep(seq(n_cores), length.out = nrow(.)))) %>%
  split(.$core) %>%
  # Apply by core
  mclapply(function(core_df) {
    core_df %>%
      group_by(trait, method, k) %>%
      do({
    
        i <- .
        
        # Cluster and cut the tree
        clusters <- cutree(i$clust[[1]], k = i$k) %>%
          data_frame(environment = names(.), cluster = .)
    
        # Assign clusters to the dataset
        S2_MET_BLUEs_clusters <- left_join(S2_MET_BLUEs, clusters, "environment") %>%
          mutate(cluster = as.factor(cluster),
                 ge = interaction(line_name, environment),
                 gc = interaction(line_name, cluster)) %>%
          filter(trait == i$trait) %>%
          droplevels()
        
        # Calculate clusters, environments, reps
        n_e <- xtabs(~ line_name + environment, data = S2_MET_BLUEs_clusters) %>%
          ifelse(. > 1, 1, .) %>%
          rowSums() %>% 
          harm_mean()
        
        # Now clusters
        n_c <- xtabs(~ line_name + cluster, data = S2_MET_BLUEs_clusters) %>%
          ifelse(. > 1, 1, .) %>%
          rowSums() %>% 
          harm_mean()
        
        # Now replicates
        n_r <- xtabs(~ line_name + environment, data = S2_MET_BLUEs_clusters) %>% 
          harm_mean()
        
        wts <- S2_MET_BLUEs_clusters$std_error^2
    
        # List of models
        forms <- formulas(~ value,
                          no_gc = ~ (1 | line_name) + environment + (1|ge) + (1|cluster/ge),
                          full = ~ (1 | line_name) + environment + (1|ge) + (1|gc) + (1|cluster/ge) )
        
        # lmer control
        lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore",
                                    calc.derivs = FALSE)
        
        # Fit all the models
        fits <- fit_with(data = S2_MET_BLUEs_clusters, lmer, forms, control = lmer_control,
                         weights = wts) %>%
          data_frame(mod = names(.), fit = .)
        
        # Extract variance components and logLikihoods
        fits_summs <- fits %>% 
          group_by(mod) %>%
          summarize(var_comp = map(fit, VarCorr) %>% map(as.data.frame),
                    llik = map(fit, logLik))
        
        # Calculate heritability
        ## Define some expressions
        herits <- list(
          exp_a = "line_name / (line_name + (gc / n_c) + (ge:cluster / n_e) + (Residual / (n_r * n_e)))",
          exp_w = "(line_name + gc) / (line_name + gc + (n_c * ((ge:cluster / n_e) + (Residual / (n_r * n_e)))))"
          ) %>%
          map_df(function(x) herit(object = subset(fits, mod == "full")$fit[[1]], exp = x,
                                   n_r = n_r, n_c = n_c, n_e = n_e))
    
        herit_summ <- herits %>%
          rename_all(str_replace_all, pattern = "exp", replacement = "herit")
        
        # Return data_frame
        bind_cols(data_frame(fit_summ = list(fits_summs)), herit_summ)
     
      })
    
  }, mc.cores = n_cores)


save_file <- file.path(result_dir, "S2_MET_cluster_heritability.RData")
save("cluster_herit_all", file = save_file)





