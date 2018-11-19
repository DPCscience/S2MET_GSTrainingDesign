## S2MET
## 
## Analysis of phenotypic data using environmental covariates
## 
## Author: Jeff Neyhart
## Last modified: November 9, 2018
## 

### Run on MSI
# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"
source(file.path(repo_dir, "source_MSI.R"))
library(lme4qtl)


# Number of cores
n_core <- 4
n_core <- detectCores()

# # Run the source script
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# library(lme4qtl)

# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))

# Modify the BLUEs for predictions
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% tp) %>%
  mutate_at(vars(environment:line_name), as.factor)


## Construct different covariance matrices using distance, covariates, or phenotypic correlation
env_cov_mats <- clust_method_df %>% 
  filter(population == "tp", !map_lgl(cov, is.null)) %>% 
  select(trait, model, env_cov_mat = cov) %>%
  arrange(trait, model)



## Test models
set.seed(1010)
pheno_data <- S2_MET_BLUEs_use %>% 
  filter(trait == "GrainYield") %>%
  droplevels() %>%
  filter(environment %in% sample(unique(.$environment), 5)) %>%
  # Copy columns
  mutate(line_name1 = line_name, environment1 = environment, ge = interaction(line_name, environment))

pheno_data <- S2_MET_BLUEs_use %>%
  droplevels()


wts <- test_data$std_error^2

# Create genotypic and environmal incidence matrices
Zg <- model.matrix(~ -1 + line_name, test_data)
Ze <- model.matrix(~ -1 + environment, test_data)

# Subset the K matrix
K1 <- K[tp_geno, tp_geno]



# First build a random effect model with g + e + ge
base_models <- pheno_data %>%
  group_by(trait) %>%
  do(fit_base = relmatLmer(value ~ 1 + (1|line_name) + (1|environment) + (1|line_name:environment), data = ., weights = .$wts))


## Split the traits and models by core
env_cov_mats_split <- env_cov_mats_use %>%
  assign_cores(n_core) %>%
  split(.$core)


env_cov_models_out <- mclapply(X = env_cov_mats_split, FUN = functions(core_df) {
  
  # Iterate over elements in the core_df
  models_out <- vector("list", nrow(core_df))
  
  for (i in seq_along(models_out)) {
    tr <- core_df$trait[i]
    
    # Subset the pheno data
    pheno_data1 <- pheno_data %>%
      filter(trait == tr) %>%
      droplevels()
    
    
    
  }
  
  
  
}


## Now iterate over each trait and covariance matrix and fit the GxE model
models_out <- vector("list", nrow(env_cov_mats_use))








# Create the g:e covariance matrix
# Ke <- subset(env_cov_mats_use, model == "pheno_dist", env_cov_mat, drop = T)[[1]]
Ke <- subset(env_cov_mats_use, model == "MYEC_Top1EC_Fstat", env_cov_mat, drop = T)[[1]]

Kge <- (Zg %*% K1 %*% t(Zg)) * (Ze %*% Ke %*% t(Ze))
dimnames(Kge) <- list(test_data$ge, test_data$ge)





# Fit a ge term with the Kge covariance matrix
fit_ge <- relmatLmer(value ~ 1 + (1|line_name) + (1|environment) + (1|line_name:environment) + (1|ge), 
                    data = test_data, weights = wts, relmat = list(ge = Kge))

lr_test(model1 = fit_base, model2 = fit_ge)








