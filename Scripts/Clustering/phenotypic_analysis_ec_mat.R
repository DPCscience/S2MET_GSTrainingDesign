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


# Number of cores
n_core <- 4
n_core <- detectCores()



# # Local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))


library(lme4qtl)
library(sommer)

# Load the clustering results
load(file.path(result_dir, "distance_method_results.RData"))

# Modify the BLUEs for predictions
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% tp) %>%
  mutate_at(vars(environment:line_name), as.factor)


## Construct different covariance matrices using distance, covariates, or phenotypic correlation
env_cov_mats <- env_rank_df %>% 
  filter((str_detect(model, "All|IPCA") & mat_set == "Malosetti") | model %in% c("pheno_dist")) %>% 
  dplyr::select(trait, set, model, env_cov_mat = cov) %>%
  arrange(trait, model)



## Test models
set.seed(1010)
pheno_data <- S2_MET_BLUEs_use %>% 
  filter(trait == "GrainYield") %>%
  filter(line_name %in% tp_geno) %>%
  droplevels() %>%
  filter(environment %in% sample(unique(.$environment), 10)) %>%
  # Copy columns
  mutate(ge = interaction(line_name, environment, drop = TRUE), ge1 = ge)

## Nest the data and combine with the relmat data
pheno_data_tomodel <- nest(group_by(pheno_data, trait)) %>%
  inner_join(env_cov_mats, .)


## Relationsiop matrix of genotypes
Kg <- K[tp_geno, tp_geno]

## Fit a full model
full_models <- group_by(pheno_data, trait) %>%
  do({
    df <- .
    
    control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
    fit <- lmer(value ~ (1|line_name) + (1|environment) + (1|ge), data = df, control = control, weights = df$std_error^2)
    
    as.data.frame(VarCorr(fit)) %>% 
      dplyr::select(term = grp, variance = vcov)
  })


## Iterate over cov mats and fit models
# cov_models_complete <- pheno_data_tomodel %>%
#   filter(set == "complete") %>%
#   group_by(trait, set, model) %>%
#   do({
#     row <- .
    
cov_models_complete <- pheno_data_tomodel %>%
  filter(set == "complete") %>%
  assign_cores(n_core = n_core) %>%
  split(.$core) %>%
  mclapply(X = ., mc.cores = n_core, FUN = function(core_df) {
    
    results_out <- vector("list", nrow(core_df))
    
    for (i in seq_along(results_out)) {
      
    row <- core_df[i,]
    df <- unnest(row, data) %>% droplevels()
    
    # Make model matrices
    mf <- model.frame(value ~ line_name + environment + ge + ge1, df)
    y <- model.response(mf)
    X <- model.matrix(~ 1, mf)
    Zg <- model.matrix(~ -1 + line_name, mf)
    Ze <- model.matrix(~ -1 + environment, mf)
    Zge <- model.matrix(~ -1 + ge, mf)
    
    G <- Zg %*% Kg %*% t(Zg)
    Ke <- row$env_cov_mat[[1]][levels(df$environment), levels(df$environment)]
    E <- Ze %*% Ke %*% t(Ze)
    
    Kge <- G * E; dimnames(Kge) <- list(df$ge, df$ge)
    # Identity matrix for GE
    Kge1 <- diag(ncol(Kge)); dimnames(Kge1) <- dimnames(Kge)
    colnames(Zge) <- colnames(Kge)
    
    fit <- mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = Kg), e = list(Z = Ze, K = Ke), ge = list(Z = Zge, K = Kge), ge1 = list(Z = Zge, K = Kge1)))
    # # Fit just gxe
    # fit1 <- mmer(Y = y, X = X, Z = list(ge = list(Z = Zge, K = Kge), ge1 = list(Z = Zge, K = Kge1)))
    
    ## Return variance components
    results_out[[i]] <- as.data.frame(lapply(fit$var.comp, as.vector))
    
  }
  
  mutate(core_df, out = results_out)
  
  })

  # as.data.frame(lapply(fit$var.comp, as.vector))
  #   
  # })
  # 







