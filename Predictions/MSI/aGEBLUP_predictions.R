## R Script to run the additive GE-BLUP model (i.e. aGE-BLUP)

# List of packages
packages <- c("tidyverse", "stringr", "sommer", "parallel")

# Directory to save the files
MSI_proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/Predictions/MSI"  
  
# Set the directory of the R packages
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))

# Load the starting data
load(file.path(MSI_proj_dir, "aGE-BLUP_start_data_MSI.RData"))

# Find the number of cores
n_cores <- detectCores()

# Divide the df into cores
model_cores <- expand.grid(environment = env_pred, trait = traits, stringsAsFactors = FALSE) %>%
  mutate(core = sort(rep(seq(n_cores), length = nrow(.)))) %>%
  split(.$core)


# Iterate over all traits and environments
aGEBLUP2_E1_results <- mclapply(X = model_cores, FUN = function(df_core) {
  
  df_core %>%
    group_by(environment, trait) %>%
    do(pred_out = {
      
      # Training and prediction phenotypes
      pheno_train <- filter(pheno, line_name %in% tp, environment != .$environment, trait == .$trait)
      pheno_pred <- filter(pheno, line_name %in% vp, environment == .$environment, trait == .$trait)
      
      # Create the y vector
      y_train <- pheno_train %>% 
        dplyr::select(line_name, environment, value) %>% 
        unite(gen_env, line_name, environment, sep = ":") %>% 
        as.data.frame() %>%
        remove_rownames() %>% 
        column_to_rownames("gen_env") %>% 
        as.matrix()
      
      # Create the fixed effect vector (environment)
      X <- model.matrix(~ 1, data = pheno_train)
      
      # Create the incidence matrix for the genotypes
      Z_g <- model.matrix(~ -1 + line_name, data = pheno_train) %>%
        structure(dimnames = list(NULL, str_replace(colnames(.), "line_name", "")))
      
      # Same for environments
      Z_e <- model.matrix(~ -1 + environment, data = pheno_train) %>% 
        structure(dimnames = list(NULL, str_replace(colnames(.), "environment", "")))
      
      
      solve_out <- mmer(Y = y_train, X = X, Z = list(G = list(Z = Z_g, K = G), 
                                                     E = list(Z = Z_e, K = E_1)))
      
      
      # Convert to data.frames
      u_hat_df <- randef(solve_out) %>%
        map(as.matrix) %>%
        map(~data.frame(id = row.names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE)) %>% 
        map(~separate(data = ., col = id, into = c("line_name", "environment"), sep = ":", 
                      fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "left", "right"))) %>%
        map(~mutate_if(., .predicate = is.character, .funs = str_replace, pattern = "line_name|environment", 
                       replacement = ""))
      
      # Measure accuracy
      pheno_pred %>% 
        mutate_at(vars(environment, line_name), parse_character) %>% 
        left_join(., dplyr::select(u_hat_df$G, -environment), by = "line_name") %>% 
        left_join(., dplyr::select(u_hat_df$E, -line_name), by = "environment") %>% 
        mutate(pred = u_hat.x + u_hat.y) %>%
        dplyr::select(-u_hat.x, -u_hat.y)
      
    }) }, mc.cores = n_cores)

# ## Use the E_2 relationship matrix
# 
# # Iterate over all traits and environments
# aGEBLUP2_E2_results <- expand.grid(environment = env_pred, trait = traits, stringsAsFactors = FALSE) %>%
#   group_by(environment, trait) %>%
#   do(pred_out = {
#     
#     # Training and prediction phenotypes
#     pheno_train <- filter(pheno, line_name %in% tp, environment != .$environment, trait == .$trait)
#     pheno_pred <- filter(pheno, line_name %in% vp, environment == .$environment, trait == .$trait)
#     
#     # Create the y vector
#     y_train <- pheno_train %>% 
#       dplyr::select(line_name, environment, value) %>% 
#       unite(gen_env, line_name, environment, sep = ":") %>% 
#       as.data.frame() %>%
#       remove_rownames() %>% 
#       column_to_rownames("gen_env") %>% 
#       as.matrix()
#     
#     # Create the fixed effect vector (environment)
#     X <- model.matrix(~ 1, data = pheno_train)
#     
#     # Create the incidence matrix for the genotypes
#     Z_g <- model.matrix(~ -1 + line_name, data = pheno_train) %>%
#       structure(dimnames = list(NULL, str_replace(colnames(.), "line_name", "")))
#     
#     # Same for environments
#     Z_e <- model.matrix(~ -1 + environment, data = pheno_train) %>% 
#       structure(dimnames = list(NULL, str_replace(colnames(.), "environment", "")))
#     
#     
#     solve_out <- mmer(Y = y_train, X = X, Z = list(G = list(Z = Z_g, K = G), 
#                                                    E = list(Z = Z_e, K = E_2)), 
#                       n.cores = n_cores)
#     
#     
#     # Convert to data.frames
#     u_hat_df <- randef(solve_out) %>%
#       map(as.matrix) %>%
#       map(~data.frame(id = row.names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE)) %>% 
#       map(~separate(data = ., col = id, into = c("line_name", "environment"), sep = ":", 
#                     fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "left", "right"))) %>%
#       map(~mutate_if(., .predicate = is.character, .funs = str_replace, pattern = "line_name|environment", replacement = ""))
#     
#     # Measure accuracy
#     pheno_pred %>% 
#       mutate_at(vars(environment, line_name), parse_character) %>% 
#       left_join(., dplyr::select(u_hat_df$G, -environment), by = "line_name") %>% 
#       left_join(., dplyr::select(u_hat_df$E, -line_name), by = "environment") %>% 
#       mutate(pred = u_hat.x + u_hat.y) %>%
#       dplyr::select(-u_hat.x, -u_hat.y)
#     
#   })
# 

# Save the results
save_file <- file.path(MSI_proj_dir, "aGE-BLUP_results.RData")
list_to_save <- c("aGEBLUP2_E1_results")
save(list = list_to_save, file = save_file)