## R Script to run the additive GE-BLUP model (i.e. aGE-BLUP)

# List of packages
packages <- c("tidyverse", "stringr", "sommer")

# Directory to save the files
MSI_proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/Predictions/MSI"  

# Set the directory of the R packages
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))

# Load the starting data
load(file.path(MSI_proj_dir, "S2MET_MSI_prediction_material.RData"))

# Find the number of cores
n_cores <- detectCores()

# Create a df of traits
results_df <- data_frame(trait = traits, results = replicate(length(traits), data_frame(environment = env_pred), simplify = FALSE))


# Create a df of traits
i_geblup2_results <- data_frame(trait = traits, 
                                results = replicate(length(traits), data_frame(environment = env_pred), simplify = FALSE))

# Iterate over all traits
for (i in seq(nrow(i_geblup2_results))) {
  
  i_geblup2_results$results[[i]] <- i_geblup2_results$results[[i]] %>%
    group_by(environment) %>%
    do(pred_out = {
      
      # Training and prediction phenotypes
      pheno_train <- filter(pheno, line_name %in% tp, environment != .$environment, trait == i_geblup2_results$trait[i])
      pheno_pred <- filter(pheno, line_name %in% vp, environment == .$environment, trait == i_geblup2_results$trait[i])
      
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
      
      Z_ge <- model.matrix(~ -1 + line_name:environment, data = pheno_train) %>%
        structure(dimnames = list(NULL, str_replace_all(colnames(.), c("line_name" = "", "environment" = "")) %>%
                                    str_split(":") %>% 
                                    map(rev) %>% 
                                    map_chr(str_c, collapse = ":") ))
      
      
      # solve_out <- mmer(Y = y_train, X = X, Z = list(G = list(Z = Z_g, K = G), 
      #                                                E = list(Z = Z_e, K = E_1), 
      #                                                GE = list(Z = Z_ge, K = O_1)))
      
      # The #2 matrices
      solve_out <- mmer(Y = y_train, X = X, Z = list(G = list(Z = Z_g, K = G),
                                                     E = list(Z = Z_e, K = E_2),
                                                     GE = list(Z = Z_ge, K = O_2)))
      
      
      # Convert to data.frames
      u_hat_df <- randef(solve_out) %>%
        map(as.matrix) %>%
        map(~data.frame(id = row.names(.), u_hat = ., row.names = NULL, stringsAsFactors = FALSE)) %>% 
        map(~separate(data = ., col = id, into = c("environment", "line_name"), sep = ":", 
                      fill = ifelse(str_detect(.$id[1], "[A-Z]{3}[0-9]{2}"), "right", "left"))) %>%
        map(~mutate_if(., .predicate = is.character, .funs = str_replace, pattern = "line_name|environment", replacement = ""))
      
      # Measure accuracy
      pheno_pred %>% 
        mutate_at(vars(environment, line_name), parse_character) %>% 
        left_join(., dplyr::select(u_hat_df$G, -environment), by = "line_name") %>% 
        left_join(., dplyr::select(u_hat_df$E, -line_name), by = "environment") %>% 
        left_join(., u_hat_df$GE, by = c("line_name", "environment")) %>% 
        mutate(pred = u_hat.x + u_hat.y + u_hat) %>%
        dplyr::select(-u_hat.x, -u_hat.y, -u_hat)
      
    })
  
}


# Save the results
# save_file <- file.path(MSI_proj_dir, "iGEBLUP_results.RData")
save_file <- file.path(MSI_proj_dir, "iGEBLUP_results2.RData")
list_to_save <- c("i_geblup2_results")
save(list = list_to_save, file = save_file)