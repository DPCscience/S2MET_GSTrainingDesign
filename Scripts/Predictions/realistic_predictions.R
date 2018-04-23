## Predictions following a realistic breeding scenario
## 
## Author: Jeff Neyhart
## Last Updated: April 19, 2018
## 
## This script will evaluate prediction accuracy under a "realistic" breeding
## scenario, where each year is predicted using the data presumed to be available
## for that year. I also try some other methods that, while violating this constraint,
## may be close approximations of what might be present in a breeding program
## 
## 

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"
source(file.path(repo_dir, "source_MSI.R"))

# Load packages
packages <- c("lme4")
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


# # The head directory
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# 
# library(lme4)


## Number of cores
n_core <- ifelse(Sys.info()["sysname"] == "Windows", 1, detectCores())

# Load the distances matrices and results
load(file.path(result_dir, "distance_methods_results.RData"))



# Modify the BLUEs for predictions
S2_MET_BLUEs_use <- S2_MET_BLUEs




### We will consider each year individually, predicting 2016 using 2015 data, and
### 2017 using 2015 and 2016 data
### 


## Create a list with validation years and the name and prediction years
## as the elements
year_list <- list(year2016 = list(val_year = 2016, train_year = c(2015)),
                  year2017 = list(val_year = 2017, train_year = c(2015, 2016)))

# Subset the distance methods
dist_method_use <- c("great_circle_dist", "ec_multi_PCA_dist")

pred_env_rank_use <- pred_env_rank_random$all %>%
  filter(dist_method %in% dist_method_use,
         population == "all")

# Map over the years
breeding_dataset <- year_list %>%
  map(function(yrs) {
    
    # Assign the years to new objects
    train_year <- yrs$train_year
    val_year <- yrs$val_year
    
    ## Extract the training data
    ## This is presumably the only data that would be available for making decisions
    train_data <- S2_MET_BLUEs_use %>% 
      filter(year %in% train_year,
             line_name %in% tp)
    
    ## Extract the validation data
    val_data <- S2_MET_BLUEs_use %>%
      filter(year == val_year,
             line_name %in% vp_geno)
    
    
    ## Split up the validation environments per trait, then create a list that
    ## stores each validation environment along with the available prediction environments
    ## 
    val_data_nest <- val_data %>%
      group_by(trait, environment) %>% 
      nest(.key = "val_data") %>%
      rename(val_environment = environment)
    
    val_data_nest1 <- select(val_data_nest, -val_data) %>% 
      left_join(., distinct(select(train_data, trait, environment)), by = "trait") %>% 
      group_by(trait, val_environment) %>% 
      do(train_environment = .$environment) %>% 
      left_join(val_data_nest, ., by = c("trait", "val_environment"))
    
    
    # Order the training environments according to the different distance metrics:
    ## 1. Great circle distance
    ## 2. Environmental covariates
    ## 4. All environments
    
    ## Select the environments to use based on the different optimizations:
    ## 1. The closest environment
    ## 2. The 5 closest environments
    ## 3. LRT picked environments
    
    # Create an empty list
    val_df_list <- vector("list", nrow(val_data_nest1))
    
    for (i in seq(nrow(val_data_nest1))) {
      val_df <- val_data_nest1[i,]
      
      # Get the training environments
      val_df_env <- as.character(val_df$train_environment[[1]])
      
      # Extract the ranks and select only the available training environments
      val_df1 <- pred_env_rank_use %>% 
        filter(environment == val_df$val_environment, trait == val_df$trait) %>% 
        mutate(train_environment = map(env_rank, ~select(., val_df_env) %>% sort %>% unlist)) %>%
        select(-env_rank)
      
      ## Determine the optimal vector of environments
      train_data1 <- train_data %>%
        filter(environment %in% val_df_env, trait == val_df$trait) %>%
        rename(env = environment)
      
      # Iterate over the possible criteria
      criteria <- eval(formals(fun = optim_env)$criterion)
      
      ## Feed this to the optim_env function
      ## Also add a row for using all of the environments
      val_df2 <- val_df1 %>% 
        group_by(dist_method) %>%
        do({
          df <- .
          envs <- names(.$train_environment[[1]])
          opt_env <- map(criteria, ~optim_env(environments = envs, data = train_data1, 
                                              criterion = .))
          # Return the modified DF
          cbind(select(df, -train_environment), criteria) %>% 
            as_data_frame() %>% 
            mutate(train_environment = opt_env)
        }) %>%
        ungroup() %>%
        bind_rows(., mutate(val_df1[1,], dist_method = "all_env", 
                            train_environment = map(train_environment, names)))
      
      val_df_list[[i]] <- left_join(select(val_df, -train_environment), val_df2, 
                                    by = c("val_environment" = "environment", "trait"))
      
    }
    
    ## Bind rows and return
    bind_rows(val_df_list)
    
  })
      



## For each validation environment, add a random order of environments
# Number of random samples
n_samples <- 50


breeding_dataset_random <- breeding_dataset %>%
  map(~{

    # Pull out the list of all environments
    train_envs <- filter(., dist_method == "all_env") %>% 
      distinct(trait, train_environment) %>% 
      group_by(trait) %>% 
      slice(1) 

    # Create a vector of samples
    train_envs %>% 
      do(random = data_frame(
        train_environment = replicate(n = n_samples, sample(.$train_environment[[1]]),
                                      simplify = FALSE))) %>% 
      unnest()
    
  })


  
### With the dataset created, now iterate over each year and then each validation
### environment, select the best dataset to use (or all data), and run predictions.
### 

# Prep the BLUEs for subsetting
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno, vp_geno)) %>% 
  mutate_at(vars(environment, line_name), as.factor) %>%
  rename(env = environment)
  


breeding_prediction_results <- breeding_dataset %>%
  map(function(yrs) {

    ## Select relevant columns and simultaneously iterate
    predict_out <- yrs %>%
      select(trait, val_data, train_environment) %>%
      as.list() %>%
      pmap(~{ # Use ~{ to write a function (remember for pmap, the list elements are numbered by default)

        # Subset the training data
        train_data <- S2_MET_BLUEs_use %>%
          filter(trait == ..1, env %in% ..3)

        # Predict and return the results
        gblup(K = K, train = train_data, test = ..2, bootreps = 1000)

      })

    # Combine these results with selected yrs data and return
    yrs %>%
      select(-val_data) %>%
      mutate(results_out = predict_out)

  })






## Iterate over each prediction environment and then add environments, one-by-one,
## to assess prediction accuracy
## 

# First do this for the intentionally ordered training environments
# We need to only use the "all_env" rows
breeding_dataset_tomodel <- breeding_dataset %>%
  bind_rows() %>%
  filter(dist_method == "all_env") %>%
  # Add cores
  assign_cores(n_core) %>%
  split(.$core)

breeding_prediction_iterative_results <- breeding_dataset_tomodel %>%
  mclapply(X = ., FUN = function(core_df) {
    
    # Iterate over the core_df rows
    results_out <- core_df %>%
      select(trait, val_data, train_environment) %>%
      pmap(~{
        
        # Reassign data
        tr <- ..1
        val <- ..2
        te <- ..3
        
        # Accumulate environments and use them to predict
        # accumulate, used this way, will add environments one-by-one
        te_acc <- accumulate(te, ~c(.x, .y)) 
        
        map_df(te_acc, ~{
            
            # Subset the training data
            train_data <- filter(S2_MET_BLUEs_use, trait == tr, env %in% .)
            
            # Predict and return the results
            gblup(K = K, train = train_data, test = val, bootreps = 1000)$boot
            
          }) %>%
          mutate(train_envs = te_acc)
        
      })
    
    ## Add the results back to the core_df
    core_df %>% 
      mutate(results_out = results_out) %>% 
      select(-val_data, -core)
    
  }, mc.cores = n_core)


## Now do the same thing for the random environments
breeding_prediction_iterative_results_random <- breeding_dataset_random %>% 
  map(~{ # Perform the operation for each year individually

    assign_cores(., n_core) %>%
      split(.$core) %>%
      # Parallelize
      mclapply(X = ., FUN = function(core_df) {
        
        # Iterate over the core_df rows
        results_out <- core_df %>%
          select(trait, train_environment) %>%
          pmap(~{
            
            # Reassign data
            tr <- ..1
            te <- ..2
            
            # Accumulate environments and use them to predict
            # accumulate, used this way, will add environments one-by-one
            te_acc <- accumulate(te, ~c(.x, .y)) 
            
            predict_out <- map(te_acc, ~{
              
              # Subset the training data
              train_data <- filter(S2_MET_BLUEs_use, trait == tr, env %in% .)
               
              # Predict and return the results
              gblup(K = K, train = train_data, bootreps = 1000)$pgv %>%
                filter(line_name %in% vp_geno)
              
            })
            
            # Create a data_frame
            data_frame(train_env = te_acc, predict_out = predict_out)
            
          })
        
        ## Add the results back to the core_df
        core_df %>% 
          mutate(results_out = results_out) %>% 
          select(-core)
        
      }, mc.cores = n_core) %>% bind_rows()
    
  })
            

## Save
save_file <- file.path(result_dir, "realistic_breeding_predictions.RData")
save("breeding_dataset", "breeding_prediction_results", "breeding_dataset_random"
     "breeding_prediction_iterative_results_random", "breeding_prediction_iterative_results",
     file = save_file)









