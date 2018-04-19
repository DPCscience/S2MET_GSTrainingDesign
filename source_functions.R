## S2MET Functions
## 
## A script with useful functions used in the prediction analysis of the S2MET
## 



## Other/utility functions
# A function to assign cores to a data.frame
assign_cores <- function(df, n_core) {
  df$core <- sort(rep(seq(n_core), length.out = nrow(df)))
  return(df)
}

# A function to replace a character vector and convert it to a factor
as_replaced_factor <- function(x, replacement) {
  x_repl <- str_replace_all(string = x, pattern = replacement)
  factor(x_repl, levels = replacement)
}


## Run a likelihood ratio test
lr_test <- function(model1, model2) {
  
  stopifnot(class(model1) %in% c("lmerMod", "merModLmerTest"))
  stopifnot(class(model2) %in% c("lmerMod", "merModLmerTest"))
  
  model_list <- list(model1 = model1, model2 = model2)
  
  # Degrees of freedom
  df_list <- sapply(X = model_list, FUN = df.residual)
  fuller_model <- names(which.min(df_list))
  red_model <- names(which.max(df_list))
  
  ## Get the log-likelihoods and store as list
  ll_list <- sapply(X = model_list, FUN = logLik)
  
  # Calculate the likelihood ratio
  lr <- -2 * (ll_list[red_model] - ll_list[fuller_model])
  # Calculate pvalue
  df <- df_list[red_model] - df_list[fuller_model]
  p_value <- pchisq(q = lr, df = df, lower.tail = FALSE)
  
  # Export a data.frame
  data.frame(
    fuller_model = fuller_model,
    df = df,
    statistic = lr,
    p_value = p_value,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

}






## Bootstrap a correlation coefficient
boot_cor <- function(x, y, boot.reps = 1000, alpha = 0.05) {
  
  # Error
  boot.reps <- as.integer(boot.reps)
  
  # Prob must be between 0 and 1
  alpha_check <- alpha > 0 | alpha < 1
  
  if (!alpha_check)
    stop("'alpha' must be between 0 and 1.")
  
  # Define a function for the correlation
  boot.cor <- function(input.data, i) {
    rep_data <- input.data[i,]
    return(cor(rep_data[,1], rep_data[,2]))
  }
  
  
  # First calculate the base statistic
  base_cor <- suppressWarnings(cor(x, y))
  
  # If the correlation is not NA, proceed
  if (!is.na(base_cor)) {
    
    # Perform the bootstrapping
    boot_results <- boot(data = cbind(x, y), statistic = boot.cor, R = boot.reps)
    
    # Standard error
    se <- sd(boot_results$t)
    # Bias
    bias <- mean(boot_results$t) - base_cor
    
    
    # Confidence interval
    ci_upper <- quantile(boot_results$t, 1 - (alpha / 2))
    ci_lower <- quantile(boot_results$t, (alpha / 2))
    
  } else {
    
    se <- bias <- ci_lower <- ci_upper <- NA
    
  }
  
  # Assemble list and return
  data.frame(cor = base_cor, se = se, bias = bias,
             ci_lower = ci_lower, ci_upper = ci_upper, row.names = NULL)
}



## Heritability functions
## Write a function to subset a distance matrix using a set of environments
subset_env <- function(dist, envs) {
  as.matrix(dist) %>%
    .[envs, envs] %>%
    as.dist()
}

## Define a function to find the harmonic mean of a variable given an xtabs object
## 'adjust' is used for reps
harm_mean.xtabs <- function(x, variable, adjust = FALSE) {
  # Get the dimension names of the xtabs
  tab_names <- names(dimnames(x))
  # Position of the variable in the array
  variable_pos <- match(c("line_name", variable), tab_names)
  
  # Find the number of "variable" in which the jth individual was observed
  variable_sum <- apply(X = x, MARGIN = variable_pos, FUN = sum)
  variable_sum_adj <- ifelse(variable_sum > 1, 1, variable_sum)
  
  # If all tab names are called, calculate the harmonic mean of the reps
  if (adjust) {
    return(harm_mean(variable_sum_adj))
    
  } else {
    return(harm_mean(rowSums(variable_sum_adj)))
    
  }
}


# Function to calculate heritability across clusters of environments
cluster_heritability <- function(object, breakup_env = TRUE) {
  # Extract the data from the model object
  mf <- model.frame(object)
  # Get the names of the mf and remove the reponse and weights
  varnames_random <- setdiff(names(mf), attr(terms(mf), "varnames.fixed"))
  
  # If breakup_env is TRUE, but location + year are not in the names, error out
  # Also set the formula for the xtabs
  if (breakup_env) {
    stopifnot(all(c("location", "year") %in% varnames_random))
  } else {
    stopifnot("environment" %in% varnames_random)
  }
  
  # Set the formula for the xtabs
  xtabs_form <- as.formula(paste("~", paste(varnames_random, collapse = " + ")))
  plot_table <- xtabs(formula = xtabs_form, data = mf)
  
  # Calculate the harmonic mean of each of location, year, cluster, and rep
  variable_harm_mean <- lapply(setdiff(varnames_random, "line_name"), FUN = harm_mean.xtabs, x = plot_table)
  variable_harm_mean <- setNames(variable_harm_mean, setdiff(varnames_random, "line_name"))
  
  # Now for reps
  rep_harm_mean <- harm_mean.xtabs(x = plot_table, variable = names(variable_harm_mean), adjust = T)
  
  ## Set the expressions to use in calculating across- and within-cluster heritability
  if (breakup_env) {
    exp_a <- "line_name / (line_name + (line_name:cluster / n_c) + (line_name:location:cluster / n_l) + (line_name:year / n_y) + (line_name:year:cluster / (n_c * n_y)) + (line_name:location:year:cluster / (n_l * n_y)) + (Residual / (n_r * n_l * n_y)))"
    exp_w <- "(line_name + line_name:cluster) / (line_name + line_name:cluster + (n_c * ( (line_name:location:cluster / n_l) + (line_name:location:year:cluster / (n_l * n_y)) +
    (Residual / (n_r * n_l * n_y)) )) + (line_name:year / n_y) + (line_name:year:cluster / (n_c * n_y)))"
    
    # Calculate heritability
    herit_list = list(across = exp_a, within = exp_w) %>%
      map_df(~herit(object = object, exp = ., n_l = variable_harm_mean$location,
                    n_y = variable_harm_mean$year, n_c = variable_harm_mean$cluster,
                    n_r = rep_harm_mean))
    
  } else {
    exp_a <- "line_name / (line_name + (line_name:cluster / n_c) + (line_name:environment:cluster / n_e) + (Residual / (n_r * n_e)))"
    
    exp_w <- "(line_name + line_name:cluster) / (line_name +  line_name:cluster + (n_c * ((line_name:environment:cluster / n_e) + (Residual / (n_r * n_e)))))"
    
    # Calculate heritability
    herit_list = list(across = exp_a, within = exp_w) %>%
      map_df(~herit(object = object, exp = ., n_e = variable_harm_mean$environment,
                    n_c = variable_harm_mean$cluster, n_r = rep_harm_mean))
    
  }
  
  # Return the heritability
  return(herit_list) 
  
}


## Function for calculating heritability in a training set
training_heritability <- function(object) {
  # Extract the data from the model object
  mf <- model.frame(object)
  # Get the names of the mf and remove the reponse and weights
  varnames_random <- setdiff(names(mf), attr(terms(mf), "varnames.fixed"))
  
  # Set the formula for the xtabs
  xtabs_form <- as.formula(paste("~", paste(varnames_random, collapse = " + ")))
  plot_table <- xtabs(formula = xtabs_form, data = mf)

  # Calculate the harmonic mean of reps
  rep_harm_mean <- harm_mean.xtabs(x = plot_table, variable = setdiff(varnames_random, "line_name"), adjust = TRUE)
  
  # Calculate the harmonic mean of environment, if present
  if ("environment" %in% varnames_random) {
    variable_harm_mean <- lapply(setdiff(varnames_random, "line_name"), FUN = harm_mean.xtabs, x = plot_table)
    variable_harm_mean <- setNames(variable_harm_mean, setdiff(varnames_random, "line_name"))
    
    exp <- "line_name / (line_name + (line_name:environment / n_e) + (Residual / (n_e * n_r)))"
    # Calculate heritability
    heritability <- herit(object = object, exp = exp, n_e = variable_harm_mean$environment, n_r = rep_harm_mean)
    
  } else {
    exp <- "line_name / (line_name + (Residual / (n_r)))"
    # Calculate heritability
    heritability <- herit(object = object, exp = exp, n_r = rep_harm_mean)

  }
  
  # Return the heritability
  return(heritability)
  
}



## Prediction functions

## A function to model variance components
calc_variance <- function(data, random_effect = c("line_name", "environment", "line_name:environment")) {
  
  # Get the weights
  wts <- data$std_error^2
  # Set the control
  control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  
  # Write the formula
  form <- as.formula(str_c("value ~ ", str_c(str_c("(1|", random_effect), ")", collapse = " + ")))
  
  # Fit the model
  fit <- lmer(formula = form, data = data, control = control, weights = wts)
  
  # Extract the variance components
  VarCorr(fit) %>% 
    as.data.frame() %>% 
    select(var_comp = grp, variance = vcov)
}






## A generic prediction function that takes training and test data and returns
## PGVs and accuracy
gblup <- function(K, train, test, bootreps = NULL) {
  
  # Create a model.frame
  mf <- model.frame(value ~ line_name + env, weights = std_error, data = train)
  
  # Verify that the number of levels of 'line_name' is equal to the dimensions of K
  stopifnot(all(nlevels(mf$line_name) == dim(K)))
  
  
  
  # Vectors and matrices
  y <- model.response(mf)
  X <- model.matrix(~ -1 + env, mf)
  # Remove columns with all zero
  X <- X[,colSums(X) != 0, drop = FALSE]
  
  Z <- model.matrix(~ -1 + line_name, mf)
  
  fit <- mixed.solve(y = y, Z = Z, K = K, X = X)
  # Extract PGVs
  pgv <- fit$u %>% 
    data.frame(line_name = names(.), pred_value = ., row.names = NULL)
  
  # Combine the PGVs with the phenotypic observations and calculate accuracy
  comb <- left_join(test, pgv, by = "line_name")
  acc <- cor(comb$value, comb$pred_value)
  
  # Bootstrap if replicates are provided
  if (!is.null(bootreps)) {
    boot <- boot_cor(x = comb$value, y = comb$pred_value, boot.reps = bootreps)
  } else {
    boot <- NA
  }
    
  
  # Return a list
  list(accuracy = acc, pgv = comb, boot = as_data_frame(boot))
}


## Define a function that calculates the locale maximum of a environmental
## distance-formatted data.frame of prediction accuracy
local_maximum <- function(df) {
  
  
  
}


















