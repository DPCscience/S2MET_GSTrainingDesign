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



## Phenotypic distance method
dist_env <- function(x, gen.col = "gen", env.col = "env", pheno.col = "yield") {
  
  # Verify that x is a data.frame
  stopifnot(is.data.frame(x))
  
  # Convert to data.frame
  x <- as.data.frame(x)
  
  # Verify columnnames exist
  if (!all(c(gen.col, env.col, pheno.col) %in% colnames(x)))
    stop("The values of gen.col or env.col or pheno.col are not columns in the x data.frame.")
  
  # Pull out environment names
  env_names <- unique(as.character(x[[env.col]]))
  # Number of environments
  n_env <- length(env_names)
  
  ## Make a table of the observations in the dataset
  x_table <- table(x[c(env.col, gen.col)])
  # Are any greater than 1?
  if (any(x_table > 1)) stop ("Only one observation of each genotype in each environment should be included in the dataset.")
  
  # Order on environment, then genotype
  x1 <- x[order(x[[env.col]], x[[gen.col]]), , drop = FALSE]
  
  # Iterate over pairs
  D_ij <- combn(x = env_names, m = 2, FUN = function(env_pairs) {
    
    ## Subset the data for this environment pair
    env_index <- x1[[env.col]] %in% env_pairs
    x_sub <- x1[env_index, , drop = FALSE]
    
    # Find the common genotypes
    geno_count <- table(x_sub[[gen.col]])
    x_geno <- names(geno_count[geno_count == 2])
    
    # If no genotypes are in common, return NA
    if (length(x_geno) == 0) return(NA)
    
    # Subset the data again
    x_sub1 <- x_sub[x_sub[[gen.col]] %in% x_geno, , drop = FALSE]
    
    # Scale the phenotypic values by the mean and sd
    pheno_scale <- tapply(X = x_sub1[[pheno.col]], INDEX = x_sub1[[env.col]],
                          function(x) as.numeric(scale(x)), simplify = FALSE)
    
    # Find the squared differences between genotypes
    pheno_scale_diff <- tapply(X = unlist(pheno_scale, use.names = FALSE),
                               INDEX = x_sub1[[gen.col]],
                               function(x) diff(x)^2, simplify = FALSE)
    
    # Calculate the mean among these squared differences and return
    mean(unlist(pheno_scale_diff))
    
  })
  
  # Empty matrix
  D_mat <- matrix(0, nrow = n_env, ncol = n_env, dimnames = list(env_names, env_names))
  
  # Add data
  D_mat[lower.tri(D_mat)] <- D_ij
  
  # Output as distance matrix
  as.dist(D_mat)
  
} # Close the fuction



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



# Function that generates regression values for a population against an environmental index
# (i.e. stability), then performs an F test to test differences among the regression coefficients
fwr <- function(formula, data, gen = "line_name", env = "environment") {
  
  out <- data %>%
    split(.[[gen]]) %>%
    map(~{
      fit <- lm(formula, data = .)
      data_frame(g = coef(fit)[1], b = coef(fit)[2], d = sum(resid(fit)^2))
    }) %>%
    map2_df(.y = names(.), ~mutate(.x, line_name = .y)) %>%
    select(line_name, names(.))

  
  # Calculate degrees of freedom
  df_b <- n_distinct(data[[gen]]) - 1
  df_d <- df_b * (n_distinct(data[[env]]) - 2)
  
  # Environmental index sums of squares
  index <- distinct_(data, env, as.character(formula)[3])[[as.character(formula)[3]]]
  e_ss <- sum(index^2)
  
  ## Significance test
  b_ss <- sum(scale(out$b, scale = FALSE)^2) * e_ss
  d_ss <- sum(out$d)
  
  # Mean squares
  b_ms <- b_ss / df_b
  d_ms <- d_ss / df_d
  
  # F stat
  Fstat <- b_ms / d_ms
  pvalue <- pf(q = Fstat, df1 = df_b, df2 = df_d, lower.tail = FALSE)
  
  ## Output the results
  data_frame(
    regression = list(out),
    Fstat = Fstat,
    pvalue = pvalue
  )
  
  
}





### Function that computes the significance of principal components via permutation.
## This function will first run a principal component analysis using data,
## then it will permute the order of the data and for each permutation re-run
## the PCA. PCs will be retained if their proportion of variance is significant when
## compared to the distribution of variance proportion under permutation
## Note that the test is one-sided, since the variance cannot be negative.
prcomp_perm <- function(x, scale = TRUE, permutations = 1000, level = 0.95) {
  
  # Scale the data, if called
  x1 <- if (scale) scale(x) else x
  
  # Fit the original PCA
  pc_fit <- prcomp(x = x1, center = FALSE, scale. = FALSE)
  # Summarize
  pc_fit_summ <- summary(pc_fit)
  
  ## Permute the data
  data_permute <- replicate(n = permutations, t(apply(X = x1, MARGIN = 1, FUN = sample)), simplify = FALSE)
  
  # Run the PCA
  pc_perm_fit_list <- map(data_permute, ~prcomp(x = ., center = FALSE, scale. = FALSE)) %>%
    map(summary) %>%
    # Pull out the variance explained
    map(~.$importance[2,])
  
  pc_perm_var <- do.call("rbind", pc_perm_fit_list)
  
  ## Compare to the original variance to get a p-value
  pvalues <- rowMeans(apply(X = pc_perm_var, MARGIN = 1, FUN = function(x) x >= pc_fit_summ$importance[2,]))
  
  # Get a confidence interval
  var_ci <- apply(X = pc_perm_var, MARGIN = 2, FUN = quantile, probs = c((1 - level) / 2, 1 - ((1 - level) / 2)))
  
  # Return a list
  list(
    tidy = tidy(pc_fit),
    summary = data.frame(PC = colnames(pc_fit$x), var_exp = pc_fit_summ$importance[2,], as.data.frame(t(var_ci)), 
                         p_value = pvalues, row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)
  )
  
}









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
    data.frame(line_name = names(.), pred_value = ., row.names = NULL, stringsAsFactors = FALSE)
  
  # If test is missing, just return the predictions
  if (missing(test)) {
    
    comb <- pgv
    acc <- boot <- NA
    
  } else {
  
    # Combine the PGVs with the phenotypic observations and calculate accuracy
    comb <- left_join(test, pgv, by = "line_name")
    acc <- cor(comb$value, comb$pred_value)
    
    # Bootstrap if replicates are provided
    if (!is.null(bootreps)) {
      boot <- boot_cor(x = comb$value, y = comb$pred_value, boot.reps = bootreps)
    } else {
      boot <- NA
    }
    
  }
  
  # Return a list
  list(accuracy = acc, pgv = comb, boot = as_data_frame(boot))
}



## Function to calculate the mean in sliding window
window_mean <- function(x, y, window = 8) {
  
  # Size of window on either side
  side_window <- ceiling(window / 2)
  # Get the min and max
  min_x <- min(x)
  max_x <- max(x)
  
  # Create a list of indices
  x_use <- map(x, ~seq(. - side_window, . + side_window) %>% .[between(., min_x, max_x)]) %>%
    map(~match(x = ., table = x))
  
  # Calculate the mean y within that x
  map_dbl(x_use, ~mean(y[.], na.rm = TRUE))
  
}




## Return the quantiles of a random sample at level alpha
quantile1 <- function(x, alpha = 0.05) {
  qs <- quantile(x, probs = c(alpha / 2, 1 - (alpha / 2)))
  `names<-`(qs, c("lower", "upper"))
}











