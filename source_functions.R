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
gblup <- function(formula, random, K, train, test, fun = c("rrblup", "sommer"), fit.env = TRUE, bootreps = NULL, add.mu = FALSE) {
  
  if (missing(formula)) formula <- value ~ 1 + environment
  if (missing(random)) random <- ~ line_name
  
  ## separate random and fixed terms
  fixed_terms <- attr(terms(formula), "term.labels")
  random_terms <- attr(terms(random), "term.labels")

  ## Combine fixed and random to one formula
  formula1 <- as.formula(paste(paste(as.character(formula)[c(2,1,3)], collapse = " "), as.character(random)[-1], sep = " + "))
  
  # Create a model.frame
  mf <- model.frame(formula1, weights = std_error, data = train)
  
  
  ## If the number of random terms is > 1, K must be a list
  if (length(random_terms) == 1) {
    stopifnot(levels(mf[[random_terms]]) %in% colnames(K))
    
    K <- list(K)
    names(K) <- random_terms
    
  } else {
    if (!is.list(K)) stop("If the number of random terms is > 1, K must be a list of relationship matrices.")
    
    ## Test names
    name_test <- map2_lgl(.x = random_terms, .y = names(K1), ~.x == .y)
    if (!all(name_test)) stop("If K is a list, the names of the list must match the random terms")
    
  }
   
  
  
  fun <- match.arg(fun)
  
  
  
  # Vectors and matrices
  y <- model.response(mf)
  
  
  if (nlevels(mf[[fixed_terms]]) <= 1 | !fit.env) {
    X <- model.matrix(~ 1, droplevels(mf))
  
    } else {
    X <- model.matrix(formula, droplevels(mf))
    
  }
  
  ## Random effects
  Z_list <- list()
  
  for (term in random_terms) {
    Zi <- model.matrix(as.formula(paste0("~ -1 + ", term)), mf)
    colnames(Zi) <- colnames(K[[term]])
    Z_list[[term]] <- Zi
    
  }
  
  
  # Split on function
  if (fun == "rrblup") {
    # Can't run with > 1 random term
    stopifnot(length(random_terms) == 1)
    
    fit <- mixed.solve(y = y, Z = Z_list[[1]], K = K[[1]], X = X)
    
    # Extract PGVs
    pgv <- fit$u %>% 
      data.frame(line_name = names(.), pred_value = ., row.names = NULL, stringsAsFactors = FALSE)
    
    beta <- fit$beta[1]
    
    
  } else if (fun == "sommer") {
    
    # R <- solve(diag(mf$`(weights)`^2))
    # fit <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Z, K = K)), R = list(res = R), silent = TRUE)
    
    random_list <- map2(Z_list, K, ~list(Z = .x, K = .y))
    fit <- sommer::mmer(Y = y, X = X, Z = random_list, silent = TRUE)
  
    # Extract PGVs
    pgv <- fit$u.hat
    pgv <- data.frame(line_name = row.names(pgv), pred_value = pgv[,1], row.names = NULL, stringsAsFactors = FALSE)
    
    beta <- fit$beta.hat[1]
    
  }
  
  if (add.mu) pgv$pred_value <- pgv$pred_value + beta
  
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
      boot <- bootstrap(x = comb$value, y = comb$pred_value, fun = "cor", boot.reps = bootreps)
    } else {
      boot <- NA
    }
    
  }
  
  # Return a list
  list(accuracy = acc, pgv = comb, boot = as_data_frame(boot))
}



## Calculate fixed effects
geno_means <- function(formula = value ~ -1 + line_name + environment, data) {
  
  ## Get the terms in the model
  terms <- terms(formula)
  term_labels <- attr(terms, "term.labels")
  
  # Check the levels of these terms
  terms_levels <- sapply(data[term_labels], n_distinct)
  
  ## Drop terms that only have 1 level
  terms_drop <- which(terms_levels <= 1)
  
  if (length(terms_drop) > 0) {
    formula_use <- formula(drop.terms(terms, terms_drop, keep.response = TRUE))
  } else {
    formula_use <- formula
  }
  
  ## Fit
  # Contrasts list
  term_labels <- attr(terms(formula_use), "term.labels")
  contrasts <- setNames(object = replicate(length(term_labels), "contr.sum", simplify = FALSE), nm = term_labels)
  
  fit <- lm(formula = formula_use, data = data, contrasts = contrasts)
  
  # Extract the blues
  coef(fit) %>%
  data_frame(line_name = names(.), value = .) %>%
    filter(str_detect(line_name, "line_name")) %>%
    mutate(line_name = factor(str_remove_all(line_name, "line_name"), levels = c(tp_geno, vp_geno)),
           env = "mean", std_error = 0)
  
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




## Function to create training and testing sets
resample_prediction <- function(data, train.exp, test.exp) {
  
  # Add row numbers to the data
  data1 <- droplevels(mutate(data, rows = seq(nrow(data))))
  
  # Convert to expressions
  train.exp1 <- parse(text = train.exp)
  test.exp1 <- parse(text = test.exp)
  
  # Create the train and test resamples
  train_rows <- subset(data1, eval(train.exp1), rows, drop = T)
  test_rows <- subset(data1, eval(test.exp1), rows, drop = T)
  
  data_frame(train = list(resample(data = data1, idx = train_rows)),
             test = list(resample(data = data1, idx = test_rows)))
  
}



## Function to perform model-based clustering
env_mclust <- function(data, min_env = 2, test.env) {
  
  ## If test.env is present, remove them from the original dataset
  if (!missing(test.env)) {
    test_data <- data[test.env,,drop = FALSE]
    train_env <- setdiff(row.names(data), test.env)
    train_data <- data[train_env,,drop = FALSE]
    
    # Otherwise use all data
  } else {
    train_data <- data
  }
  
  # Cluster
  clusmod <- Mclust(data = train_data)
  # K
  k <- G <- clusmod$G
  
  # Do any clusters have less than the minimum number of environments?
  while (any(table(clusmod$classification) < min_env)) {
    # New G
    G <- G - 1
    # Recluster
    clusmod <- Mclust(data = train_data, verbose = FALSE, G = seq(G))
  }
  
  ## Get the classification of the environments
  classif <- data_frame(environment = row.names(clusmod$data), cluster = clusmod$classification)
  # Get the mean from each cluster
  clusmean <- t(clusmod$parameters$mean) %>% 
    `row.names<-`(., paste0("cluster", seq(nrow(.))))
  
  ## If test.env is present, use the data to determine the cluster assignment for each test environment
  ## based on the means of the assigned clusters
  if (!missing(test.env)) {
    
    # Calculate the distances for each environment
    test_dist <- as.matrix(dist(rbind(clusmean, test_data)))
    # Now figure out the minimum distance from each environment to each cluster
    test_assignment <- apply(X = test_dist[row.names(test_data),row.names(clusmean), drop = FALSE], MARGIN = 1, FUN = which.min)
    # Add this data to the classif
    classif1 <- classif %>% add_row(environment = names(test_assignment), cluster = test_assignment)
    
  } else {
    classif1 <- classif
    
  }
  
  # Return the classification
  return(classif1)
  
}
  


  
## Prediction functions for cross-validation
## Functions for each model
model1 <- function(train, test) {
  fit_m1 <- lmer(value ~ 1 + (1|line_name) + (1|environment), data = train)
  pred <- test %>% 
    mutate(pred_value = predict(object = fit_m1, newdata = ., allow.new.levels = TRUE)) %>%
    select(environment, line_name, value, pred_value)
  
  return(pred)
}

model2 <- function(train, test, Kg, Ke) {
  mf <- train %>% 
    mutate(LineName = line_name,
           line_name = factor(line_name, levels = colnames(Kg))) %>%
    model.frame(value ~ line_name + LineName + environment, data = .)
  
  # Matrices
  y <- model.response(mf)
  X <- model.matrix(~ 1, data = mf)
  Zg <- model.matrix(~ -1 + line_name, mf)
  colnames(Zg) <- colnames(Kg)
  ZE <- model.matrix(~ -1 + environment, mf)
  colnames(ZE) <- unique(mf$environment)
  
  
  
  # If Ke is missing, use diagonal
  if (missing(Ke)) Ke <- diag(ncol(ZE))
  
  # Fit
  fit <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = Kg), E = list(Z = ZE, K = Ke)), silent = TRUE)
  
  # Predictions
  pred_g <- fit$u.hat$g %>%
    as.data.frame() %>%
    rownames_to_column("line_name") %>%
    rename(g = T1)
  pred_E <- fit$u.hat$E %>%
    as.data.frame() %>%
    rownames_to_column("environment") %>%
    rename(E = T1)
  
  pred <- left_join(test, pred_g, by = "line_name") %>%
    left_join(., pred_E, by = "environment") %>%
    mutate(pred_value = rowSums(select(., g, E), na.rm = TRUE)) %>%
    select(environment, line_name, value, pred_value)
  
  
  # ZE <- model.matrix(~ 1 + environment, mf, contrasts.arg = list(environment = "contr.sum"))
  # colnames(ZE) <- c("mean", head(unique(mf$environment), -1))
  # 
  # fit_alt <- mixed.solve(y = y, Z = Zg, K = Kg, X = ZE)
  # 
  # 
  # pred_g <- fit_alt$u %>% 
  #   as.data.frame() %>%
  #   rownames_to_column("line_name") %>%
  #   rename_at(vars(-line_name), ~"g")
  # pred_E <- fit_alt$beta[-1] %>% 
  #   as.data.frame() %>% 
  #   rownames_to_column("environment") %>%
  #   rename_at(vars(-environment), ~"E") %>% 
  #   add_row(environment = tail(unique(mf$environment), 1), E = -sum(.$E))
  # 
  # pred_alt <- left_join(test, pred_g, by = "line_name") %>%
  #   left_join(., pred_E, by = "environment") %>%
  #   mutate(pred_value = g + E) %>%
  #   select(environment, line_name, value, pred_value)
  
  return(pred)
  
}



model3 <- function(train, test, Kg, Ke) {
  mf <- train %>% 
    mutate(LineName = line_name,
           line_name = factor(line_name, levels = colnames(Kg)),
           interaction = interaction(line_name, environment, sep = ":")) %>%
    model.frame(value ~ line_name + LineName + environment + interaction, data = .)
  
  # Matrices
  y <- model.response(mf)
  X <- model.matrix(~ 1, data = mf)
  Zg <- model.matrix(~ -1 + line_name, mf)
  colnames(Zg) <- colnames(Kg)
  ZG <- model.matrix(~ -1 + LineName, mf)
  ZE <- model.matrix(~ -1 + environment, mf)
  colnames(ZE) <- unique(mf$environment)
  ZgE <- model.matrix(~ -1 + interaction, mf)
  colnames(ZgE) <- levels(mf$interaction)
  
  # If Ke is missing, use diagonal
  if (missing(Ke)) Ke <- diag(ncol(ZE))
  
  # K_gE <- (Zg %*% Kg %*% t(Zg)) * (ZE %*% Ke %*% t(ZE))
  K_gE <- kronecker(X = Ke, Y = Kg, make.dimnames = TRUE)
  dimnames(K_gE) <- replicate(2, colnames(ZgE), simplify = FALSE)
  
  # Fit
  fit <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = Kg), E = list(Z = ZE, K = Ke), gE = list(Z = ZgE, K = K_gE)), silent = TRUE)
  # fit <- sommer::mmer(Y = y, X = X, Z = list(gE = list(Z = ZgE, K = K_gE)), silent = TRUE)
  
  
  # Predictions
  pred_g <- fit$u.hat$g %>%
    as.data.frame() %>%
    rownames_to_column("line_name") %>%
    rename(g = T1)
  pred_E <- fit$u.hat$E %>%
    as.data.frame() %>%
    rownames_to_column("environment") %>%
    rename(E = T1)
  pred_gE <- fit$u.hat$gE %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    separate(term, c("line_name", "environment"), sep = ":") %>%
    rename(gE = T1)
  
  pred <- left_join(test, pred_g, by = "line_name") %>%
    left_join(., pred_E, by = "environment") %>%
    left_join(., pred_gE, by = c("line_name", "environment")) %>%
    mutate(pred_value = rowSums(select(., g, E, gE), na.rm = TRUE)) %>%
    select(environment, line_name, value, pred_value)
    
  # pred <- left_join(test, pred_gE, by = c("line_name", "environment")) %>%
  #   mutate(pred_value = rowSums(select(., gE), na.rm = TRUE)) %>%
  #   select(environment, line_name, value, pred_value)

  
  return(pred)
  
}



## Model 4 has a single random effect (main effect of G)
model4 <- function(train, test, Kg) {
  
  mf <- train %>% 
    mutate(line_name = factor(line_name, levels = colnames(Kg)),
           environment = as.factor(environment),
           interaction = interaction(line_name, environment, sep = ":")) %>%
    model.frame(value ~ line_name  + environment + interaction, data = .)
  
  ## First fit model accounting for G and E (no GxE)
  fit <- lm(value ~ line_name + environment, data = mf)
  # Get the marginal means
  means <- effects::Effect("line_name", fit) %>%
    as.data.frame()
  
  mf1 <- means %>%
    select(line_name, value = fit) %>%
    mutate(line_name = factor(line_name, levels = colnames(Kg))) %>%
    model.frame(value ~ line_name, data = .)
  
  # Matrices
  y <- model.response(mf1)
  X <- model.matrix(~ 1, data = mf1)
  Zg <- model.matrix(~ -1 + line_name, mf1)
  colnames(Zg) <- colnames(Kg)


  ## Prediction
  fit_alt <- mixed.solve(y = y, Z = Zg, K = Kg, X = X)

  # Extract the predictions
  pred_g <- fit_alt$u %>%
    as.data.frame() %>%
    rownames_to_column("line_name") %>%
    rename_at(vars(-line_name), ~"g")

  # Add test data
  pred_alt <- left_join(test, pred_g, by = c("line_name")) %>%
    mutate(pred_value = g) %>%
    select(environment, line_name, value, pred_value)
  
  
  return(pred_alt)
  
}


## Model 5 has a single random effect (GxE)
model5 <- function(train, test, Kg, Ke) {
  
  mf <- train %>% 
    mutate(line_name = factor(line_name, levels = colnames(Kg)),
           environment = as.factor(environment),
           interaction = interaction(line_name, environment, sep = ":")) %>%
    model.frame(value ~ line_name  + environment + interaction, data = .)
  
  # Matrices
  y <- model.response(mf)
  X <- model.matrix(~ 1, data = mf)
  Zg <- model.matrix(~ -1 + line_name, mf)
  colnames(Zg) <- colnames(Kg)
  ZE <- model.matrix(~ -1 + environment, mf)
  colnames(ZE) <- levels(mf$environment)
  # Interaction
  ZgE <- model.matrix(~ -1 + line_name:environment, mf)
  colnames(ZgE) <- levels(mf$interaction)
  
  # GxE covariance matrix
  # If Ke is missing, use diagonal
  if (missing(Ke)) Ke <- diag(ncol(ZE))
  
  # K_gE <- (Zg %*% Kg %*% t(Zg)) * (ZE %*% Ke %*% t(ZE))
  KgE <- kronecker(X = Ke, Y = Kg, make.dimnames = TRUE)
  dimnames(KgE) <- replicate(2, colnames(ZgE), simplify = FALSE)
  
  
  ## Prediction
  fit_alt <- mixed.solve(y = y, Z = ZgE, K = KgE, X = X)
  
  # Extract the predictions
  pred_gE <- fit_alt$u %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    separate(term, c("line_name", "environment"), sep = ":") %>%
    rename_at(vars(-line_name, -environment), ~"gE")
  
  # Add test data
  pred_alt <- left_join(test, pred_gE, by = c("line_name", "environment")) %>%
    mutate(pred_value = gE) %>%
    select(environment, line_name, value, pred_value)
  
  return(pred_alt)
  
}





  
