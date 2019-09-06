## S2MET Predictions
## 
## Modeling phenotypic data for analysis
## 
## Author: Jeff Neyhart
## Last updated: 29 May 2019
## 
## 


# The head directory
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

library(lme4qtl)
library(lmerTest)

## significance level
alpha <- 0.05


## Filter the BLUEs for the environments of interest
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  filter(environment %in% tp_vp_env)




## Stage-Two analysis

## Combine data
S2_MET_BLUEs_tomodel <- bind_rows(mutate(S2_MET_BLUEs_use, population = "all"), 
                                  mutate(filter(S2_MET_BLUEs_use, line_name %in% tp), population = "tp"),
                                  mutate(filter(S2_MET_BLUEs_use, line_name %in% vp), population = "vp"))

S2_MET_BLUEs_tomodel1 <- S2_MET_BLUEs %>%
  filter(environment %in% tp_vp_env)

S2_MET_BLUEs_tomodel1 <- bind_rows(mutate(S2_MET_BLUEs_tomodel1, population = "all"), 
                                   mutate(filter(S2_MET_BLUEs_tomodel1, line_name %in% tp), population = "tp"),
                                   mutate(filter(S2_MET_BLUEs_tomodel1, line_name %in% vp), population = "vp"))

S2_MET_BLUEs_tomodel <- S2_MET_BLUEs_tomodel1




### Refit models using the realized relationship matrix

lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")

stage_two_fits_GE_relmat <- S2_MET_BLUEs_tomodel %>%
  mutate_at(vars(location, year, line_name), as.factor) %>%
  group_by(trait, population) %>%
  do({
    
    df <- .
    
    ## Create interaction object
    df1 <- df %>%
      filter(line_name %in% rownames(K)) %>%
      # filter(environment %in% sample(unique(.$environment), 2)) %>%
      droplevels() %>%
      mutate(int = interaction(environment, line_name, sep = ":"),
             wts = std_error)
    
    
    # Table of lines by environments (i.e. plots)
    plot_table <- xtabs(formula = ~ line_name + environment, data = df)
    
    ## Harmonic means
    # Locations
    harm_env <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Reps
    harm_rep <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
      harm_mean()
    
    
    # Get the weights
    wts <- df1$std_error^2
    # W <- diag(wts)
    
    ## Relationship matrices
    E <- unique(df1$environment) %>%
      {`dimnames<-`(diag(x = length(.)), list(., .))}
    
    Kge <- kronecker(E, K, make.dimnames = TRUE)

    # # Relationship matrix
    # fit <- relmatLmer(value ~ 1 + environment + (1|line_name) + (1|int), relmat = list(line_name = K, int = Kge),
    #                   data = df1, control = lmer_control, weights = wts)
    # 
    # ## Test random effects
    # fit_ranova <- ranova(fit)
    
    
    
    ## Use sommer instead
    y <- df1$value
    X <- model.matrix(~ 1 + environment, df1)
    Zg <- model.matrix(~ -1 + line_name, df1)
    Zge <- model.matrix(~ -1 + int, df1)

    fit_somm <- sommer::mmer(Y = y, X = X, Z = list(line_name = list(Z = Zg, K = K), int = list(Z = Zge, K = Kge)),
                             W = W)
    
    ## Drop REs and test for significance
    fit_somm_red_noGE <- sommer::mmer(Y = y, X = X, Z = list(line_name = list(Z = Zg, K = K)), W = W)
    # LRT
    LR_GxE <- -2 * (fit_somm_red_noGE$LL - fit_somm$LL)
    p_value_GxE <- pchisq(q = LR_GxE, df = 1, lower.tail = FALSE) / 2
    
    ## Drop REs and test for significance
    fit_somm_red_noG <- sommer::mmer(Y = y, X = X, Z = list(int = list(Z = Zge, K = Kge)), W = W)
    # LRT
    LR_G <- -2 * (fit_somm_red_noG$LL - fit_somm$LL)
    p_value_G <- pchisq(q = LR_G, df = 1, lower.tail = FALSE) / 2
    
    lrt <- data_frame(
      term = c("line_name", "interaction"),
      LRT = c(LR_G, LR_GxE),
      df = 1,
      p.value = c(p_value_G, p_value_GxE)
    )
    
    ## Calculate heritability
    exp <- "line_name / (line_name + (int / n_e) + (Residual / (n_e * n_r)))"
    
    # ## Use bootstrapping to calculate a confidence interval
    # # Generate bootstrapping samples and calculate heritability using the bootMer function
    # h2_boot <- bootMer(x = fit, nsim = boot_reps, FUN = function(x) 
    #   herit(object = x, exp = exp, n_e = harm_env, n_r = harm_rep)$heritability)
    
    line_name <- fit_somm$var.comp$line_name[[1]]
    int <- fit_somm$var.comp$int[[1]]
    Residual <- fit_somm$var.comp$units[[1]]
    n_e <- harm_env
    n_r <- harm_rep
    
    h2 <- eval(parse(text = exp))
    
    
    # 
    # # Add the bootstrapping results with a confidence interval
    # h2$heritability <- tidy(h2_boot) %>% 
    #   cbind(., t(quantile(h2_boot$t, probs = c(alpha / 2, 1 - (alpha / 2))))) %>% 
    #   rename_at(vars(4, 5), ~c("lower", "upper"))
    
    # Return data_frame
    data_frame(fit = list(fit_somm), lrt = list(lrt), h2 = list(h2), n_e = harm_env, n_r = harm_rep) 
    
  }) %>% ungroup()




