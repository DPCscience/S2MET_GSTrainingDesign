## S2MET source for MSI
## 
## A script that automatically loads the data relevant for the S2MET project


# Load packages
packages <- c("dplyr", "tidyr", "stringr", "readxl", "readr", "parallel",
              "rrBLUP", "purrr")
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-unknown-linux-gnu-library/3.2/"
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"

geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Genos"
pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")



# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))

# Load an entry file
entry_list <- read_excel(file.path(data_dir, "project_entries.xlsx"))

# Load the trial metadata
trial_info <- read_csv(file.path(data_dir, "trial_metadata.csv"))

# Vector of relevant traits
traits <- c("GrainYield", "HeadingDate", "PlantHeight")


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
  filter(Class != "Check") %>% 
  pull(Line)

# Extract the tp and vp from the G matrix
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno, vp_geno),]

# Calculate the K matrix
K <- A.mat(X = s2_imputed_mat_use, min.MAF = 0, max.missing = 1)

# Create a replacement vector for each trait that adds a space between words and
# also included the unit
trait_replace <- unique(S2_MET_BLUEs$trait) %>%
  setNames(object = c("Grain Yield", "Heading Date", "Plant Height"), nm = .)
trait_replace_unit <- unique(S2_MET_BLUEs$trait) %>%
  setNames(object = c("Grain Yield\n(kg ha^-1)", "Heading Date\n(days)", "Plant Height\n(cm)"), nm = .)

# Pull out the trials with irrigation
irrig_env <- subset(trial_info, irrigated == "yes", environment)$environment
# Filter the S2 MET BLUEs for non-irrigated trials
S2_MET_BLUEs <- S2_MET_BLUEs %>% 
  filter(!environment %in% irrig_env)



# Find environments in which just data on both the TP and VP is available
tp_vp_env <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) > 1, sum(line_name %in% vp_geno) > 1) %>% 
  distinct(environment) %>% 
  pull()

# Find environments in which just data on the TP is available
tp_only_env <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) > 1, sum(line_name %in% vp_geno) == 0) %>% 
  distinct(environment) %>% 
  pull()

# Find environments in which just data on the VP is available
vp_only_env <- S2_MET_BLUEs %>% 
  group_by(environment) %>%
  filter(sum(line_name %in% tp_geno) == 0, sum(line_name %in% vp_geno) > 1) %>% 
  distinct(environment) %>% 
  pull()




##### FUNCTIONS #####

## Below are functions that were written specifically for this project


## Other/utility functions
# A function to assign cores to a data.frame
assign_cores <- function(df, n_core) {
  df$core <- sort(rep(seq(n_core), length.out = nrow(df)))
  return(df)
}


## Heritability functions
## Write a function to subset a distance matrix using a set of environments
subset_env <- function(dist, envs) {
  as.matrix(dist) %>%
    .[envs, envs] %>%
    as.dist()
}

## Define a function to find the harmonic mean of a variable given an xtabs object
harm_mean.xtabs <- function(x, variable) {
  # Get the dimension names of the xtabs
  tab_names <- names(dimnames(x))
  # Position of the variable in the array
  variable_pos <- match(c("line_name", variable), tab_names)
  
  # Find the number of "variable" in which the jth individual was observed
  variable_sum <- apply(X = x, MARGIN = variable_pos, FUN = sum)
  variable_sum_adj <- ifelse(variable_sum > 1, 1, variable_sum)
  
  # If all tab names are called, calculate the harmonic mean of the reps
  if (length(variable_pos) == length(tab_names)) {
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
  rep_harm_mean <- harm_mean.xtabs(x = plot_table, variable = names(variable_harm_mean))
  
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


## A function to make predictions within a cluster
# assess_prediction <- function(data)

























