## Testing
## 
## 

library(agridat)

kang <- kang.maize

# Covariates
# Weather covariates for each environment.
covs <- data.frame(env=c("AL85","AL86","AL87", "BR85","BR86","BR87",
                         "BC85","BC86","BC87", "SJ85","SJ86","SJ87"),
                   maxt=c(30.7,30.2,29.7,31.5,29.4,28.5, 31.9, 30.4,31.7, 32,29.6,28.9),
                   mint=c(18.7,19.3,18.5, 19.7,18,17.2, 19.1,20.4,20.3, 20.4,19.1,17.5),
                   rain=c(.2,.34,.22, .28,.36,.61, .2,.43,.2, .36,.41,.22),
                   humid=c(82.8,91.1,85.4, 88.1,90.9,88.6, 95.4,90.4,86.7, 95.6,89.5,85))

covs1 <- mutate_at(covs, vars(-env), ~scale(., scale = FALSE) %>% as.numeric())

# Add covariates to the dataset and round the yield values
kang1 <- left_join(kang, covs1) %>%
  mutate(yield = round(yield))


## Model
fit <- lm(yield ~ env + gen, kang1, contrasts = list(env = "contr.sum", gen = "contr.sum"))

# Find environmental means
env_coef <- coef(fit)[2:n_distinct(kang1$env)]
env_eff <- data.frame(env = levels(kang1$env), effect = c(env_coef, -sum(env_coef)))


# Add and regress each genotype
kang2 <- left_join(kang1, env_eff)






