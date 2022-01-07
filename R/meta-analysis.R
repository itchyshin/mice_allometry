# meta-analysis

#devtools::install_github("itchyshin/orchard_plot", subdir = "orchaRd", force = TRUE, build_vignettes = TRUE)

library(orchaRd)
library(tidyverse)
library(here)
library(metafor)
library(brms)


# reading data
# TODO data from Laura had a typo - fm_diiff_slope_se (fixed)
dat <- read_csv(here("data/data_parameters2.csv"))
dat$obs <- 1:dim(dat)[1]


# we should adjust mean as well
# folded mean
folded_mu <-function(mean, variance){
  mu <- mean
  sigma <- sqrt(variance)
  fold_mu <- sigma*sqrt(2/pi)*exp((-mu^2)/(2*sigma^2)) + mu*(1 - 2*pnorm(-mu/sigma))
  fold_mu
} 

# folded variance
folded_v <-function(mean, variance){
  mu <- mean
  sigma <- sqrt(variance)
  fold_mu <- sigma*sqrt(2/pi)*exp((-mu^2)/(2*sigma^2)) + mu*(1 - 2*pnorm(-mu/sigma))
  fold_se <- sqrt(mu^2 + sigma^2 - fold_mu^2)
  # adding se to make bigger mean
  fold_v <-fold_se^2
  fold_v
} 


dat <- dat %>% mutate(abs_int = folded_mu(fm_diff_int, fm_diff_int_se^2), 
                      abs_slope = folded_mu(fm_diff_slope, fm_diff_slope_se^2),
                      abs_lnVR = folded_mu(lnVR, VlnVR),
                      V_abs_int = folded_v(fm_diff_int, fm_diff_int_se^2), 
                      V_abs_slope = folded_v(fm_diff_slope, fm_diff_slope_se^2),
                      V_abs_lnVR = folded_v(lnVR, VlnVR))


######################
# absolute effect size 
######################

# using transform-and-analyzie 

#############
# intercept
############

modelia <- rma.mv(yi = abs_int, 
                  V= V_abs_int, 
                  random = list(~1| Category, ~1| parameter_group, ~1|obs), 
                  data = dat)
summary(modelia)
robust(modelia, cluster  =  dat$parameter_group)

funnel(modelia)
i2_ml(modelia)

orchard_plot(modelia, mod = "Int", xlab = "Difference in standarised intercepts  (F-M)")


# meta-regression
model1a <- rma.mv(yi = abs_int, V= V_abs_int, mod = ~ Category - 1,
                  random = list(~1| parameter_group, ~1|obs), 
                  data = dat)
summary(model1a)
robust(model1a, cluster  =  dat$parameter_group)
r2_ml(model1a)

orchard_plot(model1a, mod = "Category", xlab = "Difference in standarised intercepts  (F-M)", angle = 45, cb = F)




####################
# slope difference
######################

modelsa <- rma.mv(yi = abs_slope, V= V_abs_slope, 
                  random = list(~1| Category, ~1| parameter_group, ~1|obs), 
                  data = dat)
summary(modelsa) # not sig this means sometimes male is high other times female has steaper slops
robust(modelsa, cluster  =  dat$parameter_group)

funnel(modelsa)
i2_ml(modelia)

orchard_plot(modelsa, mod = "Int", xlab = "Difference in standarised intercepts  (F-M)")

# meta-regression
model2a <- rma.mv(yi = abs_slope, V= V_abs_slope,
                  mod = ~ Category - 1,
                  random = list(~1| parameter_group, ~1|obs), 
                  data = dat)
summary(model2a)
robust(model2a, cluster  =  dat$parameter_group)
r2_ml(model2a)

orchard_plot(model2a, mod = "Category", xlab = "Difference in standarised intercepts  (F-M)", angle = 45, cb = F)

#################
# sd difference
#################

modelsda <- rma.mv(yi = abs_lnVR, V= V_abs_lnVR, 
                   random = list(~1| Category, ~1| parameter_group, ~1|obs), 
                   data = dat)
summary(modelsda)
robust(modelsda, cluster  =  dat$parameter_group)
funnel(modelsda)
i2_ml(modelsda)

orchard_plot(modelsda, mod = "Category", xlab = "Relative difference in SD (lnVR: F/M)", angle = 45, cb = F)


# meta-regression
model3a <- rma.mv(yi = abs_lnVR, V= V_abs_lnVR, mod = ~ Category - 1,
                  random = list(~1| parameter_group, ~1|obs), 
                  data = dat)
summary(model3a)
r2_ml(model3a)


orchard_plot(model3a, mod = "Category", xlab = "Relative difference in SD (lnVR: F/M)", angle = 45, cb = F)

######
# it seems like best to log it for we can

# TODO - multivariate meta-analysis
cor(log(dat$abs_int), log(dat$abs_slope))
cor(log(dat$abs_int), log(dat$abs_lnVR))
cor(log(dat$abs_slope), log(dat$abs_lnVR))


plot(log(dat$abs_int), log(dat$abs_slope))
plot(log(dat$abs_int), log(dat$abs_lnVR))
plot(log(dat$abs_slope), log(dat$abs_lnVR))

##########
#  try - brms
#########

# trivariate model

mod_lnsd <- bf(log(abs_lnVR) | se(sqrt(V_abs_int)/abs_lnVR)  ~ - 1 +  Category+ (1|q|parameter_group))
mod_lnslp <- bf(log(abs_slope) | se(sqrt(V_abs_slope)/abs_slope)  ~  - 1 +  Category + (1|q|parameter_group))
mod_lnint <- bf(log(abs_int) | se(sqrt(V_abs_int)/abs_int)  ~  - 1 +  Category + (1|q|parameter_group))

fit_3 <- brm(mod_lnsd + mod_lnslp + mod_lnint, 
                 data = dat, 
                 chains = 2, cores = 2, iter = 4000, warmup = 1000
)

summary(fit_3)

# saving the model



###########################
# non-absolute effect sizes
###########################
# intercept
modeli <- rma.mv(yi = fm_diff_int, V= fm_diff_int_se^2, random = list(~1| Category, ~1|obs), data = dat)
summary(modeli)
robust(modeli, cluster = dat$Category)

funnel(modeli, yaxis = "seinv")
i2_ml(modeli)
orchard_plot(modeli, mod = "Category", xlab = "Difference in standarised intercepts  (F-M)", angle = 45, cb = F)

# meta-regression
model1 <- rma.mv(yi = fm_diff_int, V= fm_diff_int_se^2, mod = ~ Category - 1,
                 random = list(~1| Category, ~1|obs), data = dat)
summary(model1)
r2_ml(model1)

orchard_plot(model1, mod = "Category", xlab = "Difference in standarised intercepts (F-M)", angle = 45, cb = F)

# slope difference
models <- rma.mv(yi = fm_diff_slope, V= fm_diff_slope_se^2, random = list(~1| Category, ~1|obs), data = dat)
summary(models) # not sig - this means sometimes male is high other times female has steeper slops
funnel(models,  yaxis = "seinv")        
i2_ml(models)

orchard_plot(models, mod = "Category", xlab = "Difference in standarised slopes  (F-M)", angle = 45, cb = F)
# meta-regression
model2 <- rma.mv(yi = fm_diff_slope, V= fm_diff_slope_se^2, mod = ~ Category - 1,
                 random = list(~1| Category, ~1|obs), data = dat)
summary(model2)
orchard_plot(model2, mod = "Category", xlab = "Difference in standarised slopes  (F-M)", angle = 45, cb = F)

# sd difference
modelsd <- rma.mv(yi = lnVR, V= VlnVR, random = list(~1| Category, ~1|obs), data = dat)
summary(modelsd)
funnel(modelsd,  yaxis = "seinv")
i2_ml(modelsd)

orchard_plot(modelsd, mod = "Category", xlab = "relative difference in SD (lnVR: F/M)", angle = 45, cb = F)
# meta-regression
model3 <- rma.mv(yi = lnVR, V= VlnVR, mod = ~ Category - 1,
                 random = list(~1| Category, ~1|obs), data = dat)
summary(model3)
orchard_plot(model3, mod = "Category", xlab = "relative difference in SD (lnVR: F/M)", angle = 45, cb = F)

