# some test code for Laura
# install.packages("devtools")
# install.packages("tidyverse")
# install.packages("metafor")
# install.packages("patchwork")
# install.packages("R.rsp")

#devtools::install_github("itchyshin/orchard_plot", subdir = "orchaRd", force = TRUE, build_vignettes = TRUE)

library(orchaRd)
library(tidyverse)
library(here)
library(metafor)


# reading data
# TODO data from Laura had a typo - fm_diiff_slope_se (fixed)
dat <- read_csv(here("Laura/AllScenarios_dat2.csv"))
dat$obs <- 1:dim(dat)[1]

# getting lnVR to compare SDs
dat$lnVR <- log(dat$f_sd/dat$m_sd) + 1/(2*(dat$f_n-1)) - 1/(2*(dat$m_n-1))
dat$VlnVR <-  1/(2*(dat$f_n-1)) + 1/(2*(dat$m_n-1))

# getting absolute values for contrasts

# dat$abs_int <- abs(dat$fm_diff_int)
# dat$abs_slope <- abs(dat$fm_diff_slope)
# dat$abs_lnVR <- abs(dat$lnVR)

# I should add absolute versions of them too
# function to get variance for absolute values
# var_abs <- function(x, var) {
#   abs_var <- x^2 + var - ((sqrt(2/pi)*sqrt(var)*exp((-1*x^2)/(2*var))) -
#                             x*(1-2*(pnorm(x/sqrt(var)))))^2
#   return(abs_var)
# }

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

######################
# absolute effect size 
######################
# intercept
modelia <- rma.mv(yi = abs_int, V= V_abs_int, random = list(~1| Category, ~1|obs), data = dat)
summary(modelia)
funnel(modelia)
i2_ml(modelia)

orchard_plot(modelia, mod = "Int", xlab = "Difference in standardised intercepts (z values)")


# meta-regression
model1a <- rma.mv(yi = abs_int, V= V_abs_int, mod = ~ Category - 1,
                  random = ~1|obs, data = dat)
summary(model1a)
r2_ml(model1a)

orchard_plot(model1a, mod = "Category", xlab = "Difference in standardised intercepts (z values)", angle = 45, cb = F)

# slope difference
modelsa <- rma.mv(yi = abs_slope, V= V_abs_slope, random = list(~1| Category, ~1|obs), data = dat)
summary(modelsa) # not sig this means sometimes male is high other times female has steaper slops
funnel(modelsa)
i2_ml(modelia)

orchard_plot(modelsa, mod = "Int", xlab = "Difference in standardised slopes (sem-partial correlations)")

# meta-regression
model2a <- rma.mv(yi = abs_slope, V= V_abs_slope, mod = ~ Category - 1,
                  random =  ~1|obs, data = dat)
summary(model2a)
orchard_plot(model2a, mod = "Category", xlab = "Difference in sem-partial correlations(sem-partial correlations)", angle = 45, cb = F)

# sd difference
modelsda <- rma.mv(yi = abs_lnVR, V= V_abs_lnVR, random = list(~1| Category, ~1|obs), data = dat)
summary(modelsda)
funnel(modelsda)

# meta-regression
model3a <- rma.mv(yi = abs_lnVR, V= V_abs_lnVR, mod = ~ Category - 1,
                  random = ~1|obs, data = dat)
summary(model3a)
orchard_plot(model3a, mod = "Category", xlab = "log variability ratio (ratio of resiudals SD)", angle = 45, cb = F)


######

# TODO - multivariate meta-analysis
cor(dat$abs_int, dat$abs_slope)
cor(dat$abs_int, dat$abs_lnVR)
cor(dat$abs_slope, dat$abs_lnVR)


###############################
# oiginal anlayiss - do not run
###############################

# dividing data into 4!

Adat <- dat[dat$Scen == "A", ] # 16 x 20
Bdat <- dat[dat$Scen == "B", ] # 107 x 20
Cdat <- dat[dat$Scen == "C", ] # 72 x 20
Ddat <- dat[dat$Scen == "D", ] # 102 x 20

###########
# Dataset A
#############

# intercept
modelAi <- rma.mv(yi = fm_diff_int, V= fm_diff_int_se^2, random = list(~1| Category, ~1|obs), data = Adat)
summary(modelAi)
funnel(modelAi)

# meta-regression
modelA1 <- rma.mv(yi = fm_diff_int, V= fm_diff_int_se^2, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = Adat)
summary(modelA1)
orchard_plot(modelA1, mod = "Category", xlab = "Difference in z values", angle = 45)

# slope difference
modelAs <- rma.mv(yi = fm_diff_slope, V= fm_diff_slope_se^2, random = list(~1| Category, ~1|obs), data = Adat)
summary(modelAs) # not sig this means sometimes male is high other times female has steaper slops
funnel(modelAs)        
 
# meta-regression
modelA2 <- rma.mv(yi = fm_diff_slope, V= fm_diff_slope_se^2, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = Adat)
summary(modelA2)
orchard_plot(modelA2, mod = "Category", xlab = "Difference in sem-partial correlations", angle = 45)

# sd difference
modelAsd <- rma.mv(yi = lnVR, V= VlnVR, random = list(~1| Category, ~1|obs), data = Adat)
summary(modelAsd)
funnel(modelAsd)

# meta-regression
modelA3 <- rma.mv(yi = lnVR, V= VlnVR, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = Adat)
summary(modelA3)
orchard_plot(modelA3, mod = "Category", xlab = "log variability ratio", angle = 45)

############
# Dataset B
############
# intercept
modelBi <- rma.mv(yi = fm_diff_int, V= fm_diff_int_se^2, random = list(~1| Category, ~1|obs), data = Bdat)
summary(modelBi)
funnel(modelBi)

# meta-regression
modelB1 <- rma.mv(yi = fm_diff_int, V= fm_diff_int_se^2, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = Bdat)
summary(modelB1)
orchard_plot(modelB1, mod = "Category", xlab = "Difference in z values", angle = 45)

# slope difference
modelBs <- rma.mv(yi = fm_diff_slope, V= fm_diff_slope_se^2, random = list(~1| Category, ~1|obs), data = Bdat)
summary(modelBs) # not sig this means sometimes male is high other times female has steaper slops
funnel(modelBs)        

# meta-regression
modelB2 <- rma.mv(yi = fm_diff_slope, V= fm_diff_slope_se^2, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = Bdat)
summary(modelB2)
orchard_plot(modelB2, mod = "Category", xlab = "Difference in sem-partial correlations", angle = 45)

# sd difference
modelBsd <- rma.mv(yi = lnVR, V= VlnVR, random = list(~1| Category, ~1|obs), data = Bdat)
summary(modelBsd)
funnel(modelBsd)

# meta-regression
modelB3 <- rma.mv(yi = lnVR, V= VlnVR, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = Bdat)
summary(modelB3)
orchard_plot(modelB3, mod = "Category", xlab = "log variability ratio", angle = 45)

###########
# Dataset C
#############
# intercept
modelCi <- rma.mv(yi = fm_diff_int, V= fm_diff_int_se^2, random = list(~1| Category, ~1|obs), data = Cdat)
summary(modelCi)
funnel(modelCi)

# meta-regression
modelC1 <- rma.mv(yi = fm_diff_int, V= fm_diff_int_se^2, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = Cdat)
summary(modelC1)
orchard_plot(modelC1, mod = "Category", xlab = "Difference in z values", angle = 45)

# slope difference
modelCs <- rma.mv(yi = fm_diff_slope, V= fm_diff_slope_se^2, random = list(~1| Category, ~1|obs), data = Cdat)
summary(modelCs) # not sig this means sometimes male is high other times female has steaper slops
funnel(modelCs)        

# meta-regression
modelC2 <- rma.mv(yi = fm_diff_slope, V= fm_diff_slope_se^2, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = Cdat)
summary(modelC2)
orchard_plot(modelC2, mod = "Category", xlab = "Difference in sem-partial correlations", angle = 45)

# sd difference
modelCsd <- rma.mv(yi = lnVR, V= VlnVR, random = list(~1| Category, ~1|obs), data = Cdat)
summary(modelCsd)
funnel(modelCsd)

# meta-regression
modelC3 <- rma.mv(yi = lnVR, V= VlnVR, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = Cdat)
summary(modelC3)
orchard_plot(modelC3, mod = "Category", xlab = "log variability ratio", angle = 45)

########################
#######################

# all the datasets

# intercept
modeli <- rma.mv(yi = fm_diff_int, V= fm_diff_int_se^2, random = list(~1| Category, ~1|obs), data = dat)
summary(modeli)
funnel(modeli, yaxis = "seinv")
orchard_plot(model1i, mod = "Category", xlab = "Difference in z values", angle = 45, cb = F)

# meta-regression
model1 <- rma.mv(yi = fm_diff_int, V= fm_diff_int_se^2, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = dat)
summary(model1)
orchard_plot(model1, mod = "Category", xlab = "Difference in z values", angle = 45, cb = F)

# slope difference
models <- rma.mv(yi = fm_diff_slope, V= fm_diff_slope_se^2, random = list(~1| Category, ~1|obs), data = dat)
summary(models) # not sig this means sometimes male is high other times female has steaper slops
funnel(models)        

# meta-regression
model2 <- rma.mv(yi = fm_diff_slope, V= fm_diff_slope_se^2, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = dat)
summary(model2)
orchard_plot(model2, mod = "Category", xlab = "Difference in sem-partial correlations", angle = 45, cb = F)

# sd difference
modelsd <- rma.mv(yi = lnVR, V= VlnVR, random = list(~1| Category, ~1|obs), data = dat)
summary(modelsd)
funnel(modelsd)

# meta-regression
model3 <- rma.mv(yi = lnVR, V= VlnVR, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = dat)
summary(model3)
orchard_plot(model3, mod = "Category", xlab = "log variability ratio", angle = 45, cb = F)

######################
# absolute effect size 
######################
# intercept
modelia <- rma.mv(yi = abs_int, V= V_abs_int, random = list(~1| Category, ~1|obs), data = dat)
summary(modelia)
funnel(modelia)
orchard_plot(modelia, mod = "Int", xlab = "Difference in z values")

# meta-regression
model1a <- rma.mv(yi = abs_int, V= V_abs_int, mod = ~ Category - 1,
                 random = ~1|obs, data = dat)
summary(model1a)
orchard_plot(model1a, mod = "Category", xlab = "Difference in standardised intercepts (z values)", angle = 45, cb = F)

# slope difference
modelsa <- rma.mv(yi = abs_slope, V= V_abs_slope, random = list(~1| Category, ~1|obs), data = dat)
summary(modelsa) # not sig this means sometimes male is high other times female has steaper slops
funnel(modelsa)        
orchard_plot(modelsa, mod = "Int", xlab = "Difference in standardised slopes (sem-partial correlations)")

# meta-regression
model2a <- rma.mv(yi = abs_slope, V= V_abs_slope, mod = ~ Category - 1,
                 random =  ~1|obs, data = dat)
summary(model2a)
orchard_plot(model2a, mod = "Category", xlab = "Difference in sem-partial correlations(sem-partial correlations)", angle = 45, cb = F)

# sd difference
modelsda <- rma.mv(yi = abs_lnVR, V= V_abs_lnVR, random = list(~1| Category, ~1|obs), data = dat)
summary(modelsda)
funnel(modelsda)

# meta-regression
model3a <- rma.mv(yi = abs_lnVR, V= V_abs_lnVR, mod = ~ Category - 1,
                 random = ~1|obs, data = dat)
summary(model3a)
orchard_plot(model3a, mod = "Category", xlab = "log variability ratio (ratio of resiudals SD)", angle = 45, cb = F)

#######################
# Dataset A

# intercept
modelAia <- rma.mv(yi = abs_int, V= V_abs_int, random = list(~1| Category, ~1|obs), data = Adat)
summary(modelAia)
funnel(modelAia)
caterpillars(modelAia, xlab = "Difference in z values")
orchard_plot(modelAia, mod = "int", xlab = "Difference in z values")

# meta-regression
modelA1a <- rma.mv(yi = abs_int, V= V_abs_int, mod = ~ Category - 1,
                  random = ~1|obs, data = Adat)
summary(modelA1a)
orchard_plot(modelA1a, mod = "Category", xlab = "Difference in z values", angle = 45)

# slope difference
modelAsa <- rma.mv(yi = abs_slope, V= V_abs_slope, random = list(~1| Category, ~1|obs), data = Adat)
summary(modelAsa) 
funnel(modelAsa)  
caterpillars(modelAsa, xlab = "Difference in sem-partial correlations")

# meta-regression
modelA2a <- rma.mv(yi = abs_slope, V= V_abs_slope, mod = ~ Category - 1,
                  random = ~1|obs, data = Adat)
summary(modelA2a)
orchard_plot(modelA2a, mod = "Category", xlab = "Difference in sem-partial correlations", angle = 45)

# sd difference
modelAsda <- rma.mv(yi = abs_lnVR, V= V_abs_lnVR, random = list(~1| Category, ~1|obs), data = Adat)
summary(modelAsda)
funnel(modelAsda)


# meta-regression
modelA3a <- rma.mv(yi = abs_lnVR, V= V_abs_lnVR, mod = ~ Category - 1,
                  random = list(~1| Category, ~1|obs), data = Adat)
summary(modelA3a)
orchard_plot(modelA3a, mod = "Category", xlab = "log variability ratio", angle = 45)

