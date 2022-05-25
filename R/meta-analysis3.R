# meta-analysis

# info about colour blind cour
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

#devtools::install_github("itchyshin/orchard_plot", subdir = "orchaRd", force = TRUE, build_vignettes = TRUE)

library(orchaRd)
library(tidyverse)
library(here)
library(metafor)
library(brms)
library(patchwork)

#TODO - we should really renames or provide meta-data on what all column names mean (at the moment all very confusing)

# loading functions
source(here("R/function.R"), chdir = TRUE)

# reading data
# TODO data from Laura had a typo - fm_diiff_slope_se (fixed)
dat <- read_csv(here("data/data_parameters6.csv"))
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

#####################
#
#####################

test <- rma.mv(yi = Zr, mod = ~ Category - 1,
               V= VZr, 
               random = list(~1| parameter_group, ~1|obs), 
               data = dat)
t2 <- orchard_plot2(test, mod = "Category", xlab = "Fit (correlation between observations and predicted values)", angle = 45,  point.size = point.size, transfm =  "tanh") + 
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  xlim(c(-0.5, 1.5))


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


cbpl <- c("#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
          "#CC79A7", "#56B4E9", "#AA4499", "#DDCC77")

point.size = 1.5

p1 <- orchard_plot2(modelia, mod = "Int", xlab = "Absolute difference in standardized intercepts  (F-M)", angle = 45, point.size = point.size) +
  scale_y_discrete(labels = "Overall") +
  scale_fill_manual(values = "#999999") +
  scale_colour_manual(values = "#999999") +
  xlim(c(-0.5, 1.5))

# meta-regression
model1a <- rma.mv(yi = abs_int, V= V_abs_int, mod = ~ Category - 1,
                  random = list(~1| parameter_group, ~1|obs), 
                  data = dat)
summary(model1a)
robust(model1a, cluster  =  dat$parameter_group)
r2_ml(model1a)

p2 <- orchard_plot2(model1a, mod = "Category", xlab = "Absolute difference in standardized intercepts  (F-M)", angle = 45,  point.size = point.size) + 
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  xlim(c(-0.5, 1.5))

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

p3 <- orchard_plot2(modelsa, mod = "Int", xlab = "Absolute difference in standardized slopes (F-M)", angle = 45,  point.size = point.size) +
  scale_y_discrete(labels = "") +
  scale_fill_manual(values = "#999999") +
  scale_colour_manual(values = "#999999") +
  xlim(c(-1.5, 10))

# meta-regression
model2a <- rma.mv(yi = abs_slope, V= V_abs_slope,
                  mod = ~ Category - 1,
                  random = list(~1| parameter_group, ~1|obs), 
                  data = dat)
summary(model2a)
robust(model2a, cluster  =  dat$parameter_group)
r2_ml(model2a)

p4 <- orchard_plot2(model2a, mod = "Category", xlab = "Absolute difference in standardized slopes (F-M)", angle = 45, cb = F,  point.size = point.size) + 
  scale_y_discrete(labels = rep("", 9)) +
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  xlim(c(-1.5, 10))

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

p5 <- orchard_plot2(modelsda, mod = "Category", xlab = "Absolute relative difference in SD (lnVR: F/M)", angle = 45,  point.size = point.size) +
  scale_y_discrete(labels = "") +
  scale_fill_manual(values = "#999999") +
  scale_colour_manual(values = "#999999") +
  xlim(c(-0.2, 1.9))

# meta-regression
model3a <- rma.mv(yi = abs_lnVR, V= V_abs_lnVR, mod = ~ Category - 1,
                  random = list(~1| parameter_group, ~1|obs), 
                  data = dat)
summary(model3a)
r2_ml(model3a)

p6 <- orchard_plot2(model3a, mod = "Category", xlab = "Absolute relative difference in SD (lnVR: F/M)", angle = 45, cb = F,  point.size = point.size) + 
  scale_y_discrete(labels = rep("", 9)) +
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  xlim(c(-0.2, 1.9))

###############
# Fig 3
###############

(p1 + p3 + p5) / (p2 + p4 + p6)  + plot_layout(heights = c(1, 4)) + plot_annotation(tag_levels = 'A')

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

# tri-variate model
# TODO rerun this model!!

mod_lnsd <- bf(log(abs_lnVR) | se(sqrt(V_abs_lnVR)/abs_lnVR)  ~ - 1 +  Category+ (1|q|parameter_group))
mod_lnslp <- bf(log(abs_slope) | se(sqrt(V_abs_slope)/abs_slope)  ~  - 1 +  Category + (1|q|parameter_group))
mod_lnint <- bf(log(abs_int) | se(sqrt(V_abs_int)/abs_int)  ~  - 1 +  Category + (1|q|parameter_group))

fit_3b <- brm(mod_lnsd + mod_lnslp + mod_lnint,
              data = dat,
              chains = 2, cores = 2, iter = 4000, warmup = 1000,
              backend = "cmdstanr")

summary(fit_3b)

# saving the model
saveRDS(fit_3b, file = here("data", "fit_.rds"))


# creating added precisoin

dat %>%  mutate(pre_slp_int = 1/sqrt(V_abs_int/abs_int^2 + V_abs_slope/abs_slope^2),
                pre_slp_sd =  1/sqrt(V_abs_slope/abs_slope^2 + V_abs_lnVR/abs_lnVR^2),
                pre_int_sd = 1/sqrt(V_abs_int/abs_int^2 + V_abs_lnVR/abs_lnVR^2)
) -> dat 


f1 <- ggplot(data = dat) +
  geom_point(aes(x = log(abs_slope), y = log(abs_int), col = Category, size = pre_slp_int)) + 
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  labs(x = "ln(Absolute difference in standardized slopes)" , y = "ln(Absolute difference in standardized intercepts)")+
  labs(color='Trait types', size = "Precison") +
  annotate(geom="text", x=2.5, y = -5.5, label="r = 0.56 [0.42, 0.68]", size = 3)+
  theme_bw()  +
  theme(legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=10))+
  guides(col = "none")



f2 <- ggplot(data = dat) +
  geom_point(aes(x = log(abs_slope), y = log(abs_lnVR), col = Category, size = pre_slp_sd)) + 
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  labs(x = "ln(Absolute difference in standardized slopes)" , y = "ln(Absolute relative difference in SD)") +
  labs(color='Trait types', size = "Precison") +
  annotate(geom="text", x=2.5, y = -4.8, label="r = 0.19 [0.01., 0.36]", size = 3)+
  theme_bw()   +
  theme(legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=10))+
  guides(col = "none")

f3 <- ggplot(data = dat) +
  geom_point(aes(x = log(abs_int), y = log(abs_lnVR), col = Category, size = pre_int_sd)) + 
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  labs(x = "ln(Absolute difference in standardized intercepts)" , y = "ln(Absolute relative difference in SD)") +
  labs(color='Trait types', size = "Precison") +
  annotate(geom="text", x= - 0.2, y = -4.8, label="r = 0.07 [-0.10, 0.24]", size = 3)+
  theme_bw() +
  theme(legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=10))




f3/f2/f1  + plot_annotation(tag_levels = 'A')


#################
# Figure 4
#################

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

