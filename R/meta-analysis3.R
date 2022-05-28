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
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(posterior)
library(bayesplot)
#library(ggblend)

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
                      V_abs_lnVR = folded_v(lnVR, VlnVR),
                      total_n = f_n + m_n)

cbpl <- c("#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00",
          "#CC79A7", "#56B4E9", "#AA4499", "#DDCC77")

#####################
#
#####################

modelr0 <- rma.mv(yi = Zr, 
                  V= VZr, 
                  random = list(~1| Category, ~1| parameter_group, ~1|obs), 
                  data = dat)
summary(modelr0)
robust(modelr0, cluster  =  dat$parameter_group)

#funnel(modelr0)
i2_ml(modelr0)


point.size = 2
branch.size = 3.5

t1 <- orchard_plot2(modelr0, mod = "Int", xlab = "Zr (transformed variance accounted for)", angle = 45, 
                    point.size = point.size, branch.size = branch.size, k = F, N = dat$total_n) +
  scale_y_discrete(labels = "") +
  scale_fill_manual(values = "#999999") +
  scale_colour_manual(values = "#999999") #+
  #xlim(c(-0.5, 1.5))

# meta-regression
modelr1 <- rma.mv(yi = Zr, mod = ~ Category - 1,
               V= VZr, 
               random = list(~1| parameter_group, ~1|obs), 
               data = dat)

summary(modelr1)
robust(modelr1, cluster  =  dat$parameter_group)
r2_ml(modelr1)

t2 <- orchard_plot2(modelr1, mod = "Category", xlab = "Zr (transformed variance accounted for)", angle = 45,  point.size = point.size, k = F, N = dat$total_n, branch.size = branch.size,) + 
  scale_y_discrete(labels = rep("", 9)) +
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) #+
  #xlim(c(-0.5, 1.5))

#t1/t2

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

#funnel(modelia)
i2_ml(modelia)


#point.size = 1.5

p1 <- orchard_plot2(modelia, mod = "Int", xlab = "Absolute difference in standardized intercepts  (F-M)", angle = 45, point.size = point.size, N = dat$total_n, legend.on = FALSE, branch.size = branch.size,) +
  scale_y_discrete(labels = "Overall") +
  scale_fill_manual(values = "#999999") +
  scale_colour_manual(values = "#999999") #+
  #xlim(c(-0.5, 1.5))

# meta-regression
model1a <- rma.mv(yi = abs_int, V= V_abs_int, mod = ~ Category - 1,
                  random = list(~1| parameter_group, ~1|obs), 
                  data = dat)
summary(model1a)
robust(model1a, cluster  =  dat$parameter_group)
r2_ml(model1a)

p2 <- orchard_plot2(model1a, mod = "Category", xlab = "Absolute difference in standardized intercepts  (F-M)", angle = 45,  point.size = point.size, N = dat$total_n, legend.on = FALSE, branch.size = branch.size,) + 
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) #+
  #xlim(c(-0.5, 1.5))

####################
# slope difference
######################

modelsa <- rma.mv(yi = abs_slope, V= V_abs_slope, 
                  random = list(~1| Category, ~1| parameter_group, ~1|obs), 
                  data = dat)
summary(modelsa) # not sig this means sometimes male is high other times female has steaper slops
robust(modelsa, cluster  =  dat$parameter_group)

#funnel(modelsa)
i2_ml(modelia)

p3 <- orchard_plot2(modelsa, mod = "Int", xlab = "Absolute difference in standardized slopes (F-M)", angle = 45,  point.size = point.size, k = F, N = dat$total_n, legend.on = FALSE, branch.size = branch.size,) +
  scale_y_discrete(labels = "") +
  scale_fill_manual(values = "#999999") +
  scale_colour_manual(values = "#999999") #+
  #xlim(c(-1.5, 10))

# meta-regression
model2a <- rma.mv(yi = abs_slope, V= V_abs_slope,
                  mod = ~ Category - 1,
                  random = list(~1| parameter_group, ~1|obs), 
                  data = dat)
summary(model2a)
robust(model2a, cluster  =  dat$parameter_group)
r2_ml(model2a)

p4 <- orchard_plot2(model2a, mod = "Category", xlab = "Absolute difference in standardized slopes (F-M)", angle = 45, cb = F,  point.size = point.size, k = F, N = dat$total_n, legend.on = FALSE, branch.size = branch.size,) + 
  scale_y_discrete(labels = rep("", 9)) +
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) #+
  #xlim(c(-1.5, 10))

#################
# sd difference
#################

modelsda <- rma.mv(yi = abs_lnVR, V= V_abs_lnVR, 
                   random = list(~1| Category, ~1| parameter_group, ~1|obs), 
                   data = dat)
summary(modelsda)
robust(modelsda, cluster  =  dat$parameter_group)
#funnel(modelsda)
i2_ml(modelsda)

p5 <- orchard_plot2(modelsda, mod = "Category", xlab = "Absolute relative difference in SD (lnVR: F/M)", angle = 45,  point.size = point.size, k = F, N = dat$total_n, legend.on = FALSE, branch.size = branch.size,) +
  scale_y_discrete(labels = "") +
  scale_fill_manual(values = "#999999") +
  scale_colour_manual(values = "#999999") #+
  #xlim(c(-0.2, 1.9))

# meta-regression
model3a <- rma.mv(yi = abs_lnVR, V= V_abs_lnVR, mod = ~ Category - 1,
                  random = list(~1| parameter_group, ~1|obs), 
                  data = dat)
summary(model3a)
r2_ml(model3a)

p6 <- orchard_plot2(model3a, mod = "Category", xlab = "Absolute relative difference in SD (lnVR: F/M)", angle = 45, cb = F,  point.size = point.size, k = F, N = dat$total_n, legend.on = FALSE, branch.size = branch.size,) + 
  scale_y_discrete(labels = rep("", 9)) +
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) #+
  #xlim(c(-0.2, 1.9))

###############
# Fig 3
###############

(p1 | p3 | p5 | t1) / (p2 | p4 | p6 | t2)  + plot_layout(heights = c(1, 3)) + plot_annotation(tag_levels = 'A')

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
mod_lnzr <- bf(log(Zr) | se(sqrt(VZr)/Zr)  ~  - 1 +  Category + (1|q|parameter_group))

fit_4b <- brm(mod_lnsd + mod_lnslp + mod_lnint + mod_lnzr,
              data = dat,
              chains = 2, cores = 2, iter = 4000, warmup = 1000,
              backend = "cmdstanr"
              )

summary(fit_4b)

# saving the model
saveRDS(fit_4b, file = here("data", "fit2.rds"))

fit_4b <- readRDS(here("data", "fit2.rds"))

# creating added precision
# 

# 3 main parameters

# dat %>%  mutate(pre_slp_int = 1/sqrt(V_abs_int/abs_int^2 + V_abs_slope/abs_slope^2),
#                 pre_slp_sd =  1/sqrt(V_abs_slope/abs_slope^2 + V_abs_lnVR/abs_lnVR^2),
#                 pre_int_sd = 1/sqrt(V_abs_int/abs_int^2 + V_abs_lnVR/abs_lnVR^2),
#                 pre_int_fit =  1/sqrt(V_abs_int/abs_int^2 + VZr/Zr^2),
#                 pre_slp_fit =  1/sqrt(V_abs_slope/abs_slope^2 +  VZr/Zr^2),
#                 pre_sd_fit =  1/sqrt(V_abs_lnVR/abs_lnVR^2 +  VZr/Zr^2),
#                 
# ) -> dat 


f1 <- ggplot(data = dat) +
  geom_point(aes(x = log(abs_slope), y = log(abs_int), col = Category, size = total_n)) + 
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  labs(x = "ln(Absolute difference in standardized slopes)" , y = "ln(Absolute difference in standardized intercepts)")+
  labs(color='Trait types', size = "Sample size (N)") +
  annotate(geom="text", x = - 7.8, y = -1, label="r = 0.74 [0.67, 0.81]", size = 3)+
  theme_bw()  +
  theme(legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=10))+
  guides(col = "none", size = "none") 

f2 <- ggplot(data = dat) +
  geom_point(aes(x = log(abs_slope), y = log(abs_lnVR), col = Category, size = total_n)) + 
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  labs(x = "ln(Absolute difference in standardized slopes)" , y = "ln(Absolute relative difference in SD)") +
  labs(color='Trait types', size = "Sample size (N)") +
  annotate(geom="text", x= -7.5, y = 0.5, label="r = 0.09 [-0.05., 0.24]", size = 3)+
  theme_bw()   +
  theme(legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=10))+
  guides(size = "none", col = "none") 
  #scale_size_continuous(breaks = c(200, 2000, 20000), guide = guide_legend()) +

f3 <- ggplot(data = dat) +
  geom_point(aes(x = log(abs_int), y = log(abs_lnVR), col = Category, size = total_n)) + 
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  labs(x = "ln(Absolute difference in standardized intercepts)" , y = "ln(Absolute relative difference in SD)") +
  labs(color='Trait types', size = "Sample size (N)") +
  annotate(geom="text", x= - 10, y = 0.5, label="r = 0.04 [-0.10, 0.17]", size = 3)+
  theme_bw() +
  theme(legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=10)) +
  guides(size = "none", col = "none")


f4 <- ggplot(data = dat) +
  geom_point(aes(y = log(Zr), x = log(abs_int), col = Category, size = total_n)) + 
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  labs(y = "Zr (transformed variance accounted for)", x =   "ln(Absolute difference in standardized intercepts)") +
  labs(color='Trait types',size = "Sample size (N)") +
  annotate(geom="text", x= -2.5, y = -6, label="r = 0.70 [0.62., 0.77]", size = 3)+
  theme_bw()   +
  theme(legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=10)) +
  guides(size = "none") +
  theme(legend.position= c(0.03, 0.97), legend.justification = c(0, 0.97))
  

f5 <- ggplot(data = dat) +
  geom_point(aes(y = log(Zr), x = log(abs_slope), col = Category, size = total_n)) + 
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  labs(y = "Zr (transformed variance accounted for)" , x = "ln(Absolute difference in standardized slopes)") +
  labs(color='Trait types', size = "Sample size (N)") +
  annotate(geom="text", x=  0, y = -6, label="r = 0.39 [0.26, 0.51]", size = 3)+
  theme_bw() +
  theme(legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=10))+
  guides(col = "none") +
  scale_size_continuous(breaks = c(200, 2000, 20000), guide = guide_legend()) +
  theme(legend.position = c(0.03, 0.97), legend.justification = c(0, 0.97))


f6 <- ggplot(data = dat) +
  geom_point(aes(y = log(Zr), x = log(abs_lnVR), col = Category, size = total_n)) + 
  scale_fill_manual(values = cbpl) +
  scale_colour_manual(values = cbpl) +
  labs(y = "Zr (transformed variance accounted for)", x = "ln(Absolute relative difference in SD)" )+
  labs(color='Trait types', size = "Sample size (N)") +
  annotate(geom="text", x= -0.25, y = - 6, label="r = 0.16 [0.02, 0.30]", size = 3)+
  theme_bw()  +
  theme(legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=10))+
  guides(col = "none", size = "none") #+
  #scale_size_continuous(breaks = c(200, 2000, 20000), guide = guide_legend()) 


 (f3|f2)/(f1|f4)/(f5|f6)  + plot_annotation(tag_levels = 'A')


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

