# getting SMD (Cohen's d) and lnRR
# we use metafor - escalc

library(metafor)
library(tidyverse)
library(here)

dat_list <- readRDS(here("data/dat_list2.rds"))

i <- dat_list[[2]]


extra_effect <- function (i){

i %>% group_by(sex) %>% summarise(mean = mean(data_point),
                                     sd = sd(data_point),
                                     n = n()) -> sex_specific

sex_specific
sex_specific[1, "mean"]

# get parameters
parameter_name <- tolower(i[["parameter_name"]][1])
mean_female <- sex_specific[1, "mean"]
mean_male <- sex_specific[2, "mean"]
sd_female <- sex_specific[1, "sd"]
sd_male <- sex_specific[2, "sd"]
n_female <- sex_specific[1, "n"]
n_male <- sex_specific[2, "n"]

# putting it together
paras <- data.frame(parameter_name, 
                    mean_female, mean_male,
                    sd_female, sd_male,
                    n_female, n_male)
names(paras) <- c('parameter_name',
                  'mean_female', 'mean_male',
                  'sd_female', 'sd_male',
                  'n_female', 'n_male') # variance component


paras <- escalc("SMD", 
             m1i = mean_female, m2i = mean_male, 
             sd1i = sd_female, sd2i = sd_male, 
             n1i = n_female, n2i = n_male, 
             data = paras, var.names=c("SMD","v_SMD"))
paras <- escalc("ROM", 
             m1i = mean_female, m2i = mean_male, 
             sd1i = sd_female, sd2i = sd_male, 
             n1i = n_female, n2i = n_male, 
             data = paras, var.names=c("lnRR","v_lnRR"), replace = F)

invisible(paras)

}

extra_dat <-map_dfr(dat_list, extra_effect)
