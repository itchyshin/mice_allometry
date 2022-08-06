# script to get effect size estimates comparing females and males
# TODO meta-data needs to be done!!

# the data set contains 301 traits - but 4 duplicates

# package
library(nlme)
library(broom.mixed) # great package to get things out of mixed model
library(purrr)
library(metafor)
library(devtools)
library(tidyverse)
library(tibble)
library(here)
library(performance)
library(MuMIn)

# loading ddata
dat_list <- readRDS(here("data/dat_list2.rds"))

# TODO - getting data with sub-strains

short_list <- discard(dat_list, ~ .x[["nstrain"]][[1]] == 1)

# TODO getting new parameter estimates


# function to compare two models
get_strain_p<- function(i){
  
  # centering weights separately for each 
  ln_c_weight <- scale(log(i[["weight"]]), center = TRUE, scale = TRUE)
  i[,"ln_c_weight"] <- ln_c_weight
  
  if (i[["nmeta"]][1] == 1) {
    
    # model with strain as random factor 
    model_1 <- lme(log(data_point2) ~ sex*ln_c_weight, 
                   random = list(#metadata_group = ~ ln_c_weight, 
                     strain_name = ~ ln_c_weight,
                     date_of_experiment = ~ 1),
                   #weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    
    # model without strain as random factor 
    model_2 <- lme(log(data_point2) ~ sex*ln_c_weight, 
                   random = list(#metadata_group = ~ ln_c_weight, 
                     #strain_name = ~ ln_c_weight,
                     date_of_experiment = ~ 1),
                   #weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    
    
    
  }  else {

    # model with strain as random factor 
    model_1 <- lme(log(data_point2) ~ sex*ln_c_weight, 
                   random = list(metadata_group = ~ ln_c_weight, 
                                 strain_name = ~ ln_c_weight,
                                 date_of_experiment = ~ 1),
                   #weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    # model without strain as random factor 
    model_2 <- lme(log(data_point2) ~ sex*ln_c_weight, 
                   random = list(metadata_group = ~ ln_c_weight, 
                                 #strain_name = ~ ln_c_weight,
                                 date_of_experiment = ~ 1),
                   #weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    
  }
  
  # anova 
  p_value <- anova(model_1, model_2)$p[[2]]
  delta_aic <- anova(model_1, model_2)$AIC[1] - anova(model_1, model_2)$AIC[2] 
  
  # get parameters
  parameter_name <- tolower(i[["parameter_name"]][1])
  procedure_name <- i[["procedure_name"]][1]# "procedure_name"
  
  paras <- data.frame(parameter_name, procedure_name, delta_aic, p_value)
  names(paras) <- c('parameter_name', 'procedure_name', 'delta_aic', 'p_value') # variance component
  invisible(paras)
  
}

# getting ride of traits which do not run
get_para_poss <- possibly(.f = get_parmetersN, 
                          otherwise = NULL)

#run individual tests on single matrices extracted from dat_list
testInd_dat1<- short_list[[1]] #male is 1
#res1<-get_parmetersN(testInd_dat1) # this works


# # test
ln_c_weight <- scale(log(testInd_dat1[["weight"]]), center = TRUE, scale = TRUE)
testInd_dat1[,"ln_c_weight"] <- ln_c_weight

model_1 <- lme(log(data_point2) ~ sex*ln_c_weight,
                  random = list(metadata_group = ~ ln_c_weight,
                                strain_name = ~ ln_c_weight,
                                date_of_experiment = ~ 1),
                  #weights = varIdent(form = ~1 | sex),
                  control = lmeControl(opt = "optim"),
                  data = testInd_dat1)
summary(model_1)

model_2 <- lme(log(data_point2) ~ sex*ln_c_weight,
               random = list(metadata_group = ~ ln_c_weight,
                             #strain_name = ~ ln_c_weight,
                             date_of_experiment = ~ 1),
               #weights = varIdent(form = ~1 | sex),
               control = lmeControl(opt = "optim"),
               data = testInd_dat1)
summary(model_2)

anova(model_1, model_2)


#getting all we want
test_females <- broom.mixed::tidy(model_test)
test_females


# model_test2 <- lmer(log(data_point) ~ sex*ln_c_weight,
#                   random = list(metadata_group = ~ ln_c_weight, date_of_experiment = ~ 1),
#                   #weights = varIdent(form = ~1 | sex),
#                   control = lmeControl(opt = "optim"),
#                   data = testInd_dat1)
# r2_nakagawa(model_test2)

# more tests
# 
# ln_c_weight <- groupScale(log(testInd_dat1[["weight"]]) ~ testInd_dat1[["sex"]])
# testInd_dat1[,"ln_c_weight"] <- ln_c_weight
# 
# model_f <- lme((scale(log(data_point))) ~ sex*ln_c_weight, 
#                random = ~ 1|metadata_group,
#                weights = varIdent(form = ~1 | sex),
#                control = lmeControl(opt = "optim"),
#                data = testInd_dat1)


# TODO - test (1) meta-slope, (2) strain-slope and (3) both-slope

# (1) meta-slope
# this does not make much sense

# (2) strain-slope 
# fin_dat1 <-map_dfr(dat_list2, get_para_poss)
# fin_dat1 <- data.frame(fin_dat1, row.names = NULL)
# 
# dim(fin_dat1)

# [1] 377  21

# (3) both-slope
processing <-map_dfr(dat_list2, get_para_poss)
dat <- data.frame(processing, row.names = NULL)
dim(dat)

# TODO - remove "test duration"
remove <- which(dat$parameter_name == "test duration")
dat <- dat[-remove, ]
dim(dat)

write_csv(dat, here("data/test3.csv"))

dat <- read.csv(here("data/test3.csv"))




#Fin_dat[Fin_dat$parameter_name == "light side time spent",]

#remove "body weight" duplicates 
#dat<-Fin_dat[-c(27:30),] 

# grouping for category and parameter_group
#dat_category<-read_csv(here("data/cateogry_parameter2.csv")) 

dat_category<-read_csv(here("data/cateogry_parameter3.csv")) 

dat %>% left_join(dat_category, by = ("parameter_name" = "parameter_name") ) %>% arrange(Category)  -> dat

dim(dat)

#write_csv(dat, here("data/test4.csv"))

write_csv(dat, here("data/data_parameters5.csv"))

# TODO meta-data needs to be done!!

