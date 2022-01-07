# script to get effect size estimates comparing females and males
# TODO meta-data needs to be done!!

# the data set contains 301 taits - but 4 duplicates

# package
library(nlme)
library(broom.mixed) # great package to get things out of mixed model
library(purrr)
library(metafor)
library(devtools)
library(tidyverse)
library(tibble)
library(here)

# TODO getting new parameter estimates

# custom function  for within-group cenering (or z transformation)
groupScale <- function(formula, data=NULL, center=TRUE, scale=FALSE){
  if(is.null(data)) data <- model.frame(formula)
  scaled <- rep(NA,nrow(data)) #empty vector
  for(i in unique(data[,2])){
  elements <- which(data[,2]==i)
  scaled[elements] <- scale(data[elements,1], scale=scale, center=center) 
}
return(scaled)
}



# function to get what we need from these 2 models (you can include models in this function as well)
get_parmetersN<- function(i){
  
  # centering weights separately for each 
  
  ln_c_weight <- groupScale(log(i[["weight"]]) ~ i[["sex"]])
  i[,"ln_c_weight"] <- ln_c_weight
  
  model_f <- lme((scale(log(data_point))) ~ sex*ln_c_weight, 
                 random = ~ 1|metadata_group,
                 weights = varIdent(form = ~1 | sex),
                 control = lmeControl(opt = "optim"),
                 data = i)
  
  model_m <- lme((scale(log(data_point))) ~ relevel(sex, ref = "male")*ln_c_weight, 
                 random = ~ 1|metadata_group,
                 weights = varIdent(form = ~1 | sex),
                 control = lmeControl(opt = "optim"),
                 data = i)
  # getting all we want
  females <- broom.mixed::tidy(model_f)
  males <- broom.mixed::tidy(model_m)
  # gets variance weights
  weights <- attr(model_f$modelStruct$varStruct, "weights")
  male_correction <- 1/weights[which(names(weights) == "male")[1]]
  female_correction <- 1/weights[which(names(weights) == "female")[1]]
  
  # get parameters
  parameter_name <- tolower(i[["parameter_name"]][1])
  m_n <- sum(i[["sex"]] == "male") # sample size for males 
  f_n <- sum(i[["sex"]] == "female") # N fo females
  f_intercept <- as.numeric(females[1, 4])
  f_intercept_se <- as.numeric(females[1, 5])
  f_slope <- as.numeric(females[3, 4])
  f_slope_se <- as.numeric(females[3, 5])
  m_intercept <- as.numeric(males[1, 4])
  m_intercept_se <- as.numeric(males[1, 5])
  m_slope  <- as.numeric(males[3, 4])
  m_slope_se  <- as.numeric(males[3, 5])
  fm_diff_int  <- as.numeric(males[2, 4])
  fm_diff_int_se  <- as.numeric(males[2, 5])
  fm_diff_int_p  <- as.numeric(males[2, 8])
  fm_diff_slope <- as.numeric(males[4, 4])
  fm_diff_slope_se <- as.numeric(males[4, 5])
  fm_diff_slope_p <- as.numeric(males[4, 8])
  group_sd <- as.numeric(females[5, 4])
  f_sd <- as.numeric(females[6, 4])*female_correction
  m_sd <- as.numeric(females[6, 4])*male_correction
  
  
  # putting it together
  paras <- c(parameter_name, f_n, m_n, f_intercept, f_intercept_se, f_slope, f_slope_se, 
             m_intercept, m_intercept_se, m_slope, m_slope_se, 
             fm_diff_int, fm_diff_int_se, fm_diff_int_p,
             fm_diff_slope, fm_diff_slope_se, fm_diff_slope_p,
             group_sd, f_sd, m_sd)
  names(paras) <- c('parameter_name','f_n', 'm_n','f_intercept', 'f_intercept_se', 'f_slope', 'f_slope_se', 
                    'm_intercept', 'm_intercept_se', 'm_slope', 'm_slope_se', 
                    'fm_diff_int', 'fm_diff_int_se', 'fm_diff_int_p',
                    'fm_diff_slope', 'fm_diff_slope_se', 'fm_diff_slope_p',
                    'group_sd', 'f_sd', 'm_sd')
  invisible(paras)
  
}

# loading data
dat_list <- readRDS(here("data/dat_list.rds"))

# grouping for category and parameter_group
dat_category<-read_csv(here("data/cateogry_parameter2.csv")) 

#run individual tests on single matrices extracted from dat_list
testInd_dat1<- dat_list[[1]] #male is 1
res1<-get_parmetersN(testInd_dat1) # this works

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


#run function across list of matrices
Fin_dat<-map_dfr(dat_list, get_parmetersN)

Fin_dat[Fin_dat$parameter_name == "light side time spent",]

#remove "body weight" duplicates 
dat<-Fin_dat[-c(27:30),] 

dat %>% left_join(dat_category, by = ("parameter_name" = "parameter_name") ) %>% arrange(Category)  -> dat

write_csv(dat, here("data/data_parameters.csv"))

# TODO meta-data needs to be done!!

