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

# TODO - getting data with sub-strains

# loading data
allometry <- readRDS(here("data/allometryNEW.rds"))


#STEP 1 remove rows with missing  data and NA 
allometrynew<-allometry[complete.cases(allometry),]

# getting rid of NA for data_point and weight

allometrynew2 <- allometrynew %>% 
  filter(!is.na(data_point), !is.na(weight)) %>% 
  group_by(parameter_name, sex, metadata_group, strain_name) %>%
  mutate(count = n()) %>% 
  ungroup() %>% 
  group_by(parameter_name) %>% # adjusting interval data
  mutate(min_val = min(data_point),
         data_point2 = if_else(min_val > 0, data_point, data_point + abs(min_val)),
         min_val2 = min( data_point[data_point!=min(data_point)]),
         data_point2 = if_else(min_val == 0, data_point2 + min_val2, data_point2),
         ratio_int =  if_else(min_val > 0, "ratio", "interval"),
         new_min = min(data_point2),
         nmeta = n_distinct(metadata_group),
         nstrain = n_distinct(strain_name),
         sex = as.factor(sex),
         parameter_name = if_else(parameter_name == "Latency to fall_Mean", 
                                  "Latency to fall mean"  , parameter_name)) %>% 
  ungroup() %>% 
  filter(count > 49) %>% # this can be adjusted
  filter(parameter_name != "BMC/Body weight", 
         parameter_name != "Body weight",  
         parameter_name != "Body Weight", 
         parameter_name != "Body weight after experiment" , 
         parameter_name != "Body weight before experiment",
         parameter_name != "Test duration") %>% 
  filter(!is.infinite(data_point2), !is.infinite(log(data_point2))) # removing infite and 0

allometrynew2$parameter_name[allometrynew2$parameter_name == "latency to fall_mean"] 

dim(allometry)
dim(allometrynew)
dim(allometrynew2)

# the number of traits
length(unique(allometrynew2$parameter_name))

# the number of substrains
length(unique(allometrynew2$strain_name))
# check there is no 0

sum(is.infinite(log(allometrynew2$data_point2)))

# the number of interval scale traits
allometrynew2 %>% group_by(parameter_name) %>% summarise(ratio_int = ratio_int[1]) -> sum_ri
sum(sum_ri$ratio_int == "interval")


#split dataframe by parameter to generate a list of dfs
all_list<-split(allometrynew2, allometrynew2$parameter_name)

saveRDS(all_list, file = here("data", "dat_list2.rds"))


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
  
  if(i[["nmeta"]][1] == 1 && i[["nstrain"]][1] == 1){
  
    # female model 
    model_f <- lme(log(data_point2) ~ sex*ln_c_weight, 
                   random = list(#metadata_group = ~ ln_c_weight, 
                                #strain_name = ~ ln_c_weight,
                                 date_of_experiment = ~ 1),
                   weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    
    # male model
    model_m <- lme(log(data_point2) ~ relevel(sex, ref = "male")*ln_c_weight, 
                   random = list(#metadata_group = ~ ln_c_weight, 
                                 #strain_name = ~ ln_c_weight,
                                 date_of_experiment = ~ 1),
                   weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    # neutral model
    model_n <- lme(log(data_point2) ~ sex*ln_c_weight, 
                   random = list(#metadata_group = ~ ln_c_weight, 
                                 #strain_name = ~ ln_c_weight,
                                 date_of_experiment = ~ 1),
                   #weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
      
  } else if (i[["nmeta"]][1] == 1) {
    
    # female model 
    model_f <- lme(log(data_point2) ~ sex*ln_c_weight, 
                   random = list(#metadata_group = ~ ln_c_weight, 
                     strain_name = ~ ln_c_weight,
                     date_of_experiment = ~ 1),
                   weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    
    # male model
    model_m <- lme(log(data_point2) ~ relevel(sex, ref = "male")*ln_c_weight, 
                   random = list(#metadata_group = ~ ln_c_weight, 
                     strain_name = ~ ln_c_weight,
                     date_of_experiment = ~ 1),
                   weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    # neutral model
    model_n <- lme(log(data_point2) ~ sex*ln_c_weight, 
                   random = list(#metadata_group = ~ ln_c_weight, 
                     strain_name = ~ ln_c_weight,
                     date_of_experiment = ~ 1),
                   #weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    
  } else if (i[["nstrain"]][1] == 1){
    
    # female model 
    model_f <- lme(log(data_point2) ~ sex*ln_c_weight, 
                   random = list(metadata_group = ~ ln_c_weight, 
                     #strain_name = ~ ln_c_weight,
                     date_of_experiment = ~ 1),
                   weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    
    # male model
    model_m <- lme(log(data_point2) ~ relevel(sex, ref = "male")*ln_c_weight, 
                   random = list(metadata_group = ~ ln_c_weight, 
                     #strain_name = ~ ln_c_weight,
                     date_of_experiment = ~ 1),
                   weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    # neutral model
    model_n <- lme(log(data_point2) ~ sex*ln_c_weight, 
                   random = list(metadata_group = ~ ln_c_weight, 
                     #strain_name = ~ ln_c_weight,
                     date_of_experiment = ~ 1),
                   #weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
  } else {
    # female model 
    model_f <- lme(log(data_point2) ~ sex*ln_c_weight, 
                   random = list(metadata_group = ~ ln_c_weight, 
                                 strain_name = ~ ln_c_weight,
                                 date_of_experiment = ~ 1),
                   weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    
    # male model
    model_m <- lme(log(data_point2) ~ relevel(sex, ref = "male")*ln_c_weight, 
                   random = list(metadata_group = ~ ln_c_weight, 
                                 strain_name = ~ ln_c_weight,
                                 date_of_experiment = ~ 1),
                   weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
    # neutral model
    model_n <- lme(log(data_point2) ~ sex*ln_c_weight, 
                   random = list(metadata_group = ~ ln_c_weight, 
                                 strain_name = ~ ln_c_weight,
                                 date_of_experiment = ~ 1),
                   #weights = varIdent(form = ~1 | sex),
                   control = lmeControl(opt = "optim"),
                   data = i)
  }
  # getting all we want
  females <- broom.mixed::tidy(model_f)
  males <- broom.mixed::tidy(model_m)
  # gets variance weights
  weights <- attr(model_f$modelStruct$varStruct, "weights")
  male_correction <- 1/weights[which(names(weights) == "male")[1]]
  female_correction <- 1/weights[which(names(weights) == "female")[1]]
  
  # get parameters
  parameter_name <- tolower(i[["parameter_name"]][1])
  procedure_name <- i[["procedure_name"]][1]# "procedure_name"
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
  
  # variance component
  #group_sd <- as.numeric(VarCorr(model_f)[,2][2])
  #g_slope_sd <- as.numeric(VarCorr(model_f)[,2][3])
  #batch_sd <- as.numeric(VarCorr(model_f)[,2][5])
  f_sd <- as.numeric(tail(VarCorr(model_f)[,2],1))*female_correction
  m_sd <- as.numeric(tail(VarCorr(model_f)[,2],1))*male_correction
  
  # model fit
  r_m <- sqrt(MuMIn::r.squaredGLMM(model_n)[1,1])
  r_c <- sqrt(MuMIn::r.squaredGLMM(model_n)[1,2])
  # putting it together
  paras <- data.frame(parameter_name, procedure_name, 
             f_n, m_n, f_intercept, f_intercept_se, f_slope, f_slope_se, 
             m_intercept, m_intercept_se, m_slope, m_slope_se, 
             fm_diff_int, fm_diff_int_se, fm_diff_int_p,
             fm_diff_slope, fm_diff_slope_se, fm_diff_slope_p,
             f_sd, m_sd, r_m, r_c)
  names(paras) <- c('parameter_name', 'procedure_name', 
                    'f_n', 'm_n','f_intercept', 'f_intercept_se', 'f_slope', 'f_slope_se',
                    'm_intercept', 'm_intercept_se', 'm_slope', 'm_slope_se',
                    'fm_diff_int', 'fm_diff_int_se', 'fm_diff_int_p',
                    'fm_diff_slope', 'fm_diff_slope_se', 'fm_diff_slope_p',
                    'f_sd', 'm_sd', 'r_m', 'r_c') # variance component
  invisible(paras)
  
}

# getting ride of traits which do not run
get_para_poss <- possibly(.f = get_parmetersN, 
                          otherwise = NULL)


# loading data
#dat_list <- readRDS(here("data/dat_list.rds"))
dat_list2 <- readRDS(here("data/dat_list2.rds"))

#run individual tests on single matrices extracted from dat_list
testInd_dat1<- dat_list2[[2]] #male is 1
#res1<-get_parmetersN(testInd_dat1) # this works


# # test
ln_c_weight <- groupScale(log(testInd_dat1[["weight"]]) ~ testInd_dat1[["sex"]])
testInd_dat1[,"ln_c_weight"] <- ln_c_weight

model_test <- lme(log(data_point2) ~ sex*ln_c_weight,
                  random = list(metadata_group = ~ ln_c_weight,
                                strain_name = ~ ln_c_weight,
                                date_of_experiment = ~ 1),
                  weights = varIdent(form = ~1 | sex),
                  control = lmeControl(opt = "optim"),
                  data = testInd_dat1)
summary(model_test)


#getting all we want
test_females <- broom.mixed::tidy(model_test)
test_females
VarCorr(model_test)


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

