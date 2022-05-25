# analysis
#TODO - need to get contrasts between male SD and female SD - done
# TODO - sorting out to the scenario!! - done
# TODO - questions (we do not use absolute values for this analysis - so you cannot really add and subtract; e.g. time spend in light or dark areas)
# TODO - how to present adjusted version?? 
# TODO - many of repeats can be made into functions

# package
library(purrr)
library(metafor)
library(tidyverse)
library(here)
library(poolr)
library(patchwork)

# > citation("poolr")
# 
# To cite package ‘poolr’ in publications use:
#   
#   Ozan Cinar and Wolfgang Viechtbauer (2021). poolr: Methods for Pooling P-Values
# from (Dependent) Tests. R package version 1.0-0.
# https://CRAN.R-project.org/package=poolr

# first getting p values - the contrasts between males and females for 

dat <-read_csv(here("data/data_parameters5.csv"))

#assess number of traits with sig shifts in intercept and slope

# getting lnVR to compare SDs and SD

dat %>% mutate(lnVR = log(f_sd/m_sd) + 1/(2*(f_n-3)) - 1/(2*(m_n-3)), 
               VlnVR = 1/(2*(f_n-3)) + 1/(2*(m_n-3)), 
               low_lnVR = lnVR - qnorm(0.975)*VlnVR, 
               high_lnVR = lnVR + qnorm(0.975)*VlnVR,
               t_val_sd = lnVR/sqrt(VlnVR),
               p_val_sd = 2*(1-pt(abs(t_val_sd), f_n-1 + m_n-1)),
               # r squared
               Zr = atanh(r_m),
               VZr = 1/((f_n + m_n) - 3)
               ) -> dat

write_csv(dat, here("data/data_parameters6.csv"))

# significantly different 
length(which(dat$p_val_sd <= 0.05))
hist(log(dat$p_val_sd)) # p = 0.05 ~ - 3
# 151 out of 297 - this means almost all residual SD are heteroscadaistic (it is probably due to large n in our dataset) = over 50%!!

# kinda very surprising! (very powerful) - in Susi's elife paper - we only go ~ 40% of lnVR significant so it matches

# the whole data set

dim(dat)


#16 out of 375 traits sig slope diff - scenario A
dat_slopes <-dat %>%
  filter (fm_diff_slope_p <= 0.05 & fm_diff_int_p > 0.05)

dim(dat_slopes) 

#165 out of 375 traits sig intercept diff  same slope - scenario B
dat_int<- dat %>%
  filter (fm_diff_int_p <= 0.05 & fm_diff_slope_p >0.05)
dim(dat_int) 

#81 out of 375 sig intercept and slope diff - scenario C
dat_intSlopes<-dat %>%
  filter (fm_diff_int_p <= 0.05 & fm_diff_slope_p <= 0.05)
dim(dat_intSlopes) 

#113 (of out of 375) no sig difference between intercept and slope - scenario D
dat_intslopesNS<- dat %>%
  filter (fm_diff_slope_p >0.05 & fm_diff_int_p > 0.05)
dim(dat_intslopesNS) 

# here we need to collapse p values which are related
# split data into 2 ones with replications within parameter_group

dat %>% group_by(parameter_group) %>% mutate(count = n()) -> dat
# 
dat1 <- dat[which(dat$count == 1), ]
# 
dim(dat1)

length(which(dat1$p_val_sd <= 0.05))
# 112 out of 168
# 
dat2 <- dat[-which(dat$count == 1), ]
# 
# # We will use poolr to get 
# 
test_dat <-dat2[3:6, ]
# 
# Rmat <- diag(4)
# Rmat[lower.tri(Rmat)]<-0.8
# Rmat[upper.tri(Rmat)]<-0.8
# 
# p_mod <- fisher(test_dat$p_val_sd, adjust = "empirical", R = Rmat)
# 
# p_mod$p

# function to get merged p value for SD

p_mod_sd <-function(data){
  
  len <- dim(data)[1]
  Rmat <- matrix(0.8, nrow = len, ncol = len)
  diag(Rmat) <- 1
  
  p_mod <- fisher(data$p_val_sd, adjust = "liji", R = Rmat)
  p<- p_mod$p
  return(p)
  
}

# slope

p_mod_slp <-function(data){
  
  len <- dim(data)[1]
  Rmat <- matrix(0.8, nrow = len, ncol = len)
  diag(Rmat) <- 1
  
  p_mod <- fisher(data$fm_diff_slope_p, adjust = "liji", R = Rmat)
  p<- p_mod$p
  return(p)
  
}


# intersect

p_mod_int <-function(data){
  
  len <- dim(data)[1]
  Rmat <- matrix(0.8, nrow = len, ncol = len)
  diag(Rmat) <- 1
  
  p_mod <- fisher(data$fm_diff_int_p, adjust = "liji", R = Rmat)
  p<- p_mod$p
  return(p)
  
}

# test

p_mod_sd(test_dat)
p_mod_int(test_dat)
p_mod_slp(test_dat)

# nesting data into a lot of data sets and apply p_mod function

n_dat2 <- dat2 %>% group_by(parameter_group) %>%  nest()

m_dat2 <- n_dat2  %>% mutate(merged_p_sd = map_dbl(data, p_mod_sd), 
                             merged_p_int = map_dbl(data, p_mod_int),
                             merged_p_slp = map_dbl(data, p_mod_slp)
)

m_dat2 %>% unnest(data) -> dat2

length(which(m_dat2$merged_p_sd <= 0.05))
# 49
length(which(dat1$p_val_sd <= 0.05))
# 112

# (49 + 112) = 161 out of (58 + 168) = 226

#################
# creating merged p or intercepts and slopes
###################

#first just checkin gone with dat1

dim(dat1)

######
# A
######

# 11 out of 226
#(58 + 168) 

# 9 out of 168 traits sig slope diff - scenario A
dat_slopes1 <-dat1 %>%
  filter (fm_diff_slope_p <= 0.05 & fm_diff_int_p > 0.05)

dim(dat_slopes1) 

# 2 out of 58
dat_slopes2 <-m_dat2 %>%
  filter (merged_p_slp <= 0.05 & merged_p_int > 0.05)

dim(dat_slopes2) 

######
# B 
######

# 93 out of 226

# 66 out of 168 traits sig intercept diff  same slope - scenario B
dat_int1<- dat1 %>%
  filter (fm_diff_int_p <= 0.05 & fm_diff_slope_p >0.05)
dim(dat_int1) 

# 27 out of 58
dat_int2 <-m_dat2 %>%
  filter (merged_p_int <= 0.05 & merged_p_slp > 0.05)
dim(dat_int2) 

######
# C
######

# 69 out of 226

# 52 out of 168 sig. intercept and slope diff - scenario C
dat_intSlopes1<-dat1 %>%
  filter (fm_diff_int_p <= 0.05 & fm_diff_slope_p <= 0.05)
dim(dat_intSlopes1) 

# 15 out of 58
dat_intSlopes2 <-m_dat2 %>%
  filter (merged_p_int <= 0.05 & merged_p_slp <= 0.05)
dim(dat_intSlopes2) 

######
# D
######

# 41 out of 226

# 41 no sig. difference between intercept and slope - scenario D
dat_intslopesNS1 <- dat1 %>%
  filter (fm_diff_slope_p >0.05 & fm_diff_int_p > 0.05)
dim(dat_intslopesNS1) 

# 14 out of 58
dat_intslopesNS2 <- m_dat2 %>%
  filter (merged_p_int > 0.05 & merged_p_slp > 0.05)
dim(dat_intslopesNS2) 


##############################
# TODO  - re-creating Figure 2
################################

#rbind the above scenarios into one matrix with identifier letter A,B,C,D
ScenarioA<-dat_slopes %>% add_column(Scen="A")
ScenarioB<-dat_int %>% add_column(Scen="B")
ScenarioC<-dat_intSlopes %>% add_column(Scen="C")
ScenarioD<-dat_intslopesNS %>% add_column(Scen="D")
dat_add<-bind_rows(ScenarioA,ScenarioB,ScenarioC,ScenarioD)

#write.csv(dat_add,here("data/dat_add.csv"))


#for the traits with sig slope diffs (and sig and non-sig intercepts) assess sex bias in slopes

# set colour for males and females

colours <- c("#D55E00", "#009E73") # c("#882255","#E69F00") 
colours2 <-c("#D55E00", "#7D26CD", "#009E73")

#sex bias in slope parameter under scenario A
dat_p1<-dat_slopes%>%
  group_by_at(vars(Category)) %>%
  summarise(malebias = sum(m_slope > f_slope), 
            femalebias = sum(f_slope > m_slope), 
            total= malebias + femalebias, 
            malepercent = malebias*100/total, 
            femalepercent = femalebias*100/total)  


dat_p1<-gather(as.data.frame(dat_p1), 
               key = sex, 
               value = percent, 
               malepercent:femalepercent, 
               factor_key = TRUE)


dat_p1$samplesize<-with(dat_p1, 
                        ifelse(sex == "malepercent", malebias, femalebias) )

# Adding All
dat_p1 %>%  group_by(sex) %>% summarise(malebias = sum(malebias), 
                                        femalebias= sum(femalebias),
                                        total = sum(total),
) -> part

part %>% mutate(Category = "All",
                sex = c("malepercent", "femalepercent"),
                percent = c(100*(malebias[1]/total[1]), 100*(femalebias[1]/total[1])),
                samplesize = c(malebias[1] ,  femalebias[1]))-> part


#select(Category, malebias, femalebias, total, sex, percent, samplesize)
dat_p1 <- bind_rows(dat_p1, part)



p1 <- 
  ggplot(dat_p1) +
  aes(x = Category, y = percent, fill = sex) +
  geom_col() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
  geom_text(data = subset(dat_p1, samplesize != 0), aes(label = samplesize), position = position_stack(vjust = .5), 
            color = "white", size = 3.5) +
  
  scale_fill_manual(values = colours) +
  theme_bw(base_size = 18) +
  theme(strip.text.y = element_text(angle = 270, size = 10, margin = margin(t=15, r=15, b=15, l=15)), 
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(colour = NULL,linetype = "blank", fill = "gray90"),
        text = element_text(size=14),
        panel.spacing = unit(0.5, "lines"),
        panel.border= element_blank(),
        axis.line=element_line(), 
        panel.grid.major.x = element_line(linetype = "solid", colour = "gray95"),
        panel.grid.major.y = element_line(linetype = "solid", color = "gray95"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()  ) +
  coord_flip()  +
  labs(title = "Scenario A - different slopes, \n                      same intercepts")



#sex bias in intercept parameter - scenario B
dat_p2<-dat_int%>%
  group_by_at(vars(Category)) %>%
  summarise(malebias = sum(m_intercept > f_intercept), 
            femalebias = sum(f_intercept > m_intercept), 
            total= malebias + femalebias, 
            malepercent = malebias*100/total, 
            femalepercent = femalebias*100/total)  

dat_p2<-gather(as.data.frame(dat_p2), 
               key = sex, 
               value = percent, 
               malepercent:femalepercent, 
               factor_key = TRUE)

dat_p2$samplesize<-with(dat_p2, 
                        ifelse(sex == "malepercent", malebias, femalebias) )

# addeing All
dat_p2 %>%  group_by(sex) %>% summarise(malebias = sum(malebias), 
                                        femalebias= sum(femalebias),
                                        total = sum(total),
) -> part2

part2 %>% mutate(Category = "All",
                 sex = c("malepercent", "femalepercent"),
                 percent = c(100*(malebias[1]/total[1]), 100*(femalebias[1]/total[1])),
                 samplesize = c(malebias[1] ,  femalebias[1]))-> part2


#select(Category, malebias, femalebias, total, sex, percent, samplesize)
dat_p2 <- bind_rows(dat_p2, part2)


p2 <- 
  ggplot(dat_p2) +
  aes(x = Category, y = percent, fill = sex) +
  geom_col() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
  geom_text(data = subset(dat_p2, samplesize != 0), aes(label = samplesize), position = position_stack(vjust = .5), 
            color = "white", size = 3.5) +
  
  scale_fill_manual(values = colours) + 
  theme_bw(base_size = 18) +
  theme(strip.text.y = element_text(angle = 270, size = 10, margin = margin(t=15, r=15, b=15, l=15)), 
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(colour = NULL,linetype = "blank", fill = "gray90"),
        text = element_text(size=14),
        panel.spacing = unit(0.5, "lines"),
        panel.border= element_blank(),
        axis.line=element_line(), 
        panel.grid.major.x = element_line(linetype = "solid", colour = "gray95"),
        panel.grid.major.y = element_line(linetype = "solid", color = "gray95"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()  ) +
  coord_flip() +
  labs(title = "Scenario B - same slopes, \n              different intercepts")


#sex bias in sig intercept and slope parameter - scenario C
# TODO we need to fix this needs to be resolved
dat_p3<-dat_intSlopes%>%
  group_by_at(vars(Category)) %>%
  summarise(malebias = sum(m_intercept > f_intercept & m_slope > f_slope), 
            mixed = sum(m_intercept > f_intercept & m_slope < f_slope, 
                        m_intercept < f_intercept & m_slope > f_slope),
            femalebias = sum(f_intercept > m_intercept & f_slope > m_slope), 
            total= malebias + mixed + femalebias, 
            malepercent = malebias*100/total, 
            mixedpercent = mixed*100/total,
            femalepercent = femalebias*100/total)  

dat_p3<-gather(as.data.frame(dat_p3), 
               key = sex, 
               value = percent, 
               malepercent:femalepercent, 
               factor_key = TRUE)
dat_p3$samplesize<-with(dat_p3, 
                        ifelse(sex == "malepercent", malebias, ifelse(sex == "mixedpercent", mixed, femalebias)) )


# addeing All
dat_p3 %>%  group_by(sex) %>% summarise(malebias = sum(malebias), 
                                        mixed = sum(mixed),
                                        femalebias= sum(femalebias),
                                        total = sum(total),
) -> part3

part3 %>% mutate(Category = "All",
                 sex = c("malepercent", "mixedpercent", "femalepercent"),
                 percent = c(100*(malebias[1]/total[1]), 100*(mixed[1]/total[1]), 100*(femalebias[1]/total[1])),
                 samplesize = c(malebias[1] , mixed[1], femalebias[1]))-> part3


#select(Category, malebias, femalebias, total, sex, percent, samplesize)
dat_p3 <- bind_rows(dat_p3, part3)

p3 <- 
  ggplot(dat_p3) +
  aes(x = Category, y = percent, fill = sex) +
  geom_col() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
  geom_text(data = subset(dat_p3, samplesize != 0), aes(label = samplesize), position = position_stack(vjust = .5), 
            color = "white", size = 3.5) +
  
  scale_fill_manual(values = colours2) + 
  theme_bw(base_size = 18) +
  theme(strip.text.y = element_text(angle = 270, size = 10, margin = margin(t=15, r=15, b=15, l=15)), 
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(colour = NULL,linetype = "blank", fill = "gray90"),
        text = element_text(size=14),
        panel.spacing = unit(0.5, "lines"),
        panel.border= element_blank(),
        axis.line=element_line(), 
        panel.grid.major.x = element_line(linetype = "solid", colour = "gray95"),
        panel.grid.major.y = element_line(linetype = "solid", color = "gray95"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        legend.position = "none",
        axis.title.y = element_blank())+
  #axis.title.x = element_blank()  ) +
  ylab("Percentage (%)") +
  coord_flip() +
  labs(title = "Scenario C - different slopes, \n                  different intercepts") 

#sex bias in sd 

dat_p4<-dat%>% filter(p_val_sd <= 0.05) %>% 
  group_by_at(vars(Category)) %>%
  summarise(malebias = sum(m_sd > f_sd), 
            femalebias = sum(f_sd > m_sd), 
            total= malebias + femalebias, 
            malepercent = malebias*100/total, 
            femalepercent = femalebias*100/total)


dat_p4<-gather(as.data.frame(dat_p4), 
               key = sex, 
               value = percent, 
               malepercent:femalepercent, 
               factor_key = TRUE)

dat_p4$samplesize<-with(dat_p4, 
                        ifelse(sex == "malepercent", malebias, femalebias) )


# addeing All
dat_p4 %>%  group_by(sex) %>% summarise(malebias = sum(malebias), 
                                        femalebias= sum(femalebias),
                                        total = sum(total),
) -> part4

part4 %>% mutate(Category = "All",
                 sex = c("malepercent", "femalepercent"),
                 percent = c(100*(malebias[1]/total[1]), 100*(femalebias[1]/total[1])),
                 samplesize = c(malebias[1] ,  femalebias[1]))-> part4


#select(Category, malebias, femalebias, total, sex, percent, samplesize)
dat_p4 <- bind_rows(dat_p4, part4)

p4 <- 
  ggplot(dat_p4) +
  aes(x = Category, y = percent, fill = sex) +
  geom_col() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
  geom_text(data = subset(dat_p4, samplesize != 0), aes(label = samplesize), position = position_stack(vjust = .5), 
            color = "white", size = 3.5) +
  
  scale_fill_manual(values = colours) + 
  theme_bw(base_size = 18) +
  theme(strip.text.y = element_text(angle = 270, size = 10, margin = margin(t=15, r=15, b=15, l=15)), 
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(colour = NULL,linetype = "blank", fill = "gray90"),
        text = element_text(size=14),
        panel.spacing = unit(0.5, "lines"),
        panel.border= element_blank(),
        axis.line=element_line(), 
        panel.grid.major.x = element_line(linetype = "solid", colour = "gray95"),
        panel.grid.major.y = element_line(linetype = "solid", color = "gray95"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        legend.position = "none",
        axis.title.y = element_blank())+
  #axis.title.x = element_blank()  ) +
  ylab("Percentage (%)") +
  coord_flip() +
  labs(title = "Statistically significant \n sex difference in residual SDs") 

(p1 + p2) / (p3 + p4) +   plot_annotation(tag_levels = 'A')


