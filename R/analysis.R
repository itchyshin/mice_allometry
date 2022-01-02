# analysis
#TODO - need to get contrasts between male SD and female SD
# TODO - sorting out to the senario!!



# first getting p values - contrast between males and females for 



#assess number of traits with sig shifts in intercept and slope

#16 out of 300 traits sig slope diff - scenario A
Fin_dat_slopes<-Fin_dat2a %>%
  filter (fm_diff_slope_p <0.05 & fm_diff_int_p > 0.05)

length(Fin_dat_slopes) 

#109 out of 300 traits sig intercept diff  same slope - scenario B
Fin_dat_int<- Fin_dat2a %>%
  filter (fm_diff_int_p <0.05 & fm_diff_slope_p >0.05)

#73 out of 300 sig intercept and slope diff - scenario C
Fin_dat_intSlopes<-Fin_dat2a%>%
  filter (fm_diff_int_p <0.05 & fm_diff_slope_p <0.05)

#102 no sig difference between intercept and slope - scenario D
Fin_dat_intslopesNS<-Fin_dat2a %>%
  filter (fm_diff_slope_p >0.05 & fm_diff_int_p > 0.05)

#rbind the above scenarios into one matrix with identifier letter A,B,C,D
ScenarioA<-Fin_dat_slopes %>% add_column(Scen="A")
ScenarioB<-Fin_dat_int %>% add_column(Scen="B")
ScenarioC<-Fin_dat_intSlopes %>% add_column(Scen="C")
ScenarioD<-Fin_dat_intslopesNS %>% add_column(Scen="D")
AtoD<-list(ScenarioA,ScenarioB,ScenarioC,ScenarioD)
AllScenarios_dat<-do.call(rbind,AtoD)

#remove additional body weight variables (body weight ~body weight results to be deleted)
AllScenarios_dat<-AllScenarios_dat[-c(27,28,135),] #("Body weight after experiment", "Body weight before experiment", "Body weight") #remove duplicates of body weight

write.csv(AllScenarios_dat,"AllScenarios_dat.csv")


#for the traits with sig slope diffs (and sig and non-sig intercepts) assess sex bias in slopes

#sex bias in slope parameter under scenario A
meta.plot1<-Fin_dat_slopes%>%
  group_by_at(vars(Category)) %>%
  summarise(malebias = sum(m_slope > f_slope), femalebias = sum(f_slope > m_slope), total= malebias + femalebias, 
            malepercent = malebias*100/total, femalepercent = femalebias*100/total)  
as.data.frame(meta.plot1)
meta.plot1df<-gather(meta.plot1, key = sex, value = percent, malepercent:femalepercent, factor_key = TRUE)
meta.plot1df$samplesize<-with(meta.plot1df, ifelse(sex == "malepercent", malebias, femalebias) )

#sex bias in sd of the slope under scenario A
meta.plot1a<-Fin_dat_slopes%>%
  group_by_at(vars(Category)) %>%
  summarise(malebias = sum(m_sd > f_sd), femalebias = sum(f_sd > m_sd), total= malebias + femalebias, 
            malepercent = malebias*100/total, femalepercent = femalebias*100/total) 
as.data.frame(meta.plot1a)
meta.plot1adf<-gather(meta.plot1a, key = sex, value = percent, malepercent:femalepercent, factor_key = TRUE)
meta.plot1adf$samplesize<-with(meta.plot1adf, ifelse(sex == "malepercent", malebias, femalebias) )

cor.test(meta.plot1adf$femalebias,meta.plot1df$femalebias) #not significant, but the same for malebias is!
cor.test(meta.plot1adf$malebias,meta.plot1df$malebias) #males that have larger slopes have larger variance

Sexbias_slopesOnly <- 
  ggplot(meta.plot1df) +
  aes(x = Category, y = percent, fill = sex) +
  geom_col() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
  geom_text(data = subset(meta.plot1df, samplesize != 0), aes(label = samplesize), position = position_stack(vjust = .5), 
            color = "white", size = 3.5) +
  
  scale_fill_brewer(palette = "Set1") +
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
  coord_flip()

Sexbias_slopesOnly_sd <- 
  ggplot(meta.plot1adf) +
  aes(x = Category, y = percent, fill = sex) +
  geom_col() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
  geom_text(data = subset(meta.plot1adf, samplesize != 0), aes(label = samplesize), position = position_stack(vjust = .5), 
            color = "white", size = 3.5) +
  
  scale_fill_brewer(palette = "Set1") +
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
  coord_flip()

#sex bias in intercept parameter - scenario B
meta.plot2<-Fin_dat_int%>%
  group_by_at(vars(Category)) %>%
  summarise(malebias = sum(m_intercept > f_intercept), femalebias = sum(f_intercept > m_intercept), total= malebias + femalebias, 
            malepercent = malebias*100/total, femalepercent = femalebias*100/total)  
as.data.frame(meta.plot2)
meta.plot2df<-gather(meta.plot2, key = sex, value = percent, malepercent:femalepercent, factor_key = TRUE)
meta.plot2df$samplesize<-with(meta.plot2df, ifelse(sex == "malepercent", malebias, femalebias) )

#sex bias in sd of the slope - scenario B
meta.plot2a<-Fin_dat_int%>%
  group_by_at(vars(Category)) %>%
  summarise(malebias = sum(m_sd > f_sd), femalebias = sum(f_sd > m_sd), total= malebias + femalebias, 
            malepercent = malebias*100/total, femalepercent = femalebias*100/total) 
as.data.frame(meta.plot2a)
meta.plot2adf<-gather(meta.plot2a, key = sex, value = percent, malepercent:femalepercent, factor_key = TRUE)
meta.plot2adf$samplesize<-with(meta.plot2adf, ifelse(sex == "malepercent", malebias, femalebias) )

cor.test(meta.plot2adf$femalebias,meta.plot2df$femalebias) #significant, but the same for malebias is!
cor.test(meta.plot2adf$malebias,meta.plot2df$malebias) #males that have larger slopes have larger variance

Sexbias_IntOnly <- 
  ggplot(meta.plot2df) +
  aes(x = Category, y = percent, fill = sex) +
  geom_col() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
  geom_text(data = subset(meta.plot2df, samplesize != 0), aes(label = samplesize), position = position_stack(vjust = .5), 
            color = "white", size = 3.5) +
  
  scale_fill_brewer(palette = "Set1") +
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
  coord_flip()

Sexbias_IntOnly_sd <- 
  ggplot(meta.plot2adf) +
  aes(x = Category, y = percent, fill = sex) +
  geom_col() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
  geom_text(data = subset(meta.plot2adf, samplesize != 0), aes(label = samplesize), position = position_stack(vjust = .5), 
            color = "white", size = 3.5) +
  
  scale_fill_brewer(palette = "Set1") +
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
  coord_flip()


#sex bias in sig intercept and slope parameter - scenario C
meta.plot3<-Fin_dat_intSlopes%>%
  group_by_at(vars(Category)) %>%
  summarise(malebias = sum(m_intercept > f_intercept,m_slope > f_slope), femalebias = sum(f_intercept > m_intercept, f_slope > m_slope), total= malebias + femalebias, 
            malepercent = malebias*100/total, femalepercent = femalebias*100/total)  
as.data.frame(meta.plot3)
meta.plot3df<-gather(meta.plot3, key = sex, value = percent, malepercent:femalepercent, factor_key = TRUE)
meta.plot3df$samplesize<-with(meta.plot3df, ifelse(sex == "malepercent", malebias, femalebias) )

#sex bias in sd of the slope and intercept group - scenario C
meta.plot3a<-Fin_dat_intSlopes%>%
  group_by_at(vars(Category)) %>%
  summarise(malebias = sum(m_sd > f_sd), femalebias = sum(f_sd > m_sd), total= malebias + femalebias, 
            malepercent = malebias*100/total, femalepercent = femalebias*100/total) 
as.data.frame(meta.plot3a)
meta.plot3adf<-gather(meta.plot3a, key = sex, value = percent, malepercent:femalepercent, factor_key = TRUE)
meta.plot3adf$samplesize<-with(meta.plot3adf, ifelse(sex == "malepercent", malebias, femalebias) )

cor.test(meta.plot3adf$femalebias,meta.plot3df$femalebias) #significant, but the same for malebias is!
cor.test(meta.plot3adf$malebias,meta.plot3df$malebias) #males that have larger slopes have larger variance

Sexbias_IntSlopes <- 
  ggplot(meta.plot3df) +
  aes(x = Category, y = percent, fill = sex) +
  geom_col() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
  geom_text(data = subset(meta.plot3df, samplesize != 0), aes(label = samplesize), position = position_stack(vjust = .5), 
            color = "white", size = 3.5) +
  
  scale_fill_brewer(palette = "Set1") +
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
  coord_flip()

Sexbias_IntSlopes_sd <- 
  ggplot(meta.plot3adf) +
  aes(x = Category, y = percent, fill = sex) +
  geom_col() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
  geom_text(data = subset(meta.plot3adf, samplesize != 0), aes(label = samplesize), position = position_stack(vjust = .5), 
            color = "white", size = 3.5) +
  
  scale_fill_brewer(palette = "Set1") +
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
  coord_flip()