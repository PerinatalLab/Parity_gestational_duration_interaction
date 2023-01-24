library(dplyr)
library(data.table)
library(ggplot2)
library(lme4)
library(lmerTest)

dat = fread(snakemake@input[[1]])
#dat = fread("/mnt/hdd/common/karin/mfr_150202_recored_filtered_p1_all_variables.csv")




#### Step 1 - making variables as factors and remove rows with NAs ####
# Variables as factors or numeric
dat$max_grade_mor_c = factor(dat$max_grade_mor_c, levels = c(2,1,3,0))
dat$max_grade_far_c = factor(dat$max_grade_far_c, levels = c(2,1,3,0))

dat$Parity_logreg = factor(dat$Parity_logreg, levels = c(2,1,3,4))
dat$mor_birth_country_NORDIC = factor(dat$mor_birth_country_NORDIC)

dat$MALDER2 = dat$MALDER*dat$MALDER
dat$AR2 = dat$AR*dat$AR

dat$PTB_first_born = factor(dat$PTB_first_born)
dat$prev_PTD = factor(dat$prev_PTD)

dat = dat %>% filter(parity_clean<=4)   # need it by parity and not group after parity 3

#Removing NAs
dat1 = dat %>% filter(max_grade_mor_c != 0,max_grade_far_c != 0, !is.na(Parity_logreg), !is.na(KON), !is.na(AR), !is.na(mor_birth_country_NORDIC), !is.na(MALDER))

# Adding the functions in these two scripts to the enviroment
source(snakemake@params[[1]]) # "1_cleaning_modules.R"
#source("/home/karin/Parity_Project_gd/scripts/functions/1_cleaning_modules.R")




#### Step 2 - Only including spontaneous births ####
year_matrix = NULL
dat_m2 = fun_spont1990(dat1)

dat_m2 = select(dat_m2, -c(prev_PTD, PTB_first_born))

## Only Obstetrical history that is spontaneous
# First pregnancy ending preterm
dat_m2 = dat_m2 %>% group_by(lpnr_mor) %>% arrange(parity_clean) %>%
  mutate(PTB_first_born = any(row_number() == 1 & PTB ==1)*1) %>%
  mutate(PTB_first_born = ifelse(dplyr::first(parity_clean)!=1,NA,PTB_first_born))

# Previous pregnancy delivered preterm
dat_m2 = dat_m2 %>% group_by(lpnr_mor) %>% arrange(parity_clean )%>% mutate(prev_PTD = ifelse(dplyr::lag(PTB)==1,1,0))
dat_m2 = dat_m2 %>% group_by(lpnr_mor) %>% arrange(parity_clean) %>%  mutate(diff_p = parity_clean - dplyr::lag(parity_clean)) %>% mutate(prev_PTD = ifelse(diff_p == 1,prev_PTD,NA))




#### Step 3 - looking at the interactions + not stratify by parity ####
## Whole population
dat_m21 = filter(dat_m2, parity_clean>1)
dat_m21 = group_by(dat_m21, lpnr_mor) %>% filter(n()>1) 

beta_CI = function(Model, CI_min, CI_max) {
  m =character()
  for (i in 2) {
          mm = summary(Model)$coefficients[i,1]
          err = confint(Model, level=0.95, method = "Wald")[i,]
          m = rbind(m,c(err[1],mm,err[2]))
  }
  rownames(m) = c("NotByParity")
  colnames(m) = c("CI_min","Beta","CI_max")
  return(m)
}

# First born
m1 = lm(GRDBS ~ PTB_first_born + Parity_logreg + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far), data = dat_m21)
betam1=beta_CI(m1, 0.025,0.975)
print(summary(m1))

# prev PTD
dat_m22 = dat_m21 %>% group_by(lpnr_mor) %>% filter(!is.na(prev_PTD)) %>% filter(n()>1)
m2 = lm(GRDBS ~ prev_PTD + Parity_logreg + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far), data = dat_m22)
betam2=beta_CI(m2, 0.025,0.975)
print(summary(m2))


## Interactions
beta_CI = function(Model, CI_min, CI_max) {
  m =character()
  for (i in 13:14) {
          mm = summary(Model)$coefficients[i,1]
          err = confint(Model, level=0.95, method = "Wald")[i,]
          m = rbind(m,c(err[1],mm,err[2]))
  }
  rownames(m) = c("Parity_logreg3:obshist","Parity_logreg4:obhist")
  colnames(m) = c("CI_min","Beta","CI_max")

  return(m)
}

# Firstborn 
mI1 = lm(GRDBS ~ Parity_logreg*PTB_first_born + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far), data = dat_m21)
betamI1=beta_CI(mI1, 0.025,0.975)
print(summary(mI1))

#prev ptd 
mI2 = lm(GRDBS ~ Parity_logreg*prev_PTD + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far), data = dat_m22)
betamI2=beta_CI(mI2, 0.025,0.975)
print(summary(mI2))




#### Step 4 - splitting models based on parity ####
# Firstborn preterm
beta_CI = function(Model, CI_min, CI_max) {
  m =character()
  for (i in 2) {
          mm = summary(Model)$coefficients[i,1]
          err = confint(Model, level=0.95, method = "Wald")[i,]
          m = rbind(m,c(err[1],mm,err[2]))
  }
  rownames(m) = c("FirstbornPTD")
  colnames(m) = c("CI_min","Beta","CI_max")
  return(m)
}

mp2 = lm(GRDBS ~ PTB_first_born + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far),data = filter(dat_m21, parity_clean ==2))
betamp2=beta_CI(mp2, 0.025,0.975)

mp3 = lm(GRDBS ~ PTB_first_born + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far), data = filter(dat_m21, parity_clean ==3))
betamp3=beta_CI(mp3, 0.025,0.975)

mp4 = lm(GRDBS ~ PTB_first_born + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far), data = filter(dat_m21, parity_clean >=4) )
betamp4=beta_CI(mp4, 0.025,0.975)


# Previous pregnancy preterm
beta_CI = function(Model, CI_min, CI_max) {
  m =character()
  for (i in 2) {
          mm = summary(Model)$coefficients[i,1]
          err = confint(Model, level=0.95, method = "Wald")[i,]
          m = rbind(m,c(err[1],mm,err[2]))
  }
  rownames(m) = c("PrevPTD")
  colnames(m) = c("CI_min","Beta","CI_max")
  return(m)
}

m2p2 = lm(GRDBS ~  prev_PTD + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far), data = filter(dat_m22, parity_clean ==2) )
betam2p2=beta_CI(m2p2, 0.025,0.975)

m2p3 = lm(GRDBS ~  prev_PTD +  KON + AR + AR2 + mor_birth_country_NORDIC + MALDER_c + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far), data = filter(dat_m22, parity_clean ==3))
betam2p3=beta_CI(m2p3, 0.025,0.975)

m2p4 = lm(GRDBS ~  prev_PTD +  KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + FALDER + as.numeric(max_grade_mor) + as.numeric(max_grade_far),data = filter(dat_m22, parity_clean >=4))
betam2p4=beta_CI(m2p4, 0.025,0.975)




#### Steps 5 - ploting data ####
plotdata <- data.frame(Parity = c("Parity 1","Parity 2","Parity 3","Parity 1","Parity 2","Parity 3"), 
                       Odds = c(betamp2[,2],betamp3[,2],betamp4[,2],betam2p2[,2],betam2p3[,2],betam2p4[,2]),
                       CILow = c(betamp2[,1],betamp3[,1],betamp4[,1],betam2p2[,1],betam2p3[,1],betam2p4[,1]),
                       CIHigh = c(betamp2[,3],betamp3[,3],betamp4[,3],betam2p2[,3],betam2p3[,3],betam2p4[,3]),
                       Obstetrical_history = c("PTD first born","PTD first born","PTD first born","Prev PTD","Prev PTD","Prev PTD") ) # Parity 1 == second preg, Parity 2 === third preg. etc

plotdata$Parity <- factor(plotdata$Parity, levels = c("Parity 1","Parity 2","Parity 3"))
plotdata$Odds = as.numeric(plotdata$Odds)
plotdata$CILow = as.numeric(plotdata$CILow)
plotdata$CIHigh = as.numeric(plotdata$CIHigh)

text_dat <- data.frame(Parity = c('Parity 1','Parity 2','Parity 3'), lbl = c('Parity one', 'Parity two', 'Parity three'))
p = plotdata %>%
  ggplot(aes(x = Odds, y = Obstetrical_history)) +
  geom_errorbarh(aes(xmin = CILow, xmax = CIHigh),size = .6, height = .2, color = "black") +
  geom_point(aes(shape = Obstetrical_history,fill= Obstetrical_history), size = 3, stroke=0.4) +
  geom_vline(aes(xintercept = 0), linetype = 2, color="black", size=0.6) +
  coord_cartesian(xlim = c(-14, 0.3), clip="off") +
  facet_wrap(Parity~., scales ="free", ncol=1) +
  ylab("") +
  xlab("Beta (Days)") +
  theme_bw() +
  theme(panel.spacing.y = unit(10, "points"),
        axis.text.y = element_blank(),
        axis.ticks.length.y = unit(0, "points"),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.line = element_line(),
        axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 11, color="black"),
        plot.title = element_text(size =11),
        legend.text = element_text(size=10.5),
        legend.title = element_text(size=11),
        legend.position = "bottom",
        strip.text.x = element_blank(),
        text=element_text(family="Helvetica"),
        plot.margin = margin(r = 1, l = 1, t=0.05,unit = "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color ="black", fill=NA, size=0.5)
  ) +
  scale_shape_manual(labels=c("Preterm delivery firstborn","Previous delivery preterm"),values=c("PTD first born" = 24, "Prev PTD" = 23), guide = guide_legend(title.position="top",title.hjust =0.5, title="Clinical history of preterm delivery")) +
  scale_fill_manual(values=c("PTD first born" ="#756bb1","Prev PTD"= "#bcbddc"),labels=c("Preterm delivery firstborn","Previous delivery preterm"),guide = guide_legend(title.position="top",title.hjust =0.5, title = "Clinical history of preterm delivery") )  +
  geom_text(data = text_dat, aes(x = -15, y = c(0.65,0.65,0.45), label = lbl), hjust = 0, colour = "black", angle=90,size=10/.pt) +
  scale_x_continuous(breaks = seq(-14, 0.5, by = 2))

 


#ggsave("/home/karin/Parity_Project_gd/scripts/plots/Supplementary/output/obstetric_history_only_spont.png",p, width = 174, height = 105, dpi = 1200, units = "mm", device='png')
ggsave(snakemake@output[[1]],p, width = 174, height = 105, dpi = 1200, units = "mm", device='png')

model_info_final = rbind(c("first_born", nobs(mp2), nobs(mp3), nobs(mp4)),
			 c("prev_ptd",  nobs(m2p2), nobs(m2p3), nobs(m2p4)),
			 c("whole_pop_firstborn",nobs(m1),"-","-"),
		         c("whole_pop_prev_ptd", nobs(m2),"-","-"),
		         c("interaction_firstborn", nobs(mI1),"-","-"),
		         c("interaction_prev_ptd", nobs(mI2),"-","-") )


#fwrite(as.data.frame(model_info_final),"/home/karin/Parity_Project_gd/scripts/plots/Supplementary/output/clinical_history_model_info_only_spont.csv",sep=",")
fwrite(as.data.frame(model_info_final),snakemake@output[[2]],sep=",")

#fwrite(as.data.frame(plotdata),"/home/karin/Parity_Project_gd/scripts/plots/Supplementary/output/clinical_history_plotdata_only_spont.csv", sep=",")
fwrite(as.data.frame(plotdata),snakemake@output[[3]], sep=",")
