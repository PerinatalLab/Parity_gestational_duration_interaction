library(dplyr)
library(data.table)
library(ggplot2)
library(lme4)
library(lmerTest)

dat = fread(snakemake@input[[1]])
#dat = fread("/mnt/hdd/common/karin/mfr_150202_recored_filtered_p1_all_variables.csv")




#### Variables as factors ####
dat$max_grade_mor_c = factor(dat$max_grade_mor_c, levels = c(2,1,3,0))
dat$max_grade_far_c = factor(dat$max_grade_far_c, levels = c(2,1,3,0))

dat$Parity_logreg = factor(dat$Parity_logreg, levels = c(2,1,3,4))
dat$mor_birth_country_NORDIC = factor(dat$mor_birth_country_NORDIC)

dat$MALDER2 = dat$MALDER*dat$MALDER
dat$AR2 = dat$AR*dat$AR

dat$mother_herself_PTB = factor(dat$mother_herself_PTB)




#### Finding which parity the mother was born as ####
barn = dat %>% pull(lpnr_BARN)
mor_also_barn_in_mfr =  dat[dat$lpnr_mor %in% barn,]
mor_also_barn_in_mfr = mor_also_barn_in_mfr %>% select(lpnr_mor)
mor_also_barn_in_mfr = unique(mor_also_barn_in_mfr)
mor_as_barn = inner_join(dat,mor_also_barn_in_mfr, by = c("lpnr_BARN" ="lpnr_mor")) # The pregnancies in which the mothers where born

mor_as_barn = mor_as_barn %>% select(lpnr_BARN,parity_clean) #lpnr_BARN here are barn that also are mothers in mfr
colnames(mor_as_barn) = c("lpnr_BARN","parity_clean_mor")

dat = full_join(dat, mor_as_barn, by = c("lpnr_mor" ="lpnr_BARN"))

dat$parity_clean_mor = factor(dat$parity_clean_mor) 

dat = dat %>% mutate(mother_herself_PTB_c = ifelse(mother_herself_PTB == 1 & parity_clean_mor == 1,1,0)) %>% mutate(mother_herself_PTB_c = ifelse(mother_herself_PTB ==1 & parity_clean_mor != 1,2,mother_herself_PTB_c))
dat$mother_herself_PTB_c = factor(dat$mother_herself_PTB_c)

dat = dat %>% mutate(parity_mor = ifelse(parity_clean_mor == 1,1,0))




#### Removing missings ####
dat1 = dat %>% filter(max_grade_mor_c != 0,max_grade_far_c != 0, !is.na(Parity_logreg), !is.na(KON), !is.na(AR), !is.na(mor_birth_country_NORDIC), !is.na(MALDER))




#### only spontaneous births ####
source(snakemake@params[[1]]) # "1_cleaning_modules.R"
#source("~/Parity_Project1/scripts/functions/1_cleaning_modules.R")

year_matrix = NULL
dat_m2 = fun_spont1990(dat1)




#### Models ####
## Whole population 
dat_m21 = group_by(dat_m2, lpnr_mor) %>% filter(n()>1)

beta_CI = function(Model, CI_min, CI_max) {
  m =character()
  for (i in 2) {
          mm = summary(Model)$coefficients[i,1]
          err = confint(Model, level=0.95, method = "Wald")[i+2,]
          m = rbind(m,c(err[1],mm,err[2]))
  }
  rownames(m) = c("NotByParity")
  colnames(m) = c("CI_min","Beta","CI_max")
  return(m)
}


# mother herself ptb
m4 = lmer(GRDBS ~ parity_mor + mother_herself_PTB + Parity_logreg + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m21,control = lmerControl(optimizer ="bobyqa"))
betam4=beta_CI(m4, 0.025,0.975)
print(summary(m4))

# mother preterm interactions
mI4 = lmer(GRDBS ~ mother_herself_PTB*parity_mor + Parity_logreg + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m21,control = lmerControl(optimizer ="bobyqa"))
print(summary(mI4))

mI5 = lmer(GRDBS ~ mother_herself_PTB*parity_mor*Parity_logreg + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m21,control = lmerControl(optimizer ="bobyqa"))
print(summary(mI5))


## Models by parity based on interaction in mI5
beta_CI = function(Model, CI_min, CI_max) {
  m =character()
  for (i in 2) {
          mm = summary(Model)$coefficients[i,1]
          err = confint(Model, level=0.95, method = "Wald")[i,]
          m = rbind(m,c(err[1],mm,err[2]))
  }
  rownames(m) = c("motherbornpreterm")
  colnames(m) = c("CI_min","Beta","CI_max")
  return(m)
}

# Mother born as parity 0 
m2p1 = lm(GRDBS ~ mother_herself_PTB + KON + AR + AR2  + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far),data = filter(dat_m21, parity_clean ==1 & parity_clean_mor == 1) )
betam2p1=beta_CI(m2p1, 0.025,0.975)

m2p2 = lm(GRDBS ~ mother_herself_PTB + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER_c + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far),data = filter(dat_m21, parity_clean ==2 & parity_clean_mor == 1) )
betam2p2=beta_CI(m2p2, 0.025,0.975)

m2p3 = lm(GRDBS ~ mother_herself_PTB + KON + AR + AR2 + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far), data = filter(dat_m21, parity_clean ==3 & parity_clean_mor ==1) )
betam2p3=beta_CI(m2p3, 0.025,0.975)

m2p4 = lm(GRDBS ~ mother_herself_PTB + KON + AR + AR2 + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far_c), data = filter(dat_m21, parity_clean >=4 & parity_clean_mor ==1) )
betam2p4=beta_CI(m2p4, 0.025,0.975)

# Mother as parity > 0
m3p1 = lm(GRDBS ~ mother_herself_PTB + KON + AR + AR2 +  MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far),data = filter(dat_m21, parity_clean ==1 & as.numeric(parity_clean_mor)>1) )
betam3p1=beta_CI(m3p1, 0.025,0.975)

m3p2 = lm(GRDBS ~ mother_herself_PTB + KON + AR + AR2 + MALDER + MALDER2 + as.numeric(max_grade_mor)+ as.numeric(max_grade_far),data = filter(dat_m21, parity_clean ==2 & as.numeric(parity_clean_mor)>1) )
betam3p2=beta_CI(m3p2, 0.025,0.975)

m3p3 = lm(GRDBS ~ mother_herself_PTB + KON + AR + AR2 + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far), data = filter(dat_m21, parity_clean ==3 & as.numeric(parity_clean_mor)>1) )
betam3p3=beta_CI(m3p3, 0.025,0.975)

m3p4 = lm(GRDBS ~ mother_herself_PTB + KON + AR + AR2 + MALDER + MALDER2 + as.numeric(max_grade_mor)+ as.numeric(max_grade_far), data = filter(dat_m21, parity_clean >=4 & as.numeric(parity_clean_mor)>1) )
betam3p4=beta_CI(m3p4, 0.025,0.975)




#### Plotting the data ####
plotdata <- data.frame(Parity = c("Parity zero","Parity one","Parity two","Parity \u2265 three","Parity zero","Parity one","Parity two","Parity \u2265 three"), 
                       Beta = c(betam2p1[,2],betam2p2[,2],betam2p3[,2],betam2p4[,2],betam3p1[,2],betam3p2[,2],betam3p3[,2],betam3p4[,2] ),
                       CILow = c(betam2p1[,1],betam2p2[,1],betam2p3[,1],betam2p4[,1],betam3p1[,1],betam3p2[,1],betam3p3[,1],betam3p4[,1] ),
                       CIHigh = c(betam2p1[,3],betam2p2[,3],betam2p3[,3],betam2p4[,3],betam3p1[,3],betam3p2[,3],betam3p3[,3],betam3p4[,3]),
                       Model = c("Mother herself preterm, Parity = zero","Mother herself preterm, Parity = zero","Mother herself preterm, Parity = zero","Mother herself preterm, Parity = zero","Mother herself preterm, Parity > zero","Mother herself preterm, Parity > zero","Mother herself preterm, Parity > zero","Mother herself preterm, Parity > zero"))

plotdata$Parity <- factor(plotdata$Parity, levels = c("Parity zero","Parity one","Parity two","Parity \u2265 three"))
plotdata$Beta = as.numeric(plotdata$Beta)
plotdata$CILow = as.numeric(plotdata$CILow)
plotdata$CIHigh = as.numeric(plotdata$CIHigh)

p = plotdata %>%
  ggplot(aes(x = Beta, y = Model)) +
  geom_errorbarh(aes(xmin = CILow, xmax = CIHigh),size = .6, height = 0.17, color = "black") +
  geom_point(aes(shape=Model, fill = Model), size = 3, stroke=0.4) +
  geom_vline(aes(xintercept = 0), linetype = 2, color ="black", size = 0.6) +
  coord_cartesian(xlim = c(-4.5, 0.3), clip ="off") +
  facet_grid(Parity~., switch = "y") +
  ylab("") +
  xlab("Beta (days)") +
  theme_bw() +
  theme(panel.spacing.y = unit(10, "points"),
        axis.text.y = element_blank(),
        axis.ticks.length.y = unit(0, "points"),
        strip.background.y = element_blank(),
        strip.placement = "outside",
        axis.line = element_line(),
        axis.title.x = element_text(size = 11.5),
        axis.text.x = element_text(size = 11.5, color="black"),
        strip.text = element_text(size = 11.5),
        text=element_text(family="Helvetica"),
        legend.title = element_text(size=11),
        legend.text = element_text(size=10),
        legend.position = "bottom",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_rect(colour ="black", fill=NA, size=0.50)

  ) +
  scale_shape_manual(labels=c("Maternal herself preterm, Parity = zero","Mother herself preterm, Parity > zero"),values=c("Mother herself preterm, Parity = zero"=23,"Mother herself preterm, Parity > zero"=21), guide = guide_legend(title.position="top",title.hjust =0.5, title="Family history")) +
  scale_fill_manual(values=c("Mother herself preterm, Parity = zero"= "#fec44f","Mother herself preterm, Parity > zero"="#fff7bc"),labels=c("Maternal herself preterm, Parity = zero","Mother herself preterm, Parity > zero"),guide = guide_legend(title.position="top",title.hjust =0.5, title = "Family history") ) 




#### Infromation from models in this script ####
model_info_final = rbind(c("mother_herself_PTB", nobs(mp1), nobs(mp2), nobs(mp3), nobs(mp4),"-" ),
			 c("mother_her_self_PTB_P=0",  nobs(m2p1), nobs(m2p2), nobs(m2p3), nobs(m2p4),"-" ), 
			 c("mother_her_self_PTB_P>0",nobs(m3p1), nobs(m3p2), nobs(m3p3), nobs(m3p4),"-"),
			 c("whole_pop_mother_PTB_incl_morparity", nobs(m4),"-","-","-" ,m4@Gp[2]),
			 c("interaction_mother_PTB_mor_parity*mor_ptb", nobs(mI4),"-","-","-",mI4@Gp[2]),
			 c("interaction_mother_PTB_mor_parity*mor_ptb*parity", nobs(mI5),"-","-","-",mI5@Gp[2]))




#### Saving ####
#ggsave("~/Parity_Project_gd/scripts/plots/Supplementary/output/mother_born_preterm_parity.png",p, width = 19.5, height = 14, units ="cm")
ggsave(snakemake@output[[1]],p,width = 17.4, height = 13, units ="cm")

#fwrite(model_info_final,"~/Parity_Project_gd/scripts/plots/Supplementary/output/model_info_mother_born_preterm_parity.csv",sep=",")
fwrite(as.data.frame(model_info_final), snakemake@output[[2]], sep=",")

#fwrite(plotdata, "~/Parity_Project_gd/scripts/plots/Supplementary/output/plotdata_mother_born_pmother_born_preterm_parity.csv", sep=",")
fwrite(as.data.frame(plotdata), snakemake@output[[3]], sep=",")

