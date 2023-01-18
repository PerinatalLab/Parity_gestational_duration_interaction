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

dat$MALDER_c = factor(dat$MALDER_c)




#### Removing missings ####
dat1 = dat %>% filter(max_grade_mor_c != 0,max_grade_far_c != 0, !is.na(Parity_logreg), !is.na(KON), !is.na(AR), !is.na(mor_birth_country_NORDIC), !is.na(MALDER))




#### CI function for the regression models ####
beta_CI = function(Model, CI_min, CI_max) {
  m =character()
  for (i in 2:4) {
          mm = summary(Model)$coefficients[i,1]
          err = confint(Model, level=0.95, method = "Wald")[i+2,]
          m = rbind(m,c(err[1],mm,err[2]))
  }
  rownames(m) = c("P0", "P2", ">=P3")
  P1 = c(0,0,0) #ref
  m = rbind(m,P1)
  colnames(m) = c("CI_min","Beta","CI_max")
  m = m[c("P0","P1", "P2", ">=P3"),]

  return(m)
}

# OBS need to have the parity model as the first variable in all the models for the fucntion to work.




#### Only spontaneous births #### 
source(snakemake@params[[1]]) # "1_cleaning_modules.R"
#source("~/Parity_Project1/scripts/functions/1_cleaning_modules.R") 
 
year_matrix = NULL
dat_m2 = fun_spont1990(dat1)




#### At least two preg per mother ####
dat_m2 = group_by(dat_m2, lpnr_mor) %>% filter(n()>1)




#### Models ####
## Interaction
mI = lmer(GRDBS ~ Parity_logreg*MALDER + Parity_logreg*MALDER2 + KON + AR + AR2 + mor_birth_country_NORDIC + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m2,control = lmerControl(optimizer ="bobyqa"))
print(summary(mI))

## Stratifying by maternal age
dat_m21 = dat_m2 %>% filter(MALDER_c == 1) %>% group_by(lpnr_mor)  %>% filter(n()>1)
am6 = lmer(GRDBS ~ Parity_logreg + KON + AR + AR2 + mor_birth_country_NORDIC + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m21,control = lmerControl(optimizer ="bobyqa"))
betaam6=beta_CI(am6, 0.025,0.975)

dat_m21 = dat_m2 %>% filter(MALDER_c == 0) %>% group_by(lpnr_mor)  %>% filter(n()>1)
am7 = lmer(GRDBS ~ Parity_logreg + KON + AR + AR2  + mor_birth_country_NORDIC + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m21,control = lmerControl(optimizer ="bobyqa"))
betaam7=beta_CI(am7, 0.025,0.975)

dat_m21 = dat_m2 %>% filter(MALDER_c == 2) %>% group_by(lpnr_mor)  %>% filter(n()>1)
am8 = lmer(GRDBS ~ Parity_logreg + KON + AR + AR2 + mor_birth_country_NORDIC + as.numeric(max_grade_mor)+ as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m21,control = lmerControl(optimizer ="bobyqa"))
betaam8=beta_CI(am8, 0.025,0.975)

dat_m21 = dat_m2 %>% filter(MALDER_c == 3) %>% group_by(lpnr_mor)  %>% filter(n()>1)
am9 = lmer(GRDBS ~ Parity_logreg + KON + AR + AR2 + mor_birth_country_NORDIC + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m21,control = lmerControl(optimizer ="bobyqa") )
betaam9=beta_CI(am9, 0.025,0.975)

# Information form models
model_info_final = rbind(c("am6",AIC(am6), nobs(am6),am6@Gp[2]), 
			 c("am7",AIC(am7),nobs(am7),am7@Gp[2]), 
			 c("am8",AIC(am8),nobs(am8),am8@Gp[2]),
			 c("am9",AIC(am9),nobs(am9),am9@Gp[2]) ) 




#### Plotting the data ####
plotdata <- data.frame(Parity = c(rownames(betaam6), rownames(betaam7), rownames(betaam8),rownames(betaam9) ),
                       Beta = c(betaam6[,2],betaam7[,2],betaam8[,2],betaam9[,2]),
                       CILow = c(betaam6[,1],betaam7[,1],betaam8[,1],betaam9[,1]),
                       CIHigh = c(betaam6[,3],betaam7[,3],betaam8[,3],betaam9[,3]),
                       Model = c(""),
		       Maternal_age = c("Age<20","Age<20","Age<20","Age<20","20\u2264Age<30","20\u2264Age<30","20\u2264Age<30","20\u2264Age<30","30\u2264Age<40","30\u2264Age<40","30\u2264Age<40","30\u2264Age<40","Age\u226540","Age\u226540","Age\u226540","Age\u226540")

)

plotdata$Parity <- factor(plotdata$Parity, levels = c(">=P3","P2","P1","P0"))
plotdata$Maternal_age <- factor(plotdata$Maternal_age, levels = c("Age<20","20\u2264Age<30","30\u2264Age<40","Age\u226540"))
plotdata$Beta = as.numeric(plotdata$Beta)
plotdata$CILow = as.numeric(plotdata$CILow)
plotdata$CIHigh = as.numeric(plotdata$CIHigh)

plotdata = plotdata %>%
        mutate(x_min = ifelse(Maternal_age =="<20",-18.5,-1.3 ) ) %>%
        mutate(x_max = ifelse(Maternal_age == "<20",3,0.6))

p = plotdata %>%
  ggplot(aes(x = Beta, y = Parity)) +
  geom_errorbarh(aes(xmin = CILow, xmax = CIHigh),size = .6, height = .2, color = "black") +
  geom_vline(aes(xintercept = 0), linetype = 2, colour="black", size=0.6) +
  geom_point(aes(shape = Parity,fill=Parity), size = 3,stroke = 0.4) +
  facet_grid(Model~Maternal_age, switch = "y") +
  ylab("") +
  xlab("Beta (days)") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.length.y = unit(0, "points"),
        strip.background = element_blank(),
        strip.placement = "outside",
        text=element_text(family="Helvetica", size =13, color="black"),
        axis.title.x = element_text(size=11, color="black"),
        legend.title = element_text(size=11),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
     scale_shape_manual(labels=c("\u2265 three","two","one","zero"),values=c(22, 21, 24,23), guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5)) +
  scale_fill_manual(values=c(">=P3" ="#E69F00","P2"= "#56B4E9","P1"="#009E73","P0"="#CC79A7"),labels=c("\u2265 three","two","one","zero"),guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5)) + facet_wrap(Model~Maternal_age,  scales="free_x", ncol = 4) + geom_blank(aes(x = x_min)) + geom_blank(aes(x = x_max))




#### Saving ####
#ggsave("~/Parity_Project_gd/scripts/plots/Supplementary/output/parity_maternal_age_mixed.png",p, width =17.4, height = 8, units="cm")
ggsave(snakemake@output[[1]],p, width =17.4, height = 8, units = "cm")

#fwrite(model_info_final,"~/Parity_Project_gd/scripts/plots/Supplementary/output/model_info_maternal_age_mixed.csv",sep=",")
fwrite(model_info_final,snakemake@output[[2]],sep=",")

#fwrite(plotdata, "~/Parity_Project_gd/scripts/plots/Supplementary/output/plotdata_maternal_age_mixed.csv", sep=",")
fwrite(plotdata, snakemake@output[[3]], sep=",")
