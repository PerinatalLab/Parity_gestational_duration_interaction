library(data.table)
library(dplyr)
library(lme4)
library(lmerTest)


#### Loading data ####
#dat = fread("/mnt/hdd/common/karin/mfr_150202_recored_filtered_p1_all_variables.csv")
#p_id = fread("/mnt/hdd/data/swed/ver/parents.csv") 

dat = fread(snakemake@input[[1]])
p_id = fread(snakemake@input[[2]])


p_id = p_id %>% select(LopnrBarn, LopnrMor, FoddArBioMor, LopnrFar, FoddArBioFar) # selecting columns of interest
colnames(p_id) = c("lpnr_BARN", "lpnr_mor", "ar_mor", "lpnr_far", "ar_far") # changing the col names




#### Creating variables missing ####
## Full siblings children
p_id$parents = paste(p_id$lpnr_mor, p_id$lpnr_far, sep = "_") # create a variable called parents.

p_id = p_id %>% group_by(parents) %>% mutate(ID_full_siblings = ifelse(n()>1,cur_group_id(),NA)) # based on parents connect all full siblings
cat("Number of full sibling groups based on ID_full_siblings:", length(unique(p_id$ID_full_siblings)))


## Full sisters mother
# Thus conecting the mothers siblings (not the childrens sibling). 
# All mothers can not be connect her to her siblings (Do not have that information, mothers need to be children in the registry)

p_id_barn_full_S_ID = p_id %>%  ungroup() %>% select(lpnr_BARN, ID_full_siblings) # extracting the children and their sibling ID for the whole data set
p_id_mor = p_id %>% ungroup() %>% select(lpnr_mor) # extracting the mothers from the whole data set
p_id_mor = cbind(p_id_mor,p_id_mor) # When using join future down the lpnr_mor will disappear and that's why the matrix contain the same information in two columns
colnames(p_id_mor) = c("lpnr_mor","lpnr_Mor")

p_id_join=inner_join(p_id_barn_full_S_ID, p_id_mor,by = c("lpnr_BARN" = "lpnr_Mor")) %>% select(-lpnr_BARN,) # Connecting each row where the ID exist in both the mother ID column  och child ID column to identify mother sibling. Here lpnr_MOR will disapear and lpnr_mor be left in the new data frame.
p_id_join=unique(p_id_join) # to not add rows in our dataset when using left_join in the next step
p_id = left_join(p_id,p_id_join, by="lpnr_mor") # adding the sibling ID

colnames(p_id) = c("lpnr_BARN","lpnr_mor","ar_mor","lpnr_far","ar_far","parents","ID_full_siblings", "ID_mor_full_siblings")


sibling = p_id %>% ungroup() %>% select(lpnr_mor, ID_mor_full_siblings)
sibling = unique(sibling)
dat1 = left_join(dat, sibling, by="lpnr_mor")


## Full maternal sisters giving ptd in any child before "you"
dat3 = dat1 %>% group_by(lpnr_mor) %>% arrange(BFODDAT) %>% 
mutate(cumPTDpermor =  cumsum(PTB)) 

test = dat3 %>% group_by(ID_mor_full_siblings) %>% arrange(BFODDAT) %>% 
mutate(maternal_sisters_children_born_PTB_1 =  cumsum(PTB) -cumPTDpermor,0) %>% 
mutate(maternal_sisters_children_born_PTB_1 = ifelse(min(lpnr_mor) == max(lpnr_mor) | is.na(ID_mor_full_siblings),NA,maternal_sisters_children_born_PTB_1)) %>% ungroup() %>% group_by(lpnr_mor) %>%
mutate(maternal_sisters_children_born_PTB = max(maternal_sisters_children_born_PTB_1)) %>% mutate(maternal_sisters_children_born_PTB = ifelse(maternal_sisters_children_born_PTB>=1,1,0)) %>% ungroup() %>% select(lpnr_BARN,maternal_sisters_children_born_PTB)

dat4 = left_join(dat3,test, by="lpnr_BARN")

dat = dat4




#### Variables as factors ####
dat$max_grade_mor_c = factor(dat$max_grade_mor_c, levels = c(2,1,3,0))
dat$max_grade_far_c = factor(dat$max_grade_far_c, levels = c(2,1,3,0))

dat$Parity_logreg = factor(dat$Parity_logreg, levels = c(2,1,3,4))
dat$mor_birth_country_NORDIC = factor(dat$mor_birth_country_NORDIC)

dat$MALDER2 = dat$MALDER*dat$MALDER
dat$AR2 = dat$AR*dat$AR

dat$mother_herself_PTB = factor(dat$mother_herself_PTB)
dat$maternal_sisters_children_born_PTB = factor(dat$maternal_sisters_children_born_PTB)




#### Removing NAs ####
dat1 = dat %>% ungroup()%>% filter(max_grade_mor_c != 0,max_grade_far_c != 0, !is.na(Parity_logreg), !is.na(KON), !is.na(AR), !is.na(mor_birth_country_NORDIC), !is.na(MALDER))




#### Extracting the spontaneous birth ####
#source("/home/karin/Parity_Project1/scripts/functions/1_cleaning_modules.R")
source(snakemake@params[[1]])
year_matrix = NULL
dat_m2 = fun_spont1990(dat1)




#### Running models ####
## Whole population + interaction
dat_m21 = group_by(dat_m2, lpnr_mor) %>% filter(n()>1)

## maternal sibling ptb
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


m4 = lmer(GRDBS ~ maternal_sisters_children_born_PTB + Parity_logreg + KON + AR + AR2 + MALDER2  + as.numeric(max_grade_mor) + as.numeric(max_grade_far)+ (1|lpnr_mor), data = dat_m21, control = lmerControl(optimizer ="bobyqa"))
betam4 = beta_CI(m4,0.025,0.975)
print(summary(m4))
print(betam4)


mI4 = lmer(GRDBS ~ Parity_logreg*maternal_sisters_children_born_PTB + KON + AR + AR2 + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m21, control = lmerControl(optimizer ="bobyqa"))
print(summary(mI4))


## By parity
beta_CI = function(Model, CI_min, CI_max) {
  m =character()
  for (i in 2) {
          mm = summary(Model)$coefficients[i,1]
          err = confint(Model, level=0.95, method = "Wald")[i,]
          m = rbind(m,c(err[1],mm,err[2]))
  }
  rownames(m) = c("maternal_sisters_children_born_PTB")
  colnames(m) = c("CI_min","Beta","CI_max")
  return(m)
}


m2p1 = lm(GRDBS ~ maternal_sisters_children_born_PTB + KON + AR + AR2 + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far),data = filter(dat_m21, parity_clean ==1) )
betamp1=beta_CI(m2p1, 0.025,0.975)

m2p2 = lm(GRDBS ~ maternal_sisters_children_born_PTB + KON + AR + AR2 + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far),data = filter(dat_m21, parity_clean ==2) )
betamp2=beta_CI(m2p2, 0.025,0.975)

m2p3 = lm(GRDBS ~ maternal_sisters_children_born_PTB + KON + AR + AR2 + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far),data = filter(dat_m21, parity_clean ==3) )
betamp3=beta_CI(m2p3, 0.025,0.975)

m2p4 = lm(GRDBS ~ maternal_sisters_children_born_PTB + KON + AR + AR2 + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far),data = filter(dat_m21, parity_clean >= 4) )
betamp4=beta_CI(m2p4, 0.025,0.975)



model_info_final = rbind(c("Model","AIC","Number of observations", "Number of groups"),
                         c("m4",AIC(m4), nobs(m4),m4@Gp[2] ),
                         c("mI4",AIC(mI4),nobs(m4),m4@Gp[2]),
                         c("m2p1", AIC(m2p1),nobs(m2p1),"-"),
                         c("m2p2",AIC(m2p2), nobs(m2p2),"-"),
                         c("m2p3",AIC(m2p3),nobs(m2p3),"-"),
                         c("m2p4",AIC(m2p4),nobs(m2p4),"-")
                         )

#### Plotting results ####

plotdata <- data.frame(Parity = c("Parity 0","Parity 1","Parity 2","Parity \u2265 3"),
                       Beta = c(betamp1[,2],betamp2[,2],betamp3[,2],betamp4[,2]),
                       CILow = c(betamp1[,1],betamp2[,1],betamp3[,1],betamp4[,1]),
                       CIHigh = c(betamp1[,3],betamp2[,3],betamp3[,3],betamp4[,3]),
                       Model = c("Maternal sister given birth preterm","Maternal sister given birth preterm","Maternal sister given birth preterm","Maternal sister given birth preterm"))

plotdata$Parity <- factor(plotdata$Parity, levels = c("Parity 0", "Parity 1","Parity 2","Parity \u2265 3"))
plotdata$Beta = as.numeric(plotdata$Beta)
plotdata$CILow = as.numeric(plotdata$CILow)
plotdata$CIHigh = as.numeric(plotdata$CIHigh)

text_dat <- data.frame(Parity = c('Parity 0','Parity 1','Parity 2','Parity \u2265 3'), lbl = c('Parity zero','Parity one', 'Parity two', 'Parity \u2265 three'))
text_dat$Parity <- factor(text_dat$Parity, levels = c("Parity 0", "Parity 1","Parity 2","Parity \u2265 3"))
p = plotdata %>%
  ggplot(aes(x = Beta, y = Model)) +
  geom_errorbarh(aes(xmin = CILow, xmax = CIHigh),size = .6, height = .2, color = "black") +
  geom_point(aes(shape = Model,fill= Model), size = 3, stroke=0.4) +
  geom_vline(aes(xintercept = 0), linetype = 2, color="black", size=0.6) +
  coord_cartesian(xlim = c(-3.4, 0.75), clip="off") +
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
        axis.title.x = element_text(size = 11, color="black"),
        axis.text.x = element_text(size = 11, color ="black"),
        plot.title = element_text(size =11),
        legend.text = element_text(size=10.5),
        legend.title = element_text(size=11),
        legend.position = "none",
        strip.text.x = element_blank(),
        text=element_text(family="Helvetica"),
        plot.margin = margin(r = 0.5, l = 1, t=0.1,b=0.5,unit = "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color ="black", fill=NA, size=0.5)
  ) +


  scale_shape_manual(values=24) + scale_fill_manual(values="#fe9929") +  geom_text(data = text_dat, aes(x = -3.7, y = c(1, 1,1,1), label = lbl), hjust = 0.5, colour = "black", angle=90,size=10/.pt) +
  scale_x_continuous(breaks = seq(from = -3, to = 0, by = 1)) 


#### Saving ####
#ggsave("/home/karin/Parity_Project_gd/plots/family_history_maternal_sisters_given_birth_preterm.png",p, width=174, height=105, units="mm",dpi=1200)
ggsave(snakemake@output[[1]],p, width=174, height=105, units="mm",dpi=1200,bg = "transparent", device=cairo_ps)

#fwrite(model_info_final,"~/Parity_Project_gd/scripts/plots/output_test/family_history_maternal_sisters_given_birth_preterm_model_info.csv",sep=",")
#fwrite(plotdata, "~/Parity_Project_gd/scripts/plots/output_test/family_history_maternal_sisters_given_birth_preterm_plotdata.csv", sep=",")
fwrite(model_info_final,snakemake@output[[2]],sep=",")
fwrite(plotdata,snakemake@output[[3]], sep=",")
