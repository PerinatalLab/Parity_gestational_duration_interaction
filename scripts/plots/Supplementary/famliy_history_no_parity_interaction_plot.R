library(data.table)
library(dplyr)
library(lme4)
library(lmerTest)


#### Loading data ####
dat = fread(snakemake@input[[1]])
#dat = fread("/mnt/hdd/common/karin/mfr_150202_recored_filtered_p1_all_variables.csv")

p_id = fread(snakemake@input[[2]])
#p_id = fread("/mnt/hdd/data/swed/ver/parents.csv") 

p_id = p_id %>% select(LopnrBarn, LopnrMor, FoddArBioMor, LopnrFar, FoddArBioFar) # selecting columns of interest
colnames(p_id) = c("lpnr_BARN", "lpnr_mor", "ar_mor", "lpnr_far", "ar_far") # changing the col names
#nrow(p_id) = 4138617




#### Creating variables ####
## full siblings (children)
p_id$parents = paste(p_id$lpnr_mor, p_id$lpnr_far, sep = "_") # create a variable called parents.

p_id = p_id %>% group_by(parents) %>% mutate(ID_full_siblings = ifelse(n()>1,cur_group_id(),NA)) # based on parents connect all full siblings
cat("Number of full sibling groups based on ID_full_siblings:", length(unique(p_id$ID_full_siblings)))


## Full sisters (mothers)
# Thus conecting the mothers sisters (not the childrens sibling). 
# Mor in G1: can not connect her to her siblings (Do not have that information, can only connect the data I have and highest in the tree will always be a mother och a father alone whit no siblings)

p_id = p_id %>% 
  ungroup() %>%
  select(one_of(c("lpnr_BARN", "ID_full_siblings"))) %>%
  left_join(p_id, ., by=c("lpnr_mor" = "lpnr_BARN"))

colnames(p_id) = c("lpnr_BARN","lpnr_mor","ar_mor","lpnr_far","ar_far","parents","ID_full_siblings", "ID_mor_full_siblings") # obs not sbilings, sisters

sibling = p_id %>% ungroup() %>% select(lpnr_mor, ID_mor_full_siblings)
sibling = unique(sibling)
dat1 = left_join(dat, sibling, by="lpnr_mor")

rm(p_id,sibling)


## Full maternal sisters giving ptd in any child (past+ future)
dat3 = dat1 %>% group_by(lpnr_mor) %>%
  mutate(PTDpermor =  sum(PTB)) %>%
  ungroup()

test = dat3 %>% group_by(ID_mor_full_siblings)  %>%
mutate(maternal_sisters_children_born_PTB_1 =  sum(PTB) - PTDpermor) %>%
mutate(maternal_sisters_children_born_PTB_1 = ifelse(min(lpnr_mor) == max(lpnr_mor) | is.na(ID_mor_full_siblings),NA,maternal_sisters_children_born_PTB_1)) %>%
ungroup() %>% group_by(lpnr_mor) %>%
mutate(maternal_sisters_children_born_PTB = max(maternal_sisters_children_born_PTB_1)) %>%
 ungroup() %>% select(lpnr_BARN,maternal_sisters_children_born_PTB)

dat = left_join(dat3,test, by="lpnr_BARN")


##  Maternal sister born preterm
dat2 = dat %>% group_by(lpnr_mor) %>% filter(row_number()==1) %>% ungroup() %>%  group_by(ID_mor_full_siblings) %>% mutate(maternal_siblings_PTB =  ifelse(any(mother_herself_PTB==1 & !is.na(ID_mor_full_siblings), na.rm = TRUE),sum(mother_herself_PTB, na.rm=TRUE)-mother_herself_PTB,0)) %>% mutate(maternal_siblings_PTB = ifelse(min(lpnr_mor) == max(lpnr_mor) | is.na(ID_mor_full_siblings),NA,maternal_siblings_PTB)) %>% mutate(maternal_siblings_PTB = ifelse(maternal_siblings_PTB>=1,1,0)) %>% ungroup() %>% select(lpnr_mor,maternal_siblings_PTB)  # obs not siblings, sisters

dat3 = left_join(dat,dat2, by="lpnr_mor")

dat = dat3
rm(dat1,dat2,dat3)


## Father born preterm
barn = dat %>% pull(lpnr_BARN)
far_also_barn_in_mfr =  dat[dat$lpnr_far %in% barn,]
far_also_barn_in_mfr = far_also_barn_in_mfr %>% select(lpnr_far)
far_also_barn_in_mfr = unique(far_also_barn_in_mfr)
far_as_barn = inner_join(dat,far_also_barn_in_mfr, by = c("lpnr_BARN" ="lpnr_far")) # The pregnancies in which the fathers where born

far_as_barn = far_as_barn %>% mutate(father_himself_PTB = ifelse(PTB == 1, 1,0)) %>%
  select(lpnr_BARN,father_himself_PTB) #lpnr_BARN here are barn that also are fathers in mfr

dat = full_join(dat, far_as_barn, by = c("lpnr_far" ="lpnr_BARN"))

dat$father_himself_PTB = factor(dat$father_himself_PTB)




#### Varibles as Factors ####
dat$max_grade_mor_c = factor(dat$max_grade_mor_c, levels = c(2,1,3,0))
dat$max_grade_far_c = factor(dat$max_grade_far_c, levels = c(2,1,3,0))

dat$Parity_logreg = factor(dat$Parity_logreg, levels = c(2,1,3,4))
dat$mor_birth_country_NORDIC = factor(dat$mor_birth_country_NORDIC)

dat$MALDER2 = dat$MALDER*dat$MALDER
dat$AR2 = dat$AR*dat$AR

dat$mother_herself_PTB = factor(dat$mother_herself_PTB)




#### Removing NAs ####
dat1 = dat %>% ungroup()%>% filter(max_grade_mor_c != 0,max_grade_far_c != 0, !is.na(Parity_logreg), !is.na(KON), !is.na(AR), !is.na(mor_birth_country_NORDIC), !is.na(MALDER))




#### Only spontaneous birth ####
source(snakemake@params[[1]]) # "1_cleaning_modules.R"
#source("/home/karin/Parity_Project1/scripts/functions/1_cleaning_modules.R")
year_matrix = NULL
dat_m2 = fun_spont1990(dat1)



#### Models ####
## whole population
dat_m21 = group_by(dat_m2, lpnr_mor) %>% filter(n()>1)

## one random sister per sister ID
set.seed(42)
rows <- sample(nrow(dat_m21))
dat_m2 <- dat_m21[rows, ]
dat_m2 = dat_m2 %>% group_by(ID_mor_full_siblings) %>% filter(row_number()==1) %>% pull(lpnr_mor)
dat_m22 = dat_m21[dat_m21$lpnr_mor %in% dat_m2,]


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
m1 = lmer(GRDBS ~ mother_herself_PTB + Parity_logreg + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m22,control = lmerControl(optimizer ="bobyqa"))
betam1=beta_CI(m1, 0.025,0.975)
print(summary(m1))

# father himself ptb
m2 = lmer(GRDBS ~ father_himself_PTB + Parity_logreg + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m22, control = lmerControl(optimizer ="bobyqa"))
betam2 = beta_CI(m2,0.025,0.975)
print(summary(m2))

# maternal sibling ptb
m3 = lmer(GRDBS ~ maternal_siblings_PTB + Parity_logreg + KON + AR + AR2 + MALDER2  + as.numeric(max_grade_mor) + as.numeric(max_grade_far)+ (1|lpnr_mor), data = dat_m22, control = lmerControl(optimizer ="bobyqa"))
betam3 = beta_CI(m3,0.025,0.975)
print(summary(m3))

# Maternal sister delivering ptb
m4 = lmer(GRDBS ~ as.numeric(maternal_sisters_children_born_PTB) + Parity_logreg + KON + AR + AR2 + MALDER2  + as.numeric(max_grade_mor) + as.numeric(max_grade_far)+ (1|lpnr_mor), data = dat_m22, control = lmerControl(optimizer ="bobyqa"))
betam4 = beta_CI(m4,0.025,0.975)
print(summary(m4))
print(betam4)


## Interaction
# mother preterm interaction
mI1 = lmer(GRDBS ~ Parity_logreg*mother_herself_PTB + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m22,control = lmerControl(optimizer ="bobyqa"))
print(summary(mI1))

# father preterm interaction
mI2 = lmer(GRDBS ~ Parity_logreg*father_himself_PTB + KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m22, control = lmerControl(optimizer ="bobyqa"))
print(summary(mI2))

# maternal sister preterm interaction
mI3 = lmer(GRDBS ~ Parity_logreg*maternal_siblings_PTB + KON + AR + AR2 + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m22, control = lmerControl(optimizer ="bobyqa"))
print(summary(mI3))

# Maternal sister delivering ptd
mI4 = lmer(GRDBS ~ Parity_logreg*as.numeric(maternal_sisters_children_born_PTB) + KON + AR + AR2 + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m22, control = lmerControl(optimizer ="bobyqa"))
print(summary(mI4))



#### Ploting results ####
library(ggplot2)

plotdata <- data.frame(Beta = c(betam1[,2],betam3[,2],betam4[,2],betam2[,2] ),
                       CILow = c(betam1[,1],betam3[,1],betam4[,1],betam2[,1] ),
                       CIHigh = c(betam1[,3],betam3[,3],betam4[,3],betam2[,3]),
                       Model = c("Mother herself preterm","Maternal sister born preterm","Maternal sister delivering preterm","Father himself preterm"))

plotdata$Beta = as.numeric(plotdata$Beta)
plotdata$CILow = as.numeric(plotdata$CILow)
plotdata$CIHigh = as.numeric(plotdata$CIHigh)

p = plotdata %>%
  ggplot(aes(x = Beta, y = Model)) +
  geom_errorbarh(aes(xmin = CILow, xmax = CIHigh),size = .6, height = 0.17, color = "black") +
  geom_point(aes(shape=Model, fill = Model), size = 3, stroke=0.4) +
  geom_vline(aes(xintercept = 0), linetype = 2, color ="black", size = 0.6) +
  coord_cartesian(xlim = c(-4, 0.3), clip ="off") +
  ylab("") +
  xlab("Beta (days)") +
  theme_bw() +
  theme(panel.spacing.y = unit(10, "points"),
        axis.text.y = element_blank(),
        axis.ticks.length.y = unit(0, "points"),
        strip.background.y = element_blank(),
        strip.placement = "outside",
        axis.line = element_line(),
        axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 11, color="black"),
        strip.text = element_text(size = 11),
	text=element_text(family="Helvetica"),
	legend.title = element_text(size=10.5),
	legend.text = element_text(size=10.5),
	legend.position = "bottom",
	legend.direction = "vertical",
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	panel.border = element_rect(colour ="black", fill=NA, size=0.50)
				  
  ) +
  scale_shape_manual(labels=c("Mother herself preterm","Maternal sister born preterm","Maternal sister delivering preterm","Father himself preterm"),values=c("Mother herself preterm"=24,"Maternal sister born preterm"=23,"Maternal sister delivering preterm"=22,"Father himself preterm"=21), guide = guide_legend(title.position="top",title.hjust =0.5, title="Family history")) +
  scale_fill_manual(values=c("Mother herself preterm" ="#fef0d9","Maternal sister born preterm"= "#fdcc8a","Maternal sister delivering preterm"="#fc8d59","Father himself preterm"="#d7301f"),labels=c("Mother herself preterm","Maternal sister born preterm","Maternal sister delivering preterm","Father himself preterm"),guide = guide_legend(title.position="top",title.hjust =0.5, title = "Family history")) + coord_fixed(ratio=.2) # coord_cartesian( clip ="off")




#### Info from models ####
model_info_final = rbind(c("whole_pop_mother_preterm",nobs(m1),m1@Gp[2]),
                         c("whole_pop_fater_preterm", nobs(m2),m2@Gp[2]),
			 c("whole_pop_maternal_sister_ptb",nobs(m3),m3@Gp[2]),
			 c("whole_pop_maternal_sister_delivering_ptb",nobs(m4),m4@Gp[2]),
                         c("interaction_mother_preterm", nobs(mI1),mI1@Gp[2]),
                         c("interaction_fater_preterm", nobs(mI2),mI2@Gp[2]),
                         c("interaction_maternal_sister_preterm",nobs(mI3),mI3@Gp[2]),
			 c("interaction_maternal_sister_delivering_preterm",nobs(mI3),mI3@Gp[2]) )




#### Saving ####
#ggsave("/home/karin/Parity_gestational_duration_interaction/results/supplementary/family_history_no_parity_interaction_supp_fig3.png",p,dpi=1200)
ggsave(snakemake@output[[1]],p,dpi=1200)

#fwrite(as.data.frame(model_info_final),"/home/karin/Parity_gestational_duration_interaction/results/supplementary/family_history_no_parity_interaction_model_info.csv",sep=",")
fwrite(as.data.frame(model_info_final),snakemake@output[[2]],sep=",")

#fwrite(as.data.frame(plotdata),"/home/karin/Parity_gestational_duration_interaction/results/supplementary/family_history_no_parity_interaction_plotdata.csv", sep=",")
fwrite(as.data.frame(plotdata),snakemake@output[[3]], sep=",")

