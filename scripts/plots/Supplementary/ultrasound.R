library(dplyr)
library(data.table)
library(ggplot2)
library(lme4)
library(lmerTest)

dat = fread(snakemake@input[[1]])




#### Step 1 - making variables as factors and remove rows with NAs ####
dat$max_grade_mor_c = factor(dat$max_grade_mor_c, levels = c(2,1,3,0))
dat$max_grade_far_c = factor(dat$max_grade_far_c, levels = c(2,1,3,0))

dat$Parity_logreg = factor(dat$Parity_logreg, levels = c(2,1,3,4)) # The reference group is the second pregnancy (parity 1 (Parity_logreg==2))
dat$mor_birth_country_NORDIC = factor(dat$mor_birth_country_NORDIC)

dat = dat %>% mutate(unwilling_subfertility = ifelse(unwilling_subfertility>=1 ,1,unwilling_subfertility))
dat$unwilling_subfertility = factor(dat$unwilling_subfertility)

dat$diab = factor(dat$diab)
dat$preeclamspia = factor(dat$preeclamspia)

dat$MALDER2 = dat$MALDER*dat$MALDER
dat$AR2 = dat$AR*dat$AR
dat$BMI2 = dat$BMI*dat$BMI

## Remove missings
dat1 = dat %>% filter(max_grade_mor_c != 0,max_grade_far_c != 0, !is.na(Parity_logreg), !is.na(KON), !is.na(AR), !is.na(mor_birth_country_NORDIC), !is.na(MALDER))


#### Step 2 - creating an CI function for the regressions ####
beta_CI = function(Model, CI_min, CI_max) {
  m =character()
  for (i in 2:4) {
          mm = summary(Model)$coefficients[i,1]
          err = confint(Model, level=0.95, method = "Wald")[i,]
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




#### Step 3 - Extracting the spontaneous birth ####
source(snakemake@params[[1]])
year_matrix = NULL
dat_m1 = fun_spont1990(dat1)





#### Step 4 - at least two pregnancies per mother, one random of these pregnancies for the linear models. Only pregnancies estimated with ultrasound. ####
## at least two pregnancies
dat_m2 = dat_m1 %>% group_by(lpnr_mor) %>% filter(n()>1) # atleast 2 preg per mother
table(dat_m2$GRMETOD)

print("only inclung UL estimates")
dat_m1 = dat_m1 %>% filter(GRMETOD == 1 | GRMETOD == 5 | GRMETOD == 6| GRMETOD == 7)
dat_m2 = dat_m1 %>% group_by(lpnr_mor) %>% filter(n()>1) # atleast 2 preg per mother
dat_m3 = dat_m1 %>% filter(AR >=1992) %>%  filter(!is.na(unwilling_subfertility),!is.na(BMI), !is.na(diab), !is.na(preeclamspia), !is.na(smoking)) %>% group_by(lpnr_mor) %>% filter(n()>1)

## one random preg
set.seed(42)
rows <- sample(nrow(dat_m2))
dat_m21 <- dat_m2[rows, ]
dat_m21 = dat_m21 %>% group_by(lpnr_mor) %>% filter(row_number()==1)

set.seed(42)
rows <- sample(nrow(dat_m3))
dat_m31 <- dat_m3[rows, ]
dat_m31 = dat_m31 %>% group_by(lpnr_mor) %>% filter(row_number()==1)



#### Step 5 - running the models ####
## Linear regression
m6 = lm(GRDBS ~ Parity_logreg+ KON + AR + AR2, data = dat_m31)
betam6=beta_CI(m6, 0.025,0.975)

m7 = lm(GRDBS ~ Parity_logreg+ KON + AR + AR2  +mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far), data = dat_m31)
betam7=beta_CI(m7, 0.025,0.975)

m10 = lm(GRDBS ~ Parity_logreg+ KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + unwilling_subfertility + BMI + BMI2 + smoking + diab + preeclamspia, data = dat_m31)
betam10=beta_CI(m10, 0.025,0.975)


## Linear mixed models
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

m16 = lmer(GRDBS ~ Parity_logreg+ KON + AR + AR2 +(1 | lpnr_mor), data=dat_m3, control = lmerControl(optimizer ="bobyqa"))
betam16=beta_CI(m16, 0.025,0.975)

m17 = lmer(GRDBS ~ Parity_logreg+ KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2  + as.numeric(max_grade_mor) + as.numeric(max_grade_far) +(1 | lpnr_mor), data=dat_m3,control = lmerControl(optimizer ="bobyqa"))
betam17=beta_CI(m17, 0.025,0.975)

m11 = lmer(GRDBS ~ Parity_logreg+ KON + AR + AR2 + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + unwilling_subfertility + BMI + BMI2 + smoking + diab + preeclamspia +(1 | lpnr_mor), data = dat_m3, control = lmerControl(optimizer ="bobyqa")) # Nelder_Mead"
betam11=beta_CI(m11, 0.025,0.975)


model_info_final = rbind(c("Model","AIC","Number of observations", "Number of groups"),
                         c("m6",AIC(m6), nobs(m6),"-" ),
                         c("m7",AIC(m7),nobs(m7),"-"),
                         c("m10", AIC(m10),nobs(m10),"-"),
                         c("m16",AIC(m16), nobs(m16),m16@Gp[2]),
                         c("m17",AIC(m17),nobs(m17),m16@Gp[2]),
                         c("m11",AIC(m11),nobs(m11),m11@Gp[2])
                         )

print(model_info_final)

plotdata <- data.frame(Parity = c(rownames(betam6), rownames(betam7),rownames(betam16), rownames(betam17),rownames(betam10), rownames(betam11)),
                       Beta = c(betam6[,2],betam7[,2],betam16[,2],betam17[,2],betam10[,2],betam11[,2]),
                       CILow = c(betam6[,1],betam7[,1],betam16[,1],betam17[,1],betam10[,1],betam11[,1]),
                       CIHigh = c(betam6[,3],betam7[,3],betam16[,3],betam17[,3], betam10[,3],betam11[,3]),
                       Model = c("Set 1","Set 1","Set 1","Set 1","Set 2","Set 2","Set 2","Set 2","Set 1","Set 1","Set 1","Set 1","Set 2","Set 2","Set 2","Set 2","Set 3","Set 3","Set 3","Set 3","Set 3","Set 3","Set 3","Set 3"),
                       Regression = c("Linear","Linear","Linear","Linear","Linear","Linear","Linear","Linear","Linear mixed","Linear mixed","Linear mixed","Linear mixed","Linear mixed","Linear mixed","Linear mixed","Linear mixed","Linear","Linear","Linear","Linear","Linear mixed","Linear mixed","Linear mixed","Linear mixed")
)

plotdata$Parity <- factor(plotdata$Parity, levels = c(">=P3","P2","P1","P0"))
plotdata$Beta = as.numeric(plotdata$Beta)
plotdata$CILow = as.numeric(plotdata$CILow)
plotdata$CIHigh = as.numeric(plotdata$CIHigh)
plotdata$Regression = factor(plotdata$Regression, levels =c("Linear","Linear mixed"))

p = plotdata %>%
  ggplot(aes(x = Beta, y = Parity )) +
  geom_errorbarh(aes(xmin = CILow, xmax = CIHigh),size = .6, height = .2, color = "black") +
  geom_vline(aes(xintercept = 0), linetype = 2, colour="black", size=0.6) +
  geom_point(aes(shape = Parity,fill=Parity), size = 3,stroke = 0.4) +
  facet_grid(Model~Regression, switch = "y") +
  coord_cartesian(xlim = c(-1, 0.6)) +
  ylab("") +
  xlab("Beta (days)") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(color="black"),
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
  scale_fill_manual(values=c(">=P3" ="#E69F00","P2"= "#56B4E9","P1"="#009E73","P0"="#CC79A7"),labels=c("\u2265 three","two","one","zero"),guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5)) +
  facet_wrap(Regression ~ Model, scales = "free")




#### Saving ####
#similar as main fig
ggsave(snakemake@output[[1]],p, width = 174, height = 140, units = "mm", bg = "transparent")
fwrite(as.data.frame(model_info_final),snakemake@output[[2]],sep=",")
fwrite(as.data.frame(plotdata),snakemake@output[[3]], sep=",")




#### Step 7 - Distribution p0 vs p1 ####
vardat = filter(dat_m2, parity_clean<=2) %>% group_by(lpnr_mor) %>% filter(max(n())==2) %>% select(lpnr_mor, parity_clean, GRDBS) # only including mothers with both a first and second born

GD1 = vardat %>% filter(parity_clean == 1) %>% ungroup() %>%  arrange(GRDBS) %>% mutate(exp = 1:length(GRDBS)/length(GRDBS))
GD2 = vardat %>% filter(parity_clean == 2) %>% ungroup() %>% arrange(GRDBS) %>% mutate(exp = 1:length(GRDBS)/length(GRDBS))

k = left_join(GD1,GD2, by = "exp")
p = ggplot(k, aes( x= `GRDBS.x`, y = `GRDBS.y`)) +
        geom_vline(xintercept = 259, linetype = "dashed", color ="black", size = 0.6)+
        geom_point(size=0.5) +
        geom_abline( intercept=0, slope=1,size=0.6) +
        xlab("Quantile parity zero") +
        ylab("Quantile parity one") +
        ggtitle("(b)") +
        theme_bw() +
        theme(
        plot.title = element_text(face = "bold", size=11),
        axis.ticks.length.y = unit(0, "points"),
        strip.background.y = element_blank(),
        strip.placement = "outside",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10, color ="black"),
        text=element_text(family="Helvetica", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(r = 4, l = 3.5, t=0.25,unit = "cm"),
        axis.line = element_line(size=0.01),
        panel.border = element_rect(colour = "black", fill=NA, size=1)


        )
#similar as main fig
#ggsave("/home/karin/fig2_variance.png",p, width = 174, height = 140, units = "mm", bg = "transparent")

