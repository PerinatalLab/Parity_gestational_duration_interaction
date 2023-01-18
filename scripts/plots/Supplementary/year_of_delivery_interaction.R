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

#dat = dat %>% mutate(AR_c = ifelse(AR<1997,1,0),
#                           AR_c = ifelse(AR>=1997 & AR<2005,0,AR_c),
#                           AR_c = ifelse(AR>=2005,2,AR_c))
#dat$AR_c = factor(dat$AR_c)




#### Removing missings ####
dat1 = dat %>% filter(max_grade_mor_c != 0,max_grade_far_c != 0, !is.na(Parity_logreg), !is.na(KON), !is.na(AR), !is.na(mor_birth_country_NORDIC), !is.na(MALDER))




#### Only spontaneous deliveries ####
source(snakemake@params[[1]])
source("~/Parity_Project1/scripts/functions/1_cleaning_modules.R") # "1_cleaning_modules.R"

year_matrix = NULL
dat_m2 = fun_spont1990(dat1)




#### At least 2 preg per woman ####
dat_m2 = group_by(dat_m2, lpnr_mor) %>% filter(n()>1)

mI = lmer(GRDBS ~ Parity_logreg*AR + Parity_logreg*AR2 +  KON + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m2,control = lmerControl(optimizer ="bobyqa") )
print(summary(mI))




#### Estimating gestational duration for year i +/- 3 years by parity (to investigate interaction in mI) ####
name1 = function(i) {
        #i = year
        dat_m2_temp = dat_m2[dat_m2$AR>i-3 & dat_m2$AR < i+3,]
        m1 = lmer(GRDBS ~ Parity_logreg + KON + mor_birth_country_NORDIC + MALDER + MALDER2 + as.numeric(max_grade_mor) + as.numeric(max_grade_far) + (1|lpnr_mor), data = dat_m2_temp, control = lmerControl(optimizer ="bobyqa") )
        z =summary(m1)
        h1 = as.data.frame(z$coefficients[2:4,1:2])
        
        h1$n = length(unlist(m1@flist))
        #h1$n = nobs(m1)
        h1$AR = i
        h1$Parity_logreg = rownames(h1)
	
	err = confint(m1, level=0.95, method = "Wald")[4:6,]

	h1$CI_min = err[1:3,1]
	h1$CI_max = err[1:3,2]
	
        return(h1)

}

dflist= lapply(1990:2012, name1)

df = do.call("rbind", dflist)

colnames(df)[2] = "sd"

# Plotting resutls
p = ggplot(df, aes(x = AR, y= Estimate,group= Parity_logreg))  +  geom_hline(aes(yintercept = 0), linetype = 2, colour="black", size=1) +geom_point(aes(shape = Parity_logreg,fill=Parity_logreg), size = 3,stroke = 0.4) +  ylab("Beta (days)") + geom_line() +
  xlab("Year of birth") + theme_bw() +
  theme(
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
  scale_shape_manual(labels=c("\u2265 three","two","zero"),values=c(22, 21,23), guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5,title="Parity")) +
  scale_fill_manual(values=c("#E69F00","#56B4E9","#CC79A7"),labels=c("\u2265 three","two","zero"),guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5,title="Parity")) + 
  geom_ribbon(aes(ymin = as.numeric(CI_min), ymax = as.numeric(CI_max), fill = Parity_logreg),alpha=0.3,lwd=0.4,lty="dashed", show.legend = F) +
  scale_x_continuous(limits = c(1990,2012), expand = c(0, 0))

#ggsave("~/Parity_Project_gd/scripts/plots/Supplementary/output/year_of_birth.png",p, width = 17.4, height = 15, units ="cm")
ggsave(snakemake@output[[1]],p, width = 17.4, height = 15, units ="cm")



#### Gestational duration mean for year i +/- 3 years by parity (to investigate interaction in mI) ####
name1 = function(i) {
        #i = year
        dat_m2_temp = dat_m2[dat_m2$AR>i-3 & dat_m2$AR < i+3,]

	h = dat_m2_temp %>% group_by(Parity_logreg) %>% summarize(MEAN = mean(GRDBS),SD = sd(GRDBS))
	
        h$AR = i
	h$n = nrow(dat_m2_temp)
        return(h)
}

dflist = lapply(1990:2012, name1)
df = do.call("rbind", dflist)

# CI
df$SE = df$SD / sqrt(df$n)
alpha = 0.05
df$DF = df$n - 1
df$t_score = qt(p=alpha/2, df=df$DF,lower.tail=F)
df$CI_min <- df$MEAN - df$t_score * df$SE
df$CI_max <- df$MEAN + df$t_score * df$SE

# Plotting resutls
levels(df$Parity_logreg)[levels(df$Parity_logreg) == 4] = "Parity \u2265 three"                                                                     
levels(df$Parity_logreg)[levels(df$Parity_logreg) == 3] = "Parity two"
levels(df$Parity_logreg)[levels(df$Parity_logreg) == 2] = "Parity one"
levels(df$Parity_logreg)[levels(df$Parity_logreg) == 1] = "Parity zero"
df$Parity_logreg = factor(df$Parity_logreg, levels = c("Parity \u2265 three","Parity two","Parity one","Parity zero"), ordered = TRUE)

p = ggplot(df, aes(x=AR, y=MEAN, group=Parity_logreg)) + geom_point(aes(shape = Parity_logreg,fill=Parity_logreg), size = 3,stroke = 0.4)  + geom_line() + xlab("Year of birth") + ylab("Mean(Gestational duration) (days)") + theme_bw() +  theme(
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
  scale_shape_manual(values=c(22, 21,24, 23), guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5,title="Parity")) +
  scale_fill_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7"),guide = guide_legend(reverse=T,title.position="top",title.hjust =0.5,title="Parity"))+   geom_ribbon(aes(ymin = as.numeric(CI_min), ymax = as.numeric(CI_max), fill = Parity_logreg),alpha=0.3,lwd=0.4,lty="dashed", show.legend = F) + scale_x_continuous(limits = c(1990,2012), expand = c(0, 0))


#ggsave("~/Parity_Project_gd/scripts/plots/Supplementary/output/year_of_birth_meanGD.png",p, width = 17.4, height = 15, units ="cm")
ggsave(snakemake@output[[2]],p, width = 17.4, height = 15, units ="cm")

print("Sample size:")
print(summary(df$n))

 
