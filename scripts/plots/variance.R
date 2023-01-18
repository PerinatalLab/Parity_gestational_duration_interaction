library(dplyr)
library(data.table)


#### Loading data ####
dat = fread(snakemake@input[[1]])
#dat = fread("/mnt/hdd/common/karin/Parity_Project_gd/mfr_150202_recored_filtered_p1_all_variables.csv")




#### Variables as factors ####
dat$max_grade_mor_c = factor(dat$max_grade_mor_c, levels = c(2,1,3,0))
dat$max_grade_far_c = factor(dat$max_grade_far_c, levels = c(2,1,3,0))

dat$Parity_logreg = factor(dat$Parity_logreg, levels = c(2,1,3,4))
dat$mor_birth_country_NORDIC = factor(dat$mor_birth_country_NORDIC)

dat1 = dat %>% filter(max_grade_mor_c != 0,max_grade_far_c != 0, !is.na(Parity_logreg), !is.na(KON), !is.na(AR), !is.na(mor_birth_country_NORDIC), !is.na(MALDER))




#### Only including spontaneous pregnancies ####
source(snakemake@params[[1]]) # "1_cleaning_modules.R"
#source("/home/karin/Parity_Project1/scripts/functions/1_cleaning_modules.R") 
year_matrix = NULL
dat_m4 = fun_spont1990(dat1)




#### At least two pregnancies per mother ####
dat_m4 = dat_m4 %>% group_by(lpnr_mor) %>% filter(n()>1)




#### Variance in different parity ####
vardat = dat_m4

varianceParity = rbind(
c("Parity0",var(filter(vardat, Parity_logreg==1)%>%pull(GRDBS))),
c("Parity1",var(filter(vardat, Parity_logreg==2)%>%pull(GRDBS))) ,
c("Parity2",var(filter(vardat, Parity_logreg==3)%>%pull(GRDBS))),
c("Parity>=3", var(filter(vardat, Parity_logreg==4)%>%pull(GRDBS))),
c("Parity3",var(filter(vardat, parity_clean==4)%>%pull(GRDBS))) )




#### Test for significant differences in variance ####
##F-test
vardat = filter(dat_m4, as.numeric(Parity_logreg)<=2)
vartest = var.test(GRDBS ~ Parity_logreg, vardat, alternative = "two.sided")
print("Results vartest:")
print(vartest)

## Levene's test
library(car)
vardat = vardat %>% mutate(p = ifelse(Parity_logreg == 1,"A","B"))
levtest = leveneTest(GRDBS ~ as.factor(p), vardat) 
print("Results leveneTest:")
print(levtest)

## Two-sample Kolmogorov-Smirnov test
ks.test(vardat[vardat$Parity_logreg == 1,]$GRDBS,vardat[vardat$Parity_logreg == 2,]$GRDBS)




#### Density plot for each parity group ####
library(ggridges)
vardat = dat_m4                                                          
levels(vardat$Parity_logreg)[levels(vardat$Parity_logreg) == 4] = "Parity \u2265 three"
levels(vardat$Parity_logreg)[levels(vardat$Parity_logreg) == 3] = "Parity two"
levels(vardat$Parity_logreg)[levels(vardat$Parity_logreg) == 2] = "Parity one"
levels(vardat$Parity_logreg)[levels(vardat$Parity_logreg) == 1] = "Parity zero"

vardat$Parity_logreg = factor(vardat$Parity_logreg, levels = c("Parity zero","Parity one","Parity two","Parity \u2265 three"), ordered = TRUE)

p4 <- ggplot(data=vardat, aes(x=GRDBS, y = ,factor(Parity_logreg, levels = c("Parity \u2265 three","Parity two","Parity one", "Parity zero")), group=Parity_logreg, fill = Parity_logreg)) +
    geom_vline(xintercept = 259, linetype = "dashed", color ="black", size = 0.6) +	
    geom_density_ridges(alpha=.4, lwd = 0.4) +
    xlab("Gestational duration") + ylab("") +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", size=11),
        axis.ticks.length.y = unit(0, "points"),
        strip.background.y = element_blank(),
        strip.placement = "outside",
	axis.text = element_text(size =11,colour="black"), #9.5
        axis.title.x = element_text(size = 11, color="black"), #11
        legend.position= "none",
	text=element_text(family="Helvetica", size =11), #11
	panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
	panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p2 = p4 +  scale_fill_manual(values = c("#CC79A7", "#009E73","#56B4E9","#E69F00")) +
       	annotate("text", x =186, y = c(4.15,3.15,2.15,1.15), label= as.character(round(as.numeric(varianceParity[1:4,2]),digits=2)),size=10/.pt ) + 
	annotate("text", x =170, y = c(1.15,2.15,3.15,4.15), label="\U03c3^2",parse = TRUE,size=10/.pt) + 
	annotate("text", x =176, y = c(1.15,2.15,3.15,4.15), label="=",size=10/.pt) + 
	scale_y_discrete(expand = c(0.01, 0)) + 
	scale_y_discrete(expand = expansion(add= c(0.1,1.7)))
ggsave(snakemake@output[[1]],p2,width =174, height =100, units = "mm", dpi=1200,device=cairo_ps)




#### Plotting tail distibution ####
p = ggplot(vardat, aes(x=GRDBS, group=Parity_logreg)) + geom_density(aes(colour = Parity_logreg, fill = Parity_logreg), alpha =.1, lwd=0.6) + coord_cartesian(xlim = c(154,265), ylim=c(0,0.02)) + theme_bw() + theme(
        plot.title = element_text(face = "bold", size=11),
        axis.ticks.length.y = unit(0, "points"),
        strip.background.y = element_blank(),
        strip.placement = "outside",
        axis.text = element_text(size =10,colour="black"),
        axis.title = element_text(size = 11, color="black"),
        text=element_text(family="Helvetica", size =10),
        panel.grid.major = element_blank(), legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.margin = margin(r = 1, l = 0, t =0, unit = "cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=10)
  ) + xlab("Gestational duration") + ylab("") + scale_fill_manual(values = c("#CC79A7", "#009E73","#56B4E9","#E69F00"), name="Parity", guide = guide_legend(title.position="top",title.hjust =0.5)) + scale_color_manual(values = c("#CC79A7", "#009E73","#56B4E9","#E69F00"), name= "Parity", guide = guide_legend(title.position="top",title.hjust =0.5)) + geom_vline(xintercept = 259, linetype = "dashed", color ="black", size = 0.8)

#ggsave("/home/karin/Parity_Project_gd/plots/tail_density_plot.eps",p, width =174, height =100, units = "mm", dpi=1200)
ggsave(snakemake@output[[2]],p,width =174, height =100, units = "mm", dpi=1200,device=cairo_ps)




#### QQ-plot ####
vardat = filter(dat_m4, parity_clean<=2) %>% group_by(lpnr_mor) %>% filter(max(n())==2) %>% select(lpnr_mor, parity_clean, GRDBS) # only including mothers with both a first and second born

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

ggsave(snakemake@output[[3]],p, width=174, height=150 ,units = "mm", dpi=1200, device=cairo_ps)




#### Quantiles ####
## For supplementary material
vardat = dat_m4

quant =vardat %>% ungroup() %>% group_by(Parity_logreg) %>%
	summarize(ten=quantile(GRDBS,probs=0.1),
            twentyfive=quantile(GRDBS,probs=0.25),
            fifty=quantile(GRDBS,probs=0.5),
	    seventyfive = quantile(GRDBS, probs = 0.75),
	    ninety = quantile(GRDBS, probs =0.9))


#fwrite(as.data.frame(quant), "/home/karin/Parity_Project_gd/plots/quantiles.csv", sep=",")
fwrite(as.data.frame(quant), snakemake@output[[4]], sep=",")
