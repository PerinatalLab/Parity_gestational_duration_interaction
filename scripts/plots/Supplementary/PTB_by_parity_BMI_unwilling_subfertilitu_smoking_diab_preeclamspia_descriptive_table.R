library(data.table)
library(dplyr)

dat = fread(snakemake@input[[1]])
#dat = fread("/mnt/hdd/common/karin/mfr_150202_recored_filtered_p1_all_variables.csv")


#### Variables as factors ####
dat$max_grade_mor_c = factor(dat$max_grade_mor_c, levels = c(2,1,3,0))
dat$max_grade_far_c = factor(dat$max_grade_far_c, levels = c(2,1,3,0))

dat$Parity_logreg = factor(dat$Parity_logreg, levels = c(2,1,3,4))
dat$mor_birth_country_NORDIC = factor(dat$mor_birth_country_NORDIC)

dat = dat %>% mutate(unwilling_subfertility = ifelse(unwilling_subfertility>=1 ,1,unwilling_subfertility))
dat$unwilling_subfertility = factor(dat$unwilling_subfertility)

dat$diab = factor(dat$diab)
dat$preeclamspia = factor(dat$preeclamspia)

dat$MALDER2 = dat$MALDER*dat$MALDER
dat$AR2 = dat$AR*dat$AR
dat$BMI2 = dat$BMI*dat$BMI




#### Removing missings ####
dat1 = dat %>% filter(max_grade_mor_c != 0,max_grade_far_c != 0, !is.na(Parity_logreg), !is.na(KON), !is.na(AR), !is.na(mor_birth_country_NORDIC), !is.na(MALDER), !is.na(unwilling_subfertility),!is.na(BMI), !is.na(diab), !is.na(preeclamspia), !is.na(smoking))




#### Only spontaneous births ####
source(snakemake@params[[1]])
#source("/home/karin/Parity_Project1/scripts/functions/1_cleaning_modules.R") # "1_cleaning_modules.R" 
year_matrix = NULL
dat_m2 = fun_spont1990(dat1)
dat_m2 = dat_m2 %>% filter(AR >=1992) %>% group_by(lpnr_mor) %>% filter(n()>1)




#### Descriptive table ####
dat = dat_m2
# Parity 
b=table(dat$Parity_logreg)

descriptive= rbind( c(table(dat$Parity_logreg),sum(table(dat$Parity_logreg))),
                    c(table(dat$Parity_logreg)/nrow(dat),sum(table(dat$Parity_logreg))/nrow(dat)) )
rownames(descriptive) = c("N","N (%)")


# Unwilling_subfertility
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$unwilling_subfertility)[,2], sum(table(dat$Parity_logreg, dat$unwilling_subfertility)[,2])))
rownames(descriptive)[nrow(descriptive)] = "Unwilling subfertility"
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$unwilling_subfertility)[,2]/b, sum(table(dat$Parity_logreg, dat$unwilling_subfertility)[,2])/nrow(dat)))
rownames(descriptive)[nrow(descriptive)] = "Unwilling subfertility (%)"


#BMI
a = dat %>% group_by(Parity_logreg) %>% summarize(median(as.numeric(BMI)),
                                                  lowerquantile = quantile(as.numeric(BMI), probs = 0.25),
                                                  upperquantile = quantile(as.numeric(BMI), probs =0.75))

c = dat %>% ungroup() %>% summarize(median(as.numeric(BMI)),
                                                  lowerquantile = quantile(as.numeric(BMI), probs = 0.25),
                                                  upperquantile = quantile(as.numeric(BMI), probs =0.75))
descriptive = rbind(descriptive, t(rbind(a[,2:4],c)) )


#Smoking 
a = dat %>% group_by(Parity_logreg) %>% summarize(median(as.numeric(smoking)),
                                                  lowerquantile = quantile(as.numeric(smoking), probs = 0.25),
                                                  upperquantile = quantile(as.numeric(smoking), probs =0.75))

c = dat %>% ungroup() %>% summarize(median(as.numeric(smoking)),
                                                  lowerquantile = quantile(as.numeric(smoking), probs = 0.25),
                                                  upperquantile = quantile(as.numeric(smoking), probs =0.75))
descriptive = rbind(descriptive, t(rbind(a[,2:4],c)) )


#Diabetes
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$diab)[,2], sum(table(dat$Parity_logreg, dat$diab)[,2])))
rownames(descriptive)[nrow(descriptive)] = "Diabetes"
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$diab)[,2]/b, sum(table(dat$Parity_logreg, dat$diab)[,2])/nrow(dat)))
rownames(descriptive)[nrow(descriptive)] = "Diabetes (%)"


#Preeclamspia
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$preeclamspia)[,2], sum(table(dat$Parity_logreg, dat$preeclamspia)[,2])))
rownames(descriptive)[nrow(descriptive)] = "Preeclampsia"
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$preeclamspia)[,2]/b, sum(table(dat$Parity_logreg, dat$preeclamspia)[,2])/nrow(dat)))
rownames(descriptive)[nrow(descriptive)] = "Preeclampsia (%)"




#### Saving ####
#fwrite(descriptive, "/home/karin/Parity_Project_gd/plots/Supplementary/Descriptive_table_BMI_unwilling_subfertilitu_smoking_diab_preeclamspia.csv", sep=",")
fwrite(as.data.frame(descriptive), snakemake@output[[1]], sep=",")


