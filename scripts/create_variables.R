library(dplyr)
library(data.table)


#### Loading data ####
dat = fread(snakemake@input[[1]]) # Swedish Medical Birth Register
p_id = fread(snakemake@input[[2]]) # Multi-Generation Register
edu = fread(snakemake@input[[3]]) # Education Register



#### Creating variables ####

## PTB < 37 weeks
dat = dat %>% mutate(PTB = ifelse(GRDBS< (37*7),1,0))


## Sex
# male as ref
# 1=Pojke (boy); 2=Flicka (girl)
dat = dat %>% mutate(KON = ifelse(KON == 1, "M","F"))
dat$KON = factor(dat$KON, level = c("M","F"))


## Maternal age categorized:
# < 20
# 20-29
# 30-39
# >= 40

# Maternal age categorization
dat = dat %>% mutate(MALDER_c = ifelse(MALDER<20,1,0),
                     MALDER_c = ifelse(MALDER>=20 & MALDER<=29,0,MALDER_c),
                     MALDER_c = ifelse(MALDER>=30 & MALDER<=39,2,MALDER_c),
                     MALDER_c = ifelse(MALDER>=40,3,MALDER_c))
dat$MALDER_c = as.factor(dat$MALDER_c)


## Finding father to each child in  the multi-generation register and the fathers age
f_id = p_id %>% select(LopnrBarn, LopnrFar, FoddArBioFar) # selecting columns of interest
colnames(f_id) = c("lpnr_BARN","lpnr_far", "ar_far")
dat = left_join(dat,f_id,by="lpnr_BARN") # adding fathers info from parents.csv to mfr data
rm(f_id)
dat = dat %>% mutate(FALDER = AR-ar_far) # calculating how old the father was when their child was born


## Nationality  
# As maternal citizenship and the mothers birth country
dat = dat %>% mutate(swe_citizenship =as.numeric(MNAT %in% c('SVERIGE'))) %>% mutate(mor_birth_country_NORDIC = as.numeric(MFODLAND %in% c('SVERIGE','NORGE','FINLAND','ISLAND','DANMARK')))


## First child (parity 0) born preterm
dat = dat %>% group_by(lpnr_mor) %>% arrange(parity_clean) %>%
  mutate(PTB_first_born = any(row_number() == 1 & PTB ==1)*1) %>% 
  mutate(PTB_first_born = ifelse(dplyr::first(parity_clean)!=1,NA,PTB_first_born)) #parity_clean == 1 is parity 0
 

## Previous preterm delivery
dat = dat %>% group_by(lpnr_mor) %>% arrange(parity_clean )%>% mutate(prev_PTD = ifelse(dplyr::lag(PTB)==1,1,0))
dat = dat %>% group_by(lpnr_mor) %>% arrange(parity_clean) %>%  mutate(diff_p = parity_clean - dplyr::lag(parity_clean)) %>% mutate(prev_PTD = ifelse(diff_p == 1,prev_PTD,NA))


## Mother was born preterm herself
barn = dat %>% pull(lpnr_BARN)
mor_also_barn_in_mfr =  dat[dat$lpnr_mor %in% barn,]
mor_also_barn_in_mfr = mor_also_barn_in_mfr %>% select(lpnr_mor)
mor_also_barn_in_mfr = unique(mor_also_barn_in_mfr)
mor_as_barn = inner_join(dat,mor_also_barn_in_mfr, by = c("lpnr_BARN" ="lpnr_mor")) # The pregnancies in which the mothers where born

mor_as_barn = mor_as_barn %>% mutate(mother_herself_PTB = ifelse(PTB == 1, 1,0)) %>%
  select(lpnr_BARN,mother_herself_PTB) #lpnr_BARN here are barn that also are mothers in mfr

dat = full_join(dat, mor_as_barn, by = c("lpnr_mor" ="lpnr_BARN"))
dat = select(dat, -lpnr_mor.y)


## Diabetes
dat = dat %>% mutate(diab1 = ifelse(DIABETES != 1 | is.na(DIABETES), 0,1))  #diabetes according to mfr variable Diabetes

test = dat[grepl("O24|E10|E11|E12|E13|E14|648A|250A|250B|250C|250D|250E|250F|250G|250H|250X|25000| 25001| 25002| 25003| 25004| 25005| 25006| 25007| 25008| 2500",paste(dat$MDIAG1,dat$MDIAG2,dat$MDIAG3,dat$MDIAG4,dat$MDIAG5,dat$MDIAG6,dat$MDIAG7,dat$MDIAG8,dat$MDIAG9,dat$MDIAG10,dat$MDIAG11,dat$MDIAG12)),] #ICD codes (ICD-10-SE,ICD9-SE,ICD-8) related to diabetes, extracted from maternal icd diagnosis in mfr
mor_with_diabetes = test %>% pull(sq) # rows of mothers that have diabetes according to icd codes

dat = dat %>% mutate(diab2 = ifelse(sq %in% mor_with_diabetes,1,0))

dat = dat %>% mutate(diab = ifelse(diab1==1 | diab2 ==1,1,0)) # mother will have diabetes based on icd codes and the mfr variable Diabetes
dat = select(dat, -diab1,-diab2) 


## BMI
print("s1")
dat1 = dat %>% mutate(BMI = MVIKT / (MLANGD/100)^2) # BMI
#source("/home/karin/Parity_Project1/scripts/functions/1_cleaning_modules.R")
source(snakemake@params[[1]]) #fun_mBmiQC modified to not remove "bad" BMI, just set them as NA.
print("s2")
year_matrix = NULL
dat2 = fun_mBmiQC(as.data.frame(dat1)) # setting bad BMI to NA
print("s3")
dat = dat2
rm(dat1,dat2)


## Smoking 
dat = dat %>% mutate(smoking = ifelse((ROK1 ==1 | is.na(ROK1)) & (ROK0 == 1 | is.na(ROK0)) ,0,1),
                     smoking = ifelse(ROK2 == 1 |is.na(ROK2),smoking,2)) # 0 = Not smoking, 1 =  Smoking 3 months prior to the current pregnancy or/and Smoking at admission to maternal health, 2 = Smoking in pregnancy week 30-32  


## Preeclampsia
test = dat[grepl("O14|O11|O15|642E|642F|642H|63703 |63704 | 63709| 63710|6612",paste(dat$MDIAG1,dat$MDIAG2,dat$MDIAG3,dat$MDIAG4,dat$MDIAG5,dat$MDIAG6,dat$MDIAG7,dat$MDIAG8,dat$MDIAG9,dat$MDIAG10,dat$MDIAG11,dat$MDIAG12)),] #ICD codes (ICD-10-SE,ICD9-SE,ICD-8) related to preeclampsia, extracted from maternal icd diagnosis in mfr
mor_with_preeclampsia = test %>% pull(sq) # rows of mothers that have preeclampsia according to icd codes

dat = dat %>% mutate(preeclamspia = ifelse(sq %in% mor_with_preeclampsia,1,0))


## Education 
# find the maximum edu + filtering
edu = edu %>% group_by(LopNr) %>% filter(n()==1) # can not tell which of the rows are the ture one when ID for the same person exist in several rows, are removed
edu = as.data.frame(edu)
edu_grades = edu[grep("SUN2000", names(edu))] # education based on SUN2000
edu_grades[, "max"] <- apply(edu_grades, 1, max, na.rm=TRUE) # Finding highest education for each person
edu = cbind(edu, edu_grades[,"max"])
names(edu)[names(edu) == 'edu_grades[, "max"]'] = "max_grade"

# Remove reused LopNr based on AterPnr
edu_rm = edu[grep("Ater", names(edu))]
edu_rm = edu_rm %>% mutate(remove = ifelse(rowSums(edu_rm == 1,na.rm = TRUE) > 0, F, T))
edu = edu[edu_rm$remove,]

# Remove reused LopNr based on SenPnr
edu_rm = edu[grep("Sen", names(edu))]
#nr = ncol(edu_rm)
edu_rm = edu_rm %>% mutate(remove = ifelse(rowSums(edu_rm == 0,na.rm = TRUE) >0 , F, T))
edu = edu[edu_rm$remove,]
#nrow(edu) == 5828310

#Join with mfr
edu_max = edu[grep("LopNr|max_grade", names(edu))]
d_mor = left_join(dat, edu_max, by = c("lpnr_mor" = "LopNr") )
names(d_mor)[names(d_mor) == 'max_grade'] = "max_grade_mor"

d_mor_far = left_join(d_mor, edu_max, by = c("lpnr_far" = "LopNr") )
names(d_mor_far)[names(d_mor_far) == 'max_grade'] = "max_grade_far"

d_mor_far_child = left_join(d_mor_far, edu_max, by = c("lpnr_BARN" = "LopNr") )
names(d_mor_far_child)[names(d_mor_far_child) == 'max_grade'] = "max_grade_child"

dat = d_mor_far_child
rm(edu_grades,edu_max,d_mor,d_mor_far,d_mor_far_child)

# Max_grade in categories
dat = dat %>% mutate(max_grade_mor_c  = ifelse(max_grade_mor==2 | max_grade_mor==1,1,0),                   # 9 years or less
                     max_grade_mor_c = ifelse(max_grade_mor==3 | max_grade_mor ==4,2,max_grade_mor_c),     # Gymnasial utbilding (additional 2-3 years)      
                     max_grade_mor_c = ifelse(max_grade_mor >=5,3,max_grade_mor_c)) # 0 is nas             # Eftergymnasial utbildning (shorter than 3 years, 3 years or longer, postgraduate education)

dat = dat %>% mutate(max_grade_far_c  = ifelse(max_grade_far==2 | max_grade_far==1,1,0),
                     max_grade_far_c = ifelse(max_grade_far==3 | max_grade_far ==4,2,max_grade_far_c),
                     max_grade_far_c = ifelse(max_grade_far >=5,3,max_grade_far_c)) # 0 is nas


## Parity, grouping after parity 4
dat = dat %>% mutate(Parity_logreg = ifelse(as.numeric(parity_clean)<5,parity_clean,4))


#### Saving ####
fwrite(dat, snakemake@output[[1]], sep=",")
