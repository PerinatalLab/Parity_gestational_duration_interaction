###### Description of study participants ######
library(data.table)
library(dplyr)


####Loading all data, filtering and creating all variables ####
#dat = fread("/mnt/hdd/common/karin/mfr_150202_recored_filtered_p1_all_variables.csv")
dat = fread(snakemake@input[[1]])

#p_id = fread("/mnt/hdd/data/swed/ver/parents.csv") 
p_id = fread(snakemake@input[[2]])

p_id = p_id %>% select(LopnrBarn, LopnrMor, FoddArBioMor, LopnrFar, FoddArBioFar) # selecting columns of interest
colnames(p_id) = c("lpnr_BARN", "lpnr_mor", "ar_mor", "lpnr_far", "ar_far") # changing the col names
#nrow(p_id) = 4138617


### full siblings barn
p_id$parents = paste(p_id$lpnr_mor, p_id$lpnr_far, sep = "_") # create a variable called parents.

p_id = p_id %>% group_by(parents) %>% mutate(ID_full_siblings = ifelse(n()>1,cur_group_id(),NA)) # based on parents connect all full siblings
cat("Number of full sibling groups based on ID_full_siblings:", length(unique(p_id$ID_full_siblings)))


## Full siblings MOR
# Thus conecting the mothers siblings (not the childrens sibling). 
# Mor in G1: can not connect her to her siblings (Do not have that information, can only connect the data I have and highest in the tree will always be a mother och a father alone whit no siblings)
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

rm(p_id,sibling)


##  maternal sister born PTD
dat2 = dat1 %>% group_by(lpnr_mor) %>% filter(row_number()==1) %>% ungroup() %>%  group_by(ID_mor_full_siblings) %>%mutate(maternal_siblings_PTB =  ifelse(any(mother_herself_PTB==1 & !is.na(ID_mor_full_siblings), na.rm = TRUE),sum(mother_herself_PTB, na.rm=TRUE)-mother_herself_PTB,0)) %>% mutate(maternal_siblings_PTB = ifelse(min(lpnr_mor) == max(lpnr_mor) | is.na(ID_mor_full_siblings),NA,maternal_siblings_PTB)) %>% mutate(maternal_siblings_PTB = ifelse(maternal_siblings_PTB>=1,1,0)) %>% ungroup() %>% select(lpnr_mor,maternal_siblings_PTB)

dat3 = left_join(dat1,dat2, by="lpnr_mor")


## Full maternal sisters giving ptd in any child before you

dat3 = dat3 %>% group_by(lpnr_mor) %>% arrange(BFODDAT) %>% 
mutate(cumPTDpermor =  cumsum(PTB)) 

test = dat3 %>% group_by(ID_mor_full_siblings) %>% arrange(BFODDAT) %>% 
mutate(maternal_sisters_children_born_PTB =  cumsum(PTB) -cumPTDpermor,0) %>% 
mutate(maternal_sisters_children_born_PTB = ifelse(min(lpnr_mor) == max(lpnr_mor) | is.na(ID_mor_full_siblings),NA,maternal_sisters_children_born_PTB)) %>% 
#mutate(maternal_sisters_children_born_PTB = ifelse(maternal_sisters_children_born_PTB>=1,1,0)) %>% 
ungroup() %>% select(lpnr_BARN,maternal_sisters_children_born_PTB)

dat4 = left_join(dat3,test, by="lpnr_BARN")

dat = dat4


## Prep variables
## Making variables as factors and remove rows with NAs
dat = dat4
rm(dat1,dat2,dat3)

# Father born preterm
barn = dat %>% ungroup() %>% pull(lpnr_BARN)
far_also_barn_in_mfr =  dat[dat$lpnr_far %in% barn,]
far_also_barn_in_mfr = far_also_barn_in_mfr %>% ungroup() %>% select(lpnr_far)
far_also_barn_in_mfr = unique(far_also_barn_in_mfr)
far_as_barn = inner_join(dat,far_also_barn_in_mfr, by = c("lpnr_BARN" ="lpnr_far")) # The pregnancies in which the fathers where born

far_as_barn = far_as_barn %>% ungroup %>% mutate(father_himself_PTB = ifelse(PTB == 1, 1,0)) %>%
  select(lpnr_BARN,father_himself_PTB) #lpnr_BARN here are barn that also are fathers in mfr

dat = full_join(dat, far_as_barn, by = c("lpnr_far" ="lpnr_BARN"))

# Prev ptd
dat = dat %>% group_by(lpnr_mor) %>% arrange(parity_clean) %>%  mutate(diff_p = parity_clean - dplyr::lag(parity_clean)) %>% mutate(prev_PTD = ifelse(diff_p == 1,prev_PTD,NA)) ############ Remove later when rerunning the create variable script


## Remove pregnancies before 1990
dat1 = filter(dat, AR>=1990)
nrow3 = nrow(dat1)


## Only spontaneous births
#source("/home/karin/Parity_Project1/scripts/functions/1_cleaning_modules.R") 
source(snakemake@params[[1]]) # "1_cleaning_modules.R" 
year_matrix = NULL
dat_m2 = fun_spont1990(dat1)
dat = dat_m2
nrow4 = nrow(dat)


## Remove missings
dat1 = dat %>% filter(max_grade_mor_c != 0,max_grade_far_c != 0, !is.na(Parity_logreg), !is.na(KON), !is.na(AR), !is.na(mor_birth_country_NORDIC), !is.na(MALDER))
nrow2 = nrow(dat1)

dat1 = dat1 %>% group_by(lpnr_mor) %>% filter(n()>1) # atleast 2 preg per mother
nrow5 = nrow(dat1)


## Exculded pregnancies from steps above
excluded_pregnancies = rbind(c("Exclution criteria","Number of pregnancies left in data set"),
			     c("Before1990",nrow3),
			     c("Iatrogenic deliveries",nrow4),
			     c("Missing values",nrow2),
			     c("Mothers with at least 2 preg", nrow5))

print("Excluded pregnancies in last steps of data cleaning")
print(excluded_pregnancies)




#### Description of study pariticipants ####
dat = dat1
## Parity 
b=table(dat$Parity_logreg)

descriptive= rbind( c(table(dat$Parity_logreg),sum(table(dat$Parity_logreg))),
		    c((table(dat$Parity_logreg)/nrow(dat))*100,sum(table(dat$Parity_logreg))/nrow(dat)) )
rownames(descriptive) = c("N","N (%)")



## Gestational duration
a = dat %>% group_by(Parity_logreg) %>% summarize(median(GRDBS),
						  lowerquantile = quantile(GRDBS, probs = 0.25),
						  upperquantile = quantile(GRDBS, probs =0.75))
c = dat %>% ungroup() %>% summarize(median(GRDBS),
                                                  lowerquantile = quantile(GRDBS, probs = 0.25),
                                                  upperquantile = quantile(GRDBS, probs =0.75))
descriptive = rbind(descriptive, t(rbind(a[,2:4],c)) )



## Preterm delivery
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$PTB)[,2], sum(table(dat$Parity_logreg, dat$PTB)[,2])))
rownames(descriptive)[nrow(descriptive)] = "Preterm delivery"
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$PTB)[,2]/b, sum(table(dat$Parity_logreg, dat$PTB)[,2])/nrow(dat) ))
rownames(descriptive)[nrow(descriptive)] = "Preterm delivery (%)"



## Sex female
descriptive = rbind(descriptive,
		   c(table(dat$Parity_logreg, dat$KON)[,1], sum(table(dat$Parity_logreg, dat$KON)[,1])))
rownames(descriptive)[nrow(descriptive)] = "Sex female"
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$KON)[,1]/b, sum(table(dat$Parity_logreg, dat$KON)[,1])/nrow(dat) ))
rownames(descriptive)[nrow(descriptive)] = "Sex female (%)"



## Year of delviery
a = dat %>% group_by(Parity_logreg) %>% summarize(median(AR),
                                                  lowerquantile = quantile(AR, probs = 0.25),
                                                  upperquantile = quantile(AR, probs =0.75))

c = dat %>% ungroup() %>% summarize(median(AR),
                                                  lowerquantile = quantile(AR, probs = 0.25),
                                                  upperquantile = quantile(AR, probs =0.75))
descriptive = rbind(descriptive, t(rbind(a[,2:4],c)) )



## Maternal age
a = dat %>% group_by(Parity_logreg) %>% summarize(median(MALDER),
                                                  lowerquantile = quantile(MALDER, probs = 0.25),
                                                  upperquantile = quantile(MALDER, probs =0.75))
c = dat %>% ungroup() %>% summarize(median(MALDER),
                                                  lowerquantile = quantile(MALDER, probs = 0.25),
                                                  upperquantile = quantile(MALDER, probs =0.75))
descriptive = rbind(descriptive, t(rbind(a[,2:4],c)) )



## Maternal education
descriptive = rbind(descriptive,
                   cbind(t(table(dat$Parity_logreg, dat$max_grade_mor)),table(dat$max_grade_mor)))
rownames(descriptive)[(nrow(descriptive)-6):nrow(descriptive)] = paste0("maternal eduction_n_",1:7)


c = replicate(7,table(dat$Parity_logreg))
descriptive = rbind(descriptive,
                   cbind(t(table(dat$Parity_logreg, dat$max_grade_mor)/c),table(dat$max_grade_mor)/nrow(dat) ))
rownames(descriptive)[(nrow(descriptive)-6):nrow(descriptive)] = paste0("maternal eduction_%_",1:7)



## Paternal education
descriptive = rbind(descriptive,
                   cbind(t(table(dat$Parity_logreg, dat$max_grade_far)),table(dat$max_grade_far)))
rownames(descriptive)[(nrow(descriptive)-6):nrow(descriptive)] = paste0("paternal eduction_n_",1:7)


c = replicate(7,table(dat$Parity_logreg))
descriptive = rbind(descriptive,
                   cbind(t(table(dat$Parity_logreg, dat$max_grade_far)/c),table(dat$max_grade_far)/nrow(dat) ))
rownames(descriptive)[(nrow(descriptive)-6):nrow(descriptive)] = paste0("paternal eduction_%_",1:7)



## Mor birth contry Nordic
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$mor_birth_country_NORDIC)[,2], sum(table(dat$Parity_logreg, dat$mor_birth_country_NORDIC)[,2])))
rownames(descriptive)[nrow(descriptive)] = "Nordic country of Maternal birth"
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$mor_birth_country_NORDIC)[,2]/b, sum(table(dat$Parity_logreg, dat$mor_birth_country_NORDIC)[,2])/nrow(dat) ))
rownames(descriptive)[nrow(descriptive)] = "Nordic country of Maternal birth (%)"



## First born preterm
c = table(dat$Parity_logreg, is.na(dat$PTB_first_born))

descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$PTB_first_born)[,2], sum(table(dat$Parity_logreg, dat$PTB_first_born)[,2])))
rownames(descriptive)[nrow(descriptive)] = "First born preterm"
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$PTB_first_born)[,2]/c[,1], sum(table(dat$Parity_logreg, dat$PTB_first_born)[,2])/sum(c[,1]) ))
rownames(descriptive)[nrow(descriptive)] = "First born preterm (%)"



## Previous pregnancy preterm
c = table(dat$Parity_logreg, is.na(dat$prev_PTD))

descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$prev_PTD)[,2], sum(table(dat$Parity_logreg, dat$prev_PTD)[,2])))
rownames(descriptive)[nrow(descriptive)] = "Previous preterm"
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$prev_PTD)[,2]/c[,1], sum(table(dat$Parity_logreg, dat$prev_PTD)[,2])/sum(c[,1]) ))
rownames(descriptive)[nrow(descriptive)] = "Previous preterm (%)"



## Mother preterm 
c = table(dat$Parity_logreg, is.na(dat$mother_herself_PTB))

descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$mother_herself_PTB)[,2], sum(table(dat$Parity_logreg, dat$mother_herself_PTB)[,2])))
rownames(descriptive)[nrow(descriptive)] = "Mother preterm"
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$mother_herself_PTB)[,2]/c[,1], sum(table(dat$Parity_logreg, dat$mother_herself_PTB)[,2])/sum(c[,1]) ))
rownames(descriptive)[nrow(descriptive)] = "Mother preterm (%)"



## Maternal siblings born preterm
c = table(dat$Parity_logreg, is.na(dat$maternal_siblings_PTB))

descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$maternal_siblings_PTB)[,2], sum(table(dat$Parity_logreg, dat$maternal_siblings_PTB)[,2])))
rownames(descriptive)[nrow(descriptive)] = "Maternal sibling preterm"
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$maternal_siblings_PTB)[,2]/c[,1], sum(table(dat$Parity_logreg, dat$maternal_siblings_PTB)[,2])/sum(c[,1]) ))
rownames(descriptive)[nrow(descriptive)] = "Maternal sibling preterm (%)"



## Maternal sister given birth preterm
c = table(dat$Parity_logreg, is.na(dat$maternal_sisters_children_born_PTB))

descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$maternal_sisters_children_born_PTB)[,2], sum(table(dat$Parity_logreg, dat$maternal_sisters_children_born_PTB)[,2])))
rownames(descriptive)[nrow(descriptive)] = "Maternal sibling given birth preterm"
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$maternal_sisters_children_born_PTB)[,2]/c[,1], sum(table(dat$Parity_logreg, dat$maternal_sisters_children_born_PTB)[,2])/sum(c[,1]) ))
rownames(descriptive)[nrow(descriptive)] = "Maternal sibling given birth preterm (%)"



## Father preterm 
c = table(dat$Parity_logreg, is.na(dat$father_himself_PTB))

descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$father_himself_PTB)[,2], sum(table(dat$Parity_logreg, dat$father_himself_PTB)[,2])))
rownames(descriptive)[nrow(descriptive)] = "Father preterm"
descriptive = rbind(descriptive,
                   c(table(dat$Parity_logreg, dat$father_himself_PTB)[,2]/c[,1], sum(table(dat$Parity_logreg, dat$father_himself_PTB)[,2])/sum(c[,1]) ))
rownames(descriptive)[nrow(descriptive)] = "Father preterm (%)"




#### Saving ####
fwrite(as.data.frame(descriptive), snakemake@output[[1]], sep=",")
#fwrite(as.data.frame(descriptive), "/home/karin/Parity_Project_gd/plots/descriptive_table.csv", sep=",")
