# FILTERING


#### Adding the functions in these two scripts to the enviroment ####
source(snakemake@params[[1]]) # "1_cleaning_modules.R"
source(snakemake@params[[2]]) # "2_renumber_parity_to_parityF.R"




#### Cleaning the data by Removing maternal with no ID, Removing mother-kids duplications, Remove Stillbirths, Remove rows with no GestAge and Only keeping singeltons. Also cleaning parity ####

## Loading data
library(data.table)
mfr_raw = fread(snakemake@input[[1]]) # load data
print("Number of pregnancies in the data")
print(nrow(mfr_raw)) #Number of pregnancies in the data


###  Removing maternal with no ID -> Removing mother-kids duplications -> Remove Stillbirths -> Remove rows with no GestAge -> Only keeping singeltons -> cleaning parity
year_matrix = NULL  # For the function fun_momID
cleaning = function(dat){
  
  dat = fun_momID(dat)  # Removing mothers without IDs
  colnames(dat)[grep("lpnr_barn",colnames(dat))]  = "lpnr_BARN"  # Change the name the column with child id 
  print("Removed mothers without IDs, pregnancies left:")
  print(row(dat))

  dat = fun_momkidID(dat)  # Removing mother-child duplications 
  print("Removed mother-child duplications, pregnancies left:")
  print(nrow(dat))

  dat = fun_deadborn(dat)  # Removing stillbirths
  print("Removed stillbirths, pregnancies left:")
  print(nrow(dat))

  dat = fun_GAmiss(dat)  # Removing pregnancies with no gestational age 
  print("Removed pregnancies with no gestational age, pregnancies left:")
  print(nrow(dat))

  dat = fun_multipregs(dat)   # Only keeping singeltons 
  print("Only keeping singeltons, pregnancies left:")
  print(nrow(dat)) 

  nnn = recalculateParity(dat, thr_d_low = 31, thr_d_upp = 200)  # Function resulting in the parity (number of pregnancies) for each woman (cleaning parity)
  mfr = cbind(dat, data.frame(parity_clean = nnn))  # Adding the parity to the "input data" that has been filtered
  
  dat = mfr %>% filter(BORDF2==1)
  mor = dat %>% group_by(lpnr_mor) %>% filter(parity_clean==1) %>% filter(n()>1) %>% pull(lpnr_mor)
  dat = dat[!(dat$lpnr_mor %in% mor),]
  print("Only keeping singeltons, pregnancies left:")
  print(nrow(dat))

  print("Done: Removed maternal with no ID -> Removed mother-kids duplications -> Removed Stillbirths -> Removed rows with no GestAge -> Only keeping singeltons -> cleaning parity")

dat
}

mfr_new = cleaning(mfr_raw)  # Cleaning data with the function cleaning
print("Number of pregnancies in the data after the filering function:")
print(nrow(mfr_new))




####
## Filtering GD?
# No using GRDBS. 
# GRDBS range between 154 to 321 days. The number of GD > 310 days is 2774, and more common in the earlier years of the register. Maybe there are some reason that the GD was bigger than 310 days, for now they are included.
####




#### Cleaning for outliers based on deviating child weight and gestational duration correlation with a Generalised Additive Models for Location Scale and Shape ####
library(dplyr)
library(gamlss)
b=subset(mfr_new,select=c("BVIKTBS","GRDBS")) # extracting weigth and gestational duration from data
b=na.omit(b) # removing incompleacte vases of weigth and gestational duration

newy=b$BVIKTBS
newx=b$GRDBS

mBCTcs<- gamlss(BVIKTBS~cs(GRDBS),family=BCT, data=b) # A generalised additive model utilizing a Box-Cox-t distribution and cubic spline (cs) smoothing term
saveRDS(mBCTcs, file = "~/Parity_Project1/mBCTcs.rds")
h=centiles.pred(mBCTcs, xname="GRDBS",xvalues=newx, type="s",dev = c( -5, -4.5, -4,-3.5, -3, -2, -1, 0, 1, 2, 3,3.5, 4,4.5, 5)) # looks good!

mfr_new$sq = seq(nrow(mfr_new)) # to restore the original order in the end
lower<-h[,3]
upper<-h[,15]
dat2<- arrange(mfr_new,GRDBS) %>% filter(!is.na(BVIKTBS), !is.na(GRDBS)) %>% mutate(rows_to_remove = ifelse(((BVIKTBS<lower)|
                 (BVIKTBS>upper)),T,F)) # to remove = T, removing any values above 4.5 or below -4.5 standard deviations (upper and lower)
dat3 <- filter(dat2,rows_to_remove==F)
#ggplot(dat3, aes(x=GRDBS,y=BVIKTBS)) + geom_point() # control removal

dat = dat3[order(dat3$sq),] # to restore order
rm(mfr_new,mBCTcs,dat2,dat3,b)

print("Pregnancies left after cleaning for outliers based on deviating child weight and gestational duration correlation:")
print(nrow(dat))




#### IVF ####
## Filter IVF cases and making a variable telling if this mother had a problem becomming pregnant = unwilling subfertility ##
# Variable
dat = dat %>% mutate(unwilling_subfertility = ifelse(is.na(OFRIBARN) | OFRIBARN == 0,0, OFRIBARN)) # 0 = no unwilling subfertility, variable >= 1 --> years in unwilling subfertility

# Filtering
# Variables below regards action on infertility --> pregancies with a 1 (=yes) are removed
# - OFRIIATG
# - OFRIABEF
# - OFRISTIM
# - OFRIKIRU
# - OFRIICSI
# - OFRIANN
dat = dat %>% filter(is.na(OFRIIATG), is.na(OFRIABEF), is.na(OFRISTIM), is.na(OFRIKIRU), is.na(OFRIICSI), is.na(OFRIANN))

print("Pregnancies left after removing pregnancies with action on infertility:")
print(nrow(dat))



#### Saving ####
fwrite(dat, snakemake@output[[1]], sep=",")



