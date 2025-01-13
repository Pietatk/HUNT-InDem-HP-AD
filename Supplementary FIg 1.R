library(Hmisc)
library(foreign)
library(haven)
library(tidyverse)
#load data files
mydata<-  read_sav("~/Desktop/2024 02 12 Skjellegrind, data 113250/2024-02-12_113250_Data.sav") 
mydata1<-read_sav("~/Desktop/2024 08 20 Skjellegrind, fødseldato 113250/2024-08-19_113250_Data.sav")
mydata3<- read_sav("~/Desktop/2024 07 12 Skjellegrind, registerstatus 113250/2024-07-12_113250_Data.sav")
mydata2<- read_sav("~/Desktop/2024 04 23 Skjellegrind, tillegg 113250/2024-04-23_113250_Data.sav")

#remove @
library(janitor)
mydata<-clean_names(mydata)
mydata1<-clean_names(mydata1)
mydata2<-clean_names(mydata2)
mydata3<-clean_names(mydata3)
hpylori2 <- reduce(list(mydata, mydata3, mydata2), merge, by = "pid_113250", all = TRUE)
#select birth year
mydata1<-mydata1[, c("pid_113250","birth_month")]
hpylori2 <- merge(hpylori2,mydata1,by = "pid_113250", all = TRUE)
#renaming variables
hpylori2 <- hpylori2 %>% 
  rename(helicobacterpylori_cat=se_hp_ab_stat_nt2blm) #hpylori_cat (positivity)
hpylori2 <- hpylori2 %>% 
  rename(smoking=smo_stat_nt2blq1) #smoking
hpylori2 <- hpylori2 %>% 
  rename(maritalstatus=marit_stat_nt2blq1) #maritalstatus
hpylori2 <- hpylori2 %>% 
  rename(bmi=bmi_nt2blm)
hpylori2 <- hpylori2 %>% 
  rename(education=educ_nt2blq1) #education 
hpylori2 <- hpylori2 %>% 
  rename(alcohol=alc_tot_unit_w_nt2blq1) #alcohol
hpylori2 <- hpylori2 %>% 
  rename(helicobacterpylori_cont=se_hp_ab_nt2blm) #hpylori_cont (level)
hpylori2 <- hpylori2 %>% 
  rename(dementia=diag_dem_nt4eld1mix) #dementia
hpylori2 <- hpylori2 %>% 
  rename(cognitiveimpairement=diag_cog_nt4eld1mix) #cognitive impairment
hpylori2 <- hpylori2 %>% 
  rename(diabetes=dia_ev_nt2blq1,) #diabetes
hpylori2 <- hpylori2 %>% 
  rename(cvd=car_inf_ev_nt2blq1) #cvd
hpylori2 <- hpylori2 %>% 
  rename(age2=part_ag_nt2blq1) #age HUNT 2
hpylori2 <- hpylori2 %>% 
  rename(age=part_ag_nt4blm) #age HUNT 4 
hpylori2 <- hpylori2 %>% 
  rename(register=register_status.y) #register status
hpylori2 <- hpylori2 %>% 
  rename(crp=sem_crp_nt2blm) #crp
#DFS 
#invited to hunt 2
hpylori2$invhunt2 <- ifelse(hpylori2$part_nt2blm %in% c(0, 1), 1, NA)
table(hpylori2$invhunt2,useNA = "ifany") #invited n=72884
table(hpylori2$part_nt2blm,useNA = "ifany") #0=7867 1=65017
# Filter dataset based on birth year invited and participated 
birthyear <- hpylori2 %>%
  filter(birth_year.x > 1914 & birth_year.x < 1950 & invhunt2 == 1)

birthyear <- hpylori2 %>%
  filter(birth_year.x > 1914 & birth_year.x < 1950 & part_nt2blm == 1)

# Filter for available HP data
hpylori2 <- birthyear %>% filter(!is.na(helicobacterpylori_cat))

# Exclude those who moved based on register status
#hpylori2 <- hpylori2 %>% filter(register != 3) 
#invited to hunt 4
hpylori2$invhunt4 <- ifelse(hpylori2$part_nt4blm %in% c(0, 1), 1, NA)
table(hpylori2$invhunt4,useNA = "ifany") #invited n=2718
hpylori2$register <- ifelse(hpylori2$register == 5 & hpylori2$invhunt4 == 1, 1,hpylori2$register)
hpylori2$register <- ifelse(is.na(hpylori2$register), 0, hpylori2$register) #0=dead 1=alive #3 moved
table(hpylori2$register, useNA = "ifany")

###################
table(hpylori2$part_nt4blm,useNA = "ifany") #participated HUNT 4, 1(participated)=1950, 0(declined)=768, NA=2462 (2363 dead +9 moved+90 alive)
table(hpylori2$part_nt4eld1mix,useNA = "ifany") #participated HUNT 4 70+, 1=1382, 0=322, NA=3476 (2363 dead+9moved+1104 alive)

#Moved 
#9 moved internationally, i.e., register=3 
#sum(hpylori2$register == 1 & !is.na(hpylori2$obs_end_month)) #alive but moved from Trøndelag

#HUNT 4 
sum(hpylori2$part_nt4 == 0 & hpylori2$register == 0,na.rm = TRUE) #replied no to participate in HUNT 4 and dead=0
sum(hpylori2$part_nt4blm == 0 & hpylori2$register == 1,na.rm = TRUE) #replied no to participate and alive=768
sum(is.na(hpylori2$part_nt4blm) & hpylori2$register == 1) #NA and alive=90
sum(is.na(hpylori2$part_nt4blm) & hpylori2$register == 0) #NA and dead=2363

# HUNT 4 70+
sum(hpylori2$part_nt4eld1mix == 0 & hpylori2$register == 0,na.rm = TRUE) #replied no to participate in HUNT 4 and dead=0
sum(hpylori2$part_nt4eld1mix == 0 & hpylori2$register == 1,na.rm = TRUE) #replied no to participate and alive=322
sum(is.na(hpylori2$part_nt4eld1mix) & hpylori2$register == 1) #NA and alive=1104
sum(is.na(hpylori2$part_nt4eld1mix) & hpylori2$register == 0) #NA and dead=2363

#above 70
sum(is.na(hpylori2$part_nt4blm) & hpylori2$register == 1 & hpylori2$part_ag_nt4blm > 70) #NA HUNT 4 and >70=0
sum(is.na(hpylori2$part_nt4eld1mix) & hpylori2$register == 1 & hpylori2$part_ag_nt4blm > 70) #NA HUNT 4 70+ and >70=0
sum(hpylori2$part_nt4blm == 0 & hpylori2$register == 1,na.rm = TRUE & hpylori2$age > 70, na.rm = TRUE) #said no to HUNT 4 and >70=768
sum(hpylori2$part_nt4eld1mix == 0 & hpylori2$register == 1,na.rm = TRUE & hpylori2$age > 70, na.rm = TRUE) #said no to HUNT 4 70+ and >70=322

#below 70
sum(hpylori2$part_nt4blm == 0 & hpylori2$register == 1,na.rm = TRUE & hpylori2$age < 70, na.rm = TRUE) #said no to HUNT 4 and >70=768
sum(hpylori2$part_nt4eld1mix == 0 & hpylori2$register == 1,na.rm = TRUE & hpylori2$age < 70, na.rm = TRUE) #said no to HUNT 4 70+ and >70=322

#participated in both HUNT 4 and HUNT 4 70+
sum(hpylori2$part_nt4blm == 1 & hpylori2$part_nt4eld1mix == 1,na.rm = TRUE) #1382
#participated in only HUNT 4 and not HUNT 4 70+
sum(hpylori2$part_nt4blm == 1 & hpylori2$part_nt4eld1mix == 0,na.rm = TRUE) #322

#participated in only HUNT 4 and NA HUNT 4 70+
sum(is.na(hpylori2$part_nt4eld1mix) & hpylori2$part_nt4blm == 1) #0


