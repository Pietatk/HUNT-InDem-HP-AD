library(haven) #packages
library(tidyverse)
library(dplyr)
mydata<-  read_sav("~/Desktop/2024 02 12 Skjellegrind, data 113250/2024-02-12_113250_Data.sav") 
mydata1<- read_sav("~/Desktop/2024 04 23 Skjellegrind, tillegg 113250/2024-04-23_113250_Data.sav")
mydata<-merge(mydata,mydata1)
#remove @
library(janitor)
hpylori2<-clean_names(mydata)
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
  rename(crp=sem_crp_nt2blm) #crp
hpylori2 <- hpylori2 %>% 
  rename(age2=part_ag_nt2blq1) #age HUNT 2
hpylori2 <- hpylori2 %>% 
  rename(age=part_ag_nt4blm)
#factor categorical variables 
hpylori2$helicobacterpylori_cat<-factor(hpylori2$helicobacterpylori_cat,levels = c("0","1"),labels = c("Negative","Positive"))
hpylori2$sex<-factor(hpylori2$sex,levels = c("0","1"),labels = c("Female","Male"))
hpylori2$smoking<-factor(hpylori2$smoking,levels = c("0","1","2"),labels = c("Never","Previously","Currently"))
hpylori2$education<-factor(hpylori2$education,levels = c("1","2","3","4","5"),
                           labels = c("Primary","Secondary","highschool","tertiary","university"))
hpylori2$maritalstatus<-factor(hpylori2$maritalstatus,levels = c("1","2","3","4","5"),
                               labels = c("Unmarried","Married","Widow","Divorced","Separated"))
hpylori2$dementia<-factor(hpylori2$dementia,levels = c("0","1","2","3","4","5","6","7","9"),exclude = "9",
                          labels = c("None","AD","VaD","DLB","FTD","Mixed","Other","Unspecified"))
hpylori2$AD=hpylori2$dementia
hpylori2$cognitiveimpairement<-factor(hpylori2$cognitiveimpairement,levels = c("0","1","2","3","4","9"),exclude = "9",
                                      labels = c("None","aMCI","MCI","Dementia","Other"))
hpylori2$diabetes<-factor(hpylori2$diabetes,levels = c("0","1"),labels = c("No","Yes"))
hpylori2$cvd<-factor(hpylori2$cvd,levels = c("0","1"),labels = c("No","Yes"))
#re-categorize variables 
library(forcats)
hpylori2$education<-fct_collapse(hpylori2$education,"Primary"=c("Primary"),
                                 "Highschool"=c("Secondary","highschool"),"University"=c("tertiary","university")) #education
hpylori2$maritalstatus<-fct_collapse(hpylori2$maritalstatus,"Unmarried"=c("Unmarried"),
                                     "Married"=c("Married"),"Widow_divorced_separated"=c("Widow","Divorced","Separated")) #marital status
hpylori2$dementia<-fct_collapse(hpylori2$dementia,"No Dementia"=c("None"),
                                "Dementia"=c("AD","VaD","DLB","FTD","Mixed","Other","Unspecified")) 
hpylori2$cognitiveimpairement<-fct_collapse(hpylori2$cognitiveimpairement,"No Cognitive Impairment"=c("None"),
                                            "Cognitive Impairment"=c("aMCI","MCI","Dementia","Other")) 
hpylori2$AD<-fct_collapse(hpylori2$AD,"AD"=c("AD"),
                          "No AD"=c("None","VaD","DLB","FTD","Mixed","Other","Unspecified")) 
#participants inclusion 
#born after 1915 and before 1950
birthyear<- hpylori2[hpylori2$birth_year > 1914 & hpylori2$birth_year < 1950, ]
#avialble HP data 
hpylori2<- birthyear %>%filter(!is.na(helicobacterpylori_cat))
#availble HUNT 4 70+
hpylori2<- hpylori2 %>%filter(part_nt4eld1mix==1)
#select variables 
hpylori2 <- hpylori2 %>% select(helicobacterpylori_cat, smoking, sex, crp, bmi, age2, age, education, maritalstatus, alcohol, dementia, cognitiveimpairement, diabetes, cvd, AD)

#table 1
library(tableone)
  
#variable list
listvar<-c("helicobacterpylori_cat","helicobacterpylori_cont","crp","smoking", "sex","bmi","age2","age","education","maritalstatus","alcohol","dementia","E4","cognitiveimpairement","diabetes","cvd","AD","demfreesurv")
catvar <- c("smoking", "sex","education","maritalstatus","dementia","cognitiveimpairement","diabetes","cvd","AD","demfreesurv","E4")

#table 1 overall
table1<-CreateTableOne(data=hpylori2)
print(table1, showAllLevels = TRUE, nonnormal = c("age2","age","crp","helicobacterpylori_cont"),formatOptions = list(big.mark = ","))

#table 1 stratified
table1 <- CreateTableOne(listvar,hpylori2,catvar,strata=c("helicobacterpylori_cat"))
print(table1, showAllLevels = TRUE, nonnormal = c("age2","age","crp","helicobacterpylori_cont"),formatOptions = list(big.mark = ","))

#write.csv(table1,file="hpylori2.csv")
#missing
colSums(is.na(hpylori2))
  