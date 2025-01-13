library(Hmisc)
library(foreign)
library(haven)
library(tidyverse)
Mydata  <- read.delim("~/code/Hpylori/transfer_338598_files_6dfeac78/dose_PID113250_varSubset.txt")
## Håvard wrote this script to recode SNP data to ApoE genotype
Mydata$rs429358 <- Mydata[,2]
Mydata$rs7412 <- Mydata[,3]

ApoE <- data.frame(Mydata[1],Mydata$rs429358,Mydata$rs7412)

names(ApoE)[names(ApoE) == "Mydata.rs429358"] <- "rs429358"
names(ApoE)[names(ApoE) == "Mydata.rs7412"] <- "rs7412"

ApoE$r_rs429358 <- round(ApoE$rs429358)
ApoE$r_rs7412 <- round(ApoE$rs7412)

## recodes the rs7412 (a) and rs429358 (b) values to ApoE genotype
ApoErecode <- function(a,b) {
  if (a == 2 &  b == 0)     { return(1)  }
  else if (a == 1 & b == 0) { return(2)  }
  else if (a == 1 & b == 1) { return(3)  }
  else if (a == 0 & b == 0) { return(4)  }
  else if (a == 0 & b == 1) { return(5)  }
  else if (a == 0 & b == 2) { return(6)  }
  else if (a == 2 & b == 2) { return(7)  }
  else if (a == 2 & b == 1) { return(8)  }
  else if (a == 1 & b == 2) { return(9)  }
}

for (i in 1:length(ApoE$r_rs429358)) {
  ApoE$gtype_[i] <- ApoErecode(ApoE$r_rs7412[i],ApoE$r_rs429358[i])
}
ApoE$gtype__ <- factor(ApoE$gtype_, levels=c(1:9), labels=c("ApoE ε2/ε2", "ApoE ε2/ε3","ApoE ε2/ε4","ApoE ε3/ε3","ApoE ε3/ε4", "ApoE ε4/ε4", "ApoE ε1/ε1", "ApoE ε1/ε2", "ApoE ε1/ε4"))
ApoE$gtype <- factor(ApoE$gtype_, levels=c(1:6), labels=c("ApoE ε2/ε2", "ApoE ε2/ε3","ApoE ε2/ε4","ApoE ε3/ε3","ApoE ε3/ε4", "ApoE ε4/ε4"))

print(table(ApoE$gtype))

## Recodes a variable for E4 positivity
ApoE$E4 <- ifelse(ApoE$gtype_==1 | ApoE$gtype_==2 | ApoE$gtype_==4, 0, ifelse(ApoE$gtype_==3 | ApoE$gtype_==5 | ApoE$gtype_==6, 1, NA))

#Pieta save E4 as factor
ApoE$E4 <- as.factor(ApoE$E4)
#select the E4 variable in data frame 
mydata2<- ApoE %>% select(PID,E4)

mydata<-  read_sav("~/Desktop/2024 02 12 Skjellegrind, data 113250/2024-02-12_113250_Data.sav") 
#remove @
library(janitor)
mydata<-clean_names(mydata)
mydata <- mydata %>% rename(PID=pid_113250)
mydata1<- read_sav("~/Desktop/2024 04 23 Skjellegrind, tillegg 113250/2024-04-23_113250_Data.sav")
mydata1<-clean_names(mydata1)
mydata1 <- mydata1 %>% rename(PID=pid_113250)
#merge
merged_data <- merge(mydata, mydata1, by = "PID", all = TRUE)
hpylori2 <- merge(merged_data, mydata2, by = "PID", all = TRUE)
#regstat variable
mydata3<- read_sav("~/Desktop/2024 07 12 Skjellegrind, registerstatus 113250/2024-07-12_113250_Data.sav")
mydata3<-clean_names(mydata3)
mydata3 <- mydata3 %>% rename(PID=pid_113250)
hpylori2 <- merge(hpylori2, mydata3, by = "PID", all = TRUE)
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
hpylori2 <- hpylori2 %>% 
  rename(crp4=sem_crp_nt4blm) #crp HUNT 4
#DFS 
# Filter dataset based on birth year
birthyear <- hpylori2 %>% filter(birth_year.x > 1914 & birth_year.x < 1950)

# Filter for available HP data
hpylori2 <- birthyear %>% filter(!is.na(helicobacterpylori_cat))

# Exclude those who moved based on register status
#hpylori2 <- hpylori2 %>% filter(register != 3) 


#DFS
#invited to hunt 4
hpylori2$invhunt4 <- ifelse(hpylori2$part_nt4blm %in% c(0, 1), 1, NA)
table(hpylori2$invhunt4,useNA = "ifany") #invited n=2718
hpylori2$register <- ifelse(hpylori2$register == 5 & hpylori2$invhunt4 == 1, 1,hpylori2$register)
hpylori2$register <- ifelse(is.na(hpylori2$register), 0, hpylori2$register) #0=dead 
table(hpylori2$register, useNA = "ifany")

#factor categorical variables 
hpylori2$helicobacterpylori_cat<-factor(hpylori2$helicobacterpylori_cat,levels = c("0","1"),labels = c("Negative","Positive"))
hpylori2$sex<-factor(hpylori2$sex.x,levels = c("0","1"),labels = c("Female","Male"))
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
hpylori2$register<-factor(hpylori2$register,levels = c("0","1"),labels = c("Dead","Alive"))

# Initialize demfreesurv with NA for all rows
hpylori2$demfreesurv <- NA
# Set demfreesurv to FALSE for those with cognitive impairment of "Dementia"
hpylori2$demfreesurv[hpylori2$cognitiveimpairement == "Dementia"] <- FALSE
# Set demfreesurv to FALSE for those who are dead
hpylori2$demfreesurv[hpylori2$register == "Dead"] <- FALSE
# Set demfreesurv to TRUE for those with no cognitive impairment
hpylori2$demfreesurv[hpylori2$dementia == "None"] <- TRUE
# Convert demfreesurv to a factor
hpylori2$demfreesurv <- as.factor(hpylori2$demfreesurv)
# Print the count of each level in demfreesurv
table(hpylori2$demfreesurv)

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
                          "No AD"=c("None","VaD","DLB","FTD","Other","Unspecified","Mixed")) 
#select variables 
hpylori2 <- hpylori2 %>% select(helicobacterpylori_cat,helicobacterpylori_cont,smoking, sex, crp, crp4,bmi, age2, age, education, maritalstatus, alcohol, dementia, E4,cognitiveimpairement, diabetes, cvd, AD,demfreesurv)

#table 1
library(tableone)

#variable list
listvar<-c("helicobacterpylori_cat","helicobacterpylori_cont","crp","crp4","smoking", "sex","bmi","age2","age","education","maritalstatus","alcohol","dementia","E4","cognitiveimpairement","diabetes","cvd","AD","demfreesurv")
catvar <- c("smoking", "sex","education","maritalstatus","dementia","cognitiveimpairement","diabetes","cvd","AD","demfreesurv","E4")

#table 1 overall
table1<-CreateTableOne(data=hpylori2)
print(table1, showAllLevels = TRUE, nonnormal = c("age","age2","crp","crp4","helicobacterpylori_cont"),formatOptions = list(big.mark = ","))

#table 1 stratified
table1 <- CreateTableOne(listvar,hpylori2,catvar,strata=c("helicobacterpylori_cat"))
print(table1, showAllLevels = TRUE, nonnormal = c("age","age2","crp","crp4","helicobacterpylori_cont"),formatOptions = list(big.mark = ","))

#write.csv(table1,file="hpylori2.csv")

#missing
colSums(is.na(hpylori2))

