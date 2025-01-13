#linear regression HP and CRP
library(Hmisc)
library(foreign)
library(haven)
library(tidyverse)
library(ggplot2)
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

#factor categorical variables 
hpylori2$helicobacterpylori_cat<-factor(hpylori2$helicobacterpylori_cat,levels = c("0","1"),labels = c("Negative","Positive"))
hpylori2$sex<-factor(hpylori2$sex.x,levels = c("0","1"),labels = c("Female","Male"))
hpylori2$smoking<-factor(hpylori2$smoking,levels = c("0","1","2"),labels = c("Never","Previously","Currently"))
hpylori2$education<-factor(hpylori2$education,levels = c("1","2","3","4","5"),
                           labels = c("Primary","Secondary","highschool","tertiary","university"))
hpylori2$maritalstatus<-factor(hpylori2$maritalstatus,levels = c("1","2","3","4","5"),
                               labels = c("Unmarried","Married","Widow","Divorced","Separated"))
#re-categorize variables 
library(forcats)
hpylori2$education<-fct_collapse(hpylori2$education,"Primary"=c("Primary"),
                                 "Highschool"=c("Secondary","highschool"),"University"=c("tertiary","university")) #education
hpylori2$maritalstatus<-fct_collapse(hpylori2$maritalstatus,"Unmarried"=c("Unmarried"),
                                     "Married"=c("Married"),"Widow_divorced_separated"=c("Widow","Divorced","Separated")) #marital status


#Complete cases (run them separately)
hpylori2<-hpylori2 %>% filter(complete.cases(crp,sex,age,education,maritalstatus,cvd,diabetes,smoking,alcohol,bmi))
hpylori2<-hpylori2 %>% filter(complete.cases(crp4,sex,age,education,maritalstatus,cvd,diabetes,smoking,alcohol,bmi))
library(RNOmni)
hpylori2$crp4<-RankNorm(as.numeric(hpylori2$crp4))
hpylori2$crp<-RankNorm(as.numeric(hpylori2$crp))
hpylori2$helicobacterpylori_cont<-RankNorm(as.numeric(hpylori2$helicobacterpylori_cont))
hist(hpylori2$crp4)

#scatter plot
hpylori2 %>%
  ggplot(aes(x = helicobacterpylori_cont, y = crp4)) +
  geom_point(aes(color = "red")) +  # All points in red for helicobacterpylori_cont
  geom_point(aes(x = crp4, y = crp4, color = "blue")) +  # All points in blue for crp
  scale_color_identity(guide = "legend", labels = c("Helicobacter Pylori", "CRP")) +
  labs(x = "Helicobacter Pylori (Cont)", y = "CRP", title = "Scatter Plot of Helicobacter Pylori vs. CRP") +
  theme_minimal()
#correlation 
cor.test(hpylori2$helicobacterpylori_cont,hpylori2$crp4)
cor.test(hpylori2$helicobacterpylori_cont,hpylori2$crp)

library(sjPlot)
#linear regression 
#CRP HUNT 2 (not included in article)
#HP~ CRP Seropositivity
m0<- lm(crp~helicobacterpylori_cat, data = hpylori2, family = binomial) 
m1<-lm(crp~helicobacterpylori_cat+sex+age,family = binomial,data = hpylori2) 
m2<-lm(crp~helicobacterpylori_cat+sex+age+education,family = binomial,data = hpylori2) 
m3<-lm(crp~helicobacterpylori_cat+sex+age+education+maritalstatus+smoking+bmi,family = binomial,data = hpylori2) 
m4<-lm(crp~helicobacterpylori_cat+sex+age+education+maritalstatus+smoking+bmi+diabetes+cvd+alcohol,family = binomial,data = hpylori2) 
tab_model(m0,m1,m2,m3,m4,p.style = "scientific",digits.p = 2,dv.labels = c("Crude", "Model 1","Model 2","Model 3","Model 4"),show.intercept = FALSE,string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 Seropositivity and CRP",
          pred.labels = c("HP (seropositivity","Sex (Male)","Age (Years)","Education (Highschool)", "Education (University)","Marital Status (Married)","MaritalStatus (Previously Married)",
                          "Smoking (Previously)", "Smoking (Currently)","BMI (kg/m2)","Diabetes","CVD","Alcohol (cl)"))
#CRP~HP Titers
m0<- lm(crp~helicobacterpylori_cont, data = hpylori2) 
m1<-lm(crp~helicobacterpylori_cont+sex+age,data = hpylori2) 
m2<-glm(crp~helicobacterpylori_cont+sex+age+education,data = hpylori2) 
m3<-glm(crp~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi,data = hpylori2) 
m4<-glm(crp~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi+diabetes+cvd+alcohol,data = hpylori2) 
tab_model(m0,m1,m2,m3,m4,p.style = "scientific",digits.p = 2,dv.labels = c("Crude", "Model 1","Model 2","Model 3","Model 4"),show.intercept = FALSE,string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "CRP and HP Titers",
          pred.labels = c("HP (titers)","Sex (Male)","Age (Years)","Education (Highschool)", "Education (University)","Marital Status (Married)","MaritalStatus (Previously Married)",
                          "Smoking (Previously)", "Smoking (Currently)","BMI (kg/m2)","Diabetes","CVD","Alcohol (cl)"))
#CRP HUNT 4
#HP~ CRP Seropositivity
m0<- lm(crp4~helicobacterpylori_cat, data = hpylori2) 
m1<-lm(crp4~helicobacterpylori_cat+sex+age,data = hpylori2) 
m2<-lm(crp4~helicobacterpylori_cat+sex+age+education,data = hpylori2) 
m3<-lm(crp4~helicobacterpylori_cat+sex+age+education+maritalstatus+smoking+bmi,data = hpylori2) 
m4<-lm(crp4~helicobacterpylori_cat+sex+age+education+maritalstatus+smoking+bmi+diabetes+cvd+alcohol,data = hpylori2) 
tab_model(m0,m1,m2,m3,m4,p.style = "scientific",digits.p = 2,dv.labels = c("Crude", "Model 1","Model 2","Model 3","Model 4"),show.intercept = FALSE,string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 Seropositivity and CRP",
          pred.labels = c("HP (seropositivity","Sex (Male)","Age (Years)","Education (Highschool)", "Education (University)","Marital Status (Married)","MaritalStatus (Previously Married)",
                          "Smoking (Previously)", "Smoking (Currently)","BMI (kg/m2)","Diabetes","CVD","Alcohol (cl)"))
#CRP~HP Titers
m0<- lm(crp4~helicobacterpylori_cont, data = hpylori2) 
m1<-lm(crp4~helicobacterpylori_cont+sex+age,data = hpylori2) 
m2<-glm(crp4~helicobacterpylori_cont+sex+age+education,data = hpylori2) 
m3<-glm(crp4~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi,data = hpylori2) 
m4<-glm(crp4~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi+diabetes+cvd+alcohol,data = hpylori2) 
tab_model(m0,m1,m2,m3,m4,p.style = "scientific",digits.p = 2,dv.labels = c("Crude", "Model 1","Model 2","Model 3","Model 4"),show.intercept = FALSE,string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "CRP and HP Titers",
          pred.labels = c("HP (titers)","Sex (Male)","Age (Years)","Education (Highschool)", "Education (University)","Marital Status (Married)","MaritalStatus (Previously Married)",
                          "Smoking (Previously)", "Smoking (Currently)","BMI (kg/m2)","Diabetes","CVD","Alcohol (cl)"))


