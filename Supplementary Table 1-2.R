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

mydata<-  read_sav("~/Desktop/2024 02 12 Skjellegrind, data 113250/2024-02-12_113250_Data.sav") #all variable
#remove @
library(janitor)
mydata<-clean_names(mydata)
mydata <- mydata %>% rename(PID=pid_113250)
mydata1<- read_sav("~/Desktop/2024 04 23 Skjellegrind, tillegg 113250/2024-04-23_113250_Data.sav") #age at HUNT 2
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

#PGS 
mydata4<-read_sav("~/Desktop/2024 11 06 Skjellegrind, GRS 113250/2024-11-06_113250_Data.sav")
mydata4<-clean_names(mydata4)
mydata4 <- mydata4 %>% rename(PID=pid_113250)
mydata5<-read_sav("~/Desktop/2024 11 28 Skjellegrind, tilleggsvariabel 113250/2024-11-29_113250_Data.sav")
mydata5<-clean_names(mydata5)
mydata5 <- mydata5 %>% rename(PID=pid_113250)
mydata4<- merge(mydata5, mydata4, by = "PID", all = TRUE)

mydata6<- read_sav("~/Desktop/2024 12 04 Skjellegrind, tillegg 113250/2024-12-04_113250_Data.sav")
mydata6<-clean_names(mydata6)
mydata6 <- mydata6 %>% rename(PID=pid_113250)
mydata4<- merge(mydata6, mydata4, by = "PID", all = TRUE)


#renaming variables in PGS dataframe 
mydata4<- mydata4 %>% rename(PGS1=grs002249ad_ave_113861)
mydata4<- mydata4 %>% rename(PGS2=grs002280ad_ave_113861)
mydata4<- mydata4 %>% rename(PGS3=grs004034ad_ave_113861)
mydata4<- mydata4 %>% rename(PGS4=grs004228ad_ave_113861)
mydata4<- mydata4 %>% rename(PGS5=grs004229ad_ave_113861)


#PGS tertiles
#top bottom 10%
t1 <- quantile(mydata4$PGS1, probs = c(0.0, 0.1, 0.9, 1.0), na.rm = TRUE)
mydata4$PGS1_ter <- cut(mydata4$PGS1,t1,include.lowest=TRUE,labels=paste("Group", 1:3)) 
table(mydata4$PGS1_ter)
t2<-quantile(mydata4$PGS2,probs=c(0.0,0.1,0.9,1.0),na.rm = TRUE)
mydata4$PGS2_ter <- cut(mydata4$PGS2, t2, include.lowest=TRUE,labels=paste("Group", 1:3))
t3<-quantile(mydata4$PGS3,probs=c(0.0,0.1,0.9,1.0),na.rm = TRUE)
mydata4$PGS3_ter <- cut(mydata4$PGS3, t3, include.lowest=TRUE,labels=paste("Group", 1:3))
t4<-quantile(mydata4$PGS4,probs=c(0.0,0.1,0.9,1.0),na.rm = TRUE)
mydata4$PGS4_ter <- cut(mydata4$PGS4, t4, include.lowest=TRUE,labels=paste("Group", 1:3))
t5<-quantile(mydata4$PGS5,probs=c(0.0,0.1,0.9,1.0),na.rm = TRUE)
mydata4$PGS5_ter <- cut(mydata4$PGS5, t5, include.lowest=TRUE,labels=paste("Group", 1:3))

#top bottom 5%
# Define quantiles at 0%, 5%, 95%, and 100%
t1 <- quantile(mydata4$PGS1, probs = c(0.0, 0.05, 0.95, 1.0), na.rm = TRUE)
mydata4$PGS1_ter2 <- cut(mydata4$PGS1, t1, include.lowest = TRUE, labels = paste("Group", 1:3))
table(mydata4$PGS1_ter2)
t2<-quantile(mydata4$PGS2,probs=c(0.0, 0.05, 0.95, 1.0),na.rm = TRUE)
mydata4$PGS2_ter2 <- cut(mydata4$PGS2, t2, include.lowest=TRUE,labels=paste("Group", 1:3))
t3<-quantile(mydata4$PGS3,probs=c(0.0, 0.05, 0.95, 1.0),na.rm = TRUE)
mydata4$PGS3_ter2 <- cut(mydata4$PGS3, t3, include.lowest=TRUE,labels=paste("Group", 1:3))
t4<-quantile(mydata4$PGS4,probs=c(0.0, 0.05, 0.95, 1.0),na.rm = TRUE)
mydata4$PGS4_ter2 <- cut(mydata4$PGS4, t4, include.lowest=TRUE,labels=paste("Group", 1:3))
t5<-quantile(mydata4$PGS5,probs=c(0.0, 0.05, 0.95, 1.0),na.rm = TRUE)
mydata4$PGS5_ter2 <- cut(mydata4$PGS5, t5, include.lowest=TRUE,labels=paste("Group", 1:3))

#top bottom 20%
# Define quantiles at 0%, 20%, 80%, and 100%
t1 <- quantile(mydata4$PGS1, probs = c(0.0, 0.2, 0.8, 1.0), na.rm = TRUE)
mydata4$PGS1_ter3 <- cut(mydata4$PGS1, t1, include.lowest = TRUE, labels = paste("Group", 1:3))
table(mydata4$PGS1_ter3)
t2<-quantile(mydata4$PGS2,probs=c(0.0, 0.2, 0.8, 1.0),na.rm = TRUE)
mydata4$PGS2_ter3 <- cut(mydata4$PGS2, t2, include.lowest=TRUE,labels=paste("Group", 1:3))
t3<-quantile(mydata4$PGS3,probs=c(0.0, 0.2, 0.8, 1.0),na.rm = TRUE)
mydata4$PGS3_ter3 <- cut(mydata4$PGS3, t3, include.lowest=TRUE,labels=paste("Group", 1:3))
t4<-quantile(mydata4$PGS4,probs=c(0.0, 0.2, 0.8, 1.0),na.rm = TRUE)
mydata4$PGS4_ter3 <- cut(mydata4$PGS4, t4, include.lowest=TRUE,labels=paste("Group", 1:3))
t5<-quantile(mydata4$PGS5,probs=c(0.0, 0.2, 0.8, 1.0),na.rm = TRUE)
mydata4$PGS5_ter3 <- cut(mydata4$PGS5, t5, include.lowest=TRUE,labels=paste("Group", 1:3))

#merging all datafiles 
# Remove duplicate columns from hpylori2 before merging
hpylori2 <- hpylori2[, !names(hpylori2) %in% c("sex", "birth_year")]

# Now merge without duplicate columns
hpylori2 <- merge(hpylori2, mydata4[, c("PID", "PGS1", "PGS2","PGS3","PGS4","PGS5", "PGS1_ter", "PGS2_ter","PGS1_ter2", "PGS2_ter2","PGS1_ter3","PGS2_ter3","PGS3_ter","PGS3_ter2","PGS3_ter3","PGS4_ter","PGS4_ter2","PGS4_ter3","PGS5_ter","PGS5_ter2","PGS5_ter3")], 
                  by = "PID", all = TRUE)

summary(hpylori2$PGS1)

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
# Filter dataset based on birth year
birthyear <- hpylori2 %>% filter(birth_year.x > 1914 & birth_year.x < 1950)

# Filter for available HP data
hpylori2 <- birthyear %>% filter(!is.na(helicobacterpylori_cat))

# Exclude those who moved based on register status
hpylori2 <- hpylori2 %>% filter(register != 3) 

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
table(hpylori2$PGS1_ter)
table(hpylori2$PGS1_ter2)

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
#Complete cases
hpylori2<-hpylori2 %>% filter(complete.cases(sex,age,education,demfreesurv,helicobacterpylori_cat,helicobacterpylori_cont))
##############################################
#Supplementary table 1: Logistic Regression###
##############################################
library(sjPlot)
#Dementia
m0<- glm(dementia ~ helicobacterpylori_cat, data = hpylori2, family = binomial) 
m1<-glm(dementia~helicobacterpylori_cat+sex+age,family = binomial,data = hpylori2) 
m2<-glm(dementia~helicobacterpylori_cat+sex+age+education,family = binomial,data = hpylori2) 
#m3<-glm(dementia~helicobacterpylori_cat+sex+age+education+maritalstatus+smoking+bmi,family = binomial,data = hpylori2) 
#m4<-glm(dementia~helicobacterpylori_cat+sex+age+education+maritalstatus+smoking+bmi+diabetes+cvd+alcohol,family = binomial,data = hpylori2) 
tab_model(m0,m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Crude", "Model 1","Model 2"),show.intercept = FALSE,string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Sex (Male)","Age (Years)","Education (Highschool)", "Education (University)"))

#CI
m0<- glm(cognitiveimpairement ~ helicobacterpylori_cat, data = hpylori2, family = binomial) 
m1<-glm(cognitiveimpairement~helicobacterpylori_cat+sex+age,family = binomial,data = hpylori2) 
m2<-glm(cognitiveimpairement~helicobacterpylori_cat+sex+age+education,family = binomial,data = hpylori2) 
#m3<-glm(cognitiveimpairement~helicobacterpylori_cat+sex+age+education+maritalstatus+smoking+bmi,family = binomial,data = hpylori2) 
#m4<-glm(cognitiveimpairement~helicobacterpylori_cat+sex+age+education+maritalstatus+smoking+bmi+diabetes+cvd+alcohol,family = binomial,data = hpylori2) 
tab_model(m0,m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Crude", "Model 1","Model 2"),show.intercept = FALSE,string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and CI",
          pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Sex (Male)","Age (Years)","Education (Highschool)", "Education (University)"))

#AD
m0<- glm(AD ~ helicobacterpylori_cat, data = hpylori2, family = binomial) 
m1<-glm(AD~helicobacterpylori_cat+sex+age,family = binomial,data = hpylori2) 
m2<-glm(AD~helicobacterpylori_cat+sex+age+education,family = binomial,data = hpylori2) 
#m3<-glm(AD~helicobacterpylori_cat+sex+age+education+maritalstatus+smoking+bmi,family = binomial,data = hpylori2) 
#m4<-glm(AD~helicobacterpylori_cat+sex+age+education+maritalstatus+smoking+bmi+diabetes+cvd+alcohol,family = binomial,data = hpylori2) 
tab_model(m0,m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Crude", "Model 1","Model 2"),show.intercept = FALSE,string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Sex (Male)","Age (Years)","Education (Highschool)", "Education (University)"))
#DFS
m0<- glm(demfreesurv ~ helicobacterpylori_cat, data = hpylori2, family = binomial) 
m1<-glm(demfreesurv~helicobacterpylori_cat+sex+age,family = binomial,data = hpylori2) 
m2<-glm(demfreesurv~helicobacterpylori_cat+sex+age+education,family = binomial,data = hpylori2) 
#m3<-glm(demfreesurv~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi,family = binomial,data = hpylori2) 
#m4<-glm(demfreesurv~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi+diabetes+cvd+alcohol,family = binomial,data = hpylori2) 
tab_model(m0,m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Crude", "Model 1","Model 2"),show.intercept = FALSE,string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia Free Survival",
          pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Sex (Male)","Age (Years)","Education (Highschool)", "Education (University)"))
#HP titers 
#normalize 
library(RNOmni)
hpylori2$helicobacterpylori_cont<-RankNorm(as.numeric(hpylori2$helicobacterpylori_cont))
#Dementia
m0<- glm(dementia ~ helicobacterpylori_cont, data = hpylori2, family = binomial) 
m1<-glm(dementia~helicobacterpylori_cont+sex+age,family = binomial,data = hpylori2) 
m2<-glm(dementia~helicobacterpylori_cont+sex+age+education,family = binomial,data = hpylori2) 
#m3<-glm(dementia~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi,family = binomial,data = hpylori2) 
#m4<-glm(dementia~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi+diabetes+cvd+alcohol,family = binomial,data = hpylori2) 
tab_model(m0,m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Crude", "Model 1","Model 2"),show.intercept = FALSE,string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          pred.labels = c("Serum Helicobacter pylori (Titers)","Sex (Male)","Age (Years)","Education (Highschool)", "Education (University)"))


#CI
m0<- glm(cognitiveimpairement ~ helicobacterpylori_cont, data = hpylori2, family = binomial) 
m1<-glm(cognitiveimpairement~helicobacterpylori_cont+sex+age,family = binomial,data = hpylori2) 
m2<-glm(cognitiveimpairement~helicobacterpylori_cont+sex+age+education,family = binomial,data = hpylori2) 
#m3<-glm(cognitiveimpairement~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi,family = binomial,data = hpylori2) 
#m4<-glm(cognitiveimpairement~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi+diabetes+cvd+alcohol,family = binomial,data = hpylori2) 
tab_model(m0,m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Crude", "Model 1","Model 2"),show.intercept = FALSE,string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and CI",
          pred.labels = c("Serum Helicobacter pylori (Titers)","Sex (Male)","Age (Years)","Education (Highschool)", "Education (University)"))


#AD
m0<- glm(AD ~ helicobacterpylori_cont, data = hpylori2, family = binomial) 
m1<-glm(AD~helicobacterpylori_cont+sex+age,family = binomial,data = hpylori2) 
m2<-glm(AD~helicobacterpylori_cont+sex+age+education,family = binomial,data = hpylori2) 
#m3<-glm(AD~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi,family = binomial,data = hpylori2) 
#m4<-glm(AD~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi+diabetes+cvd+alcohol,family = binomial,data = hpylori2) 
tab_model(m0,m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Crude", "Model 1","Model 2"),show.intercept = FALSE,string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          pred.labels = c("Serum Helicobacter pylori (Titers)","Sex (Male)","Age (Years)","Education (Highschool)", "Education (University)"))
m0<- glm(demfreesurv~ helicobacterpylori_cont, data = hpylori2, family = binomial) 
m1<-glm(demfreesurv~helicobacterpylori_cont+sex+age,family = binomial,data = hpylori2) 
m2<-glm(demfreesurv~helicobacterpylori_cont+sex+age+education,family = binomial,data = hpylori2) 
#m3<-glm(demfreesurv~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi,family = binomial,data = hpylori2) 
#m4<-glm(demfreesurv~helicobacterpylori_cont+sex+age+education+maritalstatus+smoking+bmi+diabetes+cvd+alcohol,family = binomial,data = hpylori2) 
tab_model(m0,m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Crude", "Model 1","Model 2"),show.intercept = FALSE,string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia Free Survival",
          pred.labels = c("Serum Helicobacter pylori (Titers)","Sex (Male)","Age (Years)","Education (Highschool)", "Education (University)"))
#########################################
#Supplementary Table 2: Stratifications##
#########################################
##########Seropositivity##############
# Sex dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age,
         family = binomial, data = subset(hpylori2, sex == "Male"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age,
         family = binomial, data = subset(hpylori2, sex == "Female"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Male","Female"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)"
          ))
#sex CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age,
         family = binomial, data = subset(hpylori2, sex == "Male"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age,
         family = binomial, data = subset(hpylori2, sex == "Female"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Male","Female"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)"))

#sex AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age,
         family = binomial, data = subset(hpylori2, sex == "Male"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age,
         family = binomial, data = subset(hpylori2, sex == "Female"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Male","Female"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)"))
#sex DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age + maritalstatus + smoking + bmi + diabetes + cvd + alcohol,
         family = binomial, data = subset(hpylori2, sex == "Male"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age + maritalstatus + smoking + bmi + diabetes + cvd + alcohol,
         family = binomial, data = subset(hpylori2, sex == "Female"))
tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Male","Female"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)","Marital Status (Married)","MaritalStatus (Previously Married)",
                                                  "Smoking (Previously)", "Smoking (Currently)","BMI (kg/m2)","Diabetes","CVD","Alcohol (cl)"))

#Education
#edu dem
m1<- glm(dementia~ helicobacterpylori_cat +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Primary"))
m2<- glm(dementia~ helicobacterpylori_cat +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Highschool"))
m3<- glm(dementia~ helicobacterpylori_cat +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "University"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Primary","Highschool","University"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Sex(Male)","Age (Years)"))

#Edu cognitive impairment 
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Primary"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Highschool"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "University"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Primary","Highschool","University"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Sex(Male)","Age (Years)"))
#Edu AD 
m1<- glm(AD~ helicobacterpylori_cat +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Primary"))
m2<- glm(AD~ helicobacterpylori_cat +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Highschool"))
m3<- glm(AD~ helicobacterpylori_cat +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "University"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Primary","Highschool","University"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Sex(Male)","Age (Years)"))
#Edu demfreesurv 
m1<- glm(demfreesurv~ helicobacterpylori_cat +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Primary"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Highschool"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "University"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Primary","Highschool","University"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Sex(Male)","Age (Years)"))
#Age categorized
median_age <- median(hpylori2$age2)
median_age
hpylori2$age2<- ifelse(hpylori2$age2 > median_age, 1, 0)
hpylori2$age2<-factor(hpylori2$age2, levels = c(0, 1), labels = c("Below55", "Above55"))
table(hpylori2$age2)

#age dementia 
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           sex,
         family = binomial, data = subset(hpylori2, age2 == "Above55"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           sex,
         family = binomial, data = subset(hpylori2, age2 == "Below55"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above55","Below55"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Sex (Male)"))
#age CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           sex,
         family = binomial, data = subset(hpylori2, age2 == "Above55"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           sex,
         family = binomial, data = subset(hpylori2, age2 == "Below55"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above55","Below55"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Sex (Male)"))
#AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           sex,
         family = binomial, data = subset(hpylori2, age2 == "Above55"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           sex,
         family = binomial, data = subset(hpylori2, age2 == "Below55"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above55","Below55"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Sex (Male)"))
#DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ sex + maritalstatus + smoking + bmi + diabetes + cvd + alcohol,family = binomial, data = subset(hpylori2, age2 == "Above55"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ sex + maritalstatus + smoking + bmi + diabetes + cvd + alcohol,family = binomial, data = subset(hpylori2, age2 == "Below55"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above55","Below55"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Sex (Male)","Marital Status (Married)","MaritalStatus (Previously Married)",
                                                  "Smoking (Previously)", "Smoking (Currently)","BMI (kg/m2)","Diabetes","CVD","Alcohol (cl)"))
#CRP categorized
hpylori2$crp2<- ifelse(hpylori2$crp > 5,1,0)
hpylori2$crp2<-factor(hpylori2$crp2,levels = c("0","1"),labels = c("Below5","Above5"))
table(hpylori2$crp2)
#CRP Seropositivity
#CRP dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2 == "Above5"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2 == "Below5"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above 5","Below 5"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#CRP CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2 == "Above5"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat  +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2 == "Below5"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above","Below"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))

#CRP AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2 == "Above5"))
m2<- glm(AD~ helicobacterpylori_cat  +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2 == "Below5"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above","Below"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#CRP DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2 == "Above5"))
m2<- glm(demfreesurv~ helicobacterpylori_cat  +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2 == "Below5"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above","Below"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#E4 seropositivity
#E4 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "0"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "1"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Neg","Pos"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#E4 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "0"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat  +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "1"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Neg","Pos"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#E4 AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "0"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "1"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Neg","Pos"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#E4 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "0"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "1"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Neg","Pos"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#E4 seropositivity
#E4 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "0"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "1"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Neg","Pos"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#E4 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "0"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat  +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "1"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Neg","Pos"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#E4 AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "0"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "1"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Neg","Pos"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#############################
############PGS002249###################
#PGS1 seropositivity- 10%####
#############################
#PGS1 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))

m1<- glm(dementia~ helicobacterpylori_cat*PGS1+education+ 
           age+sex,
         family = binomial, data = hpylori2)
summary(m1)
m2<- glm(dementia~ helicobacterpylori_cat*PGS2+education+ 
           age+sex,
         family = binomial, data = hpylori2)
summary(m2)
m3<- glm(dementia~ helicobacterpylori_cont*PGS1+education+ 
           age+sex,
         family = binomial, data = hpylori2)
summary(m3)
m4<- glm(dementia~ helicobacterpylori_cont*PGS2+education+ 
           age+sex,
         family = binomial, data = hpylori2)
summary(m4)
#PGS1 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#########################
#80
#PGS2 Seropositivity 10%
#PGS2 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))

#PGS1 5%
#PGS1 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
plot(hpylori2$helicobacterpylori_cont,m3$fitted.values)
summary(m3)

predict(m3,hpylori2)
plot(hpylori2$helicobacterpylori_cont,predict(m3,hpylori2))

table(hpylori2$helicobacterpylori_cat,hpylori2$PGS1_ter,hpylori2$cognitiveimpairement)
#PGS1 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 Seropositivity
#PGS2 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 Seropositivity 20%
#PGS1 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))

#PGS1 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
############################
#PGS3=PGS004034#############
#Seropositivity PGS 3 20%
##############################
#PGS3 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))

#PGS3 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#############################
###PGS3 seropositivity 10%###
##############################
#PGS3 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))


#PGS3 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#############################
###PGS3 seropositivity 5%###
##############################
#PGS3 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))


#PGS3 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))

############################
#PGS4=PGS004228#############
#Seropositivity PGS 4 20%
##############################
#PGS4 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter3== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter3== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS4 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter3== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter3== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))

#PGS4 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter3== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter3== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS4 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter3== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter3== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#############################
###PGS4 seropositivity 10%###
##############################
#PGS4 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))


#PGS4 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS4 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS4 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#############################
###PGS4 seropositivity 5%###
##############################
#PGS3 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter2== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter2== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))


#PGS4 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter2== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter2== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS4 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter2== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter2== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS4 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter2== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter2== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS4_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
############################
#PGS5=PGS004229 (NO APOE)#############
#Seropositivity PGS 4 20%
##############################
#PGS5 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter3== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter3== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS5 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter3== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter3== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))

#PGS5 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter3== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter3== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS5 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter3== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter3== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#############################
###PGS5 seropositivity 10%###
##############################
#PGS5 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))


#PGS5 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS4 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS5 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#############################
###PGS5 seropositivity 5%###
##############################
#PGS5 dementia
m1<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter2== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter2== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))


#PGS5 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter2== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter2== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS5 AD
m1<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter2== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter2== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS5 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter2== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter2== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS5_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#############################
###HP Titers#################
#############################
# Sex dementia
m1<- glm(dementia~ helicobacterpylori_cont +education+ 
           age,
         family = binomial, data = subset(hpylori2, sex == "Male"))
m2<- glm(dementia~ helicobacterpylori_cont +education+ 
           age,
         family = binomial, data = subset(hpylori2, sex == "Female"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Male","Female"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Age (Years)"
          ))
#sex CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age,
         family = binomial, data = subset(hpylori2, sex == "Male"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age,
         family = binomial, data = subset(hpylori2, sex == "Female"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Male","Female"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Age (Years)"))

#sex AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age,
         family = binomial, data = subset(hpylori2, sex == "Male"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age,
         family = binomial, data = subset(hpylori2, sex == "Female"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Male","Female"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Age (Years)"))
#sex DFS
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age + maritalstatus + smoking + bmi + diabetes + cvd + alcohol,
         family = binomial, data = subset(hpylori2, sex == "Male"))
m2<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age + maritalstatus + smoking + bmi + diabetes + cvd + alcohol,
         family = binomial, data = subset(hpylori2, sex == "Female"))
tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Male","Female"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)","Marital Status (Married)","MaritalStatus (Previously Married)",
                                                  "Smoking (Previously)", "Smoking (Currently)","BMI (kg/m2)","Diabetes","CVD","Alcohol (cl)"))
#Education
#edu dem
m1<- glm(dementia~ helicobacterpylori_cont +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Primary"))
m2<- glm(dementia~ helicobacterpylori_cont +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Highschool"))
m3<- glm(dementia~ helicobacterpylori_cont +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "University"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Primary","Highschool","University"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Sex(Male)","Age (Years)"))

#Edu cognitive impairment 
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Primary"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Highschool"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cont +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "University"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Primary","Highschool","University"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Sex(Male)","Age (Years)"))
#Edu AD 
m1<- glm(AD~ helicobacterpylori_cont +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Primary"))
m2<- glm(AD~ helicobacterpylori_cont +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Highschool"))
m3<- glm(AD~ helicobacterpylori_cont +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "University"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Primary","Highschool","University"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Sex(Male)","Age (Years)"))
#Edu demfreesurv 
m1<- glm(demfreesurv~ helicobacterpylori_cont +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Primary"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "Highschool"))
m3<- glm(demfreesurv~ helicobacterpylori_cont +sex+ 
           age,
         family = binomial, data = subset(hpylori2, education == "University"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Primary","Highschool","University"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Sex(Male)","Age (Years)"))
###age dementia#### 
m1<- glm(dementia~ helicobacterpylori_cont+education+ 
           sex,
         family = binomial, data = subset(hpylori2, age2 == "Above55"))
m2<- glm(dementia~ helicobacterpylori_cont +education+ 
           sex,
         family = binomial, data = subset(hpylori2, age2 == "Below55"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above55","Below55"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Sex (Male)"))
#age CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           sex,
         family = binomial, data = subset(hpylori2, age2 == "Above55"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           sex,
         family = binomial, data = subset(hpylori2, age2 == "Below55"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above55","Below55"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Sex (Male)"))
# Age AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           sex,
         family = binomial, data = subset(hpylori2, age2 == "Above55"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           sex,
         family = binomial, data = subset(hpylori2, age2 == "Below55"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above55","Below55"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Sex (Male)"))

#Age DFS
m1<- glm(demfreesurv~ helicobacterpylori_cont +education+ sex + maritalstatus + smoking + bmi + diabetes + cvd + alcohol,family = binomial, data = subset(hpylori2, age2 == "Above55"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +education+ sex + maritalstatus + smoking + bmi + diabetes + cvd + alcohol,family = binomial, data = subset(hpylori2, age2 == "Below55"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above55","Below55"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Sex (Male)","Marital Status (Married)","MaritalStatus (Previously Married)",
                                                  "Smoking (Previously)", "Smoking (Currently)","BMI (kg/m2)","Diabetes","CVD","Alcohol (cl)"))
#E4 titers 
#E4 dementia
m1<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "0"))
m2<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "1"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Neg","Pos"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#E4 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "0"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont+education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "1"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Neg","Pos"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#E4 AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "0"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "1"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Neg","Pos"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
m1<- glm(demfreesurv~ helicobacterpylori_cat +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "0"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, E4 == "1"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Neg","Pos"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#CRP titers 
#CRP Dementia
m1<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2== "Above5"))
m2<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2 == "Below5"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above","Below"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#CRP CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2== "Above5"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2 == "Below5"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above","Below"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#CRP AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2== "Above5"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2 == "Below5"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above","Below"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#CRP DFS
m1<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2== "Above5"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, crp2 == "Below5"))

tab_model(m1,m2,p.style = "scientific",digits.p = 2,dv.labels = c("Above","Below"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
##################
#PGS1 titers 10%##
##################
#PGS1 dementia
m1<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 titers
#PGS2 dementia
m1<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#################
#PGS1 titers 5%##
##################
#PGS1 dementia
m1<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 titers
#PGS2 dementia
m1<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS2 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS2_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
##################
#PGS1 titers 20%##
##################
#PGS1 dementia
m1<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS1 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS1_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
##################
#PGS3 titers 20%##
##################
#PGS1 dementia
m1<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter3== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#############################
###PGS3 titers 10%###
##############################
#PGS3 dementia
m1<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))


#PGS3 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (Seropositivity)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#############################
###PGS3 titers 5%###
##############################
#PGS3 dementia
m1<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 1"))
m2<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 2"))
m3<- glm(dementia~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and Dementia",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))


#PGS3 CI
m1<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 1"))
m2<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 2"))
m3<- glm(cognitiveimpairement~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and cognitiveimpairement",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 AD
m1<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 1"))
m2<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 2"))
m3<- glm(AD~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and AD",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))
#PGS3 DFS
m1<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 1"))
m2<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 2"))
m3<- glm(demfreesurv~ helicobacterpylori_cont +education+ 
           age+sex,
         family = binomial, data = subset(hpylori2, PGS3_ter2== "Group 3"))

tab_model(m1,m2,m3,p.style = "scientific",digits.p = 2,dv.labels = c("Group1","Group2","Group=3"),string.pred = "Coeffcients",string.ci = "CI (95%)",string.p = "P-Value", title = "HP-2 and demfreesurv",
          show.intercept = FALSE, pred.labels = c("Serum Helicobacter pylori (titers)","Education (Highschool)", "Education (University)","Age (Years)", "Sex(Male)"
          ))




