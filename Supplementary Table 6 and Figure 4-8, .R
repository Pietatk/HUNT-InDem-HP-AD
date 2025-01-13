###age on a time scale
#####################
##create data file####
######################
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

#select birth month variable
mydata1<-mydata1[, c("pid_113250","birth_month")]
hpylori2 <- merge(hpylori2,mydata1,by = "pid_113250", all = TRUE)

#APOE data 
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
mydata4<- ApoE %>% select(PID,E4)
mydata4 <- mydata4 %>% rename(pid_113250=PID)
hpylori2 <- merge(hpylori2,mydata4,by = "pid_113250", all = TRUE) #merge
hpylori2 <- hpylori2 %>% rename(PID=pid_113250)#%>%mutate(PID=as.character(PID))
hpylori2%>%filter(PID==1132500005366) 

#PGS 
mydata4<-read_sav("~/Desktop/2024 11 06 Skjellegrind, GRS 113250/2024-11-06_113250_Data.sav")
mydata4<-clean_names(mydata4)
mydata4 <- mydata4 %>% rename(PID=pid_113250)#%>%mutate(PID=as.character(PID)) 
head(mydata4)
mydata5<-read_sav("~/Desktop/2024 11 28 Skjellegrind, tilleggsvariabel 113250/2024-11-29_113250_Data.sav")
mydata5<-clean_names(mydata5)
mydata5 <- mydata5 %>% rename(PID=pid_113250)
mydata4<- merge(mydata5, mydata4, by = "PID", all = TRUE)
mydata6<- read_sav("~/Desktop/2024 12 04 Skjellegrind, tillegg 113250/2024-12-04_113250_Data.sav")
mydata6<-clean_names(mydata6)
mydata6 <- mydata6 %>% rename(PID=pid_113250)
mydata4<- merge(mydata6, mydata4, by = "PID", all = TRUE)

mydata4%>%filter(is.na((grs002249ad_ave_113861)))

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


#merging all datafiles 
# Remove duplicate columns from hpylori2 before merging
hpylori2 <- hpylori2[, !names(hpylori2) %in% c("sex", "birth_year")]

# Now merge without duplicate columns
hpylori2 <- merge(hpylori2, mydata4[, c("PID","PGS1", "PGS2", "PGS3","PGS4","PGS5","PGS1_ter", "PGS2_ter","PGS3_ter","PGS4_ter","PGS5_ter")], 
                  by = "PID", all = TRUE)


###########################
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
#hpylori2 <- hpylori2 %>% 
  #rename(crp=sem_crp_nt2blm) #crp
#hpylori2 <- hpylori2 %>% 
  #rename(sex=sex.x) #sex

#factor categorical variables 
hpylori2$helicobacterpylori_cat<-factor(hpylori2$helicobacterpylori_cat,levels = c("0","1"),labels = c("Negative","Positive"))
hpylori2$sex.x<-factor(hpylori2$sex.x,levels = c("0","1"),labels = c("Female","Male"))
hpylori2$register<-factor(hpylori2$register,levels = c("1","3","5"),labels = c("Alive","Moved","Dead"))
hpylori2$smoking<-factor(hpylori2$smoking,levels = c("0","1","2"),labels = c("Never","Previously","Currently"))
hpylori2$education<-factor(hpylori2$education,levels = c("1","2","3","4","5"),
                           labels = c("Primary","Secondary","highschool","tertiary","university"))
hpylori2$maritalstatus<-factor(hpylori2$maritalstatus,levels = c("1","2","3","4","5"),
                               labels = c("Unmarried","Married","Widow","Divorced","Separated"))
hpylori2$diabetes<-factor(hpylori2$diabetes,levels = c("0","1"),labels = c("No","Yes"))
hpylori2$cvd<-factor(hpylori2$cvd,levels = c("0","1"),labels = c("No","Yes"))  

#re-categorize variables 
library(forcats)
hpylori2$education<-fct_collapse(hpylori2$education,"Primary"=c("Primary"),
                                 "Highschool"=c("Secondary","highschool"),"University"=c("tertiary","university")) #education
hpylori2$maritalstatus<-fct_collapse(hpylori2$maritalstatus,"Unmarried"=c("Unmarried"),
                                     "Married"=c("Married"),"Widow_divorced_separated"=c("Widow","Divorced","Separated")) #marital status

#########################
###filter dataset########
##########################
# Load necessary libraries
library(dplyr)
library(lubridate)

# Filter dataset based on birth year
birthyear <- hpylori2 %>% filter(birth_year.x > 1914 & birth_year.x < 1950)
# Filter for available HP data
hpylori2 <- birthyear %>% filter(!is.na(helicobacterpylori_cat))
# Exclude those who moved based on register status =9
#hpylori2 <- hpylori2 %>% filter(register != "Moved") #not running this since moved=censored in status variable 
# Remove rows where register = 1 and obs_end_month is not NA =58
#hpylori2 <- hpylori2 %>%
  #filter(!(register == "Alive" & !is.na(obs_end_month)))

summary(hpylori2$PGS1)

#############################
#####Survival Analysis Prep#######
#############################

# Define status: 0 = censored (NA), 1 = event occurred (date of death present)
hpylori2 <- hpylori2 %>%mutate(status = ifelse(is.na(obs_end_month), 0, 1)
)

# Age at study start (age at HUNT2)
hpylori2 <- hpylori2 %>% 
  rename(age_at_study_start=age2)

#Age at censored
#set NAs to cutoff date in obs_end_month
cutoff_date <- ymd('2024-07-15')
hpylori2$obs_end_month <- as.Date(hpylori2$obs_end_month)
hpylori2$obs_end_month[is.na(hpylori2$obs_end_month)] <- cutoff_date

# Calculate the age at (censored) obs_end_month
hpylori2 <- hpylori2 %>%
  mutate(age_in_days = as.numeric(difftime(obs_end_month, birth_month, units = "days")))
hpylori2 <- hpylori2 %>%
  mutate(age_at_obs_end = round(age_in_days / 365, digits=1))

#Double check 
table(hpylori2$status)
table(hpylori2$obs_end_month)
table(hpylori2$age_at_obs_end)
table(hpylori2$age_at_study_start)

#subset and check
hpylori_subset <- hpylori2 %>% 
  select(status, obs_end_month, register,age_at_obs_end)

print(hpylori_subset)

#rename birthyear
hpylori2 <- hpylori2 %>% 
  rename(birth_year=birth_year.x)

#Cut birthyear 
# Create 10-year age bins
hpylori2$age_bin10 <- cut(hpylori2$birth_year, 
                          breaks = seq(from = min(hpylori2$birth_year, na.rm = TRUE) - (min(hpylori2$birth_year, na.rm = TRUE) %% 10), 
                                       to = max(hpylori2$birth_year, na.rm = TRUE) + 10, 
                                       by = 10),
                          include.lowest = TRUE, 
                          right = FALSE, 
                          labels = paste(seq(from = min(hpylori2$birth_year, na.rm = TRUE) - (min(hpylori2$birth_year, na.rm = TRUE) %% 10), 
                                             to = max(hpylori2$birth_year, na.rm = TRUE), 
                                             by = 10), 
                                         seq(from = min(hpylori2$birth_year, na.rm = TRUE) + 9 - (min(hpylori2$birth_year, na.rm = TRUE) %% 10), 
                                             to = max(hpylori2$birth_year, na.rm = TRUE) + 9, 
                                             by = 10),
                                         sep = "-"))
hpylori2$age_bin25 <- cut(hpylori2$birth_year, 
                          breaks = seq(from = min(hpylori2$birth_year, na.rm = TRUE) - (min(hpylori2$birth_year, na.rm = TRUE) %% 25), 
                                       to = max(hpylori2$birth_year, na.rm = TRUE) + 25, 
                                       by = 25),
                          include.lowest = TRUE, 
                          right = FALSE, 
                          labels = paste(seq(from = min(hpylori2$birth_year, na.rm = TRUE) - (min(hpylori2$birth_year, na.rm = TRUE) %% 5),
                                             to = max(hpylori2$birth_year, na.rm = TRUE), 
                                             by = 25), 
                                         seq(from = min(hpylori2$birth_year, na.rm = TRUE) + 24 - (min(hpylori2$birth_year, na.rm = TRUE) %% 5),
                                             to = max(hpylori2$birth_year, na.rm = TRUE) + 24, 
                                             by = 25),
                                         sep = "-"))

table(hpylori2$age_bin10)
# Load necessary libraries
library(survival)
library(survminer)
library(ggplot2)
library(broom)
library(gtsummary)
library(broom.helpers)
######################
#####KM######
######################
# Overall
surv_fit <- survfit(Surv(age_at_study_start,age_at_obs_end, status) ~ 1, data = hpylori2)
surv_fit
######################
###Plot: Over All### (not used in the article)
#######################
summary_fit <- summary(surv_fit)
if (is.null(dim(summary_fit$table))) {
  # No strata case
  median_survival <- summary_fit$table["median"]
  num_censored <- sum(summary_fit$n.censor)
  median_survival_times <- data.frame(
    strata = "All",
    median_time = median_survival
  )} 
# Set the x-axis limits starting from 45 years
x_min <- min(hpylori2$age_at_study_start, na.rm = TRUE)
x_max <- max(hpylori2$age_at_obs_end, na.rm = TRUE)  # Maximum observed age in the dataset

# Create the survival plot with the additional annotations
plot <- ggsurvplot(
  surv_fit,                      # Survival fit object
  xlab = "Age (Years)",         # X-axis label
  ylab = "Survival Probability", # Y-axis label
  xlim = c(x_min, x_max),      # Define x-axis limits
  break.time.by = 5,             # Break X axis in time intervals by 4.
  conf.int = TRUE,               # Show confidence interval
  surv.median.line = "hv",       # Add the median survival pointer.
  ggtheme = theme_classic2() + theme( # Enhanced theme customization
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ),
  title = "Kaplan-Meier Survival Curve", # Plot title
  legend.title = "",     # Legend title
  pval = FALSE,                          # Add p-value
  ncensor.plot = TRUE,                   # Plot the number of censored subjects at time t
  risk.table = TRUE,                     # Add risk table
  fontsize = 3,                          # Used in risk table
  risk.table.height = 0.25,              # Adjust risk table height
  risk.table.y.text.col = TRUE,          # Color risk table text
  risk.table.y.text = FALSE,             # Remove risk table row labels
  risk.table.title = "Number at risk"    # Risk table title
)


# Add annotations for median survival time
for (i in 1:nrow(median_survival_times)) {
  plot$plot <- plot$plot + 
    geom_text(data = median_survival_times, aes(x = median_time[i], y = 0.05, 
                                                label = sprintf("%.1f", median_time[i])),
              color = "black", size = 3, vjust = -0.5, hjust = 1.2, fontface = "italic") 
}

# Print the plot
print(plot)
######################
#Supplementary Figure 4: Kaplan–Meier curve estimating median survival age for overall HP seropositivity.
######################
# Fit the survival model
surv_fit <- survfit(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cat, data = hpylori2)

# Print the survival fit to check the results
print(surv_fit)

# Extract summary of the fit
summary_fit <- summary(surv_fit)
print(summary_fit)

# Extract median survival times
median_survival <- summary(surv_fit)$table[, "median"]
median_survival_times <- data.frame(
  strata = names(median_survival),
  median_time = median_survival
)
# Set the x-axis limits starting from 45 years
x_min <- min(hpylori2$age_at_study_start, na.rm = TRUE)
x_max <- max(hpylori2$age_at_obs_end, na.rm = TRUE)  # Maximum observed age in the dataset

# Create the survival plot
plot <- ggsurvplot(
  surv_fit,                     # Survival fit object
  xlab = "Age (Years)",         # X-axis label
  ylab = "Survival Probability",# Y-axis label
  xlim = c(x_min, x_max),      # Define x-axis limits
  break.time.by = 5,            # Break X axis in intervals of 5 years
  conf.int = TRUE,              # Show confidence intervals
  conf.int.style = "step",      # Style of confidence intervals
  palette = c("darkgreen", "pink"), # Custom color palette
  surv.median.line = "hv",      # Add the median survival pointer
  legend.labs = c("HP Negative", "HP Positive"), # Change legend labels
  ggtheme = theme_classic() + theme(  # Enhanced theme customization
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ),
  title = "Kaplan-Meier Survival Curve", # Plot title
  legend.title = "",                # Legend title
  pval = FALSE,                     # Do not add p-value
  ncensor.plot = TRUE,             # Plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,             # Reduce the height for better visibility
  risk.table = TRUE,               # Add risk table
  fontsize = 4,                    # Font size for risk table
  risk.table.height = 0.2,        # Height of the risk table
  risk.table.y.text.col = TRUE,    # Color risk table text
  risk.table.y.text = FALSE,       # Remove risk table row labels
  risk.table.title = "Number at risk" # Risk table title
)
# Customize the y-axis scale and appearance for the censor plot
if (!is.null(plot$ncensor.plot)) {
  # Extracting the censored counts and their corresponding time points
  censor_data <- data.frame(
    time = surv_fit$time,
    censored = surv_fit$n.censor
  )}
  
# Create the censor plot
  plot$ncensor.plot <- plot$ncensor.plot +
    scale_y_continuous(
      limits = c(0, 50),  # Adjust limits for the y-axis
      breaks = seq(0, 50, by = 10)  # Adjust breaks for y-axis
    ) +
    theme(
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.y = element_text(size = 14, color = "black"), # Customize y-axis text
      axis.text.x = element_text(size = 14, color = "black"), # Customize x-axis text
      axis.line = element_line(color = "black"),
      panel.grid.major = element_line(color = "lightgray", linetype = "dashed"),
      panel.grid.minor = element_blank()
    )


# Add median survival times as annotations
library(ggrepel)
plot$plot <- plot$plot +
  geom_label_repel(data = median_survival_times, 
                   aes(x = median_time, y = 0.05, label = sprintf("%.1f", median_time)),
                   color = "black", size = 4, fontface = "italic", 
                   fill = "white",   # White background fill for the labels
                   label.padding = unit(0.2, "lines"),  # Padding around the label text
                   box.padding = unit(0.3, "lines"),    # Padding to avoid overlap with other elements
                   point.padding = unit(0.5, "lines")) +# Padding between the label and the point
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(y = "Survival (%)", title = "Kaplan-Meier Survival Curve with Median Survival Times") +
  theme(legend.position = "top")  # Remove legend if it is not needed

# Adjust y-axis limits of the risk table to ensure visibility
plot$tables$`risk.table` <- plot$tables$`risk.table` +
  coord_cartesian(ylim = c(0, 1)) +  # Adjust according to your data
  theme(legend.position = "none") # Hide legend if not needed in risk table

# Print the plot
print(plot)

#################
#Supplementary Figure 5: Kaplan–Meier curve estimating median survival age for HP seropositivity stratified by sex. 
#################
# Fit the survival model
surv_fit <- survfit(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cat+sex, data = hpylori2)

# Print the survival fit to check the results
print(surv_fit)

# Extract summary of the fit
summary_fit <- summary(surv_fit)
print(summary_fit)

# Extract median survival times
median_survival <- summary(surv_fit)$table[, "median"]
median_survival_times <- data.frame(
  strata = names(median_survival),
  median_time = median_survival
)
# Set the x-axis limits starting from 45 years
x_min <- min(hpylori2$age_at_study_start, na.rm = TRUE)
x_max <- max(hpylori2$age_at_obs_end, na.rm = TRUE)  # Maximum observed age in the dataset

# Create the survival plot
plot <- ggsurvplot(
  surv_fit,                     # Survival fit object
  xlab = "Age (Years)",         # X-axis label
  ylab = "Survival Probability",# Y-axis label
  xlim = c(x_min, x_max),      # Define x-axis limits
  break.time.by = 5,            # Break X axis in intervals of 5 years
  conf.int = TRUE,              # Show confidence intervals
  conf.int.style = "step",      # Style of confidence intervals
  palette = c("darkgreen", "pink","yellow","blue"), # Custom color palette
  surv.median.line = "hv",      # Add the median survival pointer
  legend.labs = c("HP Negative Female","HP Negative Male","HP Positive Female","HP Positive Male"), # Change legend labels
  ggtheme = theme_classic() + theme(  # Enhanced theme customization
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ),
  title = "Kaplan-Meier Survival Curve", # Plot title
  legend.title = "",                # Legend title
  pval = FALSE,                     # Do not add p-value
  ncensor.plot = FALSE,             # Plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,             # Reduce the height for better visibility
  risk.table = TRUE,               # Add risk table
  fontsize = 4,                    # Font size for risk table
  risk.table.height = 0.2,        # Height of the risk table
  risk.table.y.text.col = TRUE,    # Color risk table text
  risk.table.y.text = FALSE,       # Remove risk table row labels
  risk.table.title = "Number at risk" # Risk table title
)

# Customize the y-axis scale for the censor plot
if (!is.null(plot$ncensor.plot)) {
  plot$ncensor.plot <- plot$ncensor.plot +
    scale_y_continuous(
      limits = c(0, 150),  # Set limits for the y-axis, adjust as needed
      breaks = seq(0, 200, by = 50)  # Define breaks for the y-axis
    ) +
    theme(
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 10)
    )
}

# Add median survival times as annotations
library(ggrepel)
plot$plot <- plot$plot +
  geom_label_repel(data = median_survival_times, 
                   aes(x = median_time, y = 0.05, label = sprintf("%.1f", median_time)),
                   color = "black", size = 4, fontface = "italic", 
                   fill = "white",   # White background fill for the labels
                   label.padding = unit(0.2, "lines"),  # Padding around the label text
                   box.padding = unit(0.3, "lines"),    # Padding to avoid overlap with other elements
                   point.padding = unit(0.5, "lines")) +# Padding between the label and the point
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(y = "Survival (%)", title = "Kaplan-Meier Survival Curve with Median Survival Times") +
  theme(legend.position = "top")  # Remove legend if it is not needed

# Print the plot
print(plot)

############
#Supplementary Figure 6: Kaplan–Meier curve estimating median survival age for HP seropositivity stratified by birth cohort.
table(hpylori2$age_bin10)
# Fit the survival model
surv_fit <- survfit(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cat+age_bin25, data = hpylori2)

# Print the survival fit to check the results
print(surv_fit)

# Extract summary of the fit
summary_fit <- summary(surv_fit)
print(summary_fit)

# Extract median survival times
median_survival <- summary(surv_fit)$table[, "median"]
median_survival_times <- data.frame(
  strata = names(median_survival),
  median_time = median_survival
)
# Set the x-axis limits starting from 45 years
x_min <- min(hpylori2$age_at_study_start, na.rm = TRUE)
x_max <- max(hpylori2$age_at_obs_end, na.rm = TRUE)  # Maximum observed age in the dataset

# Create the survival plot
plot <- ggsurvplot(
  surv_fit,                     # Survival fit object
  xlab = "Age (Years)",         # X-axis label
  ylab = "Survival Probability",# Y-axis label
  xlim = c(x_min, x_max),      # Define x-axis limits
  break.time.by = 5,            # Break X axis in intervals of 5 years
  conf.int = TRUE,              # Show confidence intervals
  conf.int.style = "step",      # Style of confidence intervals
  palette = c("darkgreen", "pink","yellow","blue"), # Custom color palette
  surv.median.line = "hv",      # Add the median survival pointer
  legend.labs = c("HP Negative 1915-1939","HP Negative 1940-1964","HP Positive 1915-1939","HP Positive 1940-1964"), # Change legend labels
  ggtheme = theme_classic() + theme(  # Enhanced theme customization
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ),
  title = "Kaplan-Meier Survival Curve", # Plot title
  legend.title = "",                # Legend title
  pval = FALSE,                     # Do not add p-value
  ncensor.plot = FALSE,             # Plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,             # Reduce the height for better visibility
  risk.table = TRUE,               # Add risk table
  fontsize = 4,                    # Font size for risk table
  risk.table.height = 0.2,        # Height of the risk table
  risk.table.y.text.col = TRUE,    # Color risk table text
  risk.table.y.text = FALSE,       # Remove risk table row labels
  risk.table.title = "Number at risk" # Risk table title
)

# Customize the y-axis scale for the censor plot
if (!is.null(plot$ncensor.plot)) {
  plot$ncensor.plot <- plot$ncensor.plot +
    scale_y_continuous(
      limits = c(0, 150),  # Set limits for the y-axis, adjust as needed
      breaks = seq(0, 200, by = 50)  # Define breaks for the y-axis
    ) +
    theme(
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 10)
    )
}

# Add median survival times as annotations
library(ggrepel)
plot$plot <- plot$plot +
  geom_label_repel(data = median_survival_times, 
                   aes(x = median_time, y = 0.05, label = sprintf("%.1f", median_time)),
                   color = "black", size = 4, fontface = "italic", 
                   fill = "white",   # White background fill for the labels
                   label.padding = unit(0.2, "lines"),  # Padding around the label text
                   box.padding = unit(0.3, "lines"),    # Padding to avoid overlap with other elements
                   point.padding = unit(0.5, "lines")) +# Padding between the label and the point
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(y = "Survival (%)", title = "Kaplan-Meier Survival Curve with Median Survival Times") +
  theme(legend.position = "top")  # Remove legend if it is not needed

# Print the plot
print(plot)

############
#Supplementary Figure 7: Kaplan–Meier curve estimating median survival age for HP seropositivity stratified by APOE carrier status. 
table(hpylori2$age_bin10)
# Fit the survival model
surv_fit <- survfit(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cat+E4, data = hpylori2)

# Print the survival fit to check the results
print(surv_fit)

# Extract summary of the fit
summary_fit <- summary(surv_fit)
print(summary_fit)

# Extract median survival times
median_survival <- summary(surv_fit)$table[, "median"]
median_survival_times <- data.frame(
  strata = names(median_survival),
  median_time = median_survival
)
# Set the x-axis limits starting from 45 years
x_min <- min(hpylori2$age_at_study_start, na.rm = TRUE)
x_max <- max(hpylori2$age_at_obs_end, na.rm = TRUE)  # Maximum observed age in the dataset

# Create the survival plot
plot <- ggsurvplot(
  surv_fit,                     # Survival fit object
  xlab = "Age (Years)",         # X-axis label
  ylab = "Survival Probability",# Y-axis label
  xlim = c(x_min, x_max),      # Define x-axis limits
  break.time.by = 5,            # Break X axis in intervals of 5 years
  conf.int = TRUE,              # Show confidence intervals
  conf.int.style = "step",      # Style of confidence intervals
  palette = c("darkgreen", "pink","yellow","blue"), # Custom color palette
  surv.median.line = "hv",      # Add the median survival pointer
  legend.labs = c("HP and APOE Negative","HP Negative and APOE Positive","HP Positive APOE Negative","HP and APOE Positive"), # Change legend labels
  ggtheme = theme_classic() + theme(  # Enhanced theme customization
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ),
  title = "Kaplan-Meier Survival Curve", # Plot title
  legend.title = "",                # Legend title
  pval = FALSE,                     # Do not add p-value
  ncensor.plot = FALSE,             # Plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,             # Reduce the height for better visibility
  risk.table = TRUE,               # Add risk table
  fontsize = 4,                    # Font size for risk table
  risk.table.height = 0.2,        # Height of the risk table
  risk.table.y.text.col = TRUE,    # Color risk table text
  risk.table.y.text = FALSE,       # Remove risk table row labels
  risk.table.title = "Number at risk" # Risk table title
)

# Customize the y-axis scale for the censor plot
if (!is.null(plot$ncensor.plot)) {
  plot$ncensor.plot <- plot$ncensor.plot +
    scale_y_continuous(
      limits = c(0, 150),  # Set limits for the y-axis, adjust as needed
      breaks = seq(0, 200, by = 50)  # Define breaks for the y-axis
    ) +
    theme(
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 10)
    )
}

# Add median survival times as annotations
library(ggrepel)
plot$plot <- plot$plot +
  geom_label_repel(data = median_survival_times, 
                   aes(x = median_time, y = 0.05, label = sprintf("%.1f", median_time)),
                   color = "black", size = 4, fontface = "italic", 
                   fill = "white",   # White background fill for the labels
                   label.padding = unit(0.2, "lines"),  # Padding around the label text
                   box.padding = unit(0.3, "lines"),    # Padding to avoid overlap with other elements
                   point.padding = unit(0.5, "lines")) +# Padding between the label and the point
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(y = "Survival (%)", title = "Kaplan-Meier Survival Curve with Median Survival Times") +
  theme(legend.position = "top")  # Remove legend if it is not needed

# Print the plot
print(plot)

###########
#Supplementary Figure 8: Kaplan–Meier curve estimating median survival age for HP seropositivity stratified by CRP levels. 
##########
hpylori2$crp2<- ifelse(hpylori2$crp > 5,1,0)
hpylori2$crp2<-factor(hpylori2$crp2,levels = c("0","1"),labels = c("Below5","Above5"))
table(hpylori2$crp2)

table(hpylori2$age_bin10)
# Fit the survival model
surv_fit <- survfit(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cat+crp2, data = hpylori2)

# Print the survival fit to check the results
print(surv_fit)

# Extract summary of the fit
summary_fit <- summary(surv_fit)
print(summary_fit)

# Extract median survival times
median_survival <- summary(surv_fit)$table[, "median"]
median_survival_times <- data.frame(
  strata = names(median_survival),
  median_time = median_survival
)
# Set the x-axis limits starting from 45 years
x_min <- min(hpylori2$age_at_study_start, na.rm = TRUE)
x_max <- max(hpylori2$age_at_obs_end, na.rm = TRUE)  # Maximum observed age in the dataset

# Create the survival plot
plot <- ggsurvplot(
  surv_fit,                     # Survival fit object
  xlab = "Age (Years)",         # X-axis label
  ylab = "Survival Probability",# Y-axis label
  xlim = c(x_min, x_max),      # Define x-axis limits
  break.time.by = 5,            # Break X axis in intervals of 5 years
  conf.int = TRUE,              # Show confidence intervals
  conf.int.style = "step",      # Style of confidence intervals
  palette = c("darkgreen", "pink","yellow","blue"), # Custom color palette
  surv.median.line = "hv",      # Add the median survival pointer
  legend.labs = c("HP Negative and CRP below 5 mg/l","HP Negative and CRP above 5 mg/l","HP Positive and CRP below 5 mg/l","HP Positive and CRP above 5 mg/l"), # Change legend labels
  ggtheme = theme_classic() + theme(  # Enhanced theme customization
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  ),
  title = "Kaplan-Meier Survival Curve", # Plot title
  legend.title = "",                # Legend title
  pval = FALSE,                     # Do not add p-value
  ncensor.plot = FALSE,             # Plot the number of censored subjects at time t
  ncensor.plot.height = 0.25,             # Reduce the height for better visibility
  risk.table = TRUE,               # Add risk table
  fontsize = 4,                    # Font size for risk table
  risk.table.height = 0.2,        # Height of the risk table
  risk.table.y.text.col = TRUE,    # Color risk table text
  risk.table.y.text = FALSE,       # Remove risk table row labels
  risk.table.title = "Number at risk" # Risk table title
)

# Customize the y-axis scale for the censor plot
if (!is.null(plot$ncensor.plot)) {
  plot$ncensor.plot <- plot$ncensor.plot +
    scale_y_continuous(
      limits = c(0, 150),  # Set limits for the y-axis, adjust as needed
      breaks = seq(0, 200, by = 50)  # Define breaks for the y-axis
    ) +
    theme(
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 10)
    )
}

# Add median survival times as annotations
library(ggrepel)
plot$plot <- plot$plot +
  geom_label_repel(data = median_survival_times, 
                   aes(x = median_time, y = 0.05, label = sprintf("%.1f", median_time)),
                   color = "black", size = 4, fontface = "italic", 
                   fill = "white",   # White background fill for the labels
                   label.padding = unit(0.2, "lines"),  # Padding around the label text
                   box.padding = unit(0.3, "lines"),    # Padding to avoid overlap with other elements
                   point.padding = unit(0.5, "lines")) +# Padding between the label and the point
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(y = "Survival (%)", title = "Kaplan-Meier Survival Curve with Median Survival Times") +
  theme(legend.position = "top")  # Remove legend if it is not needed

# Print the plot
print(plot)



#################################
#Supplementary Table 6: Hazard Ratios of all-cause Mortality by HP Seropositivity and titers 
#################################

####################
#HP seropositivity##
####################
#model 1
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cat+sex.x+education, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)
sum(is.na(hpylori2$education)) #371 missing education 

#model 2: model 1+birth cohort (10 year interval)
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cat + sex.x + education + age_bin10, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 3: model 2+smoking, marital, bmi,cvd, diabetes,crp
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cat + sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 4: crp and APOE
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cat + sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp+E4, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#############
#HP titers###
#############
#normalize 
library(RNOmni)
hpylori2$helicobacterpylori_cont<-RankNorm(as.numeric(hpylori2$helicobacterpylori_cont))
#model 1
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cont+sex.x+education, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 2: model 1+birth cohort (10 year interval)
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cont + sex.x + education + age_bin10, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 3: model 2+smoking, marital, bmi,cvd, diabetes,crp
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cont + sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 4: crp and APOE
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ helicobacterpylori_cont + sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp+E4, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

##############
####PGS#######
##############
#PGS002249 
#model 1
library(RNOmni)
summary(hpylori2$PGS1)
hist(hpylori2$PGS1)
hist(hpylori2$PGS3)

cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ RankNorm(as.numeric(PGS1))+sex.x+education, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)
sum(is.na(hpylori2$education)) #371 missing education 

#model 2: model 1+birth cohort (10 year interval)
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS1+ sex.x + education + age_bin10, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 3: model 2+smoking, marital, bmi,cvd, diabetes,crp
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS1+ sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 4: crp and APOE
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS1+ sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp+E4, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#PGS1 (PGS002249)  10% group
#model 1
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS1_ter+sex.x+education, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)
sum(is.na(hpylori2$education)) #371 missing education 

#model 2: model 1+birth cohort (10 year interval)
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS1_ter+ sex.x + education + age_bin10, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 3: model 2+smoking, marital, bmi,cvd, diabetes,crp
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS1_ter+ sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 4: crp and APOE
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS1_ter+ sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp+E4, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#PGS3 (PGS004034)  10% group
#model 1
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS3_ter+sex.x+education, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)
sum(is.na(hpylori2$education)) #371 missing education 

#model 2: model 1+birth cohort (10 year interval)
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS3_ter+ sex.x + education + age_bin10, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 3: model 2+smoking, marital, bmi,cvd, diabetes,crp
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS3_ter+ sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 4: crp and APOE
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS3_ter+ sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp+E4, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#PGS4 (PGS004228)  10% group
#model 1
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS4_ter+sex.x+education, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)


#model 2: model 1+birth cohort (10 year interval)
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS4_ter+ sex.x + education + age_bin10, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 3: model 2+smoking, marital, bmi,cvd, diabetes,crp
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS4_ter+ sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 4: crp and APOE
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS4_ter+ sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp+E4, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#PGS5 (PGS004229 (NO APOE))  10% group
#model 1
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS5_ter+sex.x+education, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)


#model 2: model 1+birth cohort (10 year interval)
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS5_ter+ sex.x + education + age_bin10, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 3: model 2+smoking, marital, bmi,cvd, diabetes,crp
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS5_ter+ sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)

#model 4: crp and APOE
cox_fit <- coxph(Surv(age_at_study_start, age_at_obs_end, status) ~ PGS5_ter+ sex.x + education + age_bin10+maritalstatus+bmi+smoking+cvd+diabetes+crp+E4, data = hpylori2)
summary(cox_fit)
tbl_regression(cox_fit, exponentiate = TRUE)




