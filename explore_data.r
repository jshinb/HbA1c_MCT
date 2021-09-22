## ----load_packages---------------------------------------------------------------------------------------------------------------------------
# ====  UKB | Participant Characteristic Table =======================================================
# This script will create a merged data set and a participant characteristic table
# https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html
options(stringsAsFactors = F)
library(stringr)
library(data.table)
library(tidyverse)
library(tableone)
library(dplyr)
library(grid)
library(gratia)
library(gridExtra)
library(GenABEL)
# for tables
library(arsenal)
library(here)
library(mgsub)
library(psych)
wd = '~/OneDrive - SickKids/ukbb_insulin_resistance'
load(file.path(wd,"idata.Rd"))

# Participant characteristic table ------------------------------------------------------
tab1 = CreateTableOne(data = subset(idata,select=-eid),strata = 'Sex')
print(tab1,nonnormal = c('Age_MRI','BMI_recruitment','BMI_adjAge',
                         'SBP',"SBP_adjAge",
                         'DBP',"DBP_adjAge",
                         'HbA1c2','HbA1c_adjAge',
                         'WMH','WMH_adjAge',
                         "ICV","ICV_adjAge",
                         'logWMH','logWMH_adjAge',
                         "MCT","MCT_adjAge"))

# Exploratory plots ---------------------------------------------------------------------
## Age at baseline vs HbA1c
p1 = idata %>% ggplot(aes(x=age_base,y=HbA1c2,color=Sex)) +
  geom_point(alpha=0.3) + 
  geom_smooth() + 
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue"))+
  xlab(NULL) + ylab(NULL)+ theme(legend.position = "top")

p2 = idata %>% ggplot(aes(x=age_base,y=HbA1c2,color=Sex)) +
  geom_smooth() + 
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue"))+
  xlab(NULL) + ylab(NULL)+ theme(legend.position = "top")

grid.arrange(p1,p2,nrow=1)
rm(p1,p2);gc(reset=T)

## Age at MRI vs HbA1c
p1 = idata %>% ggplot(aes(x=Age_MRI,y=HbA1c2,color=Sex)) +
  geom_point(alpha=0.3) + 
  geom_smooth() + 
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue"))+
  xlab(NULL) + ylab(NULL)+ theme(legend.position = "top")

p2 = idata %>% ggplot(aes(x=Age_MRI,y=HbA1c2,color=Sex)) +
  geom_smooth() + 
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue"))+
  xlab(NULL) + ylab(NULL)+ theme(legend.position = "top")

grid.arrange(p1,p2,nrow=1)
png(paste("Plots/rntHbA1c_by_age.png",sep=""),
    width=6,height=4.5, units = "in", res=300)
print(p1+ theme(legend.position = "none"))
dev.off()

png(paste("Plots/rntHbA1c_by_age_smooth.png",sep=""),
    width=6,height=4.5, units = "in", res=300)
print(p2+ theme(legend.position = "none"))
dev.off()
rm(p1,p2);gc(reset=T)
