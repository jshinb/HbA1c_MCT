## ----load_packages---------------------------------------------------------------------
# 
# This script will explore data and creates scatter and smooth plots:
# 
#---------------------------------------------------------------------------------------#

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

## ----brain_volume_ICV_by_age, fig.width=8.5, fig.height=4.5----------------------------------------------------------------------------------
p1 = idata %>% ggplot(aes(x=Age_MRI,y=BrainVolume,color=Sex)) +
  geom_point(alpha=0.3) + 
  geom_smooth() + 
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue"))+
  xlab(NULL) + ylab(NULL) + theme(legend.position = 'top')
p2 = idata %>% ggplot(aes(x=Age_MRI,y=ICV,color=Sex)) +
  geom_point(alpha=0.3) + 
  geom_smooth() + 
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue"))+
  xlab(NULL) + ylab(NULL) + theme(legend.position = 'top')

grid.arrange(p1,p2,nrow=1)

png(paste("Plots/BrainVolume_by_age.png",sep=""),
    width=6,height=4.5, units = "in", res=300)
print(p1+ theme(legend.position = "none") + xlab("Age_MRI") + ylab("Brain volume (mm3)"))
dev.off()

png(paste("Plots/ICV_by_age.png",sep=""),
    width=6,height=4.5, units = "in", res=300)
print(p2+ theme(legend.position = "none") + xlab("Age_MRI") + ylab("Intracranial volume (mm3)"))
dev.off()
rm(p1,p2);gc(reset=T)

## ----MCT_by_HbA1c, fig.width=8.5, fig.height=4.5---------------------------------------------------------------------------------------------
p1 = idata %>% ggplot(aes(x=rntHbA1c,y=MCT,color=Sex)) +
  geom_point(alpha=0.2) + 
  geom_smooth() + 
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue"))+
  xlab(NULL) + ylab(NULL) + theme(legend.position = "top")
p2 = idata %>% ggplot(aes(x=rntHbA1c,y=MCT,color=Sex)) +
  geom_smooth() + 
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue"))+
  xlab(NULL) + ylab(NULL) + theme(legend.position = "top")

grid.arrange(p1,p2,nrow=1)

png(paste("Plots/MCT_by_HbA1c.png",sep=""),
    width=6,height=4.5, units = "in", res=300)
print(p1+ theme(legend.position = "none"))
dev.off()

png(paste("Plots/MCT_by_HbA1c_wo_points.png",sep=""),
    width=6,height=4.5, units = "in", res=300)
print(p2+ theme(legend.position = "none") )
dev.off()
rm(p1,p2);gc(reset=TRUE)
