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

wd = "~/OneDrive - SickKids/ukbb_insulin_resistance/"
setwd(wd)
load('data/idata_2021-06-30.RData')#idata
load('data/idataH_2021-06-30.RData')

dd.cov = fread('../ukbb/data/ukb41763.csv.gz')
dd.cov.sub = subset(dd.cov, select=c('eid','74-0.0','26521-2.0','26514-2.0'))
names(dd.cov.sub)[-1] = c('fasting_time_hours_0','ICV_mm3','BrainVolume_mm3')
describe(dd.cov.sub)
# icv_brainvolume = fread('data/ukbb_icv_brainvolume.txt.gz')#39971#same as above
idata = merge(idata,dd.cov.sub,all.x = T,sort=F)
idata = idata %>% mutate(ICV=ICV_mm3/10^3,
                         BrainVolume = BrainVolume_mm3/10^3,
                         Time=as.numeric(as.character(Time)))
# Tables ----------------------------------------------------------------------
vnames = c('Age_MRI',
           'BMI_recruitment',
           'SBP',
           'DBP',
           'Smoking',
           'EduLevel',
           'Diabetes_TypeII',
           'HbA1c2',
           'fasting_time_hours_0',
           'Time',
           #brain
           'MRI_Site',
           'MCT',
           'BrainVolume',
           "ICV")
strata='HbA1c.Level3'

#female -----------------------------------------------------------------------
tab2 = CreateTableOne(data = subset(idata,Sex=="Female",select=c(vnames,strata)),
vars = c(vnames),strata = strata)
ofile="~/Downloads/Table2.txt"
cat("---------------------------- Female ----------------------------\n",file=ofile)
capture.output(summary(tab2),file=ofile,append=T)
capture.output(print(tab2),file=ofile,append=T)
capture.output(print(tab2,nonnormal = c("Age_MRI", "BMI_recruitment", "SBP", "DBP",
"HbA1c2", "fasting_time_hours_0","Time","MCT","BrainVolume", "ICV"))), 
file=ofile,append=T)

#male -------------------------------------------------------------------------
strata='HbA1c.Level3'
tabM = CreateTableOne(data = subset(idata,Sex=="Male",select=c(vnames,strata)),
                      vars = c(vnames),strata = strata)
cat("\n\n",file=ofile)
cat("----------------------------- Male -----------------------------\n",file=ofile)
capture.output(summary(tabM),file=ofile,append=T)
capture.output(print(tabM),file=ofile,append=T)
capture.output(print(tabM,nonnormal = c("Age_MRI", "BMI_recruitment", "SBP", "DBP",
                                        "HbA1c2", "fasting_time_hours_0", "Time",
                                        "MCT", "BrainVolume", "ICV")), file=ofile,append=T)
setdiff(vnames,names(idata))

