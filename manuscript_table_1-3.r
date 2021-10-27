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
idata = idata %>% mutate(ICV=ICV_mm3/10^3,BrainVolume = BrainVolume_mm3/10^3)

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
idata = idata %>% mutate(Time=as.numeric(as.character(Time)))
#Table 1-----------------------------------------------------------------
strata = 'Sex'
tab1 = CreateTableOne(data = subset(idata,select=c(vnames,strata,'HbA1c.Level3')),
                      vars=c(vnames,'HbA1c.Level3'),strata = strata)
print(tab1,nonnormal = c("Age_MRI", "BMI_recruitment", "SBP", "DBP", 
                         "HbA1c2", "fasting_time_hours_0", #"Time", 
                         "MCT", "BrainVolume", "ICV"))
summary(tab1)

#Table 2-----------------------------------------------------------------
strata='HbA1c.Level3'
tab2 = CreateTableOne(data = subset(idata,select=c('Sex',vnames,strata)),
                      vars = c(vnames,'Sex'),strata = strata)
print(tab2,nonnormal = c("Age_MRI", "BMI_recruitment", "SBP", "DBP", 
                         "HbA1c2", #"fasting_time_hours_0", "Time", 
                         "MCT", "BrainVolume", "ICV"))
summary(tab2)

## permutation test: 
t.test.p = c()
for(i in 1:10000){
  idata.sub = idata.sub %>% mutate(HbA1c.p=sample(HbA1c.Level3))
  if(i==1){
    test.res = t.test(Time~HbA1c.Level3,data=idata.sub)
  }else{
    test.res = t.test(Time~HbA1c.p,data=idata.sub)
  }
  t.test.p = c(t.test.p,test.res$statistic)
  
}
#Table 3-----------------------------------------------------------------
strata='Sex'
tab3 = CreateTableOne(data = subset(idata,HbA1c.Level3=="High",select=c(vnames,strata)),
                      vars = vnames,strata = strata)
print(tab3,nonnormal = c("Age_MRI", "BMI_recruitment", "SBP", "DBP", 
                         "HbA1c2", "fasting_time_hours_0", "Time", 
                         "MCT", "BrainVolume", "ICV"))
summary(tab3)