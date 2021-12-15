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
#------------------------------------------------------------------------------
names(idata)

library(relaimpo)
names(idata)[str_detect(names(idata),"top")]

fit.lm = lm(top.2.MCT ~  Sex + Sex*(Age.c)  + Sex*I(Age.c^2) + Sex*Time + Sex*rntHbA1c + Sex*MRI_Site, 
            data=idataH %>% mutate(Age.c = scale(Age_MRI)))

# relaimpo::calc.relimp(fit.lm): cant' do this
summary(fit.lm)$coef
summary(fit.lm)$r.squared

#Relative importance metrics: 
#   lmg
# MRI_Site   0.0013702251
# Sex        0.0055898528
# Age.c      0.1536276274
# Time       0.0012576271
# rntHbA1c   0.0150953493
# I(Age.c^2) 0.0075787434
# Sex:Age.c  0.0009573632
# Sex:Time   0.0006403516
  
bt2 <- boot.relimp(fit.lm,b=1000,type='lmg')
booteval.relimp(bt2,lev=0.95,nodiff=TRUE)
fit.lm

fit.lmF = lm(top.2.MCT ~ Age.c  + I(Age.c^2) + Time + rntHbA1c + MRI_Site, 
            data=idataH %>% mutate(Age.c = scale(Age_MRI)) %>% subset(Sex=="Female"))
nobs(fit.lmF)
btF <- boot.relimp(fit.lmF,b=1000,type='lmg')
booteval.relimp(btF,lev=0.95,nodiff=TRUE)

fit.lmM = lm(top.2.MCT ~ Age.c  + I(Age.c^2) + Time + rntHbA1c + MRI_Site, 
             data=idataH %>% mutate(Age.c = scale(Age_MRI)) %>% subset(Sex=="Male"))
summary(fit.lmM)
btM <- boot.relimp(fit.lmM,b=1000,type='lmg')
booteval.relimp(btM,lev=0.95,nodiff=TRUE)


fit.lmFull = lm(top.2.MCT ~  Sex + Sex*(Age.c)  + Sex*I(Age.c^2) + Sex*Time + Sex*rntHbA1c + Sex*MRI_Site, 
            data=idataH %>% mutate(Age.c = scale(Age_MRI)))
fit.lm = lm(top.2.MCT ~  Sex + Sex*(Age.c)  + I(Age.c^2) + Sex*Time + rntHbA1c + MRI_Site, 
            data=idataH %>% mutate(Age.c = scale(Age_MRI)))
anova(fit.lmFull,fit.lm)#p=0.80 (F-test, anova)
bt <- boot.relimp(fit.lm,b=1000,type='lmg')
booteval.relimp(bt,lev=0.95,nodiff=TRUE)
citation("relaimpo")
