## ----load_packages---------------------------------------------------------------------
# This script will create a merged data set 
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

## ----load_data-------------------------------------------------------------------------------------------------------------------------------
setwd(wd)
load("data/ukbb_combined_MRI_2021-02-18.RData")#32589 #with HbA1c.Level

## MRI
MRI_HbA1c = fread('data/ukbb_MCT_HbA1c_2021-05-18.tsv.gz')

## Covariates
Covariates = fread('data/ukbb_covdat_2021-05-18.tsv.gz')

## ICV and brain volume
other_brianpheno = fread('data/ukbb_other_brain_phenos_2021-05-18.tsv.gz')

## T2D
dd.Diabetes = fread('SummaryDiabetesDiagnosis.csv')

## stroke/dementia
stroke_dementia = fread('data/ukbbStroke_Dementia_at_anytime.txt')

load("data/ukbb_combined_MRI_2021-02-18.RData")#32589 #with HbA1c.Level
icv_brainvolume = fread('data/ukbb_icv_brainvolume.txt.gz')#39971

## WMH
dd.WMH = fread('../ukbb/data/ukb41763.csv.gz') %>% 
  dplyr::select(eid,"25781-2.0")
dd.WMH =  dd.WMH %>% dplyr::rename("WMH" = "25781-2.0")
dd.WMH = subset(dd.WMH, !is.na(WMH))
# hist(dd.WMH$WMH)

## diagetes
dd.Diabetes = fread('SummaryDiabetesDiagnosis.csv')

## covdat for all
covdat_all = fread('data/covdat_all.tsv.gz')

## ROI-MRI data
d_MRI2 = fread("Replication_ukb_MRI_2020-11-26.tsv")#17690
d_MRI1 = fread('ThicknessData.csv')

## genotype PC for population stratification:22009-0.1 - 22009-0.5
ukb37194 = fread('~/OneDrive - SickKids/ukbb/data/ukb37194.csv.gz')

## ----data_manipulation-----------------------------------------------------------------------------------------------------------------------
tmp = merge(tmp,icv_brainvolume,sort=F,all.x=T)
tmp=subset(tmp,!is.na(HbA1c.Level))
tmp$Sex  = factor(tmp$Sex, levels=c("Male","Female"))
tmp$HbA1c.Level = factor(tmp$HbA1c.Level, levels=c("High","Normal"))
## brain volume
tmp = tmp %>% 
  mutate(abs.scaled.BV = abs(scale(BrainVolume)),abs.scaled.ICV = abs(scale(ICV)))
tmp$ICV[tmp$abs.scaled.ICV>4] <- NA
tmp$BrainVolume[tmp$abs.scaled.BV>4] <- NA
## WMH
tmp = merge(tmp,dd.WMH,all.x=T)

## diabetes
tmp = merge(tmp,dd.Diabetes,all.x=T)

idata=subset(tmp,(!eid %in% stroke_dementia$eid & !is.na(HbA1c2)))#30579
idata = merge(idata,subset(covdat_all,select=c(eid,age_base,SBP0.1,DBP0.1,SBP0.2,DBP0.2,SBP1.1,DBP1.1,SBP1.2,DBP1.2)),sort=F)
idata$SBP0 = apply(data.frame(subset(idata,select=c(SBP0.1,SBP0.2))),1,mean,na.rm=T)
idata$DBP0 = apply(data.frame(subset(idata,select=c(DBP0.1,DBP0.2))),1,mean,na.rm=T)
idata$SBP1 = apply(data.frame(subset(idata,select=c(SBP1.1,SBP1.2))),1,mean,na.rm=T)
idata$DBP1 = apply(data.frame(subset(idata,select=c(DBP1.1,DBP1.2))),1,mean,na.rm=T)
idata = idata %>% mutate(Age.c = Age_MRI - mean(Age_MRI),
                         Time = Age_MRI - age_base,
                         Time.c = Time - mean(Time))
save(idata,file=file.path(wd,"idata.Rd"))
