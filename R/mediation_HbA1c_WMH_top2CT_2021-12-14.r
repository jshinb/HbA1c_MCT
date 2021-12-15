## mediation analysis
options(stringsAsFactors = F)

library(medflex)
library(data.table)
library(tidyverse)
library(stringr)
library(psych)

names(idata)
#-----------------------#
# X = rntHbA1c
# M = logWMH
# Y = top.2.MCT
#-----------------------#
dat = idataH

dd.WMH = fread('../ukbb/data/ukb41763.csv.gz')
dd.WMH =  dd.WMH %>% dplyr::select(eid,"25781-2.0") %>% dplyr::rename("WMH" = "25781-2.0")
dd.WMH = subset(dd.WMH, !is.na(WMH))

# 0. Must scale the 'X"
dat = dat %>% mutate(rntHbA1c.scaled = (rntHbA1c-mean(rntHbA1c,na.rm=T))/sd(rntHbA1c,na.rm=T))
dat = merge(dat, dd.WMH)
dat = dat %>% 
  mutate(age = scale(Age_MRI)[,1],
         logWMH=log(1+WMH), 
         sex=ifelse(Sex=="Male",0,1),
         MRI_Site = factor(MRI_Site)) %>%
  mutate(a2 = age^2)

# 1. fit: M ~ X + C: weighting-based approach - did not work
# Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
# contrasts can be applied only to factors with 2 or more levels

medFit = glm(log(1+WMH) ~ rntHbA1c.scaled + age + a2 + sex + Time + MRI_Site, 
             family = gaussian(), 
             data=dat)

summary(medFit)
expData <- neWeight(medFit)
w <- weights(expData)## for further inspection of the empirical distribution
set.seed(20211214)
neMod1 <- neModel(top.2.MCT ~ rntHbA1c.scaled0 + rntHbA1c.scaled1 + age + a2 + sex + Time + MRI_Site,
                 family=gaussian(), expData=expData, nBoot = 5000, progress = F)
summary(neMod1)

# 2. imputation-based approach avoiding reliance on a model for the mediator 
# distribution: 
# the function both expands the data along hypothetical exposure values 
# (for each observation unit i) and imputes nested counterfactual outcomes in 
# this expanded dataset in a single run. Imputed counterfactual outcomes are 
# predictions from the imputation model that can be specified either externally 
# as a fitted model object.

dirFit <- glm(top.2.MCT ~ rntHbA1c.scaled + age +  a2+ sex + Time + MRI_Site,
              family = gaussian(), 
              data=dat)
impFit <- glm(top.2.MCT ~ rntHbA1c.scaled + logWMH + age +  a2+ sex + Time + MRI_Site,
              family = gaussian(), 
              data=dat)
summary(impFit)
expData <- neImpute(impFit)

# 3. fitting natural effect model (bootstrapping-based delta-method based B=5000)
nB=5000
neMod1 <- neModel(top.2.MCT ~ rntHbA1c.scaled0 + rntHbA1c.scaled1 + age + a2 + sex + Time + MRI_Site,
                  family=gaussian(), expData=expData, nBoot = nB, progress = T)
summary(neMod1)

# 4. results:
rm(effdecomp1)
effdecomp1 <- neEffdecomp(neMod1)
print(summary(effdecomp1))

NDE = neMod1$bootRes$t0["rntHbA1c.scaled0"]
NIE = neMod1$bootRes$t0["rntHbA1c.scaled1"]
propMed <- NIE/(NDE+NIE)
cat('\nproportion of mediated = ',propMed,'\n',sep="" )
# proportion of mediated = 0.08543659

colnames(neMod1$bootRes$t) <- names(neMod1$bootRes$t0)
head(neMod1$bootRes$t)

NDE.B = neMod1$bootRes$t[,"rntHbA1c.scaled0"]
NIE.B = neMod1$bootRes$t[,"rntHbA1c.scaled1"]
propMed.B = NIE.B/(NDE.B+NIE.B)
CI95.Boot = quantile(propMed.B,prob=c(0.025,(1-0.025)))
cat("Bootstrapping 95% CI for the porportion mediated\n",CI95.Boot)
# Bootstrapping 95% CI for the porportion mediated
# 0.05464952 0.1266088
hist(propMed.B)
abline(v=c(propMed,CI95.Boot),col="red",lty=2)

lht <- neLht(neMod1, linfct=c('rntHbA1c.scaled1 =0','rntHbA1c.scaled0=0'))
confint(lht)
?neLht
citation("medflex")
