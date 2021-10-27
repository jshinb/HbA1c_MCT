#---------------------------------------------------------------------------------------#
# 
# This script will adjust MCT for additional covariates: EduLevel, BMI, Smoking, SBP and 
# all of them.
# 
# An updated data file will be created: 
# "~/OneDrive - SickKids/ukbb_insulin_resistance/idata_wi_MCT_adjAll.Rd"
# 
#---------------------------------------------------------------------------------------#
options(stringsAsFactors = F)
x = c('stringr','tidyverse','dplyr',
      'here','data.table','psych',#table
      'GenABEL','mgcv',#rntransform
      'tableone','arsenal',#table
      'ggplot2','patchwork','pheatmap')#visualization
lapply(x,require,character.only=T);rm(x)

wd = '~/OneDrive - SickKids/ukbb_insulin_resistance'
load(file.path(wd,"idata_wi_MCT_adjBase.Rd"))
setwd(wd)

# EduLevel---------------------------------------------------------------------
idata$MCT_adjBaseEduLevel <- NA
fitF=lm(MCT_adjBase~EduLevel,data=idata,subset = Sex=="Female",
        na.action = na.exclude)
summary(fitF)
idata$MCT_adjBaseEduLevel[idata$Sex=="Female"] <- resid(fitF) + fitF$coefficients[["(Intercept)"]]
rm(fitF)

fitM=lm(MCT_adjBase~EduLevel,data=idata,subset = Sex=="Male",
        na.action = na.exclude)
summary(fitM)
idata$MCT_adjBaseEduLevel[idata$Sex=="Male"] <- resid(fitM) + fitM$coefficients[["(Intercept)"]]
rm(fitM)

## plotting
p = ggplot(subset(idata,!is.na(EduLevel)), 
           aes(x=Sex, y=rntHbA1c,color = EduLevel, by=EduLevel)) +
  geom_boxplot() 

p.MCT = ggplot(subset(idata,!is.na(EduLevel)), 
               aes(x=Sex, y=MCT_adjBase,color = EduLevel, by=EduLevel)) +
  geom_boxplot()

p.MCT.adj = ggplot(subset(idata,!is.na(EduLevel)), 
                   aes(x=Sex, y=MCT_adjBaseEduLevel,color = EduLevel, by=EduLevel)) +
  geom_boxplot()
grid.arrange(p,p.MCT,p.MCT.adj)


# BMI--------------------------------------------------------------------------
## adjust MCT for log-BMI using gam
idata$logBMI = log(idata$BMI_recruitment)
idata <- idata %>% mutate(logBMI = logBMI-mean(logBMI,na.rm=T))

idata.F = na.omit(subset(idata,select=c(eid,MCT_adjBase,Sex,logBMI),Sex=="Female"))
idata.M = na.omit(subset(idata,select=c(eid,MCT_adjBase,Sex,logBMI),Sex=="Male"))

modF_BaseBMI <- gam(MCT_adjBase ~ s(logBMI),
                    data = na.omit(subset(idata,select=c(MCT_adjBase,Sex,logBMI))), subset=Sex=="Female")
modM_BaseBMI <- gam(MCT_adjBase ~ s(logBMI),
                    data = na.omit(subset(idata,select=c(MCT_adjBase,Sex,logBMI))), subset=Sex=="Male")
idata.F$MCT_adjBaseBMI = residuals.gam(modF_BaseBMI,type = 'response') + coef(modF_BaseBMI)[['(Intercept)']]
idata.M$MCT_adjBaseBMI = residuals.gam(modM_BaseBMI, type = 'response')+ coef(modM_BaseBMI)[['(Intercept)']]

# modF_BaseBMI <- lm(MCT_adjBase ~ logBMI + I(logBMI^2), 
#                    data = na.omit(subset(idata,select=c(MCT_adjBase,Sex,logBMI))), subset=Sex=="Female")
# modM_BaseBMI <- lm(MCT_adjBase ~ logBMI + I(logBMI^2), 
#                    data = na.omit(subset(idata,select=c(MCT_adjBase,Sex,logBMI))), subset=Sex=="Male")
# idata.F$MCT_adjBaseBMI = resid(modF_BaseBMI) + coef(modF_BaseBMI)[['(Intercept)']]
# idata.M$MCT_adjBaseBMI = resid(modM_BaseBMI) + coef(modM_BaseBMI)[['(Intercept)']]

idata.MF = rbind(idata.M,idata.F)

if(any(str_detect(names(idata),"MCT_adjBaseBMI")) )
  idata = subset(idata,select=c(-MCT_adjBaseBMI))

idata = merge(
  idata,
  subset(idata.MF,select=c(eid,MCT_adjBaseBMI)),
  all.x=T,sort=F)

## checking
fit = lm(MCT_adjBaseBMI ~ log(BMI_recruitment) + Sex*log(BMI_recruitment),
         data=idata,na.action=na.exclude)
summary(fit)# everything is insignificant

p.BMI = idata %>% 
  ggplot(aes(x=logBMI,y=MCT_adjBaseBMI,color=Sex,fill=Sex))+
  geom_smooth()+
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue")) +
  theme(legend.position = "top")

p = idata %>% ggplot(aes(x=BMI_recruitment,y=rntHbA1c,color=Sex,fill=Sex))+
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  theme(legend.position = 'none') + 
  coord_trans(x="log") + 
  geom_vline(xintercept = c(18.5,25,30),linetype=2)

p1 = p + geom_point(alpha=0.3) + geom_smooth() 
p2 = p + geom_smooth()

p.MCT = idata %>% ggplot(aes(x=BMI_recruitment,y=MCT_adjBase,color=Sex,fill=Sex))+
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  theme(legend.position = 'none') +
  coord_trans(x="log") + 
  geom_vline(xintercept = c(18.5,25,30),linetype=2)

p1.MCT = p.MCT + geom_point(alpha=0.3) + geom_smooth()
p2.MCT = p.MCT + geom_smooth()

grid.arrange(p1,p2,p1.MCT,p2.MCT,nrow=2)

png(paste("Plots/logBMI_MCT_rntHbA1c.png",sep=""),
    width=8,height=3.25, units = "in", res=300)
grid.arrange(p1.MCT,p2.MCT,nrow=1)
dev.off()
rm(p.BMI,p1.MCT,p2.MCT,p1,p2)

# Smoking----------------------------------------------------------------------
idata$MCT_adjBaseSmoking <- NA
fitF=lm(MCT_adjBase~Smoking,data=idata,subset = Sex=="Female",
        na.action = na.exclude)
summary(fitF)
idata$MCT_adjBaseSmoking[idata$Sex=="Female"] <- resid(fitF) + fitF$coefficients[["(Intercept)"]]
rm(fitF)

fitM=lm(MCT_adjBase~Smoking,data=idata,subset = Sex=="Male",
        na.action = na.exclude)
summary(fitM)
idata$MCT_adjBaseSmoking[idata$Sex=="Male"] <- resid(fitM) + fitM$coefficients[["(Intercept)"]]
rm(fitM)

p = ggplot(subset(idata,!is.na(Smoking)), 
           aes(x=Sex, y=rntHbA1c,color=Smoking, by=Smoking)) +
  geom_boxplot() 

p.MCT = ggplot(subset(idata,!is.na(Smoking)), 
               aes(x=Sex, y=MCT_adjBase,color=Smoking, by=Smoking)) +
  geom_boxplot()
p.MCT.adj = ggplot(subset(idata,!is.na(Smoking)), 
                   aes(x=Sex, y=MCT_adjBaseSmoking,color=Smoking, by=Smoking)) +
  geom_boxplot()

grid.arrange(p,p.MCT,p.MCT.adj)
rm(p,p.MCT,p.MCT.adj)

# SBP0-------------------------------------------------------------------------
idata = idata %>% mutate(SBP.c = SBP0-mean(SBP0,na.rm=T))
idata.F = na.omit(subset(idata,select=c(eid,MCT_adjBase,Sex,SBP0,SBP.c),Sex=="Female"))
idata.M = na.omit(subset(idata,select=c(eid,MCT_adjBase,Sex,SBP0,SBP.c),Sex=="Male"))

modF_BaseSBP <- gam(MCT_adjBase ~ s(SBP.c), data = idata.F)
modM_BaseSBP <- gam(MCT_adjBase ~ s(SBP.c), data = idata.M)
gam.check(modF_BaseSBP)
gam.check(modM_BaseSBP)

idata.F$MCT_adjBaseSBP = residuals(modF_BaseSBP) + coef(modF_BaseSBP)[['(Intercept)']]
idata.M$MCT_adjBaseSBP = residuals(modM_BaseSBP)+ coef(modM_BaseSBP)[['(Intercept)']]
idata.MF = rbind(idata.M,idata.F)

if(any(str_detect(names(idata),"MCT_adjBaseSBP")) )
  idata = subset(idata,select=c(-MCT_adjBaseSBP))

idata = merge(
  idata,
  subset(idata.MF,select=c(eid,MCT_adjBaseSBP)),
  all.x=T,sort=F)

## checking
fit = lm(MCT_adjBaseSBP ~ Sex*SBP, data=idata,na.action=na.exclude)
summary(fit)# everything is insignificant

## plot
p = idata %>% ggplot(aes(x=SBP,y=rntHbA1c,color=Sex,fill=Sex))+
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  theme(legend.position = 'none')

p1 = p + geom_point(alpha=0.3) + geom_smooth()
p2 = p + geom_smooth()

p.MCT = idata %>% ggplot(aes(x=SBP,y=MCT_adjBase,color=Sex,fill=Sex))+
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  theme(legend.position = 'none')

p1.MCT = p.MCT + geom_point(alpha=0.3) + geom_smooth()
p2.MCT = p.MCT + geom_smooth()

idata %>% ggplot(aes(x=SBP,y=MCT_adjBaseSBP,color=Sex,fill=Sex))+
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  theme(legend.position = 'none') + geom_smooth() #+ geom_point(alpha=0.1)

idata %>% ggplot(aes(x=rntHbA1c,y=MCT_adjBaseSBP,color=Sex,fill=Sex))+
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  theme(legend.position = 'none') + geom_point(alpha=0.1) + geom_smooth()

grid.arrange(p1,p2,p1.MCT,p2.MCT,nrow=2)

png(paste("Plots/SBP_MCT_rntHbA1c.png",sep=""),
    width=8,height=3.25, units = "in", res=300)
grid.arrange(p1.MCT,p2.MCT,nrow=1)
dev.off()

# SBP1-------------------------------------------------------------------------
idata = idata %>% mutate(SBP1.c = SBP1-mean(SBP1,na.rm=T))
idata.F = na.omit(subset(idata,select=c(eid,MCT_adjBase,Sex,SBP1,SBP1.c),Sex=="Female"))
idata.M = na.omit(subset(idata,select=c(eid,MCT_adjBase,Sex,SBP1,SBP1.c),Sex=="Male"))

modF_BaseSBP1 <- gam(MCT_adjBase ~ s(SBP1.c), data = idata.F)
modM_BaseSBP1 <- gam(MCT_adjBase ~ s(SBP1.c), data = idata.M)

idata.F$MCT_adjBaseSBP1 = residuals(modF_BaseSBP1) + coef(modF_BaseSBP1)[['(Intercept)']]
idata.M$MCT_adjBaseSBP1 = residuals(modM_BaseSBP1)+ coef(modM_BaseSBP1)[['(Intercept)']]
idata.MF = rbind(idata.M,idata.F)

if(any(str_detect(names(idata),"MCT_adjBaseSBP1")) )
  idata = subset(idata,select=c(-MCT_adjBaseSBP1))

idata = merge(
  idata,
  subset(idata.MF,select=c(eid,MCT_adjBaseSBP1)),
  all.x=T,sort=F)

## checking
fit = lm(MCT_adjBaseSBP1 ~ Sex*SBP1, data=idata,na.action=na.exclude)
summary(fit)# everything is insignificant

## plot
p = idata %>% ggplot(aes(x=SBP1,y=rntHbA1c,color=Sex,fill=Sex))+
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  theme(legend.position = 'none')

p + geom_point(alpha=0.3) + geom_smooth()
p + geom_smooth()

p = idata %>% ggplot(aes(x=SBP1,y=MCT_adjBase,color=Sex,fill=Sex))+
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  theme(legend.position = 'none')
grid.arrange( (p + geom_point(alpha=0.3) + geom_smooth()),
              (p + geom_smooth()), ncol=2 )

p1.MCT = idata %>% ggplot(aes(x=SBP1,y=MCT_adjBase,color=Sex,fill=Sex))+
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  theme(legend.position = 'none') + geom_smooth() #+ geom_point(alpha=0.1)

p2.MCT = idata %>% ggplot(aes(x=rntHbA1c,y=MCT_adjBaseSBP1,color=Sex,fill=Sex))+
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  theme(legend.position = 'none') + geom_point(alpha=0.1) + geom_smooth()
p1.MCT+p2.MCT

png(paste("Plots/SBP1_MCT_rntHbA1c.png",sep=""),
    width=8,height=3.25, units = "in", res=300)
print(p1.MCT+p2.MCT)
dev.off()
rm(p1.MCT,p2.MCT,p1,p2);gc(reset=T)

# AllCov-----------------------------------------------------------------------
idata.F = na.omit(subset(idata,select=c(eid,MCT,Age.c,Time.c,MRI_Site,SBP.c,MCT_adjBase,Sex,EduLevel,logBMI,Smoking,SBP),Sex=="Female"))
idata.M = na.omit(subset(idata,select=c(eid,MCT,Age.c,Time.c,MRI_Site,SBP.c,MCT_adjBase,Sex,EduLevel,logBMI,Smoking,SBP),Sex=="Male"))
modF_adjall <- gam(MCT ~ s(Age.c) + s(Time.c) + MRI_Site + EduLevel + s(logBMI) + Smoking + s(SBP.c), data = idata.F)
modM_adjall <- gam(MCT ~ s(Age.c) + s(Time.c) + MRI_Site + EduLevel + s(logBMI) + Smoking + s(SBP.c), data = idata.M)

idata.F$MCT_adjall = residuals(modF_adjall) + coef(modF_adjall)[['(Intercept)']]
idata.M$MCT_adjall = residuals(modM_adjall)+ coef(modM_adjall)[['(Intercept)']]
idata.MF = rbind(idata.M,idata.F)

if(any(str_detect(names(idata),"MCT_adjall")) )
  idata = subset(idata,select=c(-MCT_adjall))

idata = merge(
  idata,
  subset(idata.MF,select=c(eid,MCT_adjall)),
  all.x=T,sort=F)

## checking
fit = lm(MCT_adjall ~ Sex * (Age.c + logBMI + SBP.c + MRI_Site + Smoking + EduLevel), data=idata,na.action=na.exclude)
summary(fit)# everything is insignificant

## plot
p = idata %>% ggplot(aes(x=SBP,y=MCT_adjall,color=Sex,fill=Sex))+
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  theme(legend.position = 'none')

p1 = p + geom_point(alpha=0.3) + geom_smooth()
p2 = p + geom_smooth()

p.MCT = idata %>% ggplot(aes(x=SBP1,y=MCT_adjall,color=Sex,fill=Sex))+
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  theme(legend.position = 'none')

p1 = p.MCT + geom_point(alpha=0.3) + geom_smooth()
p2 = p.MCT + geom_smooth()
p1+p2;rm(p1,p2)

save(idata,file=file.path(wd,"idata_wi_MCT_adjAll.Rd"))
