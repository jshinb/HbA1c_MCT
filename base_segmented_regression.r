#---------------------------------------------------------------------------------------#
# 
# This script will adjust MCT for the base-model covariates (age, MRI site and time 
# between blood samplong and MRI imaging), using generalize additive modelling and fit a 
# segmented linear regression model with the adjusted MCT in males and females separately.
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
load(file.path(wd,"idata.Rd"))

## ----MCTadjTime_Age_MRI_Site-----------------------------------------------------------
# REF:https://pubmed.ncbi.nlm.nih.gov/15580601/
setwd(wd)
idata.F = na.omit(subset(idata,select=c(eid,MCT,Time.c,MRI_Site,Age.c,Sex),Sex=="Female"))
idata.M = na.omit(subset(idata,select=c(eid,MCT,Time.c,MRI_Site,Age.c,Sex),Sex=="Male"))

library(mgcv)
modF_BaseMod <- gam(MCT ~ s(Age.c) + s(Time.c) + MRI_Site,
                    data = na.omit(subset(idata,select=c(MCT,Time.c,Sex,MRI_Site,Age.c))), subset=Sex=="Female")
modM_BaseMod<- gam(MCT ~ s(Age.c) + s(Time.c)  + MRI_Site,
                   data = na.omit(subset(idata,select=c(MCT,Time.c,Sex,MRI_Site,Age.c))), subset=Sex=="Male")
gam.check(modF_BaseMod)
gam.check(modM_BaseMod)

# 
idata.F$MCT_adjBase = residuals(modF_BaseMod) + coef(modF_BaseMod)[['(Intercept)']]
idata.M$MCT_adjBase = residuals(modM_BaseMod)+ coef(modM_BaseMod)[['(Intercept)']]
idata.MF = rbind(idata.M,idata.F)

if(any(str_detect(names(idata),"MCT_adjBase")) )
  idata = subset(idata,select=c(-MCT_adjBase))

idata = merge(idata,subset(idata.MF,select=c(eid,MCT_adjBase)),
              all.x=T,sort=F)

## checking
fit = lm(MCT_adjBase ~ Age_MRI + I(Age.c^2) + Time + MRI_Site + Sex*(Age_MRI + I(Age.c^2) + Time + MRI_Site),
         data=idata,na.action=na.exclude)
summary(fit)# everything is insignificant

idata %>% 
  ggplot(aes(x=Age_MRI,y=MCT_adjBase,color=Sex,fill=Sex))+
  geom_smooth()+
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue")) +
  theme(legend.position = "top")

idata %>% 
  ggplot(aes(x=Time,y=MCT_adjBase,color=Sex,fill=Sex))+
  geom_smooth()+
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue")) +
  theme(legend.position = "top") + 
  xlab("Time between blood draw and MRI (years)")

# png(paste("Plots/MCTadjBase_by_Time_and_Age.png",sep=""),
#     width=6,height=4.5, units = "in", res=300)
# grid.arrange(p.Time,p.age,nrow=1)
# dev.off()

##
p1 = idata %>% ggplot(aes(x=rntHbA1c,y=MCT_adjBase,color=Sex)) +
  geom_point(alpha=0.2) + 
  geom_smooth() + 
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue"))+
  xlab(NULL) + ylab(NULL) + 
  theme(legend.position = "top")
p2 = idata %>% ggplot(aes(x=rntHbA1c,y=MCT_adjBase,color=Sex)) +
  geom_smooth() + 
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue"))+
  xlab(NULL) + ylab(NULL) + theme(legend.position = "top")  

png(paste("Plots/MCTadjBase_by_HbA1c.png",sep=""),
    width=6,height=4.5, units = "in", res=300)
print(p1+ theme(legend.position = "none"))
dev.off()

png(paste("Plots/MCTadjBase_by_HbA1c_wo_points.png",sep=""),
    width=6,height=4.5, units = "in", res=300)
print(p2+ theme(legend.position = "none"))
dev.off()
rm(p.age,p.Time,p1,p2);gc(reset=T)

## ----segmented_regression--------------------------------------------------------------
library(segmented)
idata$MCT_scaled = scale(idata$MCT_adjBase)

## male
fitM = lm(MCT_adjBase~rntHbA1c,data=idata,Sex=="Male")
summary(fitM)$coef
my.segM <- segmented(fitM,seg.Z = ~ rntHbA1c)#37.3
my.segM$psi
summary(my.segM)
slope(my.segM)
my.fitted <- fitted(my.segM)
idata.M = na.omit(subset(idata,select = c(eid,Sex,rntHbA1c,MCT_scaled), Sex=="Male"))
my.model.M <- data.frame(rntHbA1c = idata.M$rntHbA1c, MCT = my.fitted)

## female
idataF = subset(idata,Sex=="Female")

psiF = NULL
fitF = lm(MCT_adjBase~rntHbA1c,data=idataF)
summary(fitF)$coef
my.segF <- segmented(fitF,seg.Z = ~ rntHbA1c, psi = 0)#40.1
psiF = rbind(psiF,my.segF$psi)
# for(bi in 1:999){
# inds = sample(1:nrow(idataF),replace=T)
# fitF = lm(MCT_adjBase~rntHbA1c,data=idataF[inds,])
# summary(fitF)$coef
# my.segF <- segmented(fitF,seg.Z = ~ rntHbA1c)#40.1
# psiF = rbind(psiF,my.segF$psi)
# }
summary(my.segF)
slope(my.segF)
my.fitted <- fitted(my.segF)
idata.F = na.omit(subset(idata,select = c(eid,Sex,rntHbA1c,MCT_scaled), Sex=="Female"))
my.model.F <- data.frame(rntHbA1c = idata.F$rntHbA1c, MCT = my.fitted)

## plot the fitted model
pA = ggplot(my.model.F, aes(x = rntHbA1c, y = MCT)) + 
  geom_line(color="darkred")+
  geom_line(data=my.model.M,aes(x = rntHbA1c, y = MCT),color="darkblue") 
print(pA)
rm(pA);gc(reset=T)

## ----MCT_HbA1c_vs_HbA1c_tests----------------------------------------------------------
(cutoffM = my.segM$psi[1,"Est."])#0.7715445 (wo Age2) #0.7714078 (wi Age2)
d = subset(idata,Sex=="Male")
d = d[order(d$rntHbA1c),]
cat("cutoff-Male (mmol/mol)\n")
plot(d$rntHbA1c,abs(d$rntHbA1c-cutoffM))
# cutoffM = d$rntHbA1c[which.min(abs(d$rntHbA1c-cutoffM))];
(HbA1c.cutoffM = d$HbA1c2[which.min(abs(d$rntHbA1c-cutoffM))])#37.3

(cutoffF = my.segF$psi[1,"Est."])#1.406941 (wo Age2), 1.425325 (wi s(Age)) #1.43121/1.431371
d = subset(idata,Sex=="Female")
d = d[order(d$rntHbA1c),]
plot(d$rntHbA1c,abs(d$rntHbA1c-cutoffF))
# (cutoffF = d$rntHbA1c[which.min(abs(d$rntHbA1c-cutoffF))])
cat("cutoff-Female (mmol/mol)\n")
(HbA1c.cutoffF = d$HbA1c2[which.min(abs(d$rntHbA1c-cutoffF))])#40.1,#40

idata$HbA1c.Level3= NA
ind = !is.na(idata$rntHbA1c)
idata$HbA1c.Level3[ind & (idata$Sex=="Male" & idata$rntHbA1c>cutoffM)] <- "High"
idata$HbA1c.Level3[ind & (idata$Sex=="Male" & idata$rntHbA1c<=cutoffM)] <- "Normal"
idata$HbA1c.Level3[ind & (idata$Sex=="Female" & idata$rntHbA1c>cutoffF)] <- "High"
idata$HbA1c.Level3[ind & (idata$Sex=="Female" & idata$rntHbA1c<=cutoffF)] <- "Normal"
print(table(idata$HbA1c.Level3,idata$Sex))
#         Male Female
# High    3374   1070
# Normal 11158  14977
p.normal = subset(idata,HbA1c.Level3=="Normal") %>%
  ggplot(aes(x=rntHbA1c,y=MCT_adjBase,color=Sex,fill=Sex)) + 
  # geom_point(alpha=0.2)+
  geom_smooth()+
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue"))+
  theme(legend.position="top")+
  coord_cartesian(ylim = c(2.35,2.485))+
  ylab(NULL)

p.high = subset(idata,HbA1c.Level3=="High") %>%
  ggplot(aes(x=rntHbA1c,y=MCT_adjBase,color=Sex,fill=Sex)) + 
  geom_smooth() +
  scale_color_manual(values=c("Female"="darkred","Male"="darkblue")) + 
  theme(legend.position="top") +
  coord_cartesian(ylim = c(2.35,2.485)) +
  ylab(NULL)

png('Plots/HbA1c_MCTadjBase_Normal_High.png',
    width=8,height=4.5,units='in',res=300)
grid.arrange(p.normal,p.high,nrow=1)
dev.off()

##
fit.High = lm(MCT_scaled~rntHbA1c*Sex,data=idata,subset=HbA1c.Level3=="High")
summary(fit.High)
relaimpo::calc.relimp(fit.High)#TotalR2 = 2.52%: R2 = 0.79% (HbA1c), R2_sex =1.7%, joint: 0.034%
my.model.F$Sex = "Female"
my.model.M$Sex = "Male"

p.MCTadj.rntHbA1c <- idata %>% 
  ggplot(aes(x=rntHbA1c,y=MCT_adjBase,color=Sex,fill=Sex))+
  # geom_point(alpha=0.2)+
  geom_smooth(linetype=3,size=0.5) +
  scale_color_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  scale_fill_manual(values = c("Female"= 'darkred','Male'="darkblue")) + 
  geom_vline(xintercept = cutoffM, color='darkblue')+ #-0.1467856
  geom_vline(xintercept = cutoffF, color='darkred') + #0.2019919
  theme(legend.position="top") +
  geom_line(data=my.model.F,aes(x = rntHbA1c, y = MCT)) +
  geom_line(data=my.model.M,aes(x = rntHbA1c, y = MCT))

png('Plots/fitted_HbA1c_MCTadjBase_Normal_High.png',
    width=5,height=4.5,units='in',res=300)
p.MCTadj.rntHbA1c
dev.off()
rm(p.high,p.MCTadj.rntHbA1c,p.normal)

save(idata, file=file.path(wd,"idata_wi_MCT_adjBase.Rd"))

#------------------------------------------------------------------------------
idata = idata %>% mutate(x2 = as.numeric(HbA1c.Level3=="High"))
idata$x2.star = idata$x2
idata$x2.star[idata$Sex=="Male"] = (idata$rntHbA1c[idata$Sex=="Male"]-cutoffM)*idata$x2[idata$Sex=="Male"]
idata$x2.star[idata$Sex=="Female"] = (idata$rntHbA1c[idata$Sex=="Female"]-cutoffF)*idata$x2[idata$Sex=="Female"]
fitM = lm(MCT_adjBase ~ rntHbA1c + x2.star, data=idata, subset=Sex=="Male")
fitF = lm(MCT_adjBase ~ rntHbA1c + x2.star, data=idata, subset=Sex=="Female")
fitI = lm(MCT_adjBase ~ Sex*(rntHbA1c + x2.star), data=idata)
summary(fitM)
