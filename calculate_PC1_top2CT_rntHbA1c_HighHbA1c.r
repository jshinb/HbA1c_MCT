#------------------------------------------------------------------------------#
# Created: Sep 30, 2021
# 
# This script will 
# 1. calculate PC1 of top2CT and rntHbA1c in High-HbA1c group
# 2. generate PLINK input files (that will be run in SciNet)
# 
#------------------------------------------------------------------------------#

options(stringsAsFactors = F)
x = c('stringr','tidyverse','dplyr',
      'here','data.table','psych',#table
      'GenABEL','mgcv',#rntransform
      'tableone','arsenal',#table
      'ggplot2','patchwork','pheatmap','viridis')#visualization
lapply(x,require,character.only=T);rm(x)

# wd and load data ------------------------------------------------------------
wd = '~/OneDrive - SickKids/ukbb_insulin_resistance'
setwd(wd)

## data with adjusted MCT
load(file.path(wd,"idata_wi_roiCT_adjBase.Rd"))

## genotype PC for population stratification:22009-0.1 - 22009-0.5
ukb37194 = fread('~/OneDrive - SickKids/ukbb/data/ukb37194.csv.gz')
genoPC = subset(ukb37194,eid %in% idata$eid[idata$HbA1c.Level3=="High"], select=c('eid',paste("22009-0.",1:5,sep='')))##1070(F) + 3374(M)
names(genoPC)[-1] = paste('genoPC',1:5,sep='')

# PCA of average of top 2 regions: superior temporal and superior frontal -----
require(FactoMineR)
idata = idata %>% mutate(top.2.MCT = apply(data.frame(subset(idata,select=c(superiorfrontal,superiortemporal))),1,mean,na.rm=F))
idata = merge(idata,genoPC)
idata = idata %>% mutate(genoPC1=scale(genoPC1),genoPC2=scale(genoPC2),genoPC3=scale(genoPC3),genoPC4=scale(genoPC4),genoPC5=scale(genoPC5))
idata.F = na.omit(subset(idata,select=c(eid,top.2.MCT,Time.c,MRI_Site,Age.c,Sex,HbA1c.Level3,genoPC1,genoPC2,genoPC3,genoPC4,genoPC5),Sex=="Female"))
idata.M = na.omit(subset(idata,select=c(eid,top.2.MCT,Time.c,MRI_Site,Age.c,Sex,HbA1c.Level3,genoPC1,genoPC2,genoPC3,genoPC4,genoPC5),Sex=="Male"))#3369

tmp.High = subset(idata,HbA1c.Level3=="High")#4444
base::prop.table(table(tmp.High$HbA1c.Level3, tmp.High$Sex))#76%(M), 24%(F)

library(mgcv)
modF_BaseModH <- gam(top.2.MCT ~ s(Age.c) + s(Time.c,k=5) + MRI_Site + s(genoPC1)+ s(genoPC2)+ s(genoPC3)+ s(genoPC4)+ s(genoPC5),
                     data = idata.F, subset=HbA1c.Level3 =="High")
modM_BaseModH <- gam(top.2.MCT ~ s(Age.c) + s(Time.c,k=5)  + MRI_Site + s(genoPC1)+ s(genoPC2)+ s(genoPC3)+ s(genoPC4)+ s(genoPC5),
                     data = idata.M, subset=HbA1c.Level3 =="High")
gam.check(modF_BaseModH)
gam.check(modM_BaseModH)

idata.F = subset(idata.F,HbA1c.Level3=="High")
idata.M = subset(idata.M,HbA1c.Level3=="High")
idata.F$top2.MCT_adjBase = residuals(modF_BaseModH) + coef(modF_BaseModH)[['(Intercept)']]
idata.M$top2.MCT_adjBase = residuals(modM_BaseModH)+ coef(modM_BaseModH)[['(Intercept)']]
idata.MF = rbind(idata.M,idata.F)#4439

if( any(str_detect(names(idata),"top2.MCT_adjBase")) )
  idata = subset(idata,select=c(-top2.MCT_adjBase))

idataH = merge(idata,subset(idata.MF,select=c(eid,top2.MCT_adjBase)),sort=F)
summary(lm(top2.MCT_adjBase~Sex*Age.c,data=idataH))

# correlation matrix plot -----------------------------------------------------
library(corrplot)
roi1=which(names(idata)=="bankssts");roi34=which(names(idata)=="insula")
roi =  names(idata)[roi1:roi34]
cormat = cor(subset(idataH,select=c(roi,'top.2.MCT','top2.MCT_adjBase','Age.c')),use='p')
corp1 = corrplot(cormat,method='shade',order='hclust',hclust.method = 'ward.D2',
                 tl.cex = 0.75,tl.col = 'grey20')
cormatF = cor(subset(idataH,Sex=="Female",
                     select=c(roi,'top.2.MCT','top2.MCT_adjBase','Age.c')),use='p')
cormatF = cormatF[rownames(corp1$corr),rownames(corp1$corr)]
corrplot(cormatF,method='shade',#order='hclust',hclust.method = 'ward.D2',
         tl.cex = 0.75,tl.col = 'grey20')
cormatM = cor(subset(idataH,Sex=="Male",select=c(roi,'top.2.MCT','top2.MCT_adjBase','Age.c')),use='p')
cormatM = cormatM[rownames(corp1$corr),rownames(corp1$corr)]
corrplot(cormatM,method='shade',#order='hclust',hclust.method = 'ward.D2',
         tl.cex = 0.75,tl.col = 'grey20')
# -----------------------------------------------------------------------------
tmp.High = idataH;rm(idataH);gc(reset=T)
(print(table(tmp.High$HbA1c.Level3, tmp.High$Sex)))
(base::prop.table(table(tmp.High$HbA1c.Level3, tmp.High$Sex)))
#      Male Female
# High 3369 1070
# High 0.76 0.24
pca.dat = subset(tmp.High,select=c(eid,Sex,top2.MCT_adjBase,HbA1c2,rntHbA1c))
pca.dat = pca.dat %>% mutate(rntHbA1c.high = rntransform(HbA1c2))

p1=pca.dat %>% ggplot(aes(x=HbA1c2,color=Sex,fill=Sex))+
  geom_histogram(alpha=0.5) + theme(legend.position=c(0.95,0.95),legend.justification = c('right','top'))
p2=pca.dat %>% ggplot(aes(x=rntHbA1c.high,color=Sex,fill=Sex))+
  geom_histogram(alpha=0.5)+ theme(legend.position="none")
p3=pca.dat %>% ggplot(aes(x=rntHbA1c,color=Sex,fill=Sex))+
  geom_histogram(alpha=0.5)+ theme(legend.position="none")
p1+p2+p3

rownames(pca.dat) = pca.dat$eid
PC.A = PCA(na.omit(subset(pca.dat,select=c(top2.MCT_adjBase,rntHbA1c))),graph = F)
PC.F = PCA(na.omit(subset(pca.dat,Sex=="Female",select=c(top2.MCT_adjBase,rntHbA1c))),
           graph = F)
PC.M = PCA(na.omit(subset(pca.dat,Sex=="Male",select=c(top2.MCT_adjBase,rntHbA1c))),
           graph = F)
plot(PC.A, choix="var")
plot(PC.F, choix="var")
plot(PC.M,choix='var')

pc1.loadings = data.frame( rbind(
"Sex-Combined"=PC.A$var$coord[,1],
"Male"=PC.M$var$coord[,1],
"Female"=PC.F$var$coord[,1]) )

#              top2.MCT_adjBase   rntHbA1c
# Sex-Combined       -0.7460768  0.7460768
# Male                0.7559875 -0.7559875
# Female              0.7536307 -0.7536307
PC.A$ind$coord[,1] = sign(PC.A$var$coord[,1])[1]*PC.A$ind$coord[,1]; PC.A$var$coord[,1]=sign(PC.A$var$coord[1,1])*PC.A$var$coord[,1]
PC.M$ind$coord[,1] = sign(PC.M$var$coord[,1])[1]*PC.M$ind$coord[,1]; PC.M$var$coord[,1]=sign(PC.M$var$coord[1,1])*PC.M$var$coord[,1]
PC.F$ind$coord[,1] = sign(PC.F$var$coord[,1])[1]*PC.F$ind$coord[,1]; PC.F$var$coord[,1]=sign(PC.F$var$coord[1,1])*PC.F$var$coord[,1]

PC1.sex = rbind(data.frame(eid = rownames(PC.M$ind$coord), PC1.sex = PC.M$ind$coord[,1]),#12802 
                data.frame(eid = rownames(PC.F$ind$coord), PC1.sex = PC.F$ind$coord[,1]))#4338
PC1.A = data.frame(eid = rownames(PC.A$ind$coord), PC1 = PC.A$ind$coord[,1])#8464

PC1.data = merge(PC1.A,PC1.sex)
PC1.data = merge(subset(tmp.High,select=c(eid,Sex)),PC1.data)
with(PC1.data,plot(PC1,PC1.sex,pch=20));abline(0,1,col="red")
with(subset(PC1.data,Sex=="Male"),points(PC1,PC1.sex,col="orange"))
with(subset(PC1.data,Sex=="Male"),cor(PC1,PC1.sex))#0.999988
with(subset(PC1.data,Sex=="Female"),cor(PC1,PC1.sex))#0.9887885

PC1.data = PC1.data %>% dplyr::rename(PC1.all=PC1) %>% dplyr::rename(PC1=PC1.sex)
PC1.data %>% ggplot(aes(x=PC1.all,y=PC1,color=Sex)) + geom_point() + geom_abline(slope=1,intercept = 0)

tmp.PC1.data = merge(dplyr::select(PC1.data,-Sex),tmp.High,sort=F)
cormat = cor(dplyr::select(tmp.PC1.data, roi,top.2.MCT,top2.MCT_adjBase,Age.c,PC1,PC1.all),use='p')
corp1 = corrplot(cormat,method='shade',order='hclust',hclust.method = 'ward.D2',
                 tl.cex = 0.75,tl.col = 'grey20')
#
cormatF = cor(dplyr::select(subset(tmp.PC1.data,Sex=="Female"),
                            roi,top.2.MCT,top2.MCT_adjBase,Age.c,PC1,PC1.all),use='p')
cormatF = cormatF[rownames(corp1$corr),rownames(corp1$corr)]
corrplot(cormatF,method='shade',#order='hclust',hclust.method = 'ward.D2',
         tl.cex = 0.75,tl.col = 'grey20')
#
cormatM = cor(dplyr::select(subset(tmp.PC1.data,Sex=="Male"),
                            roi,top.2.MCT,top2.MCT_adjBase,Age.c,PC1,PC1.all),use='p')
cormatM = cormatM[rownames(corp1$corr),rownames(corp1$corr)]
corrplot(cormatM,method='shade',#order='hclust',hclust.method = 'ward.D2',
         tl.cex = 0.75,tl.col = 'grey20')
# -----------------------------------------------------------------------------
## write out files
to_write = dplyr::select(tmp.PC1.data,eid,Sex,top2.MCT_adjBase,PC1) %>% mutate(rntPC1=rntransform(PC1))
to_write$Sex_Male0 = 0
to_write$Sex_Male0[to_write$Sex=="Female"] = 1

#6. Write pheno and covfiles for PLINK analysis =========================================
# original: untransformed
# sensitivity: rnt-transformed PC1
pheno_out = to_write %>% mutate(FID = eid) %>% 
  dplyr::rename(IID = eid) %>%
  dplyr::select(FID,IID,top2.MCT_adjBase,PC1,Sex_Male0,Sex)

describe(pheno_out)
# original
# write_delim(subset(pheno_out,select=c(FID,IID,top2.MCT_adjBase,PC1)),
#             file='pheno_PC1_HighHbA1c_top2_CT_adjBase_all.txt',delim=" ")
# 
# write_delim(subset(pheno_out,Sex=="Male",select=c(FID,IID,top2.MCT_adjBase,PC1)),
#             file='pheno_PC1_HighHbA1c_top2_CT_adjBase_male.txt',delim=" ")
# 
# write_delim(subset(pheno_out,Sex=="Female",select=c(FID,IID,top2.MCT_adjBase,PC1)),
#             file='pheno_PC1_HighHbA1c_top2_CT_adjBase_female.txt',delim=" ")

# sensitivity - with rntPC1
pheno_out = pheno_out %>% mutate(PC1.org = PC1)
pheno_out = pheno_out %>% mutate(PC1 = rntransform(PC1.org))
p1 = ggplot(pheno_out,aes(x=PC1.org)) + geom_histogram()
p2 = ggplot(pheno_out,aes(x=PC1)) + geom_histogram()
p1+p2

write_delim(subset(pheno_out,select=c(FID,IID,PC1)),
            file='pheno_PC1_HighHbA1c_top2_CT_adjBase_all.txt',delim=" ")

write_delim(subset(pheno_out,Sex=="Male",select=c(FID,IID,PC1)),
            file='pheno_PC1_HighHbA1c_top2_CT_adjBase_male.txt',delim=" ")

write_delim(subset(pheno_out,Sex=="Female",select=c(FID,IID,PC1)),
            file='pheno_PC1_HighHbA1c_top2_CT_adjBase_female.txt',delim=" ")

write_delim(subset(pheno_out,select=c(FID,IID,Sex_Male0)),
            file='cov_PC1_PC1_HighHbA1c_top2_CT_adjBase_all.txt',delim=" ")

#6. Write fam for PLINK analysis ==============================================
cleaned_fam = fread('~/OneDrive - SickKids/ukbb/data/ukb37194_cleaned_famfile.csv')
head(cleaned_fam)
write_delim(subset(cleaned_fam,V1 %in% pheno_out$FID),
            'famfile_PC1_HighHbA1c_adjBase_all.fam',delim=" ")
