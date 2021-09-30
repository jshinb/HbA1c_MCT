#------------------------------------------------------------------------------#
# Created: Sep 23, 2021
# 
# This script will fit segmented regression models roiMCTadj ~ HbA1c to the data.
# 
# segmented linear regression model with the adjusted MCT in males and females separately.
#
# 1. Should I do it in the sex-combined sample?
# 2. I could not reproduce the 'alladj' results - Why?
# Assuming that the slopes differ at the same break point identified in MCT
# 
#------------------------------------------------------------------------------#
options(stringsAsFactors = F)
x = c('stringr','tidyverse','dplyr',
      'here','data.table','psych',#table
      'GenABEL','mgcv',#rntransform
      'tableone','arsenal',#table
      'ggplot2','patchwork','pheatmap','viridis')#visualization
lapply(x,require,character.only=T);rm(x)

# functions -------------------------------------------------------------------
get_results = function(fitC,fitM,fitF,fitI){
  vnames=c('rntHbA1c2','x2.star')
  var.statC = sum(vcov(fitC)[vnames,vnames])
  slopeHC =sum(summary(fitC)$coef[vnames,"Estimate"])
  wald.testC  = (slopeHC/sqrt(var.statC))^2
  pC = pchisq(wald.testC,df = 1,lower.tail = F)
  
  var.statM = sum(vcov(fitM)[-1,-1])
  slopeHM =sum(summary(fitM)$coef[vnames,"Estimate"])
  wald.testM  = (slopeHM/sqrt(var.statM))^2
  pM = pchisq(wald.testM,df = 1,lower.tail = F)
  
  var.statF = sum(vcov(fitF)[-1,-1])
  slopeHF =sum(summary(fitF)$coef[vnames,"Estimate"])
  wald.testF  = (slopeHF/sqrt(var.statF))^2
  pF = pchisq(wald.testF,df = 1,lower.tail = F) #
  
  iterm_name = c('SexFemale:rntHbA1c2','SexFemale:x2.star')
  iterm_index = which(colnames(vcov(fitI)) %in% iterm_name)
  slopeI = sum(summary(fitI)$coef[iterm_name,"Estimate"])
  var.statI = sum(vcov(fitI)[iterm_index,iterm_index])
  wald.testI  = (slopeI/sqrt(var.statI))^2
  pI = pchisq(wald.testI,df = 1,lower.tail = F) #
  
  resC = data.table(Term = c("Intercept","beta.N","beta.diff.HN"),summary(fitC)$coef[,-3]) %>% rename(beta="Estimate", SE="Std. Error",P="Pr(>|t|)")
  resM = data.table(Term = c("Intercept","beta.N","beta.diff.HN"),summary(fitM)$coef[,-3]) %>% rename(beta="Estimate", SE="Std. Error",P="Pr(>|t|)")
  resF = data.table(Term = c("Intercept","beta.N","beta.diff.HN"),summary(fitF)$coef[,-3])%>% rename(beta="Estimate",SE="Std. Error",P="Pr(>|t|)")
  
  resC = rbind(resC,data.table(Term="beta.H",beta=slopeHC,SE=sqrt(var.statC),P=pC))
  resM = rbind(resM,data.table(Term="beta.H",beta=slopeHM,SE=sqrt(var.statM),P=pM))
  resF = rbind(resF,data.table(Term="beta.H",beta=slopeHF,SE=sqrt(var.statF),P=pF))

  resC = data.table(resC,Sex="Sex-Combined")
  resM = data.table(resM,Sex="Male")
  resF = data.table(resF,Sex="Female")
  
  resIH = data.table(Term="beta.IH",beta=slopeI,SE=sqrt(var.statI),P=pI,Sex="interaction")
  resIN = data.table(Term="beta.IN",
                     beta=summary(fitI)$coef["SexFemale:rntHbA1c2","Estimate"],
                     SE=summary(fitI)$coef["SexFemale:rntHbA1c2","Std. Error"],
                     P=summary(fitI)$coef["SexFemale:rntHbA1c2","Pr(>|t|)"],
                     Sex="interaction")
  res = rbind(resC,resM,resF,resIN,resIH)
  res
}

fit_function = function(yname,data,adj){
  analdat = data
  analdat$y = analdat[[yname]]
  fitC = lm(y ~ rntHbA1c2 + x2.star, data=analdat) #does this make sense?
  fitM = lm(y ~ rntHbA1c2 + x2.star, data=analdat, subset=Sex=="Male")
  fitF = lm(y ~ rntHbA1c2 + x2.star, data=analdat, subset=Sex=="Female")
  fitI = lm(y ~ Sex*(rntHbA1c2 + x2.star), data=analdat)
  
  res = data.table(get_results(fitC=fitC,fitM=fitM,fitF=fitF,fitI=fitI),adjustment=adj)
  res
} 

# wd and load data ------------------------------------------------------------
wd = '~/OneDrive - SickKids/ukbb_insulin_resistance'
setwd(wd)

## data with adjusted MCT
load(file.path(wd,"idata_wi_MCT_adjAll.Rd"))
cutoffM = 0.7714078
cutoffF = 1.425325
cutoff = (cutoffM + cutoffF)/2
cutdiff=cutoffF - cutoff

idata = idata %>% mutate(rntHbA1c2 = ifelse(Sex=="Male", rntHbA1c+cutdiff,rntHbA1c-cutdiff))
idata$x2 = as.numeric(idata$HbA1c.Level3=="High")
idata$x2.star = idata$x2
idata$x2.star[idata$Sex=="Male"] = (idata$rntHbA1c2[idata$Sex=="Male"]-cutoff)*idata$x2[idata$Sex=="Male"]
idata$x2.star[idata$Sex=="Female"] = (idata$rntHbA1c2[idata$Sex=="Female"]-cutoff)*idata$x2[idata$Sex=="Female"]

## ROI-MRI data
d_MRI2 = fread("Replication_ukb_MRI_2020-11-26.tsv")#17690
d_MRI1 = fread('ThicknessData.csv')

## ----roi_ctx-----------------------------------------------------------------
group.levels = c("Sex-Combined","Male","Female")
d_MRI2 = d_MRI2 %>% rename(eid=ID)
str(d_MRI2)

cols_extract = c('eid',names(d_MRI1)[str_detect(names(d_MRI1),"lh_")|str_detect(names(d_MRI1),"rh_")])
d_MRI = unique(rbind(subset(d_MRI1,!eid%in%d_MRI2$eid, select=cols_extract),subset(d_MRI2,select=cols_extract)))#38,664

roi = str_split(names(d_MRI)[str_detect(names(d_MRI),"lh_")],"_",simplify=T)[,2]
mean.ctx = subset(d_MRI,select='eid')
for(roii in roi){
  rois = paste(c('lh','rh'),roii,'thickness',sep="_")
  roi.ctx = apply(subset(data.frame(d_MRI),select=rois),1,mean,na.rm=T)
  roi.ctx[ !is.na(roi.ctx) & abs(scale(roi.ctx))>4 ] <- NA #outlier
  print(sum(!is.na(roi.ctx) & abs(scale(roi.ctx))>4))#8
  mean.ctx = data.table(mean.ctx, roi.ctx )
}
names(mean.ctx)[-1] <- roi

idata = subset(idata,select=names(idata)[!names(idata) %in% roi])
idata = subset(idata, select = names(idata)[!names(idata) %in% paste(roi,"1",sep='.')])
idata = merge(idata,mean.ctx,by='eid',sort=F)

##adjustment - Base model
M.H.eid = idata$eid[idata$Sex =="Male" & idata$HbA1c.Level3 =="High"]#
M.N.eid = idata$eid[idata$Sex =="Male" & idata$HbA1c.Level3 =="Normal"]
F.H.eid = idata$eid[idata$Sex =="Female" & idata$HbA1c.Level3 =="High"]
F.N.eid = idata$eid[idata$Sex =="Female" & idata$HbA1c.Level3 =="Normal"]
repi=1
select.eid = c(M.N.eid,
               F.N.eid,
               M.H.eid,
               F.H.eid)

rm(roii,fit_BaseF,fit_BaseM)
res_Base_allroi = NULL
for(i in 1:34){
  roii=roi[i]
  analdat = subset(idata,eid%in% select.eid,
                   select=c('eid',roii,'Time.c','Sex','Age.c',"MRI_Site","HbA1c.Level3")) %>%  dplyr::rename('y'=roii)
  tab = table(analdat$Sex,analdat$HbA1c.Level3)
  base::prop.table(tab)
  analdat = na.omit(analdat)
  
  idata.F = subset(analdat,Sex=="Female")
  idata.M = subset(analdat,Sex=="Male")
  
  modF_BaseMod <- gam(y ~ s(Time.c) + s(Age.c) + MRI_Site, 
                      data = idata.F)
  modM_BaseMod <- gam(y ~ s(Time.c) + s(Age.c) + MRI_Site, 
                      data = idata.M)
  
  idata.F = data.table(subset(idata.F,select="eid"),
                       y.resid=residuals(modF_BaseMod) + coef(modF_BaseMod)[['(Intercept)']])
  idata.M = data.table(subset(idata.M,select="eid"),
                       y.resid=residuals(modM_BaseMod) + coef(modM_BaseMod)[['(Intercept)']])
  
  idata.MF = rbind(idata.M,idata.F) 
  
  merge(idata,idata.MF,all.x=T) %>% 
    ggplot(aes(x=Time.c,y=y.resid,color=Sex,fill=Sex)) +
    geom_smooth() + 
    scale_color_manual(values=c("Female"="darkred","Male"="darkblue")) +
    theme(legend.position = "top")
  
  ## piecewise-linear regression
  ## Base
  resid.data = merge(idata,idata.MF,all.x=T)
  resid.data$y.resid = scale(resid.data$y.resid)
  fit_BaseM = lm(y.resid ~ rntHbA1c2 + x2.star, data=resid.data, subset=Sex=="Male")
  fit_BaseF = lm(y.resid ~ rntHbA1c2 + x2.star, data=resid.data, subset=Sex=="Female")
  fit_BaseI = lm(y.resid ~ Sex*(rntHbA1c2 + x2.star), data=resid.data)
  fit_BaseC = lm(y.resid ~ Sex + rntHbA1c2 + x2.star, data=resid.data)
  res_Base = data.table(get_results(fitC=fit_BaseC,fitM=fit_BaseM,fitF=fit_BaseF,fitI=fit_BaseI),
                        adjustment="base")
  res_Base$roi = roii
  res_Base_allroi = rbind(res_Base_allroi,res_Base)
}

res_High = subset(res_Base_allroi,Term %in% c("beta.H","beta.IH"))
res_High.C = subset(res_High,Sex=="Sex-Combined")
res_High.C = res_High.C[order(res_High.C$beta),]
roi.level = res_High.C$roi
res_High$roi = factor(res_High$roi,levels=(roi.level))
res_High = res_High[order(res_High$roi),]

res_High.F = subset(res_High,Sex=="Female")
res_High.M = subset(res_High,Sex=="Male")
res_High.I = subset(res_High,Sex=="interaction")

res_Base_Horizontal = merge(res_High.C,res_High.M,by='roi',sort=F)
res_Base_Horizontal = merge(res_Base_Horizontal,res_High.F,by='roi',sort=F)
res_Base_Horizontal$roi = factor(res_Base_Horizontal$roi,levels=roi.level)
res_Base_Horizontal = res_Base_Horizontal[order(res_Base_Horizontal$roi),]

var(res_High.C$beta)#0.004161753
var(res_High.M$beta)#4.637529e-05
var(res_High.F$beta)#1.494967e-04

cor.MF = cor(res_Base_Horizontal$beta.x,res_Base_Horizontal$beta.y)
res_High$Sex = factor(res_High$Sex, levels=group.levels)
res_High$roi = factor(res_High$roi,levels=rev(roi.level))
res_High %>% ggplot(aes(x=roi,y=beta,color=Sex,by=Sex)) +
  geom_point() + 
  geom_path() + 
  theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),legend.position = 'top') + 
  scale_color_manual(values=c("Sex-Combined" = "black","Male"="darkblue","Female"="darkred"))

# FDR-P: sex-stratified/interaction separately
res_High$FDR.P = NA
res_High$FDR.P[!is.na(res_High$Sex) & res_High$Sex=="Sex-Combined"] = p.adjust(res_High$P[!is.na(res_High$Sex) & res_High$Sex=="Sex-Combined"],method = 'fdr')
res_High$FDR.P[!is.na(res_High$Sex) & res_High$Sex!="Sex-Combined"] = p.adjust(res_High$P[!is.na(res_High$Sex) & res_High$Sex!="Sex-Combined"],method = 'fdr')
res_High$FDR.P[is.na(res_High$Sex)] = p.adjust(res_High$P[is.na(res_High$Sex)],method = 'fdr')
head(res_High)

## plot
library(viridis)
res = res_High
res$CI95L = res$beta - 1.96*res$SE
res$CI95U = res$beta + 1.96*res$SE
print(range(res$CI95L,res$CI95U));
siglevel=0.05
nonsignif_label=unique(res$signif[res$P >= siglevel])
res$Sex = factor(res$Sex,levels=rev(group.levels))
res$signif.FDR = ""
res$signif.FDR[res$FDR.P <0.05] <- "*"
res$signif = paste("FDR-P<",siglevel,sep='')
res$signif[res$FDR.P >= siglevel] <-  paste("FDR-P>=",siglevel,sep='')
signif_label="FDR-P<0.05"

# plot: create plot data frame ------------------------------------------------
plotd = res
plotd$roi <- factor(plotd$roi,levels=rev(roi.level))
plotd$p.label = as.character(format(plotd$P,scientific = T,digits = 3))
plotd$p.label = paste("p=",plotd$p.label,sep='')
plotd <- plotd[order(plotd$roi),]

# plot: create initial ggplot object ----------------------------------------------------
yintercept = 0
pos <- position_dodge(width=0)
ggp <- subset(plotd,!is.na(Sex)) %>% ggplot(aes(y=beta, x=roi, ymin=CI95L, ymax=CI95U,
                                                group=Sex, color=Sex, linetype=Sex,
                                                shape=signif,fill=signif)) +
  geom_point(size=6, stroke=1, position=pos) +
  geom_path(aes(group=Sex), position = pos) +
  geom_errorbar(position=pos, width=0.25,size=1) + 
  geom_hline(yintercept=yintercept, size=2) +
  scale_color_manual(values=c("Sex-Combined"="black","Male"=scales::alpha('darkblue',0.5),"Female"=scales::alpha("darkred",0.5))) +
  scale_shape_manual(values=c('FDR-P<0.05'=16,'FDR-P>=0.05'=21)) +
  scale_fill_manual(values=c('FDR-P<0.05'='black','FDR-P>=0.05'='white')) +
  #scale_y_continuous( limits = ylim, breaks = ylim.breaks) +
  coord_flip() +
  theme_bw() +
  #ggtitle(myggtitle) +
  labs(y="Effect estimates in SD-units (95% CI)") +
  theme(#panel.grid.major.y = element_line(colour = c("gray90", "gray90"), size=0.4),
    panel.grid.major.x = element_line(colour = "gray60", linetype="dashed", size=0.4),
    panel.grid.minor=element_blank(),
    panel.border = element_rect(color=NA),
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.8, "lines"),
    plot.title = element_text(size=12, hjust=0.5),
    strip.text.x = element_text(size=16, hjust=0.5, angle=0, face="bold"),
    strip.text.y = element_text(size=12, hjust=0, angle=30, face="bold"),#not being used
    strip.background = element_rect(fill=NA, color=NA),
    axis.line = element_line(color="black"),
    axis.ticks = element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_text(size=24),
    axis.title.x=element_text(size=20),
    axis.text.x=element_text(size=20),
    legend.position="top",
    legend.title=element_blank(),
    # legend.title=element_text(size=16),
    legend.text = element_text(size=16)
  )

# plot: annotation ----------------------------------------------------------------------
geom.text.size = 8
theme.size = (14/5) * geom.text.size

## sex-combined
ann_text <- subset(plotd, Sex=="Sex-Combined")
ann_text <- ann_text[order(ann_text$roi),]
p <- ggp + 
  geom_text(data = ann_text,label = ann_text$p.label,
            mapping = aes(x = roi, y = 0.38, label = label),
            colour="grey35",size=geom.text.size,  
            vjust = 0.5, hjust = 0.5,position = pos)
rm(ann_text)

## male
ann_text <- subset(plotd, Sex=="Male")
ann_text <- ann_text[order(ann_text$roi),]
p <- p + 
  geom_text(data = ann_text,label = ann_text$p.label,
            mapping = aes(x = roi, y = 0.58, label = label),
            colour="grey35",size=geom.text.size,  
            vjust = 0.5, hjust = 0.5,position = pos)
rm(ann_text)

## female
ann_text <- subset(plotd, Sex=="Female")
ann_text <- ann_text[order(ann_text$roi),]
p <- p + 
  geom_text(data = ann_text,label = ann_text$p.label,
            mapping = aes(x = roi, y = 0.78, label = label),
            colour="grey35",size=geom.text.size,  
            vjust = 0.5, hjust = 0.5,position = pos) 
rm(ann_text)

## interaction
ann_text <- subset(plotd, is.na(Sex))
ann_text <- ann_text[order(ann_text$roi),]
p <- p + 
  geom_text(data = ann_text,label = ann_text$p.label,
            mapping = aes(x = roi, y = 0.98, label = label),
            colour="grey35",size=geom.text.size,  
            vjust = 0.5, hjust = 0.5,position = pos)
rm(ann_text)

ggp1 = p + 
  scale_y_continuous(limits =c(-0.515,1.0), 
                     breaks = 0.1*c(-5,-4,-3,-2,-1,1,2,3)) +
  guides(fill = guide_legend(override.aes = list(shape = 21,size=4)),
         color = guide_legend(override.aes = list(size=4,linetype=1))) + 
  theme(legend.key.size = unit(1.5, "cm"))

# plot: create png file -----------------------------------------------------------------
png(paste('Plots/HbA1cLevels_roi_base_',Sys.Date(),'.png',sep=''),
    width=(12*1.75),height=(12*1.5),units='in',res=300)
print(ggp1+theme(strip.text.y.left = element_text(angle = 0,size=24)))
dev.off()

# table: write out a table ----------------------------------------------------
profiles_to_write = subset(res_Base_Horizontal,select=c(roi,beta.x,beta.y,beta)) %>% 
  dplyr::rename(beta.combined=beta.x,beta.male=beta.y,beta.female=beta)
(cor(subset(profiles_to_write,select=-roi)))#0.7283397 (wo Age2) #r=0.73069792 (wi Age2)
write_tsv(profiles_to_write,
          paste("~/OneDrive - SickKids/ukbb_insulin_resistance/Results_Tables/res_Base_HighHbA1c_",Sys.Date(),".tsv",sep=""))

dont.run <- function(){
  # test ----------------------------------------------------------------------
  roii = 'superiortemporal'
  analdat = subset(idata,eid%in% select.eid,select=c('eid',roii,'Time.c','Sex','Age.c',"MRI_Site","HbA1c.Level3")) %>%  
    dplyr::rename('y'=roii)
  analdat = subset(idata,eid%in% select.eid,
                   select=c('eid',roii,'Time.c','Sex','Age.c',"MRI_Site","HbA1c.Level3")) %>%  dplyr::rename('y'=roii)
  tab = table(analdat$Sex,analdat$HbA1c.Level3)
  base::prop.table(tab)
  analdat = na.omit(analdat)
  
  idata.F = subset(analdat,Sex=="Female")
  idata.M = subset(analdat,Sex=="Male")
  
  modF_BaseMod <- gam(y ~ s(Time.c) + s(Age.c) + MRI_Site, 
                      data = idata.F)
  modM_BaseMod <- gam(y ~ s(Time.c) + s(Age.c) + MRI_Site, 
                      data = idata.M)
  
  idata.F = data.table(subset(idata.F,select="eid"),
                       y.resid=residuals(modF_BaseMod) + coef(modF_BaseMod)[['(Intercept)']])
  idata.M = data.table(subset(idata.M,select="eid"),
                       y.resid=residuals(modM_BaseMod) + coef(modM_BaseMod)[['(Intercept)']])
  
  idata.MF = rbind(idata.M,idata.F) 
  
  merge(idata,idata.MF,all.x=T) %>% 
    ggplot(aes(x=Time.c,y=y.resid,color=Sex,fill=Sex)) +
    geom_smooth() + 
    scale_color_manual(values=c("Female"="darkred","Male"="darkblue")) +
    theme(legend.position = "top")
  
  ## piecewise-linear regression
  ## Base
  resid.data = merge(idata,idata.MF,all.x=T)
  resid.data$y.resid = scale(resid.data$y.resid)
  fit_BaseM = lm(y.resid ~ rntHbA1c2 + x2.star, data=resid.data, subset=Sex=="Male")
  fit_BaseF = lm(y.resid ~ rntHbA1c2 + x2.star, data=resid.data, subset=Sex=="Female")
  fit_BaseI = lm(y.resid ~ Sex*(rntHbA1c2 + x2.star), data=resid.data)
  fit_BaseC = lm(y.resid ~ Sex + rntHbA1c2 + x2.star, data=resid.data)
  
  ylim=c(-1,0.25)
  p1 = resid.data %>% ggplot(aes(x=rntHbA1c,y=y.resid,color=Sex)) + geom_smooth(method="gam") + 
    geom_vline(xintercept = c(cutoffM,cutoff,cutoffF), color=c("darkblue","black","darkred"), linetype="dotted")+
    ylab(NULL)  + theme(legend.position=c(0.05,0.05),legend.justification = c('left','bottom')) + 
    xlim(range(resid.data$rntHbA1c,resid.data$rntHbA1c2))
  p2 = resid.data %>% ggplot(aes(x=rntHbA1c2,y=y.resid,color=Sex)) + geom_smooth(method="gam") + 
    geom_vline(xintercept = c(cutoffM,cutoff,cutoffF), color=c("darkblue","black","darkred"), linetype="dotted")+
    ylab(NULL) + theme(legend.position='none') + xlim(range(resid.data$rntHbA1c,resid.data$rntHbA1c2))
  p3 = resid.data %>% ggplot(aes(x=rntHbA1c2,y=y.resid)) + geom_smooth(method="gam") + 
    geom_vline(xintercept = c(cutoffM,cutoff,cutoffF), color=c("darkblue","black","darkred"), linetype="dotted") +
    ylab(NULL) + xlim(range(resid.data$rntHbA1c,resid.data$rntHbA1c2))
  design.layout='
  12
  33
  '
  p.patched = (p1+(p2+ylab(NULL)))+p3 + plot_layout(design=design.layout) + coord_cartesian(ylim=ylim)
  p.patched + plot_annotation(
    title = paste0('Adjusted mean cortical thickness vs. HbA1c'),
    subtitle = roii,
    caption = 'compare the sex-differential vs. common breakpoints',
    tag_levels = 'A'
  )
}