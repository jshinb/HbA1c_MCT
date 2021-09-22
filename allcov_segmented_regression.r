#---------------------------------------------------------------------------------------#
# Created: Sep 22, 2021
# 
# This script will fit segmented regression models to the data
# 
# segmented linear regression model with the adjusted MCT in males and females separately.
#
# 1. Should I do it in the sex-combined sample?
# 2. I could not reproduce the 'alladj' results - Why?
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

get_results = function(fitM,fitF,fitI){
  var.statM = sum(vcov(fitM)[-1,-1])
  slopeHM =sum(summary(fitM)$coef[c('rntHbA1c','x2.star'),"Estimate"])
  wald.testM  = (slopeHM/sqrt(var.statM))^2
  pM = pchisq(wald.testM,df = 1,lower.tail = F)
  
  var.statF = sum(vcov(fitF)[-1,-1])
  slopeHF =sum(summary(fitF)$coef[c('rntHbA1c','x2.star'),"Estimate"])
  wald.testF  = (slopeHF/sqrt(var.statF))^2
  pF = pchisq(wald.testF,df = 1,lower.tail = F) #
  
  iterm_name = c('SexFemale:rntHbA1c','SexFemale:x2.star')
  iterm_index = which(colnames(vcov(fitI)) %in% iterm_name)
  slopeI = sum(summary(fitI)$coef[iterm_name,"Estimate"])
  var.statI = sum(vcov(fitI)[iterm_index,iterm_index])
  wald.testI  = (slopeI/sqrt(var.statI))^2
  pI = pchisq(wald.testI,df = 1,lower.tail = F) #
  
  resM = data.table(Term = c("Intercept","beta.N","beta.diff.HN"),summary(fitM)$coef[,-3]) %>% rename(beta="Estimate", SE="Std. Error",P="Pr(>|t|)")
  resF = data.table(Term = c("Intercept","beta.N","beta.diff.HN"),summary(fitF)$coef[,-3])%>% rename(beta="Estimate",SE="Std. Error",P="Pr(>|t|)")
  resM = rbind(resM,data.table(Term="beta.H",beta=slopeHM,SE=sqrt(var.statM),P=pM))
  resF = rbind(resF,data.table(Term="beta.H",beta=slopeHF,SE=sqrt(var.statF),P=pF))
  resM = data.table(resM,Sex="Male")
  resF = data.table(resF,Sex="Female")
  resIH = data.table(Term="beta.IH",beta=slopeI,SE=sqrt(var.statI),P=pI,Sex="interaction")
  resIN = data.table(Term="beta.IN",
                     beta=summary(fitI)$coef["SexFemale:rntHbA1c","Estimate"],
                     SE=summary(fitI)$coef["SexFemale:rntHbA1c","Std. Error"],
                     P=summary(fitI)$coef["SexFemale:rntHbA1c","Pr(>|t|)"],
                     Sex="interaction")
  res = rbind(resM,resF,resIN,resIH)
  res
}

fit_function = function(yname,data,adj){
  analdat = data
  analdat$y = analdat[[yname]]
  fitM = lm(y ~ rntHbA1c + x2.star, data=analdat, subset=Sex=="Male")
  fitF = lm(y ~ rntHbA1c + x2.star, data=analdat, subset=Sex=="Female")
  fitI = lm(y ~ Sex*(rntHbA1c + x2.star), data=analdat)
  
  res = data.table(get_results(fitF=fitF,fitM=fitM,fitI=fitI),adjustment=adj)
  res
} 

wd = '~/OneDrive - SickKids/ukbb_insulin_resistance'
load(file.path(wd,"idata_wi_MCT_adjAll.Rd"))

## ----MCT_HbA1c_test_in_HighNormal,fig.height=9, fig.width=8----------------------------------------------------------------------------------
#https://online.stat.psu.edu/stat501/lesson/8/8.8
setwd(wd)
cutoffM = 0.7714078
cutoffF = 1.425325
idata$x2 = as.numeric(idata$HbA1c.Level3=="High")
idata$x2.star = idata$x2
idata$x2.star[idata$Sex=="Male"] = (idata$rntHbA1c[idata$Sex=="Male"]-cutoffM)*idata$x2[idata$Sex=="Male"]
idata$x2.star[idata$Sex=="Female"] = (idata$rntHbA1c[idata$Sex=="Female"]-cutoffF)*idata$x2[idata$Sex=="Female"]

ynames = c("MCT_adjBase", "MCT_adjBaseEduLevel", "MCT_adjBaseBMI", "MCT_adjBaseSmoking", 
           "MCT_adjBaseSBP","MCT_adjall")
adjs = c("base", "edu", "bmi", "smoking", "sbp", "alladj")
data.frame(ynames,adjs)
res_MCT = NULL
for(i in 1:length(adjs)){
  res_MCT = rbind(res_MCT,fit_function(yname=ynames[i],data=idata,adj=adjs[i])) 
}
res = res_MCT

# generate columns for plotting ---------------------------------------------------------
siglevel=0.05
color.levels=c("Male_All", "Male_High", "Male_Normal",
               "Female_All", "Female_High", "Female_Normal", 
               "interaction_All","interaction_High", "interaction_Normal")
adj.levels = c('base','edu','bmi','smoking','sbp','alladj')
hba1c.level=c( "Normal","High")

res = res %>% mutate(CI95L = beta-1.96*SE, 
                     CI95U = beta+1.96*SE,
                     HbA1c.Level = NA,
                     signif = paste("P<",siglevel,sep=''))
res$HbA1c.Level[str_detect(res$Term,"beta.H")|str_detect(res$Term,"beta.IH")] <- "High"
res$HbA1c.Level[str_detect(res$Term,"beta.N")|str_detect(res$Term,"beta.IN")] <- "Normal"
res$signif[res$P >= siglevel] <-  paste("P>=",siglevel,sep='')

signif_label=unique(res$signif[res$P < siglevel])
nonsignif_label=unique(res$signif[res$P >= siglevel])

res = res %>% mutate(Sex_HbA1c.Level = factor(paste(Sex,HbA1c.Level,sep="_"),levels=color.levels),
                     adjustment = factor(adjustment,levels=rev(adj.levels)),
                     x = factor(adjustment, levels = adj.levels))

res = res %>% mutate(HbA1c.Level = factor(HbA1c.Level,levels=rev(hba1c.level)),
                         group = paste(Sex,adjustment,sep="_"),
                         p.label = paste("p=",as.character(format(P,scientific = T,digits = 3)),sep=''),
                         Sex_signif = paste(Sex,signif,sep="_"))
res = res[order(res$adjustment,res$HbA1c.Level),]

# plot ----------------------------------------------------------------------------------
library(viridis)
# color.vals = c("darkblue","#1F78B4","#A6CEE3",
#                "darkred","#E31A1C","#FB9A99",
#                "grey","black","darkgrey")
color.vals = c("darkblue","darkblue","darkblue",
               "darkred","darkred","darkred",
               "grey","black","darkgrey")
names(color.vals) = color.levels

pointcols <- c('black',viridis(9)[1:7])#c("black","lightseagreen","grey","pink",'lightblue','lightblue')
shape_values = c(21,22,23,24,25,21,22,23)
names(shape_values) <- names(pointcols) <- adj.levels

pos <- position_dodge(width=1)
myggtitle="HbA1c"

dd <- subset(res,Sex!='interaction'& HbA1c.Level %in% c('Normal',"High"))
dd <- dd[order(dd$adjustment),]
dd <- dd[order(dd$HbA1c.Level,dd$adjustment),]

yintercept=0
print(range(dd$CI95L,dd$CI95U))
ylim = c(-0.04,0.0675); ylim.breaks = 0.01*seq(from=-5,to=3,by=1)
pos <- position_dodge(width=1)

ggp <- ggplot(dd,aes(y=beta, x=adjustment, ymin=CI95L, ymax=CI95U,
                     shape=x, group=Sex, linetype=Sex,
                     color=Sex,
                     fill=Sex_signif,
                     label=p.label)) +
  geom_point(size=7, stroke=1.5,position=pos) +  
  geom_errorbar(position=pos, width=0.25,size=0.75) + 
  geom_hline(yintercept=yintercept, size=1) +
  geom_path(aes(group=Sex),position=pos) +
  facet_grid(rows=vars(HbA1c.Level),scales = "free_y", switch='y') +
  scale_color_manual(name=NULL,values=c("Male"="darkblue","Female"="darkred")) +
  scale_shape_manual(values=shape_values) +
  scale_fill_manual(name=NULL,values=c('Male_P<0.05'="darkblue",'Male_P>=0.05'="white",
                                       'Female_P<0.05'="darkred",'Female_P>=0.05'="white")) +
  guides(color = guide_legend(override.aes = list(shpae=23,linetype=0,size=3)),
         fill = guide_legend(override.aes = list(shpae=21,size=3))) +
  scale_y_continuous( limits = ylim, breaks = ylim.breaks) +
  coord_flip() +
  theme_bw() +
  #ggtitle(myggtitle) +
  labs(y="Effect estimates in SD-units (95% CI)") +
  theme(#panel.grid.major.y = element_line(colour = c("gray90", "white"), size=6),
    panel.grid.major.x = element_line(colour = "gray60", linetype="dashed", size=0.4),
    panel.grid.major.y = element_blank(),
    panel.grid.minor=element_blank(),
    panel.border = element_rect(color=NA),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(0.8, "lines"),
    plot.title = element_text(size=12, hjust=0.5),
    strip.text.y = element_text(size=20, hjust=0.5, angle=0, face="bold"),
    strip.background = element_rect(fill=NA, color=NA),
    axis.line = element_line(color="black"),
    axis.ticks = element_blank(),
    axis.title.y=element_blank(),
    axis.title.x=element_text(size=16),
    axis.text.x=element_text(size=16,angle =0, vjust = 0, hjust=0.5),
    axis.text.y=element_text(size=12),#element_blank(),
    legend.position="top",
    # legend.title=element_blank()
    legend.title=element_text(size=16)
  )

## annotation
geom.text.size = 6
theme.size = (14/5) * geom.text.size
## 
ann_text.m <- subset(res, Sex=="Male" & !is.na(HbA1c.Level))
ann_text.f <- subset(res, Sex=="Female" & !is.na(HbA1c.Level))
ann_text.i <- subset(res, Sex=="interaction" & !is.na(HbA1c.Level))

## male
ann_text = ann_text.m %>% mutate(HbA1c.Level=factor(HbA1c.Level,levels=rev(hba1c.level)))
ann_text <- ann_text[order(ann_text$HbA1c.Level,ann_text$adjustment),]
p <- ggp + 
  geom_text(data = ann_text,label = ann_text$p.label,
            mapping = aes(x = adjustment, y = 0.0375, label = label),
            colour="grey50",size=geom.text.size,  
            vjust = 0.5, hjust = 0.5,position = pos)
rm(ann_text)

#female
ann_text = ann_text.f %>% mutate(HbA1c.Level=factor(HbA1c.Level,levels=rev(hba1c.level)))
ann_text <- ann_text[order(ann_text$HbA1c.Level,ann_text$adjustment),]
p <- p + 
  geom_text(data = ann_text,label = ann_text$p.label,
            mapping = aes(x = adjustment, y = 0.0520, label = label),
            colour="grey50",size=geom.text.size,  
            vjust = 0.5, hjust = 0.5,position = pos) 
rm(ann_text)

# interaction
ann_text <- ann_text.i %>% mutate(HbA1c.Level=factor(HbA1c.Level,levels=rev(hba1c.level)))
ann_text <- ann_text[order(ann_text$HbA1c.Level,ann_text$adjustment),]
ann_text$Sex_signif = str_replace(ann_text$Sex_signif,"interaction","Male")

p <- p + 
  geom_text(data = ann_text,label = ann_text$p.label,
            mapping = aes(x = adjustment, y = 0.0675, label = label),
            colour="grey50",size=geom.text.size,  
            vjust = 0.5, hjust = 0.5,position = pos) 

pdf(paste('Plots/tmp_MCT_rntHbA1c_in_different_sex_notAll_',Sys.Date(),'.pdf',sep=''),
    width=16,height=9,onefile = T)
print(p)
dev.off()
