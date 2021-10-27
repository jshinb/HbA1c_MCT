#------------------------------------------------------------------------------#
# Created: Oct 1, 2021
# 
# This script will 
# 1. create manhattan and QQ plots for the GWAS of PC1 of top2CT and HbA1c 
# in high HbA1c group
# 
#------------------------------------------------------------------------------#

options(stringsAsFactors = F)
x = c('stringr','tidyverse','dplyr',
      'here','data.table','psych',#table
      'GenABEL','mgcv',#rntransform
      'tableone','arsenal',#table
      'ggplot2','patchwork','pheatmap','viridis',#visualization
      'R.utils','qqman','gap','GWASTools')# manhattan and qqplot: do I use GWASTools?
lapply(x,require,character.only=T);rm(x)

source('~/OneDrive - SickKids/CHARGE_MetabolomicsGWAS/scripts/manhattan.plot-MacBook.R')
source('~/OneDrive - SickKids/CHARGE_MetabolomicsGWAS/scripts/annotateSNPRegions.r')

clean_PLINK_GWAS_file <- function(gwas_file,
                                  gwasres_dir,
                                  wd,
                                  pheno,group,adj)
{
  require(stringr)
  main = paste(pheno,'_',adj," (",group,")",sep="")
  
  setwd(gwasres_dir)
  # dd = fread(file.path(gwasres_dir,gwas_file))
  dd = fread(gwas_file)
  setwd(wd)
  if(group=="SNPxSex"){
    dd = subset(dd,TEST=="ADDxSex_Male0")
  }
  if(any(names(dd)=="V1")){
    dd = dplyr::select(dd,-V1)
  }
  if(any(str_detect(names(dd),"TEST"))){
    dd = dplyr::select(dd,-TEST)
  }
  if(any(str_detect(names(dd),"T_STAT"))){
    dd = dplyr::select(dd,-T_STAT)
  }
  names(dd) <- str_remove(names(dd),"#")
  dd$MarkerName = paste(dd$CHROM,dd$POS,sep=":")
  write_delim(dd,file.path(gwasres_dir,paste('cleaned_',resfile,sep="")),delim="\t")
}

create_manhattan_plot <- function(gwasres_dir,
                                  cleaned_gwas_file,
                                  pheno,group,adj,
                                  ymax=NULL,
                                  width=11,
                                  height=3,
                                  main,
                                  colvalues=c("#69b3a2", "#404080"),
                                  man_plot_dir)
{
  dd = fread(file.path(gwasres_dir,cleaned_gwas_file))
  manhattan.p = manhattan.plot(chr=dd$CHROM,
                               pos=dd$POS,
                               pvalue = dd$P,
                               direction=1,#sign(dd$BETA),
                               main=main,
                               col=colvalues,
                               sig.level=5e-8,
                               ymax=ymax)
  
  pfile=paste(man_plot_dir,"/mh.",paste("ukbb",main,sep="_"),"_",Sys.Date(),".png",sep="")
  png(pfile, width=width, height=height, units="in", res=300)
  print(manhattan.p)
  dev.off()
  ret = list(gwas.dd = dd, N=mean(dd$OBS_CT),pheno=pheno,group=group, 
             minP=min(dd$P), minpSNP = dd$ID[which.min(dd$P)])
  ret
}
#------------------------------------------------------------------------------
wd = '~/OneDrive - SickKids/ukbb_insulin_resistance'
manplot_dir = file.path(wd,'manhattan_plots')
setwd(wd)
dir.create(manplot_dir)

#------------------------------------------------------------------------------#
gwasres_dir='~/OneDrive - SickKids/GWAS_results/ukbb_th_crp_hba1c/HighHbA1c_top2CT_adjBase_20211001/'
setwd(gwasres_dir)
dir()

resfiles=c("HighHbA1c_adjBase_all_PC1.glm.linear.txt.gz")
groups = str_split(resfiles,"_",simplify = T)[,1]
phenos = str_remove(str_split(resfiles,"_",simplify = T)[,4],".glm.linear.txt.gz")
adj="age_sex_genoPC"

resfile.df = data.frame(resfile=resfiles,group=groups,pheno=phenos)
for(i in 1:nrow(resfile.df)){
  print(i)
  resfile=resfile.df$resfile[i]
  group=resfile.df$group[i]
  pheno=resfile.df$pheno[i]
  clean_PLINK_GWAS_file(gwas_file=resfile,group = group,pheno = pheno,gwasres_dir = gwasres_dir,wd = wd, adj=adj)
}

#
fs=c("cleaned_HighHbA1c_adjBase_all_PC1.glm.linear.txt.gz")
minp = c()
for(i in 1:length(fs)){
  print(i)
  dd = fread(file.path(gwasres_dir,fs[i]))
  minp = c(minp,min(dd$P,na.rm=T))
}

fs = data.frame(file=fs,minP = minp)
if(any(dd$P==0)){
  dd$P[dd$P==0] = 4.940656e-324
  if(minp==0) minp=4.940656e-324
  write_delim(dd,file.path(gwasres_dir,fs$file),delim="\t")
}

##
# manhattan and qqplots
##
ymax=max(10,max(-log10(minp))*1.15)
for(i in 1:nrow(resfile.df)){
  print(i)
  resfile=resfile.df$resfile[i]
  group=resfile.df$group[i]
  pheno=resfile.df$pheno[i]
  
  cat("---------------------------\n")
  cat(pheno,group,adj,resfile,sep='\n')
  cat("---------------------------\n")
  
  res_gwas = create_manhattan_plot(
    gwasres_dir=gwasres_dir,
    cleaned_gwas_file = paste('cleaned',resfile,sep='_'),
    pheno=pheno,
    group=group,
    adj="adjBase",
    ymax=ymax,
    main=paste(pheno,group,adj,sep="_"),
    width=11,height=3,
    man_plot_dir=manplot_dir)
  
  pval = na.omit(res_gwas$gwas.dd$P)
  if(str_count(resfile,pheno)>1){
    qqfile=str_replace(str_remove(resfile,pheno),'txt.gz',paste0("_",Sys.Date(),'_qqplot.png'))
    qqfile=str_replace(qqfile,"__","_")
  }else{
    qqfile=str_replace(resfile,'.txt.gz',paste0("_",Sys.Date(),'_qqplot.png'))
  }
  cat(qqfile,"\n")
  
  gc(reset=T)
  png(file.path(manplot_dir,qqfile),width=8,height=8,
      units="in",res=300)
  r = gcontrol2(pval)
  title(main = paste(pheno,'_',group,'_',adj,' (lambda=',round(r$lambda,4),")",sep=""))
  dev.off()
  
  print(r$lambda)
}