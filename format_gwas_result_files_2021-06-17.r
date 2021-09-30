# module load gcc/8.3.0 r/3.5.3;R
# rfile='/home/t/tpaus/jshinb/ukbb_gwas_scripts/format_gwas_result_files_2021-06-17.r'; source(rfile)
options(stringsAsFactors = F)
library(data.table)

##
# arguments
args=(commandArgs(TRUE))
if(length(args) > 0){
  print(args)
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  } #for(i in 1:length(args))
  mypars = unlist(strsplit(dirnamegroup,"-"))
  dirname <- mypars[1]
  group <- mypars[2]
  pheno <- mypars[3]
}else{#if(length(args) > 0)
  pheno='PC1.resid'
  group='all'
  dirname=paste('ukbb_gwas_',pheno,'.',group,sep="")
} #else

##
phenos=c('top2.MCT_adjBase','PC1')
adj='adjBase'
resdir=paste('/scratch/t/tpaus/jshinb/PLINK/HighHbA1c_top2CT_adjBase_2021-06-16/',adj,'/',sep="")
setwd(resdir)
for( group in c('all','male','female','SNPxSex') ){
#for( group in c('SNPxSex') ){
  cat(group,"\n",sep="")
  nmarker_all = NULL
  for( pheno in phenos ){
    cat(pheno,"\n",sep="")
    gwas=NULL
    nmarker = c()
    for(i in 1:22){
      cat(i,",")
      resf = paste('pheno_',adj,'_',group,'_',i,'.',pheno,'.glm.linear',sep='')
      gwasi = subset(fread(resf),!is.na(P))	
      nmarker = c(nmarker,nrow(gwasi))
      gwas = rbind(gwas,gwasi)
    }#for(i in 1:22)
    cat("\n")
    gwasfname=paste(resf,'.txt',sep='')
    gwasfname=gsub(paste('pheno_',adj,"_",group,'_',i,'.',sep=''),
                   paste('HighHbA1c_',adj,"_",group,'_',sep=''),gwasfname)
    write.table(gwas,gwasfname,sep='\t',quote=F,col.names=T,row.names=T)
    system(paste('gzip',gwasfname))
  }# for(outcome)
  nmarker_all = cbind(nmarker_all,nmarker)
  df.nmarker_all = data.frame(chr=1:length(nmarker_all),nMarker=nmarker_all)
  write.table(df.nmarker_all,paste(pheno,group,'nmarker_all.tsv',sep='_'),sep="\t",row.names=F, col.names=T, quote=F)
} #for( group )
