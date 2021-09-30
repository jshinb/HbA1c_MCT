#---------------------------------------------------------------------------------------#
# Created: Sep 30, 2021
# 
# This script will run VH in the sex-combined sample.
#
#---------------------------------------------------------------------------------------#
options(stringsAsFactors = F)
x = c('stringr','tidyverse','dplyr',
      'here','data.table','psych',#table
      'GenABEL','mgcv',#rntransform
      'tableone','arsenal',#table
      'ggplot2','patchwork','pheatmap','viridis')#visualization
lapply(x,require,character.only=T);rm(x)

options(stringsAsFactors = F)
require(readr)
require(data.table)
require(tidyverse)
require(dplyr)

perm.test.gene.specific=FALSE
setwd("~/OneDrive - SickKids/4752955")
#rm(list=ls())
#-------------------------------------------------------------------------------#
# 
# Input parameters to set (n=6)
# 1. hemisphereParameter: 'lh' (left hemisphere) or 'rh' (right hemisphere), 
# 2. ReferenceGeneListFile: Path to the file contatining Gene Symbols, correlation 
#                           between BrainSpan and Allen expression levels and Zeisel 
#                           Cell types for the genes 
#    remaining after 2-step QC. 
# 3. wd: working directory where the downloaded files are loacated.
# 4. nreps: number of re-sampling replications for association test
# 5. fig_cex.main: for graphical display: scale of the title for each panel showing 
#                  the empirical density of cell-type-specific expression-phenotype  
#                  correlation coefficients.
# 6. fig_ylim: NULL or c(min,max), where min and max correspnd to minimum and maximum
#              values of y-axis. If NULL, each panel uses its local minimum and maximum
#			   estimated density values for each cell type; consequently, the y-axis scale 
#              can vary by cell type. Otherwise, all the panel will use the same y-axis 
#              scale specified by c(min,max).
#
# **** Example: ****
# hemisphereParameter = "lh"
# UserProfileFile = 'Profile_males_corcoef_thickness_age.csv'
# wd = "~b/Downloads/4752955/" 
# nreps = 10000
# siglevel = 0.05
# fig_cex.main = 1  
# fig_ylim = NULL # or c(0,1.6) # set a global y-axis scale
#
# Notes
# UserPhenoFile must be clneaned and formatted before loading:
# * The file must be saved as a csv.
# * The file must contain a column vector contatining 34 cortical-region-specific values 
#   with rownames containing the 34 regions.
# * It is assumed that these values are based on QC-ed data set 
#   (without the MRI-QC-failed, outlier measurements).  
# 
# Outputs (n=2): two files containing the folloiwng plot and table are created in 
#                'output' directory 
# 1. 3-row-by-3-column Figure showing empirical distributions for expression-phenotype 
#    correlation coefficients and the correponding (1-'siglevel')*100% critical values 
#    under null hypothesis of no association cell-type phenotype profile obtained from 
#    re-sampling distributions with 'nreps' replications
# 2. Re-sampling based test results, where the test statistic is mean expression-phenotype
#    coefficients (cell-type-specific)
#-------------------------------------------------------------------------------#

#------------------------- USER INPUTS #1 - #6 --------------------------#
# Users can complete this part and source the file (see the **example** above)
# To source file, Open R (or R studio) and type the following command line:
# source('CreateCorrelationCoefficients.r', echo=TRUE)
#
hemisphereParameter = "lh"
wd = "~/OneDrive - SickKids/4752955"
nreps = 100000
siglevel = 0.05
fig_cex.main = 1.5
fig_ylim = c(0,4)
roi_colname = 'roi'
#------------------------------------------------------------------------#

## Load necessary libraries 
## 
if("stringr" %in% rownames(installed.packages()) == FALSE) {install.packages("stringr")}
if("scales" %in% rownames(installed.packages()) == FALSE) {install.packages("scales")}

library(ggplot2); library(dplyr); library(data.table)
library(R.utils); library(readr); library(stringr)
library(readxl); library(scales)

## Set working directory to where the files are downloaded
setwd(wd) 
dir.create('output')

## Load and format rownames (i.e., cortical region names) for the Allen Gene expression data
geneExpressionMatrix <- read.table("AllenHBA_DK_ExpressionMatrix.tsv")
regionMapping <- read.table("DKRegionStatistics.tsv")
regionMapping <- subset(regionMapping, Hemisphere == hemisphereParameter)
rownames(regionMapping) <- gsub("-", ".", rownames(regionMapping))

## Load gene symbols and other information for the reference panel 
refpanel <- read.csv('Reference_Consistent_Genes_ObtainedBy2StageFiltering.tsv',
                     sep="\t",stringsAsFactors = F)

## Match the column names between expression and cortical phenotype data
## Load and format profile data (from the user)
## 6 TAG signif for mean cortical thickness 'RESULTS_TAG_thickness_assoc_adj_age_sex.xlsx'
#------------------------------------------------------------------------------------#
# thinning profile
# thinning=T
UserProfileFile='~/OneDrive - SickKids/ukbb_insulin_resistance/Results_Tables/res_Base_HighHbA1c_2021-07-22.tsv'
yname="beta.combined"
ctx.pheno.profile = subset(fread(UserProfileFile,data.table = F),select=c("roi",yname));
which.r = 'r_pearson'#or 'r_spearman'
thinning=T
if(thinning){
  ctx.pheno.profile[,yname] = (-1)*ctx.pheno.profile[,yname]
}else{
  ctx.pheno.profile[,yname] = ctx.pheno.profile[,yname]
}
rownames(ctx.pheno.profile) <- ctx.pheno.profile$roi

expressionForGene <- geneExpressionMatrix[refpanel$GeneSymbol,rownames(regionMapping)]
expressionForGene <- t(expressionForGene)

expressionForGene.df <- data.frame(expressionForGene,stringsAsFactors = F)
names(expressionForGene.df) <- colnames(expressionForGene)
expressionForGene <- expressionForGene.df
rm(expressionForGene.df)

expressionForGene_tmp = data.frame(rownames(expressionForGene),expressionForGene)
names(expressionForGene_tmp)[1] = roi_colname
rownames(expressionForGene) = str_remove(rownames(expressionForGene),paste('ctx',hemisphereParameter,"",sep="."))
matching.ind <- rownames(ctx.pheno.profile) %in% rownames(expressionForGene)
print(sum(matching.ind))

# #some or all regions match between phenotype and gene-expression files
# ctx.pheno.profile <- ctx.pheno.profile[matching.ind,,drop=F] #nrow <= 34
# # all regions are in the phenotype
# 
# if(nrow(ctx.pheno.profile) < 34){
#   missing.ind <- !rownames(expressionForGene) %in% rownames(ctx.pheno.profile)
#   missing.regions <- rownames(expressionForGene)[missing.ind]
#   warning("Missing phenotype values exist for the following regions: ",
#           paste0(str_replace(missing.regions,paste0("ctx[.]",hemisphereParameter,"[.]"),""),collapse = ","))
#   missing.pheno.vals <- data.frame(x=rep(NA,sum(missing.ind)))
#   row.names(missing.pheno.vals) <- missing.regions
#   ## augment the original data and the missing values
#   ctx.pheno.profile <- rbind(ctx.pheno.profile,missing.pheno.vals)
# }
# 
# if(sum(rownames(expressionForGene) == rownames(ctx.pheno.profile)) < 34){
#   o <- match(rownames(expressionForGene),rownames(ctx.pheno.profile))
#   ctx.pheno.profile <- ctx.pheno.profile[o,,drop=F] 
#   if(sum(rownames(expressionForGene) == rownames(ctx.pheno.profile))){
#     cat("The rows are matched with respect to the cortical regions.\n")
#   }
# }else if( sum(rownames(expressionForGene) == rownames(ctx.pheno.profile)) == 34 ){
#   cat("The rows between Allen expression and cortical phenotype profile data are already matched with respect to the cortical regions - no need to do anything.\n")
# }

## Calculate correlation coefficients between gene-expression and cortical-phenotype profiles
## Extract gene expression levels for the 2511 genes in refpanel
tmp_expressionForGene <- merge(expressionForGene,ctx.pheno.profile,by=0)
rownames(tmp_expressionForGene) <- tmp_expressionForGene$Row.names
names(tmp_expressionForGene) <- str_replace(names(tmp_expressionForGene),"Row.names",roi_colname)
expressionForGene = tmp_expressionForGene[,names(expressionForGene)[names(expressionForGene)!=roi_colname]]
ctx.pheno.profile = tmp_expressionForGene[,names(ctx.pheno.profile)[names(ctx.pheno.profile)!=roi_colname],drop=F]
rm(tmp_expressionForGene);gc(reset=T)

getP <- function(mylist){
  p = mylist$p.value
  p
}
print(identical(rownames(ctx.pheno.profile),rownames(expressionForGene)))
tmp = data.frame(ctx.pheno.profile,expressionForGene)
cor_test_p <- apply(tmp[,-1],2,cor.test,y=tmp[,1],method = 'p',use="p")
cor_test_s <- apply(tmp[,-1],2,cor.test,y=tmp[,1],method = 's',use="p")

cor_ReferencePanel <- NULL
cor_ReferencePanel <- rbind(cor_ReferencePanel,cor(tmp[,1],tmp[,-1],method = 'p',use="p"))
cor_ReferencePanel <- rbind(cor_ReferencePanel,sapply(cor_test_p,getP))
cor_ReferencePanel <- rbind(cor_ReferencePanel,cor(tmp[,1],tmp[,-1],method = 's',use="p"))
cor_ReferencePanel <- rbind(cor_ReferencePanel,sapply(cor_test_s,getP))
cor_ReferencePanel <- data.frame(colnames(cor_ReferencePanel),t(cor_ReferencePanel))  
names(cor_ReferencePanel) <- c('GeneSymbol','r_pearson','p_pearson','r_spearman','p_spearman')
cor_ReferencePanel$p_pearson.FDR <- p.adjust(cor_ReferencePanel$p_pearson,method = 'fdr')
cor_ReferencePanel$p_spearman.FDR <- p.adjust(cor_ReferencePanel$p_spearman,method = 'fdr')

cor_ReferencePanel$GeneSymbol = str_replace(cor_ReferencePanel$GeneSymbol,"HLA.DRA","HLA-DRA")
cor_ReferencePanel$GeneSymbol = str_replace(cor_ReferencePanel$GeneSymbol,"LY86.AS1","LY86-AS1")
cor_ReferencePanel$GeneSymbol = str_replace(cor_ReferencePanel$GeneSymbol,"MCM3AP.AS1","MCM3AP-AS1")
print(setdiff(cor_ReferencePanel$GeneSymbol,refpanel$GeneSymbol))

cor_ReferencePanel.temp <- merge(refpanel,cor_ReferencePanel,sort=F)
cor_ReferencePanel <- cor_ReferencePanel.temp
rm(cor_ReferencePanel.temp)

plotting=T
UserProfileFile = unlist(str_split(UserProfileFile,"/"))
UserProfileFile = UserProfileFile[length(UserProfileFile)]
file.ending=unlist(str_split(UserProfileFile,"[.]"))[length(unlist(str_split(UserProfileFile,"[.]")))]
for(metaboi in 1:ncol(ctx.pheno.profile)){
  profile_name = names(ctx.pheno.profile)[metaboi]
  tab_filename = paste('output/VH_Table_with',nreps,'replications',profile_name,
                       UserProfileFile,sep="_")
  tab_filename = str_replace(tab_filename,file.ending,"tsv")
  
  CorPlotFile = str_replace(tab_filename,'Table','Plot')
  CorPlotFile = str_replace(CorPlotFile,'tsv','png')
  if(plotting){
    png(CorPlotFile,width=10,height=8.5,pointsize=17,bg="white",units="in",res=250)
    par(mfrow=c(3,3),mar=c( 2.6, 3.6, 2.6, 1.1),oma=c(2,0,0,0),lwd=3)
  }
  celltypes = sort(unique(cor_ReferencePanel$CellType)[!is.na(unique(cor_ReferencePanel$CellType))])
  celltypes = celltypes[c(2,9,5,1,6,8,4,7,3)]
  r_resampling <- vector(mode = "list",length=length(celltypes))
  r.T.celli <- r.T.celli.p <- NULL
  rm(i)
  for(celli in 1:9){
    celltype = celltypes[celli]
    gene_ind = !is.na(cor_ReferencePanel$CellType) & cor_ReferencePanel$CellType==celltype
    celltype_genes = cor_ReferencePanel$GeneSymbol[gene_ind] 
    
    # print the range of donor to median correlation
    range.r = round(range(geneExpressionMatrix[celltype_genes,"Average.donor.correlation.to.median"]),2)
    cat(paste(celltype," - Average correlation between donors and median expression: (",range.r[1],",",range.r[2],")\n",
              sep=""))
    size.i=sum(gene_ind)
    
    ## perform re-sampling test (2-sided)
    r_resampling[[celli]] <- matrix(NA,nrow=size.i,ncol=nreps)
    x=1:length(cor_ReferencePanel$GeneSymbol)
    for(j in 1:nreps){
      row.ind = sample(x,size = size.i,replace = F)
      r_resampling[[celli]][,j] <- cor_ReferencePanel[row.ind,which.r]
    }
    names(r_resampling)[celli] <- celltype
    
    ## test statistics (T): mean correlation coefficient
    avgr = mean(cor_ReferencePanel[cor_ReferencePanel$GeneSymbol %in% celltype_genes,which.r])
    avgr_reps = apply(r_resampling[[celli]],2,mean)
    mean_under_null = mean(avgr_reps)
    avgr_reps = c(avgr,avgr_reps)# augment the observed statistic
    ## adjusted for the 'overall' mean
    avgr_reps = avgr_reps-mean_under_null
    
    p = sum(abs(avgr_reps) >= abs(avgr_reps[1]))/length(avgr_reps)
    r.T.celli <- c(r.T.celli,avgr_reps[1])
    r.T.celli.p <- c(r.T.celli.p, p)
    
    ## plotting empirical density of expression-phenotype correlation coefficients for each cell type
    if(plotting){
      ntest=1
      unAdjCI = quantile((avgr_reps),probs = c((siglevel/2)/ntest,(1-(siglevel/2)/ntest)))
      res_r=cor_ReferencePanel[gene_ind,which.r]
      d_res=density(res_r)
      xlim=range(density(cor_ReferencePanel[,which.r])$x)
      plot(d_res,main=celltype,las=1,cex.main=fig_cex.main,
           xlab="",yaxs="i",ylab="",xlim=xlim,ylim=fig_ylim)
      abline(v=(mean(res_r)-mean_under_null),col="black",lwd=3,lty=2) ## fixed
      myCI = unAdjCI#-mean_under_null
      if(is.null(fig_ylim)){
        polygon(c(myCI,rev(myCI)),c(c(0,0),c(max(d_res$y),max(d_res$y))),
                col=scales::alpha("grey",0.5),border=NA)
      }else{#!is.null(ylim)
        polygon(c(myCI,rev(myCI)),c(c(0,0),c(max(fig_ylim),max(fig_ylim))),
                col=scales::alpha("grey",0.5),border=NA)
      }#else
      # for(repi in 1:nreps){
      #   lines(density(r_resampling[[celli]][,repi]),col=scales::alpha('lightblue',0.5))
      # }
      # lines(d_res,col="black")
    }#if plotting
  }#celli in 1:9 (celltype)
  if(plotting){
    title(xlab="Correlation coefficient (r)",outer=T,line=0,cex.lab=1.25)
    title(ylab="Density",outer=T,line=-1,cex=0.9,cex.lab=1.25)
    dev.off()
  }
  ## Construct the result table
  res_tab <- data.frame(Cell_Type=celltypes,
                        nGenes = sapply(r_resampling,nrow),
                        avgr = round(r.T.celli,3),
                        p = r.T.celli.p,
                        fdr.p = p.adjust(r.T.celli.p,method="fdr"), 
                        stringsAsFactors=F)
  print(res_tab)
  
  ## write out the result table to 'res_tab,tab_filename'
  write_tsv(res_tab,tab_filename)
  write_tsv(cor_ReferencePanel,paste("cor_ReferencePanel",profile_name,UserProfileFile,sep="_"))
}# for metaboi

cor_ReferencePanel$CellType <- factor(cor_ReferencePanel$CellType,levels=res_tab$Cell_Type)
if(perm.test.gene.specific){
  tmp = select(ctx.pheno.profile,beta)
  tmp.cor = cor(tmp,expressionForGene)
  for(i in 1:10000){
    tmp$beta = sample(tmp$beta)
    tmp.cor = rbind(tmp.cor,cor(tmp,expressionForGene) )
  }
  head(tmp.cor[,1:2])
  get_p_perm <- function(x){
    p_perm = sum(abs(x)>=abs(x[1]))/length(x)
    p_perm
  }
  p_perm_profile = apply(tmp.cor,2,get_p_perm)
  identical(names(p_perm_profile),cor_ReferencePanel$GeneSymbol)
  cor_ReferencePanel$p_perm = p_perm_profile
  head(cor_ReferencePanel)
  profile.name=yname
  tab = table(cor_ReferencePanel$CellType,sign(cor_ReferencePanel[,which.r]),useNA='a')
  tab_0.05 =  table(cor_ReferencePanel$CellType[cor_ReferencePanel$p_perm<0.05],sign(cor_ReferencePanel[cor_ReferencePanel$p_perm<0.05,which.r]),useNA='a')
}
write_tsv(cor_ReferencePanel,paste("cor_ReferencePanel",profile_name,UserProfileFile,sep="_"))

# S1.Pyramidal.male = subset(cor_ReferencePanel,CellType=="S1.Pyramidal" & r_pearson> 0.109)
# Ependymal.male = subset(cor_ReferencePanel,CellType=="Ependymal" & r_pearson< -0.060)
# S1.Pyramidal.male = S1.Pyramidal.male[order(S1.Pyramidal.male$r_pearson),]
# Ependymal.male = Ependymal.male[order(Ependymal.male$r_pearson),]
# write_tsv(S1.Pyramidal.male,'S1.Pyramidal.male.tsv')
# write_tsv(Ependymal.male,'Ependymal.male.tsv')
