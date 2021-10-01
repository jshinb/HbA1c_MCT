#------------------------------------------------------------------------------#
# Created: Oct 1, 2021
# 
# This script will show the post-GWAS analyses
# 1. LocusZoom - regional association plots 
# 2. LDSC - Genetic Correlation
# 3. PRSice - polygenic risk scores of BMI and AD
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

# 7. Run GWAS with PLINK
## ----locuszoom, eval=F-------------------------------------------------------
## resdir='~/OneDrive - SickKids/GWAS_results/ukbb_th_crp_hba1c/adjBase_top2_CT_20210616'
## setwd(resdir)
## fs='cleaned_HighHbA1c_adjBase_all_PC1.glm.linear.txt.gz
## cleaned_HighHbA1c_adjBase_male_PC1.glm.linear.txt.gz
## cleaned_HighHbA1c_adjBase_female_PC1.glm.linear.txt.gz'
## fs =unlist(str_split(fs,"\n"))
## 
## 
## get_gwas_sub = function(x,chr,b0,b1){
##   dd = fread(x)
##   dd.sub = subset(dd,CHROM==chr&(POS>=b0 & POS<b1))
##   write_tsv(dd.sub,str_replace(x,"cleaned","LocusZoom"))
## }
## 
## chr=6
## b0=32578127-10^6
## b1=32578127+10^6
## lapply(fs,get_gwas_sub,chr=chr,b0=b0,b1=b1)
## c((32578127-250000),(32578127+250000))/10^6


## ----rg_PC1------------------------------------------------------------------
#https://jinghuazhao.github.io/Omics-analysis/BMI/


## ----PRS_AD_BMI--------------------------------------------------------------
# https://librarysearch.library.utoronto.ca/discovery/fulldisplay?docid=cdi_swepub_primary_oai_prod_swepub_kib_ki_se_232691445&context=PC&vid=01UTORONTO_INST:UTORONTO&lang=en&search_scope=UTL_AND_CI&adaptor=Primo%20Central&tab=Everything&query=any,contains,A%20principal%20component%20approach%20to%20improve%20association%20testing%20with%20polygenic%20risk%20scores&offset=0
