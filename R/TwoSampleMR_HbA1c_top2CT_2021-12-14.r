# HbA1c -> CT
library(TwoSampleMR)
# Methods to be applied:
mr.methods = mr_method_list()[c(11,10,6,13,3),'obj']
cat(mr.methods,sep="\n")
# mr_ivw_fe
# mr_ivw_mre
# mr_weighted_median
# mr_weighted_mode
# mr_egger_regression

# Indcate directory with EINGMA cortical gwas studies
enigma_dir = '/Users/jshin/OneDrive - SickKids/ukbb_insulin_resistance/ENIGMA_GWAS_results_top2CT/'
dir(enigma_dir)

# List available GWASs
ao <- available_outcomes()

# List available GWASs
ao <- available_outcomes()

# Identify potential exposure (HbA1c) GWAS studies
# traits="Glycated haemoglobin (HbA1c)
# Glycated haemoglobin
# Glycated hemoglobin levels
# HbA1C"
# traits = unlist(str_split(traits,"\n"))
# subset(ao, trait %in% traits & population == "European")

# Get instruments
exposure_dat <- extract_instruments("ukb-d-30750_irnt")#UKB
subset(ao, id=="ukb-d-30750_irnt")

# Get effects of instruments on outcome
outcome_files = c('ENIGMA3_mixed_se_wo_Mean_superiorfrontal_thickavg_20200522.txt.gz',#(p=0.79)
                  'ENIGMA3_mixed_se_wo_Mean_superiortemporal_thickavg_20200522.txt.gz',#(p=0.64)
                  'ENIGMA3_mixed_se_wTHICK_Mean_superiorfrontal_thickavg_20200522.txt.gz',#p=0.51
                  'ENIGMA3_mixed_se_wTHICK_Mean_superiortemporal_thickavg_20200522.txt.gz')#p=0.65

ofile=outcome_files[i]

get2SampleMR_results = function(exposure_dat,ofile,exposure,outcome){
  outcome_dat <- read_outcome_data(
    snps = exposure_dat$SNP,
    filename = file.path(enigma_dir,ofile),
    sep = " ",
    snp_col = "SNP",
    beta_col = "BETA1",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "FREQ1",
    pval_col = "P"
  )
  
  # Harmonise the exposure and outcome data
  dat <- harmonise_data(exposure_dat,outcome_dat)
  
  # Run 2 sample MR
  mr_pleiotropy_test(dat)#p=0.53, 0.67 (superiorforntal, wo and wi MeanTh); p=0.64, 0.66 (superiortemporal, wo and wi MeanTh)
  mr_results <- mr(dat, method_list=mr.methods)
  p = mr_scatter_plot(mr_results,dat)
  ret=list(p=p, mr_results=mr_results, exposure=exposure, outcome=outcome)
  ret
}

anal1 = get2SampleMR_results(exposure_dat=exposure_dat,
                             ofile=outcome_files[1],
                             exposure='HbA1c',
                             outcome='superiorfrontal_CT_wo_MeanTH')


anal2 = get2SampleMR_results(exposure_dat=exposure_dat,
                             ofile=outcome_files[2],
                             exposure = 'HbA1c',
                             outcome = 'superiortemporal_CT_wo_MeanTH')

anal3 = get2SampleMR_results(exposure_dat=exposure_dat,
                             ofile=outcome_files[3],
                             exposure='HbA1c',
                             outcome='superiorfrontal_CT_wi_MeanTH')


anal4 = get2SampleMR_results(exposure_dat=exposure_dat,
                             ofile=outcome_files[2],
                             exposure = 'HbA1c',
                             outcome = 'superiortemporal_CT_wi_MeanTH')

# results
library(gapminder)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(gt)
library(mmtable2)

write_2sampleMR_results = function(ResList,fname="~/Desktop/test.png"){
  mr_results_melt = melt(dplyr::select(ResList$mr_results,method,nsnp,b,se,pval)) 
  
  tab1 <- mr_results_melt %>%
    mmtable(cells = round(value,4)) + 
    header_top(variable) +
    header_left(method)
  
  png(fname,width=11,height=6,units="in",res=300)
  p <- ResList$p[[1]] + 
    ggtitle(paste(ResList$exposure,"->",ResList$outcome)) + 
    theme(legend.position = 'bottom')
  print(p)
  dev.off()

  tab1#800x300 -> copy
} 
write_2sampleMR_results(anal4)
