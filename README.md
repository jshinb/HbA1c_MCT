# HbA1c_MCT
Examine association between glycated hemoglogin and cortical thickness in adults

- Methods
  1. Sample characteristics
      1. Sample selection: individuals with HbA1c and MRI data (data download date: XX-XX-2020)
  2. Correlation structure among variables  
  3. Association between HbA1c and MCT - by segmented regression analysis
      1. Base model: MCT ~ HbA1c + Age + Sex + MRI site
      2. Covariates: education level (high/low), BMI, smoking, and/or systolic blood pressure
  4. Regional association between HbA1c and MCT in high-level of HbA1c
      1. It is already conditioned on HbA1c
      2. In a way, it is a stratified analysis
  5. Virtual histology
  6. GWAS of PC1 of HbA1c and MCT in high-HbA1c group
      1. Potential issue: HbA1c varies less...because we only took a subset based on cutoff values [do we want this?]
      2. The distribution of PC1 is almost truncated normally distributed because HbA1c is truncated (no transformation will give normal distribution).
      3. Script: obtain PC1 in males and females separately, look at the distributions, run GWAS
  7. Regional association plots (LocusZoom)
  8. Genomic landscape
  9. Look up genes associated with the top locus in the brain cortex (GTEx database)
