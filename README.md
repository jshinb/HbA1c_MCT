# HbA1c_MCT
Examine association between glycated hemoglobin and cortical thickness in adults

- Data cleaning and manipulation
  1. List of variables downloaded:
  2. List of variables to be created:
  
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
      2. The distribution of PC1 is skewed.
      3. Script: obtain PC1 in males and females separately, look at the distributions, run GWAS
      4. We further restricted to SNPs with minor allele frequency (MAF) > 0.1% and HWE p-value > 1e-10 in the 337,199 QC positive individuals, an INFO score > 0.8 (directly from UK Biobank), leaving 10.8 million SNPs for analysis. [http://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas]
  7. Regional association plots (LocusZoom)
  8. Genomic landscape
  9. Look up genes associated with the top locus in the brain cortex (GTEx database)
  10. Identify genes that are co-expressed with each of the eQTL-genes in brain cortex during adulthood (BrainSpan and harmonized database [by Nadine])
  11. Over-representation tests for the co-expressed genes of each eQTL-gene.
  12. Genetic correlations: AD and BMI
  13. Associations between PC1 and AD-PRS/BMI-PRS
