# 1. GCTA - input file
cd "/Users/jeanshin/OneDrive - SickKids/GWAS_results/ukbb_th_crp_hba1c/HighHbA1c_top2CT_adjBase_20211001"
f=cleaned_HighHbA1c_adjBase_all_PC1_with_a1freq.glm.linear.txt.gz
echo "SNP A1 A2 freq b se p N" > HighHbA1c_adjBase_all_PC1_for_gcta.txt
gunzip -c $f | awk '{if($3 != "ID"){print $3, $6, $4, $7, $9, $10, $11, $8}}' >> HighHbA1c_adjBase_all_PC1_for_gcta.txt
gzip HighHbA1c_adjBase_all_PC1_for_gcta.txt
scp HighHbA1c_adjBase_all_PC1_for_gcta.txt.gz jshin@data.ccm.sickkids.ca:/hpf/projects/zdenka30/GCTA_JS/HighHbA1c_top2CT_adjBase_2021-10-01/

# 2. Cluster
# 'HighHbA1c_PC1_top2CT_HbA1c_cojo_top1_2021-10-05.sh'