#!/bin/bash 
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=1:00:00
#SBATCH --job-name=adjBase
#SBATCH --output=adjBase_%j.txt
#SBATCH --mail-type=FAIL

chr=$1
resdir=$2
adj=$3

echo $chr

mkdir $resdir
mkdir ${resdir}/${adj}
cd ${resdir}/${adj}

mybfile=ukb43688_imp_chr${chr}_v3_filtered
keepfile=/home/t/tpaus/jshinb/ukbb_gwas_data/famfile_PC1_HighHbA1c_top2_CT_${adj}_all.fam
ls -l $keepfile

module load NiaEnv/2019b intel/2019u5 plink/.experimental-2.00alpha
plink2 --bgen /project/t/tpaus/tpaus/UKBB/datasets/gene_data/ukb_imp_chr${chr}_v3.bgen --sample /project/t/tpaus/tpaus/UKBB/datasets/gene_data/gene_data_backup/backup_gene_data_sample/ukb43688_imp_chr${chr}_v3_s487314.sample --maf 0.01 --mac 100 --geno 0.05 --keep $keepfile --remove /home/t/tpaus/mbernard/ukbiobank/subjects_to_remove_list.txt --hwe 1e-15 --extract /home/t/tpaus/mbernard/ukbiobank/chrom${chr}_selected_markers.txt --make-bed --out ${mybfile}

echo "female analysis"
group=female
pfile=/gpfs/fs1/home/t/tpaus/jshinb/ukbb_gwas_data/pheno_PC1_HighHbA1c_top2_CT_${adj}_${group}.txt
ls -l $pfile
plink2 --bfile ${mybfile} --glm omit-ref --pheno ${pfile} --out pheno_${adj}_${group}_${chr}

echo "male analysis"
group=male
pfile=/gpfs/fs1/home/t/tpaus/jshinb/ukbb_gwas_data/pheno_PC1_HighHbA1c_top2_CT_${adj}_${group}.txt
ls -l $pfile
plink2 --bfile ${mybfile} --glm omit-ref --pheno ${pfile} --out pheno_${adj}_${group}_${chr}

echo "interaction analysis"     
group=all
pfile=/gpfs/fs1/home/t/tpaus/jshinb/ukbb_gwas_data/pheno_PC1_HighHbA1c_top2_CT_${adj}_${group}.txt
covfile=/gpfs/fs1/home/t/tpaus/jshinb/ukbb_gwas_data/cov_PC1_PC1_HighHbA1c_top2_CT_${adj}_${group}.txt
ls -l $pfile
ls -l $covfile
plink2 --bfile ${mybfile} --glm omit-ref --pheno ${pfile}  --out pheno_${adj}_${group}_${chr}
plink2 --bfile ${mybfile} --pheno ${pfile} --covar ${covfile} --glm interaction --parameters 1-7 --tests 1, 3 --out pheno_${adj}_SNPxSex_${chr}
