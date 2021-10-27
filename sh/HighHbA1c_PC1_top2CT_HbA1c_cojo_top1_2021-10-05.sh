#------------------------------------------------------------------#
# created on 2021-OCT-05
# this script is to run conditional-joint analyses for GWAS of PC1
# of top2CT and HbA1c in the high-HbA1c group.
#------------------------------------------------------------------#

#!/bin/bash
#$ -cwd
cd /hpf/projects/zdenka30/GCTA_JS/HighHbA1c_top2CT_adjBase_2021-10-01/genedata
module load gcta

gwasf=../HighHbA1c_adjBase_all_PC1_for_gcta.txt.gz
gcta64  --mbfile ukb43688_imp_v3_filtered.list --maf 0.01 --cojo-file $gwasf --cojo-cond topSNP.list --out test_top1
#qsub -N"cojo" -l mem=64g,vmem=64g,walltime=2:00:00 HighHbA1c_PC1_top2CT_HbA1c_cojo_top1_2021-10-05.sh