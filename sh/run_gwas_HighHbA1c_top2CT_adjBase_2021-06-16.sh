#####
sfile=/home/t/tpaus/jshinb/ukbb_gwas_scripts/gwas_HighHbA1c_top2CT_adjBase_2021-06-16.sh

resdir=/scratch/t/tpaus/jshinb/PLINK/HighHbA1c_top2CT_adjBase_2021-06-16

mkdir $resdir
cd $resdir

adj=adjBase
for chr in {2..22}
do
sbatch $sfile $chr $resdir $adj
done

# runsfile=/home/t/tpaus/jshinb/ukbb_gwas_scripts/run_gwas_HighHbA1c_top2CT_adjBase_2021-06-16.sh;sh $runsfile