############## Disorder Data Extraction ################
# module load gcc/9.2.0 r/4.0.3;R
########################################################
rm(list=ls());gc()
library(tidyverse)
library(data.table)
library(stringr)
extract_variables = function(fname,fieldID,fieldName){
  if (!file.exists(fname)){stop ("file does not exist")
  }else{
    d = fread(fname)    
    d.sub = subset(d, select=c('eid',fieldID))
    names(d.sub) = c('eid',fieldName)
    d.sub
  }
}

ukbb_data_dir = '/project/t/tpaus/tpaus/UKBB/datasets'

# mri = fread('/gpfs/fs1/home/t/tpaus/jshinb/ukbb_gwas_data/freesurfer_combined.txt')
mri = fread('d_MCT_38864.tsv.gz')%>% dplyr::rename(ID=eid)
hba1c = fread('/gpfs/fs1/home/t/tpaus/jshinb/HbA1c_data.tsv')
hba1c_no_na = subset(hba1c,!is.na(V3),select=c(V1,V3))#n=502489
names(hba1c_no_na) = c("ID","HbA1c")#
mri_hba1c = merge(hba1c_no_na,mri)
print(dim(mri_hba1c))#n=36,059
mri_hba1c = mri_hba1c %>% dplyr::rename(eid=ID)

#
f='/project/t/tpaus/tpaus/UKBB/datasets/ukb_01-04-2020/41449/ukb41449.csv'
dd.genetic.eth = fread(f)
dd.genetic.eth = subset(dd.genetic.eth,select=c('eid','22006-0.0','22019-0.0','22001-0.0','22027-0.0','22005-0.0'))
names(dd.genetic.eth)[-1] = c('white.yes','aneuploidy','genetic.sex','heterozygosity','missing_call_rate')

# genotype
anal_data=subset(mri_hba1c,eid %in% dd.genetic.eth$eid[!is.na(dd.genetic.eth$genetic.sex)])
dim(anal_data)#35878
dd_anal = anal_data#mri_hba1c
dim(dd_anal)

# white
#table(dd.genetic.eth$white.yes[dd.genetic.eth$eid%in% dd_anal$eid])# 31252
#dd_anal= subset(dd_anal,eid %in% dd.genetic.eth$eid[!is.na(dd.genetic.eth$white.yes)])

# unrelated
kinship.info = fread('/gpfs/fs1/home/t/tpaus/jshinb/ukbb_gwas_data/ukb43688_rel_s488264.dat')
kinship.info = subset(kinship.info,Kinship>0)
ind2 = kinship.info$ID1 %in%kinship.info$ID2
ind3 = kinship.info$ID1 > kinship.info$ID2
ID1 = kinship.info$ID1
ID2 = kinship.info$ID2
kinship.info$ID1[ind3] <- ID2[ind3]
kinship.info$ID2[ind3] <- ID1[ind3]
kinship.info = kinship.info[order(kinship.info$ID1,kinship.info$ID2),]

kinship.anal = subset(kinship.info,ID1 %in% dd_anal$eid & kinship.info$ID2 %in% dd_anal$eid);print(dim(kinship.anal))#1101
# kinship.anal = subset(kinship.info,ID1 %in% idata$eid & kinship.info$ID2 %in% idata$eid);print(dim(kinship.anal))#1101
ids = unique(c(kinship.anal$ID1,kinship.anal$ID2));print(length(ids))#634

##
tmp.fam.ids = sort(unique(c(kinship.anal$ID1,kinship.anal$ID2)));length(tmp.fam.ids)#2064
mat = matrix(NA,nrow=length(tmp.fam.ids),ncol=length(tmp.fam.ids))
rownames(mat) <- colnames(mat) <- tmp.fam.ids
print(dim(mat))

for(i in 1:nrow(mat)){
  id1 = as.numeric(rownames(mat)[i]);#print(id1)
  for(j in c(1:nrow(mat))[-i]){
    id2 = as.numeric(rownames(mat)[j]);#print(id2)
    ind = (kinship.anal$ID1==id1 & kinship.anal$ID2==id2);#print(sum(ind))
    ind = ind | (kinship.anal$ID1==id2 & kinship.anal$ID2==id1)
    #print(sum(ind))#1
    if(sum(ind)==1){
      mat[i,j] <- mat[j,i] <- kinship.anal$Kinship[ind]
    }
  }
  cat("*",sep='')
}
nrel = c()
for(i in 1:nrow(mat)){
  nrel = c(nrel,sum(!is.na(mat[,i])))
}
diag(mat) <- 1

o = order(nrel,decreasing = T)
ids = rownames(mat)[o]
fam.ids = list()
i = 1
while(length(ids)>0){
  print(i)
  id = ids[1]
  fam.idsj = rownames(mat)[!is.na(mat[,id])]
  for(idj in fam.idsj[fam.idsj!=id]){
    print(idj)
    fam.idsj = c(fam.idsj,rownames(mat)[!is.na(mat[,idj])])
  }
  fam.idsj = unique(fam.idsj)
  mat[fam.idsj,fam.idsj]
  fam.ids[[i]] = fam.idsj
  ids = ids[!ids %in% fam.idsj]
  i = i+1
}

dup.inds = unlist(fam.ids)[duplicated(unlist(fam.ids))]
detect.dup.id = function(x,id){
  ret = any(x == id)
  ret
}

if(length(dup.inds)>0){
  for(i in 1:length(dup.inds)){
    print(which(sapply(fam.ids,detect.dup.id,id=dup.inds[i])))  
  }
  ## manually
  fam.ids_wo_dup = fam.ids#737 families
  length(fam.ids_wo_dup[[1]] )#9
  fam.ids_wo_dup[[1]] = unique(c(fam.ids_wo_dup[[1]],fam.ids_wo_dup[[737]]))
  length(fam.ids_wo_dup[[1]] )#10
  fam.ids_wo_dup[[737]] <- NA
  
  length(unlist(fam.ids_wo_dup)[!is.na(unlist(fam.ids_wo_dup))])#1511
  length(unlist(fam.ids_wo_dup[!is.na(fam.ids_wo_dup)]))#1511
  fam.ids_wo_dup = fam.ids_wo_dup[!is.na(fam.ids_wo_dup)]#736
  length(unique(unlist(fam.ids_wo_dup)))
  fam.ids = fam.ids_wo_dup
  rm(fam.ids_wo_dup)
}

include.ids = c()
for(i in 1:length(fam.ids)){
  sample.ids = fam.ids[[i]]
  include.ids = c(include.ids,sample(sample.ids,size = 1));print(length(include.ids))
}
remove.ids = unlist(fam.ids)[!unlist(fam.ids)%in% include.ids]#323
length(remove.ids)

dd_anal2 = subset(dd_anal,!eid %in% remove.ids);print(dim(dd_anal2))
dim(dd_anal2)#35555

# white
table(dd.genetic.eth$white.yes[dd.genetic.eth$eid%in% dd_anal2$eid],useNA='a')# 30962
#       1  <NA>
#   30964  4591
#35555 - 30962
# other QC:
dd_anal2= subset(dd_anal2,eid %in% dd.genetic.eth$eid[!is.na(dd.genetic.eth$white.yes)])
dd.genetic.eth = subset(dd.genetic.eth,eid %in% dd_anal2$eid)#30817
dim(dd.genetic.eth)
#
# disease
#==============================================================================#
# disease status for exclusion
#==============================================================================#
fe1=file.path(ukbb_data_dir,'ukb42388_18062020/ukb42388.csv')
var.names.e1 = c(
  "AD" = '42020-0.0',
  "dementia" = '42018-0.0',
  "stroke" = '42006-0.0')
e1 = extract_variables(fe1,fieldID=var.names.e1,fieldName=names(var.names.e1))
excl.dementia <- e1$eid[!is.na(e1$dementia)]
excl.stroke <- e1$eid[!is.na(e1$stroke)]
excl.dementia.stroke = unique(c(excl.dementia,excl.stroke))
print(length(excl.dementia.stroke))#15439
sum(excl.dementia.stroke %in% dd_anal2$eid)#326

dd_anal2 = subset(dd_anal2,!eid %in% excl.dementia.stroke)
print(dim(dd_anal2))#30638

dd.genetic.eth = subset(dd.genetic.eth,eid %in% dd_anal2$eid)#
dim(dd.genetic.eth)#

QC_rm_ind =!is.na(dd.genetic.eth$aneuploidy);print(sum(QC_rm_ind))#12
QC_rm_ind = QC_rm_ind | !is.na(dd.genetic.eth$heterozygosity);print(sum(QC_rm_ind))#59-12=47
dd_anal2 = subset(dd_anal2, !eid %in% dd.genetic.eth$eid[QC_rm_ind])#35293
print(dim(dd_anal2))#30577

dd.genetic.eth = subset(dd.genetic.eth,eid %in% dd_anal2$eid)#
table(dd.genetic.eth$white.yes,useNA="a")

load("data/idata_2022-02-03.RData")
nrow(idata)#30579
nrow(idata)-nrow(dd_anal2)
