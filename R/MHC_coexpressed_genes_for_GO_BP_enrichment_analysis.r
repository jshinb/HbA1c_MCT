options(stringsAsFactors = F)
library(stringr)
library(data.table)
library(tidyverse)
library(tableone)
library(dplyr)
library(grid)
library(gratia)
library(gridExtra)
library(GenABEL)

# for tables
library(arsenal)
library(here)
library(mgsub)
library(psych)

# functions -------------------------------------------------------------------#
get_cor_pearson = function(i) {
  pearson_r = cor(subset(t_gene_expression,select=select.cols[i]),t_gene_expression)
  pearson_r
}

get_cor_spearman_pvalue = function(x,i) {
  y=unlist(subset(t_gene_expression,select=select.cols[i]))
  pval = cor.test(x,y,method="s")$p.value
  pval
}

get_cor_pearson_pvalue = function(x,i) {
  y=unlist(subset(t_gene_expression,select=select.cols[i]))
  pval = cor.test(x,y,method="p")$p.value
  pval
}

#------------------------------------------------------------------------------#

setwd('/Users/jshin/OneDrive\ -\ SickKids/ukbb_insulin_resistance/PANTHER_HLAgenes')
gene_colums = fread('genes_matrix_csv/columns_metadata.csv')
gene_rows = fread('genes_matrix_csv/rows_metadata.csv')
head(gene_colums)
head(gene_rows)

ctx.rois = c("OFC", #orbital frontal 
             "DFC", #dorsolateral
             "MFC", #medial frontal
             "VFC", #ventrolateral prefrontal
             "M1C", #primaray motor
             "S1C", #primary somatosensory
             "M1C-S1C", 
             "IPC", #posteroinferior parietal,
             "A1C", #primary auditory 
             "STC", #posterior superior temporal
             "ITC", #inferior temporal
             "V1C") #primary visual (occipital)

cortex_ind = gene_colums$structure_acronym %in% ctx.rois
adult_ind = gene_colums$age %in% c("21 yrs", "23 yrs", "30 yrs", "36 yrs", "37 yrs","40 yrs")
sum(cortex_ind & adult_ind)#63

gene_expression = fread('genes_matrix_csv/expression_matrix.csv')
head(gene_expression[,1:3])
gene_expression = gene_expression %>% dplyr::select(-V1)
select.cols = names(gene_expression)[which(cortex_ind & adult_ind)]
gene_expression = subset(gene_expression,select=select.cols)

##
gene.names <- 'HLA-DQA2
HLA-DQB1-AS1
HLA-DQB2
HLA-DRB1
HLA-DRB5
HLA-DRB6'
gene.names<- MHC_genes <- unlist(str_split(gene.names,'\n'))

t_gene_expression = data.table(t(gene_expression) )
sd.gene.expression = apply(t_gene_expression,2,sd)
gene_rows_adult_cortex = subset(gene_rows,sd.gene.expression>0)#44,049
t_gene_expression = subset(t_gene_expression,
                           select=paste("V",gene_rows_adult_cortex$row_num,sep=""))

gene.rows.test = subset(gene_rows_adult_cortex,gene_symbol %in% gene.names)
select.cols = paste("V",gene.rows.test$row_num,sep="")
cormat.test = NULL
for(i in 1:nrow(gene.rows.test)){
  cormat.test = rbind(cormat.test,get_cor_pearson(i))
}

cormat.test = data.table(t(cormat.test))
cormat.test = data.table(ensembl_gene_id=gene_rows_adult_cortex$ensembl_gene_id,
                         gene_symbol=gene_rows_adult_cortex$gene_symbol,
                         cormat.test)

subset(cormat.test,gene_symbol %in% gene.rows.test$gene_symbol)
names(cormat.test)[-c(1:2)] = gene.rows.test$gene_symbol

# spearman correlation p-values ------------------------------------------
cormat.test.p = apply(t_gene_expression,2,get_cor_spearman_pvalue,i=1)
for(i in 2:nrow(gene.rows.test)){
  print(i)
  tmp.p = apply(t_gene_expression,2,get_cor_spearman_pvalue,i=i)
  cormat.test.p = cbind(cormat.test.p,tmp.p)
}
dim(cormat.test.p)

cormat.test.p = data.table(cormat.test.p)
names(cormat.test.p) = gene.rows.test$gene_symbol

cormat.test.p = data.table(ensembl_gene_id = gene_rows_adult_cortex$ensembl_gene_id,
                           gene_symbol = gene_rows_adult_cortex$gene_symbol,
                           cormat.test.p)

subset(cormat.test.p,gene_symbol %in% gene.rows.test$gene_symbol)
names(cormat.test.p)[-c(1:2)] = gene.rows.test$gene_symbol
spearman.cor.p = cormat.test.p
rm(cormat.test.p)

# pearson correlation p-values ------------------------------------------
cormat.test.p = apply(t_gene_expression,2,get_cor_pearson_pvalue,i=1)
for(i in 2:nrow(gene.rows.test)){
  print(i)
  tmp.p = apply(t_gene_expression,2,get_cor_pearson_pvalue,i=i)
  cormat.test.p = cbind(cormat.test.p,tmp.p)
}
dim(cormat.test.p)

cormat.test.p = data.table(cormat.test.p)
names(cormat.test.p) = gene.rows.test$gene_symbol

cormat.test.p = data.table(ensembl_gene_id = gene_rows_adult_cortex$ensembl_gene_id,
                           gene_symbol = gene_rows_adult_cortex$gene_symbol,
                           cormat.test.p)

subset(cormat.test.p,gene_symbol %in% gene.rows.test$gene_symbol)
names(cormat.test.p)[-c(1:2)] = gene.rows.test$gene_symbol
pearson.cor.p = cormat.test.p

# adjust p-values -------------------------------------------------------------
included.genes = c("HLA-DRB1", "HLA-DRB5", "HLA-DQB2", "HLA-DQA2")
p.FDR.spearman = p.adjust(unlist(subset(spearman.cor.p,select=included.genes)))
p.FDR.pearson = p.adjust(unlist(subset(pearson.cor.p,select=included.genes)))

p.FDR.spearman = matrix(p.FDR.spearman,ncol=4,byrow=F)
p.FDR.pearson = matrix(p.FDR.pearson,ncol=4,byrow=F)
colnames(p.FDR.spearman) <- colnames(p.FDR.pearson) <- included.genes
rownames(p.FDR.spearman) <- rownames(p.FDR.pearson) <- cormat.test.p$gene_symbol
p.FDR.spearman = data.table(ensembl_gene_id=cormat.test.p$ensembl_gene_id,gene_symbol = cormat.test.p$gene_symbol,p.FDR.spearman)
p.FDR.pearson = data.table(ensembl_gene_id=cormat.test.p$ensembl_gene_id,gene_symbol = cormat.test.p$gene_symbol,p.FDR.pearson)

FDR.siglevel = 0.05
#gi=1; 
for(gi in 1:4){
  genei = included.genes[gi];print(genei)
  
  ind = cormat.test[[genei]] >0;print(sum(ind))
  ind = ind & p.FDR.pearson[[genei]] <FDR.siglevel;print(sum(ind))#1917,328,742,163
  ind = ind & cormat.test$gene_symbol != genei;print(sum(ind))
  ind_pos = ind
  
  ofile_pos = paste(wd,"PANTHER_HLAgenes/pos_coexp_",genei,".tsv",sep='')
  out_pos = data.table(gene=cormat.test$ensembl_gene_id[ind_pos],logp=-log10(p.FDR.pearson[[genei]][ind_pos]))
  out_pos = out_pos[order(out_pos$logp,decreasing=T),]
  write_tsv(out_pos,file=ofile_pos)
  
  ind = cormat.test[[genei]] <0;print(sum(ind))
  ind = ind & p.FDR.pearson[[genei]] <FDR.siglevel;print(sum(ind))#916,241,455
  ind_neg = ind
  
  ofile_neg = paste(wd,"PANTHER_HLAgenes/neg_coexp_",genei,".tsv",sep='')
  out_neg = data.table(gene=cormat.test$ensembl_gene_id[ind_neg],logp=-log10(p.FDR.pearson[[genei]][ind_neg]))
  out_neg = out_neg[order(out_neg$logp,decreasing=T),]
  write_tsv(out_neg,file=ofile_neg)
}

# approach 2 ------------------------------------------------------------------
library(tidyr)
cormat_long <- gather(cormat.test, key=MHC_genes, value = r, `HLA-DRB1`:`HLA-DQA2`,factor_key=FALSE)
head(cormat_long )
hist(cormat_long$r)
print(cormat.quantile <- quantile(cormat_long$r, probs=c(0.01,0.99)))

poscor_genes = unique(cormat_long$ensembl_gene_id[cormat_long$r> cormat.quantile[2]]);print(length(poscor_genes))#1220
negcor_genes = unique(cormat_long$ensembl_gene_id[cormat_long$r< cormat.quantile[1]]);print(length(negcor_genes))#1224
negcor_genes = unique(cormat_long$ensembl_gene_id[cormat_long$r< -cormat.quantile[2]]);print(length(negcor_genes))#1224

write.table(poscor_genes,"~/Downloads/poscor_MHC_genes_v2.txt",
            col.names=F,row.names=F,quote=F)
write.table(negcor_genes,"~/Downloads/negcor_MHC_genes_v2.txt",
            col.names=F,row.names=F,quote=F)

mat_t_gene_expression = as.matrix(t_gene_expression)
# random.quantile = NULL
for(i in 1:5000){
  random.genes = sample(names(t_gene_expression),6,replace = F)
  r.random = apply(mat_t_gene_expression[,random.genes],2,cor,y=mat_t_gene_expression)
  write.header = F
  do.append = T
  if(i==1) {
    write.header = T
    do.append = F}
  write.table(rbind(quantile(r.random,probs = c(0.005,0.995))),file='random.quantile.txt',sep="\t",
                            quote=F,col.names = write.header, row.names=F,append=do.append)
  #random.quantile <- rbind( random.quantile, quantile(r.random,probs = c(0.005,0.995)) )
}
random.quantile = fread('random.quantile.txt')
genes = gene.rows.test$gene_symbol
median.random.quantile = apply(random.quantile,2,median)#
#         0.5%      99.5% 
#   -0.6244163  0.7455479 
poscor = cormat.test[[genes[1]]] > median.random.quantile[2];print(sum(poscor))
for(i in c(2:length(genes))){
  print(genes[i])
  poscor = poscor | cormat.test[[genes[i]]] > median.random.quantile[2];print(sum(poscor))
}

write.table(gene_rows_adult_cortex$ensembl_gene_id[poscor],
            paste("~/Downloads/poscor_MHC_genes_",median.random.quantile[2],".txt",sep=""),
            col.names=F,row.names=F,quote=F)
#
pos_genes = NULL
for(i in 1:length(genes)){
  o <- order(cormat.test[[genes[i]]],decreasing = T)
  pos_genes = cbind(pos_genes,cormat.test$ensembl_gene_id[o][1:441])
}
pos_genes = data.frame(pos_genes)
names(pos_genes)=genes
write.table(pos_genes, paste("~/Downloads/poscor_MHC_genes_equalNoGenes.txt"),
            sep="\t", col.names=T,row.names=F, quote=F)
#
neg_genes = NULL
for(i in 1:length(genes)){
  o <- order(cormat.test[[genes[i]]],decreasing = F)
  neg_genes = c(neg_genes,cormat.test$ensembl_gene_id[o][1:74])
}
neg_genes = unique(neg_genes)#368
write.table(neg_genes,
            paste("~/Downloads/negcor_MHC_genes_equalNoGenes.txt",sep=""),
            col.names=F,row.names=F,quote=F)

#
negcor = cormat.test[[genes[1]]] < random.quantile[1];print(sum(negcor))
for(i in c(2:length(genes))){
  negcor = negcor | cormat.test[[genes[i]]] < median(random.quantile[,1]);print(sum(negcor))
}
write.table(gene_rows_adult_cortex$ensembl_gene_id[negcor],
            paste("~/Downloads/negcor_MHC_genes_",median(random.quantile[,1]),".txt",sep=""),
            col.names=F,row.names=F,quote=F)


###
library(readxl)
poscor_genes <- read_excel("PANTHER_GO_BP_analysis_poscor_MHC_genes_2021-07-02.xlsx", 
                           sheet = "Genes")

cormat_long %>% ggplot(aes(x=r,color=MHC_genes, fill=MHC_genes)) + 
  # geom_histogram(position = 'identity',alpha=0.3) + 
  # geom_histogram(aes(y =..density..), 
  #                alpha=0.2, 
  #                position='identity') + 
  geom_density(alpha=0.2) + 
  geom_vline(xintercept = c(random.quantile,0))

pos_select.cols = gene_rows_adult_cortex$row_num[gene_rows_adult_cortex$ensembl_gene_id %in% gene_rows_adult_cortex$ensembl_gene_id[poscor]]
pos_select.cols  = paste('V',pos_select.cols,sep='')

mRNA_pos = subset(t_gene_expression,select=pos_select.cols)
mRNA_pos = data.frame(mRNA_pos)
colnames(mRNA_pos) <- gene_rows_adult_cortex$gene_symbol[poscor]
cormat_pos = cor(mRNA_pos)
dim(cormat_pos)

genes.random = c(gene.rows.test$gene_symbol,sample(colnames(cormat_pos),40,replace = F))
genes.random = unique(genes.random)
pdf('corrplot_pos_309genes.pdf',height=11,width=8)
corrplot(cormat_pos[genes.random,genes.random],order = 'hclust', hclust.method = 'ward.D2',tl.cex = 0.7)
dev.off()

write.table(table(cormat_long$MHC_genes[cormat_long$r>median.random.quantile[2]],
      cormat_long$gene_symbol[cormat_long$r>median.random.quantile[2]]),
      'pos_table.txt',quote=F,row.names = T,col.names = T, sep="\t")

MHC.genes=gene.rows.test$gene_symbol
op = par(mfrow=c(2,3))
for(gi in MHC.genes){
  hist(cormat.test[[gi]],xlim=c(-1,1),
       main=gi,xlab='')
  abline(v=median.random.quantile,lty=2,col="darkgreen")
}
title(xlab='correlation coefficient (r) of 44,049 genes',outer=T,line=-1.5)
par(op)

cormat.test.matrix = as.matrix(subset(cormat.test,gene_symbol %in% MHC.genes,select=-c(ensembl_gene_id,gene_symbol)))
rownames(cormat.test.matrix) = colnames(cormat.test.matrix)
corrplot::corrplot(cormat.test.matrix,method = 'number',type='upper',
                   #order='hclust',hclust.method = 'ward.D2',
                   tl.cex = 0.8, tl.col="grey20", tl.srt = 30)
tmp.mat = subset(cormat.test,select=-c(ensembl_gene_id,gene_symbol))
# for(gi in MHC.genes){
#   print(max(tmp.mat[[gi]]))
#   if(!is.na(max(tmp.mat[[gi]])))
#     tmp.mat[[gi]][tmp.mat[[gi]]==max(tmp.mat[[gi]])] <- NA
# }
apply(tmp.mat,2,quantile,probs=c(0.01,0.99))
