###############################################################################
#
# Created on Dec 8, 2021
# To create bar graphs on Figure 5 (working-version)
# This may change later and not be used in the end (but just in case) 
#
###############################################################################

rm(list=ls());gc()
# setting working directory ---------------------------------------------------
# after running 'MHC_coexpressed_genes_for_GO_BP_enrichment_analysis.r'
wd = '/Users/jshin/OneDrive\ -\ SickKids/ukbb_insulin_resistance/PANTHER_HLAgenes/gProfiler_res_12-8-2021'
setwd(wd)

# load packages and define any functions --------------------------------------
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

# for GO table/plot
library(ggplot2)
library(patchwork)
library(ontoProc)
library(ontologyIndex)
library(ontologyPlot)
library(GOxploreR)

create_enrichment_barplot <- function(genei,direction,go,panel.label=NULL){
  #inputs-----------------------------------------------------------------------#
  # genei = gene_names[2]; direction = "neg"; go = getGeneOnto()
  #-----------------------------------------------------------------------------#
  fname=paste("gProfiler_12-8-2021",direction,"coexp",genei,sep="_")
  fname=paste(fname,"csv",sep=".")
  if( file.exists(fname) ){
    d=subset(fread(fname),source=="GO:BP");dim(d)
    if(nrow(d)>0){
      d_terms = d$term_id
      # get the numbers of child terms for each GO-BP terms
      d_GOplot = onto_plot2(go,d_terms)
      d_GOedges = unlist(lapply(d_GOplot@edgeL,sapply,length))
      d_GOedges[d_GOedges ==0]
      names(d_GOedges) = str_remove(names(d_GOedges),"[.]edges")
      names(d_GOedges) = str_split(names(d_GOedges),":",simplify = T)[,2]
      names(d_GOedges) = paste("GO",names(d_GOedges),sep=":")
      setdiff(names(d_GOedges),d_terms)
      
      # d_specific = subset(d,term_id %in% names(d_GOedges[d_GOedges==0]))
      # d_specific[,1:5]
      # length(d_GOedges)
      # 
      # d_GOedges[setdiff(names(d_GOedges),d_terms)]
      # 
      # cat(d_terms,sep=",")
      # d_specific[,1:4]
      priotized_d_terms = prioritizedGOTerms(d_terms, organism = "Human", sp=TRUE, domain = "BP")
      d_specific = subset(d, term_id %in% priotized_d_terms$HF)
      
      # plot with d_specific --------------------------------------------------------
      terms_ordered = d_specific$term_name[order(d_specific$adjusted_p_value)]
      theme_set(theme_bw())
      gene_inputs = readxl::read_xlsx('../HLA-genes_used_for_gProfiler_2021-12-08.xlsx',
                                      sheet=paste(direction,"coexp",genei,sep="_"))
      nG = nrow(gene_inputs)
      d_ggtitle = paste(panel.label,genei,": ",str_to_sentence(direction),"tively (N = ",nG,")",sep="")
      p <- d_specific[order(d_specific$adjusted_p_value),] %>% 
        mutate(term_name=factor(term_name,levels=rev(terms_ordered))) %>%
        ggplot(aes(x=term_name, y=negative_log10_of_adjusted_p_value,fill=negative_log10_of_adjusted_p_value)) + 
        geom_bar(stat='identity', width=.5)  +
        # theme(axis.text.x = element_text(angle=-30,hjust = 0,vjust = 0),
        #       plot.margin = margin(1, 10, 0, 0.5, "cm")) + 
        xlab(NULL) + 
        ylab("-log10(adjusted p-value)") + 
        ggtitle(d_ggtitle) + 
        ylim(0,16) +
        coord_flip() + 
        scale_fill_gradient(name="-log10(adjusted p-value)",low="yellow",high="#00ba38",
                            limits=c(0,16))  
      ret = list(gene=genei, p=p, d_specific=d_specific, d=d)
    }else{
      ret = list(gene=genei, p=NA,d_specific=NA, d=d, message="no enriched GO-BP terms")#no enriched term
    }#if/else(nrow(d)>0)
  }else{
    ret = list(gene=genei, p=NA,d_specific=NA,d=NA,message='not enough genes in the test set')
  }#if/else(file.exists(fname) )
  ret
}

# -----------------------------------------------------------------------------
gene_names = c("HLA-DRB5", "HLA-DRB1", "HLA-DQA2", "HLA-DQB2")
gene_names = factor(unique(gene_names),levels=c("HLA-DRB5", "HLA-DRB1", "HLA-DQA2", "HLA-DQB2"))
gene_names = sort(gene_names)
go = getGeneOnto()
#inputs-----------------------------------------------------------------------#
#create_enrichment_barplo: genei = gene_names[2]; #direction = "pos"; #go = getGeneOnto()
#-----------------------------------------------------------------------------#
p_pos <- p_neg <- vector(mode="list",length=length(gene_names))
for(gi in 1:length(gene_names)){
  p_pos[[gi]] = create_enrichment_barplot(genei=gene_names[[gi]],direction = "pos", go=go)
  p_neg[[gi]] = create_enrichment_barplot(genei=gene_names[[gi]],direction = "neg", go=go)
}
names(p_pos) <- names(p_neg) <- gene_names
gi=2

# DR genes
h_DR = c(nrow(p_pos[["HLA-DRB5"]]$d_specific),
nrow(p_neg[["HLA-DRB5"]]$d_specific),
nrow(p_pos[["HLA-DRB1"]]$d_specific),
nrow(p_neg[["HLA-DRB1"]]$d_specific))

p_DR = p_pos[[1]]$p+theme(legend.position = "n")+ylab(NULL)+
  p_neg[[1]]$p+theme(legend.position = "n")+ylab(NULL)+
  p_pos[[2]]$p+theme(legend.position = "n")+ylab(NULL)+
  p_neg[[2]]$p+theme(legend.position = "n")+
  plot_layout(heights = h_DR/2)

png(paste("HLA-DR_most_specific_enriched_terms_",Sys.Date(),".png",sep=''),
       width=8,height=8,units='in',res=300)
print(p_DR)
dev.off()

# DQ genes
h_DQ = c(pos_DQA2 = nrow(p_pos[[3]]$d_specific),
         neg_DQA2 = 2,
         pos_DQB2 = 2,
         neg_DQB2 = nrow(p_neg[[4]]$d_specific))

p_DQ = p_pos[[3]]$p+theme(legend.position = "n")+ylab(NULL)+
  plot_spacer()+ggtitle("tmp") +
  plot_spacer()+
  p_neg[[4]]$p+theme(legend.position = "n")+
  plot_layout(heights = h_DQ/2)

png(paste("HLA-DQ_most_specific_enriched_terms_",Sys.Date(),".png",sep=''),
    width=8,height=8,units='in',res=300)
print(p_DQ)
dev.off()
