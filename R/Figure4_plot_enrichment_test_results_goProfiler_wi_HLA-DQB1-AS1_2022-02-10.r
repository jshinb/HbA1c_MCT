###############################################################################
#
# Created on Dec 8, 2021
# To create bar graphs on Figure 5 (working-version)
# This may change later and not be used in the end (but just in case) 
# Reference: https://www.nature.com/articles/s41598-020-73326-3 (for prioritizing terms)
# Modified on Jan 4, 2022
# 
###############################################################################

rm(list=ls());gc()
# setting working directory ---------------------------------------------------
# after running 'MHC_coexpressed_genes_for_GO_BP_enrichment_analysis.r'
wd = '/Users/jshin/OneDrive\ -\ SickKids/ukbb_insulin_resistance/PANTHER_HLAgenes/gProfiler_res_1-6-2022'
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
library("viridis") 
library(ggpattern)

library(ontoProc)
library(ontologyIndex)
library(ontologyPlot)
library(GOxploreR)
go = getGeneOnto()

create_enrichment_barplot <- function(genei,direction,go,panel.label=NULL,ymax=7){
  #inputs-----------------------------------------------------------------------#
  # genei = gene_names[2]; direction = "neg"; go = getGeneOnto()
  #-----------------------------------------------------------------------------#
  gene_inputs = readxl::read_xlsx('../HLA-genes_used_for_gProfiler_2022-01-06.xlsx',
                                  sheet=paste(direction,"coexp",genei,sep="_"))
  nG = nrow(gene_inputs)
  fillcolor = ifelse(genei %in% c("HLA-DQB2","HLA-DQA2"), "black", "grey")
  print(fillcolor)
  fname=paste("gProfiler_1-6-2022",direction,"coexp",genei,sep="_")
  fname=paste(fname,"csv",sep=".")
  if( file.exists(fname) ){
    d=subset(fread(fname),source=="GO:BP" & term_size <1000);dim(d)
    # d=subset(fread(fname),source=="GO:BP" & (term_size >=10 & term_size <250));dim(d)
    if(nrow(d)>0){
      d_terms = d$term_id
      # get the numbers of child terms for each GO-BP terms
      d_GOplot = onto_plot2(go,d_terms)
      d_GOedges = unlist(lapply(d_GOplot@edgeL,sapply,length))
      d_GOedges[d_GOedges ==0]
      names(d_GOedges) = str_remove(names(d_GOedges),"[.]edges")
      names(d_GOedges) = str_split(names(d_GOedges),":",simplify = T)[,2]
      # names(d_GOedges) = str_split(names(d_GOedges),":",simplify = T)[,2]
      names(d_GOedges) = paste("GO",names(d_GOedges),sep=":")
      setdiff(names(d_GOedges),d_terms)
      priotized_d_terms = prioritizedGOTerms(d_terms, organism = "Human", sp=TRUE, domain = "BP")
      # priotized_d_terms = prioritizedGOTerms(d_terms, organism = "Human", sp=FALSE, domain = "BP")
      d_specific = subset(d, term_id %in% priotized_d_terms$HF)
      # plot with d_specific --------------------------------------------------------
      long_term = "antigen processing and presentation of exogenous peptide antigen via MHC class II"
      long_term_2line = "antigen processing and presentation of\nexogenous peptide antigen via MHC class II"
      
      d_specific = d_specific %>% mutate(term_name = ifelse(term_name==long_term,long_term_2line,term_name))
      terms_ordered = d_specific$term_name[order(d_specific$adjusted_p_value)]
      theme_set(theme_bw())
      d_ggtitle = paste(panel.label,genei,": ",str_to_sentence(direction),"tively (N = ",nG,")",sep="")
      if(direction=="neg"){
        d_specific = d_specific %>% mutate(y=negative_log10_of_adjusted_p_value)
        p_ylim = c(0,ymax)
      }else{#pos
        d_specific = d_specific %>% mutate(y=negative_log10_of_adjusted_p_value)
        p_ylim = c(0,ymax)
      }
      p <- d_specific[order(d_specific$adjusted_p_value),] %>% 
        mutate(term_name=factor(term_name,levels=rev(terms_ordered))) %>%
        ggplot(aes(x=term_name, y=y)) + #,fill=negative_log10_of_adjusted_p_value
        geom_bar(stat='identity', width=.5, fill=fillcolor)  +
        # theme(axis.text.x = element_text(angle=-30,hjust = 0,vjust = 0),
        #       plot.margin = margin(1, 10, 0, 0.5, "cm")) + 
        xlab(NULL) + 
        ylab("-log10(adjusted p-value)") + 
        #ggtitle(d_ggtitle) + 
        ylim(p_ylim) +
        coord_flip() #+ 
      # scale_fill_gradient(name="-log10(adjusted p-value)",low="yellow",high="#00ba38",
      #                     limits=c(0,11))  
      ret = list(gene=genei, p=p, d_specific=d_specific, d=d,nG=nG)
    }else{
      p = plot_spacer()
      ret = list(gene=genei, p=p ,d_specific=NA, d=d, message="no enriched GO-BP terms",nG=nG)#no enriched term
    }#if/else(nrow(d)>0)
  }else{
    p <- plot_spacer()
    ret = list(gene=genei, p=p, d_specific=NA,d=NA,message='not enough genes in the test set',nG=nG)
  }#if/else(file.exists(fname) )
  ret
}
# -----------------------------------------------------------------------------
gene_names = c("HLA-DRB5", "HLA-DRB1","HLA-DQB1-AS1", "HLA-DQB2", "HLA-DQA2");L=length(gene_names)
gene_names = factor(unique(gene_names),levels=gene_names)
gene_names = sort(gene_names)
#inputs-----------------------------------------------------------------------#
#create_enrichment_barplo: genei = gene_names[2]; #direction = "pos"; #go = getGeneOnto()
#-----------------------------------------------------------------------------#
p_pos <- p_neg <- vector(mode="list",length=L)
for(gi in 1:L){
  p_pos[[gi]] = create_enrichment_barplot(genei=gene_names[[gi]],direction = "pos", go=go) 
  p_neg[[gi]] = create_enrichment_barplot(genei=gene_names[[gi]],direction = "neg", go=go)
}
names(p_pos) <- names(p_neg) <- gene_names

# DR genes
h_plot <- c()
for(gi in 1:L){
  hp = nrow(p_pos[[gi]]$d_specific)
  hp = ifelse (any(hp),hp,NA)
  hn = nrow(p_neg[[gi]]$d_specific)
  hn = ifelse (any(hn),hn,NA)
  h_plot = c(h_plot,c(hp,hn))
  print(h_plot)
}
h_plot[is.na(h_plot)] <- 2
t_margin <- b_margin <- 1
font.size = 14
ymax=7
p_ylim_pos = c(ymax,0)
p_ylim_neg = c(0,-ymax)

p1 = p_pos[[1]]$p+
  theme(legend.position = "n")+
  ylab(NULL)+theme(axis.ticks.y=element_blank())+ 
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=font.size)) + 
  theme(plot.margin = unit(c(t_margin, 2.5, b_margin, 0), "pt")) 

p2 = p_neg[[1]]$p+
  theme(legend.position = "n")+
  ylab(NULL)+
  scale_x_discrete(position = "top") + 
  theme(axis.ticks.y=element_blank()) +
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=font.size))+ 
  theme(plot.margin = unit(c(t_margin, 0, b_margin, 2.5), "pt"))

p3 = p_pos[[2]]$p+theme(legend.position = "n")+
  ylab(NULL)+
  theme(axis.ticks.y=element_blank())+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=font.size))+ 
  theme(plot.margin = unit(c(t_margin, 2.5, b_margin, 0), "pt"))

p4 = p_neg[[2]]$p+theme(legend.position = "n")+
  ylab(NULL)+
  scale_x_discrete(position = "top")+
  theme(axis.ticks.y=element_blank())+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=font.size))+ 
  theme(plot.margin = unit(c(t_margin, 0, b_margin, 2.5), "pt"))  

p5 = p_pos[[3]]$p+
  theme(legend.position = "n")+
  ylab(NULL)+
  theme(axis.ticks.y=element_blank())+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=font.size))  + 
  theme(plot.margin = unit(c(t_margin, 2.5, b_margin, 0), "pt"))

p6 = p_neg[[3]]$p

p7 = p_pos[[4]]$p

p8 = p_neg[[4]]$p+
  theme(legend.position = "n")+
  ylab(NULL)+
  scale_x_discrete(position = "top")+
  theme(axis.ticks.y=element_blank())+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=font.size))+ 
  theme(plot.margin = unit(c(t_margin, 0, b_margin, 2.5), "pt"))  
p8_axis = p_neg[[4]]$p+
  theme(legend.position = "n")+
  ylab(NULL)+
  scale_x_discrete(position = "top")+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.y = element_text(size=font.size),axis.text.x = element_text(size=font.size))+ 
  theme(plot.margin = unit(c(t_margin, 0, b_margin, 2.5), "pt")) 

p9 = p_pos[[5]]$p+
  theme(legend.position = "n")+
  ylab(NULL)+
  theme(axis.ticks.y=element_blank())+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=font.size))+ 
  theme(plot.margin = unit(c(t_margin, 2.5, b_margin, 0), "pt"))  

p9_axis = p_pos[[5]]$p+
  theme(legend.position = "n")+
  ylab(NULL)+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.y = element_text(size=font.size),axis.text.x = element_text(size=font.size))+ 
  theme(plot.margin = unit(c(t_margin, 2.5, b_margin, 0), "pt")) 

p10 = p_neg[[5]]$p

height = c(7,15,2,7,23)
p <- ((p1 + p2) + ggtitle(""))/ (p3 + p4+ ggtitle(""))/ (p5+ ggtitle("") + p6) / (p7 + p8+ ggtitle("")) / (p9+ ggtitle("") + p10+ ggtitle("")) + 
  # plot_layout(widths = c(1, 1), heights = height/2)
  plot_layout(heights = height/2)

p_axis <- ((p1 + p2) + ggtitle(""))/ (p3 + p4+ ggtitle(""))/ (p5+ ggtitle("") + p6) / (p7 + p8_axis+ ggtitle("")) / (p9_axis+ ggtitle("") + p10+ ggtitle("")) + 
  # plot_layout(widths = c(1, 1), heights = height/2)
  plot_layout(heights = height/2)

png(paste("HLA_genes_most_specific_enriched_terms_",Sys.Date(),".png",sep=''),
    width=14,height=16,units='in',res=300)
print(p)
dev.off()

png(paste("HLA_genes_most_specific_enriched_terms_",Sys.Date(),"_wi_axis.png",sep=''),
    width=14,height=16,units='in',res=300)
print(p_axis)
dev.off()

nG <- c()
for(gi in 1:L){
  # nG <- rbind(nG,c(p_pos[[gi]]$nG,p_neg[[gi]]$nG))
  nG <- c(nG,c(p_pos[[gi]]$nG,p_neg[[gi]]$nG))
}
cat(paste("n = ",nG,sep=''),sep='\n')

#
vnames=c("source", "term_name", "term_id", "adjusted_p_value", "negative_log10_of_adjusted_p_value", 
         "term_size", "query_size", "intersection_size", "effective_domain_size", 
         "intersections", "y")
text.cols = viridis(L);names(text.cols) = gene_names

d_specific_all_genes = NULL
for(i in 1:L){
  gi = gene_names[i]
  d_pos = (p_pos[[ gi ]]$d_specific)
  if(typeof(d_pos)=="logical"){
    d_pos = matrix(rep(NA,length(vnames)),nrow=1)
    colnames(d_pos) = vnames
    d_pos = data.table(d_pos)
  }
  d_pos$gene_name = gi
  d_pos$direction = "pos"
  d_pos$text.col = text.cols[i]
  
  d_neg = (p_neg[[ gi ]]$d_specific)
  if(typeof(d_neg)=="logical"){
    d_neg = matrix(rep(NA,length(vnames)),nrow=1)
    colnames(d_neg) = vnames
    d_neg = data.table(d_neg)
  }
  d_neg$gene_name = gi
  d_neg$direction = "neg"
  d_neg$text.col = scales::alpha(text.cols[i],0.8)
  
  d_tmp = rbind(d_pos,d_neg)
  d_specific_all_genes = rbind(d_specific_all_genes,d_tmp)
}

d_specific_all_genes$term_name[is.na(d_specific_all_genes$source)] <- c("*","**","***")#paste(NA,1:3,sep='')
dup.term_name = 'antigen processing and presentation of\nexogenous peptide antigen via MHC class II'
dup.term_name2 = paste(dup.term_name,"*",sep='')
d_specific_all_genes= d_specific_all_genes %>%
  mutate(y = ifelse(direction=="pos",negative_log10_of_adjusted_p_value,-negative_log10_of_adjusted_p_value)) %>%
  mutate(term_name = ifelse(term_name == dup.term_name & gene_name=="HLA-DQB1-AS1",dup.term_name2,term_name))
terms_ordered = d_specific_all_genes$term_name
d_specific_all_genes = d_specific_all_genes %>% 
  mutate(term_name2=factor(term_name,levels=rev(terms_ordered)))
d_specific_all_genes = d_specific_all_genes %>% arrange(term_name2)
d_specific_all_genes = d_specific_all_genes %>% 
  mutate(gene_direction=paste(gene_name,direction,sep="_"))
# text.col = c(rep(c("#018571","#5E3C99"),36),"#018571")#rev(d_specific_all_genes$text.col)
text.col = c(rep(c("black",grey(0.2)),36),"black")#rev(d_specific_all_genes$text.col)
p <-  d_specific_all_genes %>% 
  mutate(x=term_name2) %>%
  mutate(direction = factor(direction,levels=c('pos','neg'))) %>%
  filter(gene_name %in% gene_names[1:3]) %>%
  ggplot(aes(x=x, y=y, color=gene_name, fill=gene_direction)) + 
  geom_bar(stat='identity', width=.5)+#, fill=fillcolor)  +
  scale_color_manual(values=c("HLA-DRB5"="grey","HLA-DRB1"="grey","HLA-DQB1-AS1"="grey",
                              "HLA-DQB2"="black","HLA-DQA2"="black")) +
  scale_fill_manual(values=c("HLA-DRB5_pos"="grey","HLA-DRB1_pos"="grey","HLA-DQB1-AS1_pos"="grey",
                             "HLA-DQB2_pos"="black","HLA-DQA2_pos"="black",
                             "HLA-DRB5_neg"="white","HLA-DRB1_neg"="white","HLA-DQB1-AS1_neg"="white",
                             "HLA-DQB2_neg"="white","HLA-DQA2_neg"="white")) +
  theme(axis.text.x = element_text(angle=0,hjust = 0.5,vjust = 0.5,colour=text.col),
        axis.text.y = element_text(size=10,angle=0,hjust=1,vjust=0.5,colour=text.col),
        legend.position = 'none',
        plot.margin = margin(t=5.5, r=(5.5*0), b=5.5, l=5.5, "pt"))+
  scale_x_discrete(labels = function(x) stringr::str_remove_all(x,"[*]")) +
  #scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) + 
  # facet_grid( ~ gene_name)+
  xlab(NULL) + 
  ylab("-log10(adjusted p-value)") + 
  # ylim(-7,7) +
  scale_y_continuous(breaks = c(-4,0,4),labels=c('4','0','4'),limits = c(-7,7)) +
  coord_flip()
p1 = p

#

p <- d_specific_all_genes %>% 
  mutate(x=term_name2) %>%
  mutate(direction = factor(direction,levels=c('pos','neg'))) %>%
  filter(gene_name %in% gene_names[-c(1:3)]) %>%
  ggplot(aes(x=x, y=y, color=gene_name, fill=gene_direction)) + 
  geom_bar(stat='identity', width=.5)+#, fill=fillcolor)  +
  scale_color_manual(values=c("HLA-DRB5"="grey","HLA-DRB1"="grey","HLA-DQB1-AS1"="grey",
                              "HLA-DQB2"="black","HLA-DQA2"="black")) +
  scale_fill_manual(values=c("HLA-DRB5_pos"="grey","HLA-DRB1_pos"="grey","HLA-DQB1-AS1_pos"="grey",
                             "HLA-DQB2_pos"="black","HLA-DQA2_pos"="black",
                             "HLA-DRB5_neg"="white","HLA-DRB1_neg"="white","HLA-DQB1-AS1_neg"="white",
                             "HLA-DQB2_neg"="white","HLA-DQA2_neg"="white")) +
  # facet_grid( ~ gene_name)+
  theme(axis.text.x = element_text(angle=0,hjust = 0.5,vjust = 0.5,colour=text.col),
        axis.text.y = element_text(size=10,angle=0,hjust=1,vjust=0.5,colour=text.col),
        legend.position = 'none',
        plot.margin = margin(t=5.5, r=(5.5*0), b=5.5, l=5.5, "pt"))+
  scale_x_discrete(labels = function(x) stringr::str_remove_all(x,"[*]")) +
  #scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) 
  xlab(NULL) + 
  ylab("-log10(adjusted p-value)") + 
  # ylim(-7,7) +
  scale_y_continuous(breaks = c(-4,0,4),labels=c('4','0','4'),limits = c(-7,7)) +
  coord_flip() 
p2 = p

png(paste("HLA_genes_GTEx_",Sys.Date(),"_landscape.png",sep=''),
    height=8.5,width=15,units='in',res=300)
print(p1 + p2)
dev.off()