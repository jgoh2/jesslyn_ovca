# Returns top differentially expressed genes between two conditions. Also returns volcano plots when specified. 
# Date: 05/25/2020
# Author: Jesslyn Goh

library(tidyverse) # ggplot, dplyr, read_r (read_xlsx)
library(Seurat)
library(readxl)
library(ggplot2)
library(here)
library(dplyr)
library(EnhancedVolcano)
library(patchwork)

DEAnalysis_code <- function(seurat_object, markers = NULL, group.by, group.1, group.2, graph = FALSE, oxphosGenes = NULL, utrGenes = NULL){
  if(is.null(markers)){
  Idents(seurat_object) <- group.by
  #compute differential gene expression of all genes for ident.1 vs. ident.2 (wilcox)
  DEgenes <- FindMarkers(seurat_object, ident.1 = group.1 , ident.2 = group.2, logfc.threshold = 0)
  }
  
  if(!(is.null(markers))){
    DEgenes = markers 
  if (graph == TRUE){
  #volcano plot for p values of ALL genes and label ALL GENES
  all = EnhancedVolcano(DEgenes, lab=rownames(DEgenes), x='avg_logFC', y='p_val_adj', 
    pCutoff = 0.01, col = c("black", "black", "black", "red2"), 
    title= paste(group.1, "vs", group.2, "DE Genes (lab all Genes)"), subtitle= "FC cutoff: 1.0, p cutoff: 0.01")
  
  return (all)
  }
  
  if(!(is.null(oxphosGenes))){
    #only label OXPHOS genes 
    all_lab_oxphos = EnhancedVolcano(DEgenes, lab= rownames(DEgenes), x='avg_logFC', y='p_val_adj', 
        pCutoff = 0.01, col = c("black", "black", "black", "red2"), selectLab=oxphosGenes, 
        title=paste(group.1, "vs", group.2, "DE Genes (lab oxphos genes)"), subtitle= "FC cutoff: 1.0, p cutoff: 0.01")
    
    #only graph OXPHOS genes
    DEgenes.oxphos = DEgenes %>% rownames_to_column("Gene") %>% filter(Gene %in% oxphosGenes) #filter oxphos genes
    DEgenes.oxphos = column_to_rownames(DEgenes.oxphos, "Gene")
    just_oxphos = EnhancedVolcano(DEgenes.oxphos,lab= rownames(DEgenes.oxphos), x='avg_logFC', y='p_val_adj', 
        pCutoff = 0.01, col = c("black", "black", "black", "red2"), 
        title=paste(group.1, "vs", group.2, "Oxphos DE Genes"), subtitle= "FC cutoff: 1.0, p cutoff: 0.01")
    
    return (all_lab_oxphos + just_oxphos)
  }
  
  if(!(is.null(utrGenes))){
    #only label UTR genes 
    all_lab_utr = EnhancedVolcano(DEgenes, lab= rownames(DEgenes), x='avg_logFC', y='p_val_adj', 
        pCutoff = 0.01, col = c("black", "black", "black", "red2"), selectLab=utrGenes, 
        title=paste(group.1, "vs", group.2, "DE Genes (lab utr genes)"), subtitle= "FC cutoff: 1.0, p cutoff: 0.01")
    
    #only graph UTR genes
    DEgenes.utr = DEgenes %>% rownames_to_column("Gene") %>% filter(Gene %in% utrGenes) #filter utr genes
    DEgenes.utr = column_to_rownames(DEgenes.utr, "Gene")
    just_utr = EnhancedVolcano(DEgenes.utr,lab= rownames(DEgenes.utr), x='avg_logFC', y='p_val_adj', 
        pCutoff = 0.01, col = c("black", "black", "black", "red2"), 
        title=paste(group.1, "vs", group.2, "UTR DE Genes"), subtitle= "FC cutoff: 1.0, p cutoff: 0.01")
    
    return (all_lab_utr + just_utr)
  }
  }
  return (DEgenes)
}
