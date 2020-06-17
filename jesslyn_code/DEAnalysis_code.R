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

DEAnalysis_code <- function(seurat_object, markers = NULL, group.by, group.1, group.2, graph = FALSE, geneset = NULL){
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
    pCutoff = 0.05, FCcutoff = 0.5, col = c("black", "black", "black", "red2"), 
    title= paste(group.1, "vs", group.2, "DE Genes (label all Genes)"), subtitle= "LogFC cutoff: 0.5, p cutoff: 0.05")
  
  return (all)
    }
  
  if(!(is.null(geneset))){
    #only graph genes in specified geneset
    DEgenes.gs = DEgenes %>% rownames_to_column("Gene") %>% filter(Gene %in% geneset) #filter out genes in geneset
    DEgenes.gs = column_to_rownames(DEgenes.gs, "Gene")
    just_gs = EnhancedVolcano(DEgenes.gs,lab= rownames(DEgenes.gs), x='avg_logFC', y='p_val_adj', 
        pCutoff = 0.05, FCcutoff = 0.5, col = c("black", "black", "black", "red2"), 
        title=paste(group.1, "vs", group.2, "Geneset DE Genes"), subtitle= "LogFC cutoff: 0.5, p cutoff: 0.05")
    
    return (just_gs)
    }

  }
  return (DEgenes)
}
