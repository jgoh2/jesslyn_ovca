# Returns fgsea result for a specific geneset and treatment condition
# Date: 05/25/2020
# Author: Jesslyn Goh

library(tidyverse) # ggplot, dplyr, read_r (read_xlsx)
library(Seurat)
library(readxl)
library(ggplot2)
library(msigdbr)
library(here)
library(dplyr)
library(devtools)
#library(DESeq2)
#devtools::install_github("immunogenomics/presto")
library(presto)
library(fgsea)

GSEA_code <- function(seurat_object, group.by, fgseaGS, condition = NULL, ranks = NULL, paired = FALSE, group.1 = NULL, group.2 = NULL){

  if (is.null(ranks)){
  #if didn't input ranks, we need to compute the ranks first
    if(paired == FALSE){
      #step 1. run statistical test on all genes for each condition
      seurat.genes <- wilcoxauc(seurat_object, group.by) #across all treatment statuses 
      condition.genes <- seurat.genes %>% filter(group == condition) %>% arrange(desc(logFC)) %>% select(feature, logFC)
      #filter for a specific condition of interest
    }
  
    else{
      Idents(seurat_object) <- group.by
      seurat.genes <- FindMarkers(seurat_object, ident.1 = group.1 , ident.2 = group.2, logfc.threshold = 0)
      #returns a single result that encompasses both conditions (group1 vs. group2), so the 
      #enrichment plot will be for group.1 vs. group2 instead of only for group1 relative to the other groups
      condition.genes <- seurat.genes %>% rownames_to_column %>% arrange(desc(avg_logFC)) %>% select(rowname, avg_logFC) 
    }
  
    #step 3. create a ranked vector for specified condition using logFC
    condition.ranks <- deframe(condition.genes) #deframe converts two-column data frames to a named vector or list, using the first column as name and the second column as value
    
    return(condition.ranks)
  }
  
  else{ #if inputted ranks, use the inputted ranks to calculate gsea results
    condition.ranks <- ranks
    #step 4. run fgsea
    condition.results<- fgsea(fgseaGS, stats = condition.ranks, minSize = 15, maxSize = 500, nperm = 1000) %>% as_tibble() %>% select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(desc(NES), padj)
    
    return(condition.results)
  }

}