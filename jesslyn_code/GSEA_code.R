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
library(DESeq2)
#devtools::install_github("immunogenomics/presto")
library(presto)
library(fgsea)

GSEA_code <- function(seurat_object, group.by, fgseaGS, condition, ranks = FALSE){
#step 1. run statistical test on all genes for each condition
  seurat.genes <- wilcoxauc(seurat_object, group.by) #across treatment statuses 
  
  #step 3. create a ranked vector for specified condition
  condition.genes <- seurat.genes %>% filter(group == condition) %>% arrange(desc(auc)) %>% select(feature, auc)
  condition.ranks <- deframe(condition.genes) #deframe converts two-column data frames to a named vector or list, using the first column as name and the second column as value
  
  #step 4. run fgsea
  condition.results<- fgsea(fgseaGS, stats = condition.ranks, nperm = 1000) %>% as_tibble() %>% select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(desc(NES), padj)
  
  if(ranks == FALSE){
  return(condition.results)}
  
  else{
    return(condition.ranks)
  }
}