# Returns ranked differentially expressed genes between two conditions from FindMarkers.
# Date: 07/14/2020
# Author: Jesslyn Goh

library(tidyverse) # ggplot, dplyr, read_r (read_xlsx)
library(Seurat)
library(here)
library(dplyr)
library(GGally)
library(ggpubr)

# 1) Finds DE genes betwen two conditions 
# 2) Arrange them by padj ascending, then by logFC (or avg_diff) descending. Also computes negative log10 padj 
# 3) Rank the genes (1 being the top of the list)
# 4) Add columns that indicate whether each gene is considered DE (padj < 0.05, abs(logFC) > 0.5)
# 5) Returns the output

FindMarkers_code <- function(seurat_object, group.by, group.1 = NULL, group.2 = NULL, stattest = "wilcox", slot = "data", rank.by = "padj"){
  
  markers <- FindMarkers(seurat_object, group.by = group.by, ident.1 = group.1, ident.2 = group.2, test.use = stattest, slot = slot, logfc.threshold = 0)
  
  if (slot == "data"){
    if(rank.by == "padj"){
    markers <- markers %>% rownames_to_column %>% arrange(p_val_adj, -avg_logFC) %>% select(rowname, avg_logFC, p_val_adj)
    }
    else{
      #rank by logFC
      markers <- markers %>% rownames_to_column %>% arrange(-avg_logFC, p_val_adj) %>% select(rowname, avg_logFC, p_val_adj)
    }
    markers <- markers %>% 
      mutate("data.-log10padj" = -(log10(markers$p_val_adj))) %>% 
      mutate("data.rank" = seq(from=1, to=length(markers$rowname))) %>% 
      mutate("data.DE" = (markers$avg_logFC > 0.5 | markers$avg_logFC < -0.5) & (markers$p_val_adj < 0.05))
  }
  else{
    markers <- markers %>% rownames_to_column %>% arrange(p_val_adj, -avg_diff) %>% select(rowname, avg_diff, p_val_adj) 
    markers <- markers %>% 
      mutate("scale.-log10padj" = -(log10(markers$p_val_adj))) %>% 
      mutate("scale.data.rank" = seq(from=1, to=length(markers$rowname))) %>% 
      mutate("scale.data.DE" = (markers$avg_diff > 0.5 | markers$avg_diff < -0.5) & (markers$p_val_adj < 0.05))
  }
  return (markers)
}


