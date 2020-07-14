# R packages --------------------------------------------------------------
# pacman and BiocManager for easy loading/installation
if (!require(pacman)) install.packages("pacman"); library(pacman)
if (!require(BiocManager)) install.packages("BiocManager"); library(BiocManager)

# reproducibility tools
p_load(renv)
p_load(here)
p_load(rrtools)
p_load(workflowr)
p_load(styler)

# data analysis
p_load(tidyverse)
p_load(readxl)
p_load(plyr)
p_load(reshape2)
p_load(Matrix)
p_load(glue)
p_load(broom)
p_load(vroom)

# biology
p_load(msigdbr)
p_load(Seurat)
p_load(fgsea)

# visualization
p_load(gt)

# other
p_load(beepr)

# custom scripts
source(here::here('code','DotPlot2.R'))
source(here::here('code','DEAnalysis_code.R'))
source(here::here('code','GSEA_code.R'))
source(here::here('code','seurat_tab.R'))
