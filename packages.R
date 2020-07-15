# R packages --------------------------------------------------------------
# pacman and BiocManager for easy loading/installation
if (!require(renv)) install.packages("renv"); library(renv)
if (!require(BiocManager)) install.packages("BiocManager"); library(BiocManager)

renv::init()

# custom scripts
source(here::here('code','DotPlot2.R'))
source(here::here('code','DEAnalysis_code.R'))
source(here::here('code','GSEA_code.R'))
source(here::here('code','seurat_tab.R'))

