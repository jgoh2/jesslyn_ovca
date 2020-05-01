# read_ovca.R - Load ovca data and save in directory for Seurat reading
# Author: Mike Cuoco
# Date: 4/1/2020

# TO DO: make this into a function 

library(tidyverse)
library(Matrix)
library(vroom)
library(glue)

# Split data by 10x run
data_10x = vroom("/Volumes/jesslyn_ovca/data/Izar_2020/GSE146026_Izar_HGSOC_ascites_10x_log.tsv.gz")
data_10x = column_to_rownames(data_10x, "Cell_ID")

# Split data by 10x run
samples = data_10x[4,] %>% t() %>% as.vector()

for (i in unique(samples)) {
  
  # subset data frame
  sub = samples == i
  sub_data_10x = data_10x[,sub]
  
  dir.out = glue("~/jesslyn_ovca/data/Izar_2020GSE146026_Izar_HGSOC_ascites_10x_log/{i}")
  if (!dir.exists(dir.out)) { dir.create(dir.out, rec=T) }
  
  # save barcodes file
  substr(sub_data_10x[1,1:ncol(sub_data_10x)],13,nchar(sub_data_10x[1,1:ncol(sub_data_10x)])) %>%
    as.data.frame() %>%
    write_tsv(glue("/Volumes/jesslyn_ovca/data/GSE146026_Izar_HGSOC_ascites_10x_log/{i}/barcodes.tsv"), col_names = F)
  
  # save genes file
  cbind(rownames(sub_data_10x)[8:nrow(sub_data_10x)],rownames(sub_data_10x)[8:nrow(sub_data_10x)]) %>%
    as.data.frame() %>%
    write_tsv(glue("/Volumes/jesslyn_ovca/data/GSE146026_Izar_HGSOC_ascites_10x_log/{i}/genes.tsv"),col_names = F)
  
  # save matrix file
  as.matrix(sub_data_10x[8:nrow(sub_data_10x),1:ncol(sub_data_10x)]) %>% 
    Matrix(sparse = T) %>%
    writeMM(glue("/Volumes/jesslyn_ovca/data/GSE146026_Izar_HGSOC_ascites_10x_log/{i}/matrix.mtx"))
}

