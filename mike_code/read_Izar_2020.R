# read_Izar_2020.R - read Izar_2020 data and save as Seurat objects in RDS files

# Load packages
library(tidyverse)
library(Matrix)
library(Seurat)
library(vroom)
library(glue)
library(beepr)

# read in housekeeping gene list
hkgenes = read_lines("gene_lists/tirosh_house_keeping.txt", skip = 2)

# read and save 10x data to Seurat object 
data_10x = vroom("data/Izar_2020/GSE146026_Izar_HGSOC_ascites_10x_log.tsv.gz")
data_10x = column_to_rownames(data_10x, "Cell_ID")
data_10x[1,] = substr(data_10x[1,1:ncol(data_10x)],13,nchar(data_10x[1,1:ncol(data_10x)]))
colnames(data_10x) = data_10x[1,]
data_10x = data_10x[-1,]
system.time(seurat_obj <- CreateSeuratObject(counts = data_10x, project = "Izar 2020: cohort 1 - 10x"))

# assign meta data from rows of matrix to @metadata slot
for (i in c("patient","time","sample-ID","clst","TSNE-x","TSNE-y")){
  meta = as.vector(seurat_obj[["RNA"]][i])
  names(meta) = Cells(seurat_obj)
  seurat_obj[[i]] = meta
  print(glue("loaded {i} into metadata"))
}

# remove those same rows from matrix
seurat_obj = subset(seurat_obj, features = rownames(seurat_obj)[-(1:6)])
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
hkgenes.found = which(toupper(rownames(seurat_obj[["RNA"]]@data)) %in% hkgenes)
n.expressed.hkgenes <- seurat_obj[["RNA"]]@data[hkgenes.found, ] > 0 
seurat_obj <- AddMetaData(object = seurat_obj, metadata = Matrix::colSums(n.expressed.hkgenes), col.name = "n.exp.hkgenes")
saveRDS(seurat_obj, "data/Izar_2020/Izar_2020_10x.RDS")

# clear memory
rm(seurat_obj)
rm(data_10x)

# read and save SS2 data to Seurat object 
data_SS2 = vroom("data/Izar_2020/GSE146026_Izar_HGSOC_ascites_SS2_log.tsv.gz")
data_SS2 = column_to_rownames(data_SS2, "Cell_ID")
system.time(seurat_obj <- CreateSeuratObject(counts = data_SS2, project = "Izar 2020: cohort 2 - SS2"))

# assign meta data from rows of matrix to @metadata slot
for (i in c("tSNE1","tSNE2","Patient","Time","clst")){
  meta = as.vector(seurat_obj[["RNA"]][i])
  names(meta) = Cells(seurat_obj)
  seurat_obj[[i]] = meta
  print(glue("loaded {i} into metadata"))
}

seurat_obj = subset(seurat_obj, features = rownames(seurat_obj)[-(1:4)])
hkgenes.found = which(toupper(rownames(seurat_obj[["RNA"]]@data)) %in% hkgenes)
n.expressed.hkgenes <- seurat_obj[["RNA"]]@data[hkgenes.found, ] > 0 
seurat_obj <- AddMetaData(object = seurat_obj, metadata = Matrix::colSums(n.expressed.hkgenes), col.name = "n.exp.hkgenes")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")

saveRDS(seurat_obj, "data/Izar_2020/Izar_2020_SS2.RDS")

# read and save SS2 data to Seurat object 
data_PDX = vroom("data/Izar_2020/PDX_OvCa_processed.txt.gz")
data_PDX = column_to_rownames(data_PDX, "...1")
meta_df = colnames(data_PDX) %>% 
        strsplit(split = "_") %>% 
        unlist() %>%
        matrix(ncol = 3, byrow = T) %>%
        as.data.frame()
colnames(meta_df) = c("model_ID","mouse_ID","cell_ID")
colnames(data_PDX) = as.character(meta_df$cell_ID)
system.time(seurat_obj <- CreateSeuratObject(counts = data_PDX, project = "Izar 2020: PDX data"))

# assign meta data from rows of matrix to @metadata slot
for (i in c("mouse_ID","model_ID")){
  meta = meta_df[[i]]
  names(meta) = Cells(seurat_obj)
  seurat_obj[[i]] = meta
  print(glue("loaded {i} into metadata"))
}

seurat_obj = subset(seurat_obj, features = rownames(seurat_obj))
hkgenes.found = which(toupper(rownames(seurat_obj[["RNA"]]@data)) %in% hkgenes)
n.expressed.hkgenes <- seurat_obj[["RNA"]]@data[hkgenes.found, ] > 0 
seurat_obj <- AddMetaData(object = seurat_obj, metadata = Matrix::colSums(n.expressed.hkgenes), col.name = "n.exp.hkgenes")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")

saveRDS(seurat_obj, "data/Izar_2020/Izar_2020_PDX.RDS")

beep(4)
