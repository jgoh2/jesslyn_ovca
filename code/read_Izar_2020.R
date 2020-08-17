# read_Izar_2020.R - read Izar_2020 data and save as Seurat objects in RDS files

# Load packages
source(here::here('packages.R'))
library(vroom)
library(tidyverse)
library(reshape2)
library(readxl)
library(Matrix)
library(beepr)

# read in files for scoring, etc -----------
hkgenes = read_lines("data/gene_lists/tirosh_house_keeping.txt", skip = 2)
pat_info = read_excel("data/Izar_2020/Table_S1_FINAL_Rev3.xlsx")
clust_cell_10x = read_excel("data/Izar_2020/Table_S2.xlsx", skip = 2)[1,] %>% 
  t() %>% as_tibble() %>% filter(V1 != "Cell type") %>% rownames_to_column()
names(clust_cell_10x) = c("clst","cell.type")
clust_cell_10x$clst = as.numeric(clust_cell_10x$clst)
clust_cell_SS2 = data.frame(clst = c(1,2,3,4,5,6,7,8,9),
                            cell.type = c(rep("Malignant",6),"Fibroblast","Macrophage","Malignant"))

#read in hallmarks and cell cycle genes for scoring
# Read in ccgenes
ccgenes = read_lines("data/gene_lists/regev_lab_cell_cycle_genes.txt")
s.genes <- ccgenes[1:43]
g2m.genes <- ccgenes[44:97]

#Read in hallmarks of interest
hallmark_names = read_lines("data/gene_lists/hallmarks.txt")
hallmark.list <- vector(mode = "list", length = length(hallmark_names))
names(hallmark.list) <- hallmark_names

for(hm in hallmark_names){
  if(file.exists(glue("data/gene_lists/hallmarks/{hm}_updated.txt"))){
    file <- read_lines(glue("data/gene_lists/hallmarks/{hm}_updated.txt"), skip = 1)
    hallmark.list[[hm]] <- file
  }
  else{
    file <- read_lines(glue("data/gene_lists/extra/{hm}.txt"), skip =2)
    hallmark.list[[hm]] <- file
  }
}

# 10x data --------------------
# read and save 10x data to Seurat object 
data_10x = vroom("data/Izar_2020/GSE146026_Izar_HGSOC_ascites_10x_log.tsv.gz")
data_10x = column_to_rownames(data_10x, "Cell_ID")
data_10x[1,] = substr(data_10x[1,1:ncol(data_10x)],13,nchar(data_10x[1,1:ncol(data_10x)]))
colnames(data_10x) = data_10x[1,]
data_10x = data_10x[-1,]
system.time(seurat_obj <- CreateSeuratObject(counts = data_10x, project = "Izar 2020: cohort 1 - 10x"))

# assign meta data
for (i in c("patient","time","sample-ID","clst","TSNE-x","TSNE-y")){
  meta = as.vector(seurat_obj[["RNA"]][i])
  names(meta) = Cells(seurat_obj)
  seurat_obj[[i]] = meta
  print(glue("loaded {i} into metadata"))
}


seurat_obj = subset(seurat_obj, features = rownames(seurat_obj)[-(1:6)]) # remove those same rows from matrix
## percent mito
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-") 
## HK genes
hkgenes.found = which(toupper(rownames(seurat_obj[["RNA"]]@data)) %in% hkgenes)
n.expressed.hkgenes <- seurat_obj[["RNA"]]@data[hkgenes.found, ] > 0 
seurat_obj <- AddMetaData(object = seurat_obj, metadata = Matrix::colSums(n.expressed.hkgenes), col.name = "n.exp.hkgenes")
## cell type
seurat_obj = AddMetaData(seurat_obj, left_join(seurat_obj@meta.data, clust_cell_10x)$cell.type, col.name = "cell.type")
# remove cells with no assignment
cells_keep = rownames(seurat_obj@meta.data)[!is.na(seurat_obj@meta.data$cell.type)]
seurat_obj = subset(seurat_obj, cells = cells_keep)

## patient treatment status
pnt = pat_info %>% 
  mutate(Patient = parse_number(`Patient ID`)) %>%
  select(`Patient`, `Treatment status of sample`) %>% 
  filter(`Patient`<=6)
pnt = rbind(pnt[1:2,], c(2.1, "On-treatment"), pnt[-(1:2),])
pnt = rbind(pnt[1:6,], c(5.1, "After 1 cycle of chemotherapy"), pnt[-(1:6),])
pnt[6,2] = "Treatment-naïve"
pnt = mutate(pnt, sample.ID = as.numeric(c("1976", "3250", "3250.1", "3266", "3281", "3288", "3288.1", "3290")))
seurat_obj = AddMetaData(seurat_obj, left_join(seurat_obj@meta.data, pnt)$`Treatment status of sample`, col.name = "treatment.status")
## better sample IDs
samp = data.frame(sample.ID = as.numeric(c("1976", "3250", "3250.1", "3266", "3281", "3288", "3288.1", "3290")),
                  new.sample.ID = c("1","2.0","2.1","3","4","5.0","5.1","6"))
seurat_obj = AddMetaData(seurat_obj, left_join(seurat_obj@meta.data, samp)$new.sample.ID, col.name = "sample.ID")
## Jesslyn AddModuleScore here!
seurat_obj <- AddModuleScore(seurat_obj, features = hallmark.list, name = names(hallmark.list), nbin = 50, search = T)
seurat_obj<- CellCycleScoring(seurat_obj, g2m.features = g2m.genes, s.features = s.genes)

# Save and clear memory
saveRDS(seurat_obj, "data/Izar_2020/Izar_2020_10x.RDS")
rm(seurat_obj)
rm(data_10x)

# SS2 data ----------------
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
## Percent mito
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
## HK genes
hkgenes.found = which(toupper(rownames(seurat_obj[["RNA"]]@data)) %in% hkgenes)
n.expressed.hkgenes <- seurat_obj[["RNA"]]@data[hkgenes.found, ] > 0 
seurat_obj <- AddMetaData(object = seurat_obj, metadata = Matrix::colSums(n.expressed.hkgenes), col.name = "n.exp.hkgenes")
## treatment status
pnt = pat_info %>% 
  mutate(Patient = parse_number(`Patient ID`)) %>%
  select(`Patient`, `Treatment status of sample`)
pnt[5,2] = "Treatment-naïve"
new_meta = left_join(seurat_obj@meta.data, pnt)
new_meta$`Treatment status of sample`[new_meta$Time == 0] = "Treatment-naïve"
seurat_obj = AddMetaData(seurat_obj, new_meta$`Treatment status of sample`, col.name = "treatment.status")
## cell type
seurat_obj = AddMetaData(seurat_obj, left_join(seurat_obj@meta.data, clust_cell_SS2)$cell.type, col.name = "cell.type")

##sample ID
sample_ID <- rep("NA", ncol(seurat_obj))
for(patient in unique(seurat_obj$Patient)){
  for (time in unique(seurat_obj$Time)){
    sample_ID[seurat_obj$Patient == patient & seurat_obj$Time == time] <- glue("{patient}.{time}")
  }
}
seurat_obj <- AddMetaData(seurat_obj, sample_ID, col.name = "sample")

## Jesslyn AddModuleScore here!
seurat_obj <- AddModuleScore(seurat_obj, features = hallmark.list, name = names(hallmark.list), nbin = 50, search = T)
seurat_obj<- CellCycleScoring(seurat_obj, g2m.features = g2m.genes, s.features = s.genes)

# Save and clear memory
saveRDS(seurat_obj, "data/Izar_2020/Izar_2020_SS2.RDS")
rm(seurat_obj)
rm(data_SS2)

# PDX data ------------
# read and save PDX data to Seurat object 
data_PDX = vroom("data/Izar_2020/PDX_OvCa_processed.txt.gz")
data_PDX = column_to_rownames(data_PDX, "...1")
pnt = list("MRD" = c(496,926,504,530,493,495),
           "relapse" = c(477,500,917,918,486,523,529,485),
           "vehicle" = c(494,933,505,534,476)) %>% reshape2::melt(value.name = "mouse_ID")
meta_df = colnames(data_PDX) %>% 
        strsplit(split = "_") %>% 
        unlist() %>%
        matrix(ncol = 3, byrow = T) %>%
        as.data.frame()
colnames(meta_df) = c("model_ID","mouse_ID","cell_ID")
meta_df$mouse_ID = gsub("Mice","",meta_df$mouse_ID) %>% as.numeric()
meta_df$model_ID = gsub("Model","",meta_df$model_ID)
meta_df = left_join(meta_df, pnt, by = "mouse_ID")
names(meta_df)[4] = "treatment.status"
colnames(data_PDX) = as.character(meta_df$cell_ID)
system.time(seurat_obj <- CreateSeuratObject(counts = data_PDX, project = "Izar 2020: PDX data"))

# assign meta data from rows of matrix to @metadata slot
for (i in c("mouse_ID","model_ID","treatment.status")){
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

## Jesslyn AddModuleScore here!
seurat_obj <- AddModuleScore(seurat_obj, features = hallmark.list, name = names(hallmark.list), nbin = 50, search = T)
seurat_obj<- CellCycleScoring(seurat_obj, g2m.features = g2m.genes, s.features = s.genes)

saveRDS(seurat_obj, "data/Izar_2020/Izar_2020_PDX.RDS")
rm(seurat_obj)
rm(data_PDX)

beep(4)
