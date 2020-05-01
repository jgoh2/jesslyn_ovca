# jesslyn_ovca
Analyzing gene expression data of cells in ascites samples of patients with ovarian cancer. Data provided by Izar et al.
Our goal is to detect how treatment affects the gene expression profile of malignant cells, and how the gene expression profile of 
cycling persister cells differ from those that are not resistant to chemotherapy treatment. 

## collaborators: Mike and Jesslyn 

### tasks
- [x] Load in the Seurat objects (10x, SmartSeq, PDX) 
  - [ ] Add relevant metadata slots 
    - [x] Cell type assignment for each cluster
    - [x] Treatment info assignment for each sample ID 
    - [x] Add cell cycle state (using gene list) 
    - [x] Score cells using gene list (describe biological functions)
    - [ ] Add enrichment
- [ ] Compute patient metrics 
  - [x] Number of patients in each cohort 
  - [x] Number of cells per patient 
  - [x] Number of cells per type per patient 
  - [x] Number of cells per treatment status 
  - [x] Number of patients per treatment status
  - [x] Number of cells per treatment status per cell type
  - [x] Number of malignant cells per patient per treatment status 
- [ ] Exploratory data analysis 
  - [x] Identify highly variable genes 
  - [x] Dimensional Reduction (Run and plot PCA) 
  - [x] Run Clustering 
  - [ ] Run and Plot UMAP 
  - [x] Find Marker Genes and create a Heatmap 
- [x] More specific Exploratory data analysis 
  - [x] Subset the Seurat object into only Malignant and only Nonmalignant cells 
  - [x] Perform the same analysis as above 
  
- [ ] Deeper Analysis 
  - [ ] Why do the cells cluster in such a way?
    - Try to color the cells by cell cycle phase 
    - Try to color the cells by treatment 
  - [ ] What are the identities/functions of each marker gene?
    - Are the computed marker genes typically associated with this specific cell type?
  - [ ] How does treatment affect the level of cycling malignant cells (after treatment should see less cycling)
    - [ ] Stack Bar Plot: show what fraction of cells is in S, G1, and G2M
      - By sample (treatment status) 
      - By cell type 


