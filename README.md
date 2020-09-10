# CHARACTERIZING METABOLIC DYNAMICS OF OVARIAN CANCER CELLS FOLLOWING CARBOPLATIN TREATMENT
Developed by [Jesslyn Goh](https://github.com/jgoh2) and [Mike Cuoco](https://github.com/mikecuoco)

In this meta-analysis, we aim to identify gene expression signatures associated with  **drug tolerant ovarian cancer persister cells**. (can add a brief background sentence about persister cells here). The mechanistic characterization of persister cell survival and proliferation will expand our understanding of cancer's ability to resist traditional cytotoxic agents and inform the development of new drugs to prevent disease relapse. 

Our single-cell RNA-seq dataset is kindly provided by the authors ["A single-cell landscape of high-grade serous ovarian cancer." (Izar et al. 2020)](https://doi.org/10.1038/s41591-020-0926-0). 

For more background information about ovarian cancer and drug tolerant persister cells and our ovarian cancer dataset , please visit [our project website](https://jgoh2.github.io/jesslyn_ovca/).

## SPECIFIC AIMS
1. Identify candidate drug tolerant cancer persister cells among available single-cell RNA-seq datasets.

2. Characterize the gene expression signatures unique to the identified persister populations.



## OUR ANALYSIS WORKFLOW 
The code for our scRNA-seq analysis can be found in the `analysis/` folder and helper functions that are used in our analysis can be found in the `code/` folder. 

We split our scRNA-seq analysis into five parts: 

  1. **Load Data and Create Seurat Objects:** `code/read_Izar_2020.R` 
  
         a) Load in the raw count matrix and Create Seurat Object
         b) Assign Metadata identities including: 
              * Patient ID
              * Time
              * Sample ID
              * Treatment Status 
         c) Score cells for cell cycle and hallmark genesets 
         d) Save Seurat Object
          
  2. **Process Data:** `old/edited/02_Izar2020_SS2_Load_Plots.Rmd` and `old/edited/03_Izar2020_PDX_Load.Rmd`
  
         a) Load in SS2 Seurat Object from Part 1 and subset by patient/model. Continue analysis separately for each patient/model. 
         b) Scale and center data, and FindVariableFeatures (prepare data for dimensionality reduction)
         c) Dimensionality Reduction (PCA + UMAP)
         d) Save Seurat Objects 
      
  3. **Exploratory Data Analysis** `analysis/02.0_Izar2020_SS2_ExploratoryAnalysis.Rmd` and `analysis/03.0_Izar2020_PDX_ExploratoryAnalysis.Rmd`
  
         a) Load in Seurat Objects from Part 2. Analyze separately for each patient/model
         b) Compute summary metrics for SS2 data such as: 
              * Number of cells per patient per treatment 
              * Number of cells per treatment per cell cycle phase 
         c) Visualize how cells separate based on metadata identities via UMAP and PCA
              * Intermodel heterogeneity: How do cells separate by sample? 
              * Intramodel heterogeneity: 
                  - How do cells separate by treament status / sample? 
                  - How do cells separate by cell cycle phase?
            
  4. **DE Analysis** `analysis/02.1_Izar2020_SS2_DEAnalysis.Rmd` and `analysis/03.1_Izar2020_PDX_DEAnalysis.Rmd`
    
         a) TYPE #1 DE ANALYSIS: Visualizing and Quantifying Differentially Expression on Predefined Gene Ontology Genesets
         b) TYPE #2 DE ANALYSIS: Identifying Individual DE Genes 
         c) TYPE #3 DE ANALYSIS: Gene Ontology Functional Enrichment Analysis 
          
  5. **CELL CYCLE ANALYSIS** `analysis/02.2_Izar2020_SS2_CellCycleAnalysis.Rmd` and `analysis/03.2_Izar2020_PDX_CellCycleAnalysis.Rmd`
    
         a) Examine whether there is a correlation between treatment condition and cell cycle phase through DE Analysis of cell cycle gene-set scores
         b) Examine a correlation between cell cycle phase and the expression of our gene-sets of interest
        
Please find our analysis workflow, our code, and results on [our project website](https://jgoh2.github.io/jesslyn_ovca/).    
