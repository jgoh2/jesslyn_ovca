# CHARACTERIZATION OF METABOLIC DYNAMICS OF OVARIAN PERSISTER CELLS FOLLOWING CARBOPLATIN TREATMENT
Developed by [Jesslyn Goh](https://github.com/jgoh2) and [Mike Cuoco](https://github.com/mikecuoco)

We characterize the gene expression signatures of **High-grade Serous Ovarian Cancer cells** across different treatment conditions. 
We aim to identify differential gene or gene-set expression that are associated with the drug tolerant cancer persister cells. The identification of underlying biological mechanisms 
that enable persister cells to survive and proliferate following drug treatment will inform future work in the creation of new drugs that specifically eradicate 
these persisters cells and prevent cancer relapse. For more background information about Ovarian Cancer and Non-genetic drug tolerant persister cells, 
please visit [our project website](https://jgoh2.github.io/jesslyn_ovca/).

Our single-cell Ovarian cancer data is kindly provided by Izar and Tirosh et al., from their paper A single-cell landscape of high-grade serous ovarian cancer. 
For more background information about our Ovarian Cancer data, please visit [our project website](https://jgoh2.github.io/jesslyn_ovca/).

## OUR ANALYSIS WORKFLOW 
Our code for our scRNA-seq analysis can be found in the **analysis** folder and helper functions that are used in our analysis can be found in the **code** folder. 

We split our scRNA-seq analysis into five parts: 
  1. **Load Data and Create Seurat Objects**
  
    a. The code to this part of our analysis is in the [read_Izar_2020.R] file in the [code] folder. During this part of our analysis we: 
            i. Load in the raw count matrix and Create Seurat Object 
            ii. Assign Metadata identities including: 
                - Patient ID
                - Time
                - Sample ID
                - Treatment Status 
            iii. Score cells for **cell cycle** and **hallmark genesets** 
            iv. Save Seurat Object
          
  2. **Process Data**
  
    a. The code to this part of our analysis can be found in the [02_Izar2020_SS2_Load_Plots.Rmd] file for the SS2 data and in the 
      [03_Izar2020_PDX_Load.Rmd] file in the [old/edited] folder for the PDX data. During this part of our analysis we: 
            i. Load in SS2 Seurat Object from Part 1 and subset by patient/model. Continue analysis separately for each patient/model. 
            ii. Scale and center data, and FindVariableFeatures (prepare data for dimensionality reduction)
            iii. Dimensionality Reduction (PCA + UMAP)
            vi. Save Seurat Objects 
      
  3. **Exploratory Data Analysis**
  
    a. The code to this part of our analysis can be found in the [02.0_Izar2020_SS2_ExploratoryAnalysis.Rmd] and [03.0_Izar2020_PDX_ExploratoryAnalysis.Rmd] files in the [analysis] folder. During this part of our analysis we: 
            i. Load in Seurat Objects from Part 2. Analyze separately for each patient/model
            ii. Compute summary metrics for SS2 data such as: 
                - Number of cells per patient per treatment 
                - Number of cells per treatment per cell cycle phase 
            iii. Visualize how cells separate based on metadata identities via UMAP and PCA
                - Intermodel heterogeneity: How do cells separate by sample? 
                - Intramodel heterogeneity: 
                    * How do cells separate by treament status / sample? 
                    * How do cells separate by cell cycle phase?
            
  4. **DE Analysis**
    
    a. The code to this part of our analysis can be found in the [02.1_Izar2020_SS2_DEAnalysis.Rmd] and [03.1_Izar2020_PDX_DEAnalysis.Rmd] files in the [analysis] folder. 
        i. TYPE #1 DE ANALYSIS: Visualizing and Quantifying Differentially Expression on Predefined Gene Ontology Genesets
        ii. TYPE #2 DE ANALYSIS: Identifying Individual DE Genes 
        iii. TYPE #3 DE ANALYSIS: Gene Ontology Functional Enrichment Analysis 
          
  5. **CELL CYCLE ANALYSIS**
    
    a. The code to this part of our analysis can be found in the [02.2_Izar2020_SS2_CellCycleAnalysis.Rmd] and [03.2_Izar2020_PDX_CellCycleAnalysis.Rmd] files in the [analysis] folder. 
        i. Examine whether there is a correlation between treatment condition and cell cycle phase through DE Analysis of cell cycle gene-set scores
        ii. Examine a correlation between cell cycle phase and the expression of our gene-sets of interest
        
Please find our analysis workflow, our code, and results on [our project website](https://jgoh2.github.io/jesslyn_ovca/).    
