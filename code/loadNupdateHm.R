# Loads in hallmarks from msigbr, filter out hallmarks of interest, update them, and save updated hallmarks in folder
# Date: 06/17/2020
# Author: Jesslyn Goh

#load in a list of hallmark names that we are interested in 
hallmark_names = read_lines("gene_lists/hallmarks.txt")
#load in all hallmarks, category H stands for hallmark geneset 
#filter out genesets we're interested in, and only select for columns with the geneset name and gene symbols
hallmark = msigdbr(species = "Homo sapiens", category = "H") %>% filter(gs_name %in% hallmark_names) %>% select(gs_name, gene_symbol) 
#convert hallmark into a list
hallmark.list = vector(mode = "list")
update = NULL
for(hm in unique(hallmark$gs_name)){
  update = select(filter(hallmark, gs_name %in% hm),gene_symbol) %>% as.list()
  update[["gene_symbol"]] <- update[["gene_symbol"]] %>% UpdateSymbolList(several.ok = TRUE)
  hallmark.list = c(hallmark.list, update)
}
names(hallmark.list) <- unique(hallmark$gs_name)

#save updated list of hallmark genes 
for(hm in unique(hallmark$gs_name)){
  write_tsv(as.data.frame(hallmark.list[[hm]]), path = glue("hallmarks/{hm}_updated.txt"))
}