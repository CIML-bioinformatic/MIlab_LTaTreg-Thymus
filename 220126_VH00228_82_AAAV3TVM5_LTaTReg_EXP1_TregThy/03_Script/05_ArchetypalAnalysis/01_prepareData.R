# ##################################################
# Load the Seurat object from previous analysis step
# and prepare data for the scAANet analysis
# ##################################################


## @knitr prepare_data

sc10x.seurat = readRDS( PATH_SEURAT_ROBJECT_FILE)
sc10x_counts_df = as.data.frame( t( as.matrix( sc10x.seurat@assays$RNA@counts)))
