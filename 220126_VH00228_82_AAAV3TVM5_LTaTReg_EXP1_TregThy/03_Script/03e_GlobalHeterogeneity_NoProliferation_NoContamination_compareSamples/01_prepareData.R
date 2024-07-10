# ##################################################
# This script reads data from previous analysis
# ##################################################




# READ DATA
# ---------

## @knitr loadData



# Get the Seurat object and the expression count table
sc10x.seurat = readRDS( file.path( PATH_GLOBALHETEROGENEITY[[ SAMPLE_NAME]], FILENAME_SEURAT_OBJECT))

# Look at the number of cells per samples and cluster
Idents( sc10x.seurat) = "seurat_clusters"
cells_clusterid = Idents( sc10x.seurat)
Idents( sc10x.seurat) = "HTO_classification"
cells_sample = Idents( sc10x.seurat)
pander( table( cells_clusterid, cells_sample))

# Plot the UMAP of cells as reference
DimPlot( sc10x.seurat, group.by = "seurat_clusters") + ggtitle( "UMAP of cells with clusters")

