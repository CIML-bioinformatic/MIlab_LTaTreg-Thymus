# ##################################################
# This script reads data from previous analysis
# ##################################################




# READ DATA
# ---------s

## @knitr loadData



# Get the Seurat object and the expression count table
sc10x.seurat = readRDS( file.path( PATH_GLOBALHETEROGENEITY[[ SAMPLE_NAME]], FILENAME_SEURAT_OBJECT))

# Marker genes per sample
# cluster_marker_genes_list = list()
# for( sample in SAMPLE_SET){
#   cluster_marker_genes_list[[ sample]] = read.table( file.path( PATH_GLOBALHETEROGENEITY[[ sample]], FILENAME_CLUSTER_MARKERGENES),
#                                                      sep ="\t", header = TRUE, row.names = 1)
# }


# Look at the number of cells per samples and cluster
Idents( sc10x.seurat) = "seurat_clusters"
cells_clusterid = Idents( sc10x.seurat)
Idents( sc10x.seurat) = "HTO_classification"
cells_sample = Idents( sc10x.seurat)
pander( table( cells_clusterid, cells_sample))
