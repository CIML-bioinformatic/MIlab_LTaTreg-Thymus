# ################################################
# This script aims to read and filter data
# ################################################


## @knitr prepare_data

# READ THE DATA
# ...............

## @knitr prepare_data

# Load the Seurat object from heterogenity analysis
sc10x.seurat = readRDS( SEURAT_OBJECT_HETEROGENEITY_ANALYSIS)

# If required, keep only the cells from a single sample
if( SAMPLE_NAME != SAMPLE_ALL){
  Idents( sc10x.seurat) = "HTO_classification"
  sc10x.seurat = subset( sc10x.seurat, idents = SAMPLE_NAME)
}

# Convert Seurat Object to cell_data_set object
cds_sc10X_seurat = as.cell_data_set( sc10x.seurat)
