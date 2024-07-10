# #########################################
# This script execute a simple dataset
# merge with no integration
# #########################################




# MERGE WITH NO INTEGRATION
# ----------------------------

## @knitr mergeWithoutIntegration

# Merge the dataset into a single Seurat Object
merge_seurat_object <- merge( NormalizeData( owenetal_seurat_object), y = NormalizeData( project_seurat_object), add.cell.ids = c("Owenetal", "RUN1-THY"), project = "MergedProject")

# Apply Seurat standard pipeline
merge_seurat_object <- ScaleData( merge_seurat_object, do.scale = FALSE, verbose = FALSE)
merge_seurat_object = FindVariableFeatures( merge_seurat_object, selection.method = "vst")
merge_seurat_object <- RunPCA( merge_seurat_object, npcs = 30, verbose = FALSE)
merge_seurat_object <- RunUMAP( merge_seurat_object, reduction = "pca", dims = 1:30)
# merge_seurat_object <- FindNeighbors( merge_seurat_object, reduction = "pca", dims = 1:30)
# merge_seurat_object <- FindClusters( merge_seurat_object, resolution = 0.5)

DimPlot( merge_seurat_object, label = TRUE) + theme( legend.position = "None") + ggtitle ("Merge with no integration")

# Clean the data to free memory
rm( "merge_seurat_object")

