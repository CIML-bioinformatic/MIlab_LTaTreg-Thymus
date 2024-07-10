# #########################################
# This script execute a Seurat integration
# of the two datasets
# #########################################


## @knitr mergeWithSeuratIntegration

## Make Seurat type integration of the two datasets
## ................................................
options( future.globals.maxSize = 8000 * 1024^2)

cat("\n \n")
cat("#### Integration of the datasets")
cat("\n \n")

# Create a list with the two Seurat objects of the datasets
seurat_object_list = list( owenetal_seurat_object, project_seurat_object)

# Get the features and anchors to produce the integration
features <- SelectIntegrationFeatures( object.list = seurat_object_list)
datasets.anchors <- FindIntegrationAnchors( object.list = seurat_object_list, anchor.features = features, reference = 1)

# Create the integration Seurat object
seurat_integration_seurat_object <- IntegrateData( anchorset = datasets.anchors)

# Specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay( seurat_integration_seurat_object) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat_integration_seurat_object <- ScaleData( seurat_integration_seurat_object, do.scale = FALSE, verbose = FALSE)
seurat_integration_seurat_object <- RunPCA( seurat_integration_seurat_object, npcs = 30, verbose = FALSE)
seurat_integration_seurat_object <- RunUMAP( seurat_integration_seurat_object, reduction = "pca", dims = 1:30)
# seurat_integration_seurat_object <- FindNeighbors( seurat_integration_seurat_object, reduction = "pca", dims = 1:30)
# seurat_integration_seurat_object <- FindClusters( seurat_integration_seurat_object, resolution = 0.5)

# Change in metadata the name of the origin ID to match suitable names
origin_meta_data = seurat_integration_seurat_object@meta.data$orig.ident
origin_meta_data = gsub( "GEX", "OwenEtAl", origin_meta_data, fixed = TRUE)
origin_meta_data = gsub( "RUN1", "LTaTReg_Thy", origin_meta_data, fixed = TRUE)
seurat_integration_seurat_object = AddMetaData( seurat_integration_seurat_object, metadata = origin_meta_data, col.name = "orig.project")

# Plot the integrated object through different views
options( repr.plot.height = 6, repr.plot.width = 9)
DimPlot( seurat_integration_seurat_object, label = TRUE, split.by = "orig.project",) + theme( legend.position = "None") + ggtitle ("Merge with Seurat integration (by origin)")
options( repr.plot.height = 6, repr.plot.width = 6)
DimPlot( seurat_integration_seurat_object, label = TRUE) + theme( legend.position = "None") + ggtitle ("Merge with Seurat integration (all)")

Idents( seurat_integration_seurat_object) = "hash.ID"
for( current_sample_type in levels( Idents( seurat_integration_seurat_object))){
  cells_to_plot = names( Idents( seurat_integration_seurat_object))[ which( Idents(seurat_integration_seurat_object) == current_sample_type)]
  print( DimPlot( subset( seurat_integration_seurat_object, cells = cells_to_plot), group.by = "CellType", label = TRUE, repel = TRUE)
         + ggtitle( paste( "Cell type dispersion in", current_sample_type)) + theme( legend.position = "None")
        )
}


## Prediction of cell types on LTA Treg Thymus dataset from Owen et al dataset cell type information
## ....................................................................................................

cat("\n \n")
cat("#### Prediction of cell types on LTA Treg Thymus dataset from Owen et al dataset cell type information")
cat("\n \n")

# Transfert the Cell Type labels from Owne et al data to LTA Treg Thymus data
transfert.anchors <- FindTransferAnchors(reference = owenetal_seurat_object, query = project_seurat_object, dims = 1:30)
predictions <- TransferData( anchorset = transfert.anchors, refdata = owenetal_seurat_object$CellType, dims = 1:30)
project_seurat_object <- AddMetaData( project_seurat_object, metadata = predictions$predicted.id, col.name = "CellType")

# Add new predicted cell type also to integrated Seurat Object
Idents( owenetal_seurat_object) = "CellType"
Idents( project_seurat_object) = "CellType"
seurat_integration_seurat_object@meta.data$CellType <- c( Idents( owenetal_seurat_object)[ intersect( Cells( seurat_integration_seurat_object), Cells( owenetal_seurat_object))],
                                                           Idents( project_seurat_object)[ intersect( Cells( seurat_integration_seurat_object), Cells( project_seurat_object))])

# Plot the result in the integrated dataset
options( repr.plot.height = 6, repr.plot.width = 9)
DimPlot( seurat_integration_seurat_object, label = TRUE, group.by = "CellType", split.by = "orig.project") + theme( legend.position = "None") + ggtitle ("Merge with Seurat integration\n(split by origin, grouped by Cell Types)")

# Plot the predicted cell types to the original project UMAP 
options( repr.plot.height = 6, repr.plot.width = 6)
DimPlot( project_seurat_object, label = TRUE, group.by = "seurat_clusters") + theme( legend.position = "None") + ggtitle ("LTa Treg Thymus UMAP with Seurat Clusters")
DimPlot( project_seurat_object, label = TRUE, group.by = "CellType") + theme( legend.position = "None") + ggtitle ("LTa Treg Thymus UMAP with Predicted Cell Types")

# Make quantification of the relation between Seurats Clusters and predicted cell types
as.data.frame.matrix( table( project_seurat_object$CellType, project_seurat_object$seurat_clusters)) %>% kable( caption = "Predicted cell types against seurat clusters") %>% kable_styling()

# Make quantification of the relation between Condition and predicted cell types
DimPlot( project_seurat_object, label = TRUE, group.by = "CellType", split.by = "HTO_classification") + theme( legend.position = "None") + ggtitle ("LTa Treg Thymus UMAP with Predicted Cell Types split by Condition")
chi2_test = chisq.test( project_seurat_object$CellType, project_seurat_object$HTO_classification)

as.data.frame.matrix( chi2_test$observed) %>% kable( caption = "Predicted cell types against Condition") %>% kable_styling()
as.data.frame.matrix( chi2_test$residuals) %>% kable( caption = "Predicted cell types against Condition Pearson residuals") %>% kable_styling()

# # remove integrated object to spare memory
# rm( "seurat_integration_seurat_object")