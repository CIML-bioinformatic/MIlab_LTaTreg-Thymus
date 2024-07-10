# #########################################
# This script make Harmony type integration 
# of the two datasets
# #########################################


## @knitr mergeWithHarmonyIntegration

# Find the genes in common between the two datasets
common_genes = intersect( row.names( owenetal_seurat_object@assays$RNA@counts), row.names( project_seurat_object@assays$RNA@counts))

# Create a Seurat Object with the merge of the two datasets counts tables
harmony_integration_seurat_object <- CreateSeuratObject( counts = cbind( owenetal_seurat_object@assays$RNA@counts[ common_genes, ],
                                                                         project_seurat_object@assays$RNA@counts[ common_genes, ]),
                                                         project = "Harmony_Integration", min.cells = 3)

# Apply a classic Seurat pipeline
harmony_integration_seurat_object = NormalizeData( harmony_integration_seurat_object, verbose = FALSE)
harmony_integration_seurat_object = FindVariableFeatures( harmony_integration_seurat_object, selection.method = "vst", nfeatures = 2000)
harmony_integration_seurat_object = ScaleData( harmony_integration_seurat_object, verbose = FALSE)
harmony_integration_seurat_object = RunPCA( harmony_integration_seurat_object, npcs = 30, verbose = FALSE)
harmony_integration_seurat_object = RunUMAP( harmony_integration_seurat_object, dims = 1:30, verbose = FALSE)

# Add a meta data with the dataset of origin
harmony_integration_seurat_object@meta.data$orig.project <- factor( c(rep("OwenEtAl", ncol( owenetal_seurat_object@assays$RNA@counts)),
                                                              rep("LTaTReg_Thy", ncol( project_seurat_object@assays$RNA@counts))))

Idents( owenetal_seurat_object) = "HTO_classification"
Idents( project_seurat_object) = "HTO_classification"

harmony_integration_seurat_object@meta.data$HTO_classification <- c( na.exclude( Idents( owenetal_seurat_object)[ Cells( harmony_integration_seurat_object)]),
                                                                     na.exclude( Idents( project_seurat_object)[ Cells( harmony_integration_seurat_object)]))


Idents( owenetal_seurat_object) = "CellType"
Idents( project_seurat_object) = "seurat_clusters"

harmony_integration_seurat_object@meta.data$clusters <- c( na.exclude( Idents( owenetal_seurat_object)[ Cells( harmony_integration_seurat_object)]),
                                                                     na.exclude( Idents( project_seurat_object)[ Cells( harmony_integration_seurat_object)]))


# Plot the UMAP of the merged data before Harmony
# DimPlot(object = harmony_integration_seurat_object, reduction = "pca", pt.size = .1, group.by = "orig.project")
# DimPlot(object = harmony_integration_seurat_object, reduction = "umap", pt.size = .1, group.by = "orig.project")
# VlnPlot(object = harmony_integration_seurat_object, features = "PC_1", group.by = "orig.project")

# Run the Harmony integration
harmony_integration_seurat_object = Run_Harmony( object = harmony_integration_seurat_object,
                                                  group.by.vars = "orig.project", assay.use = "RNA", plot_convergence = TRUE, max.iter.harmony = 20)
harmony_integration_seurat_object = RunUMAP( harmony_integration_seurat_object, reduction = "harmony", dims = 1:20) 


# Plot the data after Harmony integration
# DimPlot(object = harmony_integration_seurat_object, reduction = "harmony", pt.size = .1, group.by = "orig.project")
# VlnPlot(object = harmony_integration_seurat_object, features = "harmony_1", group.by = "orig.project")
DimPlot(harmony_integration_seurat_object, reduction = "umap", group.by = "orig.project", pt.size = .1, split.by = 'orig.project')
DimPlot(harmony_integration_seurat_object, reduction = "umap", group.by = "HTO_classification", pt.size = .1, split.by = 'orig.project')
DimPlot(harmony_integration_seurat_object, reduction = "umap", group.by = "clusters", pt.size = .1, split.by = 'orig.project', label = TRUE) + theme( legend.position = "None")

