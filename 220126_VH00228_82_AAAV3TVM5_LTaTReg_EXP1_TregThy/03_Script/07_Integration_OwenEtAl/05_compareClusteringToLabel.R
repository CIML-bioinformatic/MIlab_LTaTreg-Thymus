# #########################################################################
# This script aims to compare the clustering of cells in the
# No Contamination/No proliferation analysis to the cell labels
# from the Owen et al. data
# #########################################################################


## @knitr compareClusteringAndLabel

# Get the analysis performed previously with the removal of contaminating and proliferating cells
noprofliferation_nocontamination_seurat_object = readRDS( PATH_NOCONTAMINATION_NOPROLIFERATION_SEURAT_OBJECT_RDS)

# Transfer the inferred cell type labels from Owen et al integration to the loaded analysis
Idents( seurat_integration_seurat_object) = "CellType"
noprofliferation_nocontamination_seurat_object = AddMetaData( noprofliferation_nocontamination_seurat_object,
                                                              metadata = Idents( seurat_integration_seurat_object)[ Cells( noprofliferation_nocontamination_seurat_object)],
                                                              col.name = "OwenEtAl.infered.celltype")

# Look at the distribution of clusters and cell types
DimPlot( noprofliferation_nocontamination_seurat_object, group.by = "seurat_clusters", label = TRUE) + 
  theme( legend.position = "None")
DimPlot( noprofliferation_nocontamination_seurat_object, group.by = "OwenEtAl.infered.celltype", label = TRUE) + 
  theme( legend.position = "None")
DimPlot( noprofliferation_nocontamination_seurat_object, group.by = "OwenEtAl.infered.celltype", label = TRUE) + 
  facet_wrap( .~OwenEtAl.infered.celltype)  + 
  theme( legend.position = "None")

# Make table of the relation between Seurat Clusters and predicted cell types
as.data.frame.matrix( table( noprofliferation_nocontamination_seurat_object$OwenEtAl.infered.celltype, noprofliferation_nocontamination_seurat_object$seurat_clusters)) %>% kable( caption = "Predicted cell types against seurat clusters") %>% kable_styling()

# Make quantification of the relation between Seurats Clusters and predicted cell types
chi2_test = chisq.test( noprofliferation_nocontamination_seurat_object$OwenEtAl.infered.celltype, noprofliferation_nocontamination_seurat_object$seurat_clusters)

as.data.frame.matrix( chi2_test$residuals) %>% kable( caption = "Predicted cell types against Clusters Pearson residuals") %>% kable_styling()
