# #########################################
# This script reads and filters sc10x  data
# #########################################




# READ DATA
# ---------

## @knitr loadData

# Load current project data from previous setp
project_seurat_object = readRDS( PATH_SEURAT_OBJECT_PROJECT)

# Load the Seurat object from Owen et al. 2022 study (PMID: 36038290)
owenetal_seurat_object = readRDS( PATH_SEURAT_OBJECT_OWENETAL)
# Save the Active Idents of the object to a meta.data
owenetal_seurat_object = AddMetaData( owenetal_seurat_object, metadata = Idents( owenetal_seurat_object), col.name = "CellType")

# Plot the two dataset in their original UMAP
DimPlot( owenetal_seurat_object) + ggtitle( "Owen et al dataset")
DimPlot( project_seurat_object) + ggtitle( "LTA Treg Thymus project")

# Compute the marker genes between the cell types defined by Owen et al.
owenetal_all_markers = FindAllMarkers( object          = owenetal_seurat_object,
                                      test.use        = FINDMARKERS_METHOD,
                                      only.pos        = FINDMARKERS_ONLYPOS,
                                      min.pct         = FINDMARKERS_MINPCT,
                                      logfc.threshold = FINDMARKERS_LOGFC_THR,
                                      return.thresh   = FINDMARKERS_PVAL_THR,
                                      random.seed     = SEED,
                                      verbose         = .VERBOSE);

# Keep only the most significant markers
topMarkers = by( owenetal_all_markers, owenetal_all_markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( x);
});

# Merge marker genes in a single data.frame
topMarkersDF = do.call( rbind, topMarkers);

# Save markers list as 'tsv' table
write.table( topMarkersDF,
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "OwenEtAL_MarkerGenes.tsv")),
             quote = FALSE,
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t")

