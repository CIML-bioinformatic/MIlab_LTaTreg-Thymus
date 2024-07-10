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

# Prepare fixed colors for cluster
Idents( sc10x.seurat) = "seurat_clusters"
clustersColor = hue_pal()( nlevels( Idents( sc10x.seurat)))
names( clustersColor) = levels( Idents( sc10x.seurat))

### Remove eventual NULL (empty) list elements from list of genes to be monitored
monitoredGroupEmpty = sapply( MONITORED_GENES, is.null);
if(any( monitoredGroupEmpty)) warning( paste("Following group(s) of genes to be monitored will be ignored because empty:", paste( names(monitoredGroupEmpty)[monitoredGroupEmpty], collapse=" - ")));
MONITORED_GENES = MONITORED_GENES[! monitoredGroupEmpty];

# Check whether genes in MONITORED_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
matchMonitoredGenes = match( toupper( unlist( MONITORED_GENES)), toupper( rownames( GetAssayData( sc10x.seurat))));
matchMonitoredGenes = match( ( unlist( MONITORED_GENES)), ( rownames( GetAssayData( sc10x.seurat))));
monitoredGenesNotFound = unique( unlist( MONITORED_GENES)[is.na( matchMonitoredGenes)]);
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( paste0("'", monitoredGenesNotFound, "'"), collapse=" - ")));

# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MONITORED_GENES = relist( rownames( GetAssayData( sc10x.seurat))[ matchMonitoredGenes ], skeleton = MONITORED_GENES); # Does not work with NULL list elements (removed earlier)

# Finally remove names that did not match
MONITORED_GENES = lapply( MONITORED_GENES, na.omit)
