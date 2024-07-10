# ##################################################
# This script reads data from previous analysis
# ##################################################




# READ DATA
# ---------s

## @knitr loadData

# Get the Seurat object and the expression count table
sc10x.all = readRDS( file.path( PATH_GLOBALHETEROGENEITY[[ SAMPLE_NAME]], FILENAME_SEURAT_OBJECT))

# Select the cells of the cluter group to study
Idents( sc10x.all) = "seurat_clusters"
sc10x = subset( sc10x.all, idents = CLUSTER_GROUP_LIST[[ CLUSTER_GROUP_ID]])

# Keep trace of the original cluster annotation
sc10x@meta.data$seurat_clusters = droplevels( sc10x@meta.data$seurat_clusters)
sc10x = AddMetaData( sc10x, metadata = sc10x@meta.data$seurat_clusters, col.name = "original_seurat_clusters")

# Look at the number of cells per samples and cluster
Idents( sc10x) = "seurat_clusters"
cells_clusterid = Idents( sc10x)
Idents( sc10x) = "HTO_classification"
cells_sample = Idents( sc10x)
pander( table( cells_clusterid, cells_sample))


### Remove eventual NULL (empty) list elements from list of genes to be monitored
# monitoredGroupEmpty = sapply( MONITORED_GENES, is.null);
# if(any( monitoredGroupEmpty)) warning( paste("Following group(s) of genes to be monitored will be ignored because empty:", paste( names(monitoredGroupEmpty)[monitoredGroupEmpty], collapse=" - ")));
# MONITORED_GENES = MONITORED_GENES[! monitoredGroupEmpty];

# Check whether genes in MONITORED_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
#matchMonitoredGenes = match( toupper( unlist( MONITORED_GENES)), toupper( rownames( GetAssayData( sc10x))));
# matchMonitoredGenes = match( ( unlist( MONITORED_GENES)), ( rownames( GetAssayData( sc10x))));
# monitoredGenesNotFound = unique( unlist( MONITORED_GENES)[is.na( matchMonitoredGenes)]);
# if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( paste0("'", monitoredGenesNotFound, "'"), collapse=" - ")));
# # Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
# MONITORED_GENES = relist( rownames( GetAssayData( sc10x))[ matchMonitoredGenes ], skeleton = MONITORED_GENES); # Does not work with NULL list elements (removed earlier)
# # Finally remove names that did not match
# MONITORED_GENES = lapply( MONITORED_GENES, na.omit);


### Remove eventual NULL (empty) list elements from list of genes in modules
modulesGroupEmpty = sapply( MODULES_GENES, is.null);
if(any( modulesGroupEmpty)) warning( paste("Following module(s) of genes will be ignored because empty:", paste( names(modulesGroupEmpty)[modulesGroupEmpty], collapse=" - ")));
MODULES_GENES = MODULES_GENES[! modulesGroupEmpty];

# Check whether genes in MODULES_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
#matchModulesGenes = match( toupper( unlist( MODULES_GENES)), toupper( rownames( GetAssayData( sc10x))));
matchModulesGenes = match( ( unlist( MODULES_GENES)), ( rownames( GetAssayData( sc10x))));
modulesGenesNotFound = unique( unlist( MODULES_GENES)[is.na( matchModulesGenes)]);
if(any( is.na( matchModulesGenes))) warning( paste( "Following gene(s) from modules list could not be found in experimental data and will be ignored:", paste( paste0("'", modulesGenesNotFound, "'"), collapse=" - ")));
# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MODULES_GENES = relist( rownames( GetAssayData( sc10x))[ matchModulesGenes ], skeleton = MODULES_GENES); # Does not work with NULL list elements (removed earlier)
# Finally remove names that did not match
MODULES_GENES = lapply( MODULES_GENES, na.omit);


### Transfer genes in very small modules (MODULES_GENES) to be analyzed individually (MONITORED_GENES)
# modulesToTransfer = sapply(MODULES_GENES, length) < MONITORED_GENES_SMALL_MODULE
# if(any(modulesToTransfer))
# {
#   warning( paste0( "Following Module(s) contained very few genes (<", MONITORED_GENES_SMALL_MODULE, "). These genes were transfered to 'Monitored genes' to be analyzed individually: ", paste( names(modulesToTransfer)[modulesToTransfer], collapse = " - ")));
#   # Alter name so they can be recognized in list of Monitored genes
#   names( MODULES_GENES)[modulesToTransfer] = paste0( "MOD_", names( MODULES_GENES)[modulesToTransfer]);
#   MONITORED_GENES = c( MONITORED_GENES, MODULES_GENES[modulesToTransfer]);
#   MODULES_GENES = MODULES_GENES[!modulesToTransfer];
# }

