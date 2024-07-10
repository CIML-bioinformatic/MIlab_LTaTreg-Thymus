###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "06_DynamicAnalysis_Velocyto_Archetype"
ANALYSIS_STEP_LITTERAL_DESCRIPTION = "Analysis of cell dynamics using scVelo"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Location of the velocyto loom files
VELOCYTO_LOOM_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "04_DynamicAnalysis_Velocyto", "loom", "mm10.loom")

# Seurat object from Archetypal analysis
AUTOENCODER_LAYER_INTERNAL_NODES = 5
AUTOENCODER_LAYER_IN_OUT_NODES = 512
SEURAT_OBJECT_ARCHETYPAL_ANALYSIS = file.path( PATH_EXPERIMENT_OUTPUT, "05_ArchetypalAnalysis", paste0( "archetype_K", AUTOENCODER_LAYER_INTERNAL_NODES, "_N", AUTOENCODER_LAYER_IN_OUT_NODES, "_seuratObject_final.RDS"))

## Genes monitored as modules (tsv file, one column for each group of genes)
MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, "05_Module", "Modules.csv"),
                                     sep = ",",
                                     header = TRUE,
                                     stringsAsFactors = FALSE,
                                     row.names = NULL, fill = TRUE));
MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

MONITORED_GENES = sort( unique( c( MODULES_GENES[[ "Resting_Treg"]], MODULES_GENES[[ "Early_Activated"]], MODULES_GENES[[ "Effector_Treg"]], MODULES_GENES[[ "LateActivated_AgedTreg"]])))
