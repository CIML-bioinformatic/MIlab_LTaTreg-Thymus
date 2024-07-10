###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "04_DynamicAnalysis_Velocyto"
ANALYSIS_STEP_LITTERAL_DESCRIPTION = "Analysis of cell dynamics using Velocyto"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Location of the velocyto loom files
VELOCYTO_LOOM_FILE = file.path( PATH_ANALYSIS_OUTPUT, "loom", "mm10.loom")

# Seurat object from Heterogeniety analysis
SEURAT_OBJECT_HETEROGENEITY_ANALYSIS = file.path( PATH_EXPERIMENT_OUTPUT, "03a_GlobalHeterogeneity", paste0( outputFilesPrefix, "seuratObject_final.RDS"))

## Genes monitored as modules (tsv file, one column for each group of genes)
MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, "05_Module", "Modules.csv"),
                                     sep = ",",
                                     header = TRUE,
                                     stringsAsFactors = FALSE,
                                     row.names = NULL, fill = TRUE));
MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

MONITORED_GENES = sort( unique( c( MODULES_GENES[[ "Resting_Treg"]], MODULES_GENES[[ "Early_Activated"]], MODULES_GENES[[ "Effector_Treg"]], MODULES_GENES[[ "LateActivated_AgedTreg"]])))
