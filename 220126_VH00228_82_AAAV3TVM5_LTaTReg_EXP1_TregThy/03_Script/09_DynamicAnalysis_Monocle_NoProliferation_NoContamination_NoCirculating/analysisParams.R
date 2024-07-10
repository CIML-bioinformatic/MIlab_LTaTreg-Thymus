###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "09_DynamicAnalysis_Monocle_NoProliferation_NoContamination_NoCirculating"
ANALYSIS_STEP_LITTERAL_DESCRIPTION = "Analysis of cell dynamics using Monocle 3 (No proliferation, no contamination, no circulating)"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Seurat object from Heterogeniety analysis
SEURAT_OBJECT_HETEROGENEITY_ANALYSIS = file.path( PATH_EXPERIMENT_OUTPUT, "03f_GlobalHeterogeneity_NoProliferation_NoContamination_NoCirculating", 
                                                  "All-TregThy", 
                                                  "ClusterGroup_NotProliferating_NoContamination_NoCirculating", 
                                                  "0.8", 
                                                  "LTaTReg_220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy_seuratObject_final.RDS")

## Genes monitored as modules (tsv file, one column for each group of genes)
MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, "05_Module", "Modules.csv"),
                                     sep = ",",
                                     header = TRUE,
                                     stringsAsFactors = FALSE,
                                     row.names = NULL, fill = TRUE));
MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

MONITORED_GENES = sort( unique( unlist( MODULES_GENES)))


## Genes groups to analyze (csv file, one column for each group of genes)
GROUP_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, 
                                              "02_MonitoredGenesLists",
                                              "GeneGroups.csv"),
                                   sep = ",", quote = '"',
                                   header = TRUE,
                                   stringsAsFactors = FALSE,
                                   row.names = NULL, fill = TRUE));
GROUP_GENES = Map('[', GROUP_GENES, lapply(GROUP_GENES, function(x){ which( nchar( x)>0)}))

## Monocle trajectory between custom pairs of vertex
TRAJECTORY_CUSTOM_VERTEX_DF = data.frame( c( "Y_76", "Y_44"), 
                                          c( "Y_76", "Y_16"), 
                                          c( "Y_76", "Y_57"), 
                                          c( "Y_41", "Y_44"), 
                                          c( "Y_41", "Y_16"), 
                                          c( "Y_41", "Y_57"), 
                                          c( "Y_30", "Y_44"), 
                                          c( "Y_1", "Y_16"), 
                                          c( "Y_2", "Y_16"), 
                                          c( "Y_30", "Y_16")
                              )

## List of KEGG pathways of interest to study in details
MONITORED_KEGG_PATHWAY_LIST = list( "Glycolysis_Gluconeogenesis" = "mmu00010",
                                    "Oxidative_phosphorylation" = "mmu00190"
                                  )

