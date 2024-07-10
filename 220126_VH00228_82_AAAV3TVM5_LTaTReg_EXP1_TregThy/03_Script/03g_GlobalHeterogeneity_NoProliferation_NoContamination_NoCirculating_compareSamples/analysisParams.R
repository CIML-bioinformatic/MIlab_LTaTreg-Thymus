###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "03g_GlobalHeterogeneity_NoProliferation_NoContamination_NoCirculating_compareSamples"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Compare heterogeneity between samples (conditions WT vs LTa KO) (No Proliferation, No Contamination)"


# Seurat Objects from previous analysis step
PATH_GLOBALHETEROGENEITY = list()
for( sample in SAMPLE_SET){
  PATH_GLOBALHETEROGENEITY[[ sample]] = file.path( PATH_EXPERIMENT_OUTPUT, "03f_GlobalHeterogeneity_NoProliferation_NoContamination_NoCirculating", sample, "ClusterGroup_NotProliferating_NoContamination_NoCirculating", "0.8")
}

# Filename of the Seurat Object from previous analysis
FILENAME_SEURAT_OBJECT = "LTaTReg_220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy_seuratObject_final.RDS"

# File name of cluster gene marker files
FILENAME_CLUSTER_MARKERGENES = "LTaTReg_220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy_MarkerGenes.tsv"

# Path to the Cell cycle gene lists
CELL_CYCLE_SPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "03_CellCycle", "S_phase_genes.csv")))
CELL_CYCLE_G2MPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "03_CellCycle", "G2M_phase_genes.csv")))

# Path to Heat Shock stress genes list
PATH_HS_STRESS_MARKER_GENES_TABLE_FILE = file.path( PATH_EXPERIMENT_REFERENCE, "06_HeatShock", "coregene_df-FALSE-v3.csv")

#### General

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (for Seurat3 using 'future')
NBCORES = 4;

# Number of cells above which use ggplot instead of interactive plotly
PLOT_RASTER_NBCELLS_THRESHOLD = 20000;

# Alpha error for the FindMakers analysis
MARKER_GENES_ALPHA_THRESHOLD = 0.05


#### Lists of genes of interest

# Number of genes under which genes of a module (MODULES_GENES) are transfered to be analyzed individually (MONITORED_GENES)
MONITORED_GENES_SMALL_MODULE = 5;
MODULES_CONTROL_SIZE = 100;


## Contamination-related genes (txt file, one gene name by line)
#CONTAMINATION_GENES = readLines( file.path( PATH_PROJECT_EXTERNALDATA, "ContaminationGenes.txt"))
CONTAMINATION_GENES=NULL


## List of KEGG pathways of interest to study in details

MONITORED_KEGG_PATHWAY_LIST = list( "Fatty_acid_metabolism" = "mmu01212",
                                      "Glycolysis_Gluconeogenesis" = "mmu00010",
                                      "Citrate_cycle_TCA_cycle" = "mmu00020",
                                      "Pentose_phosphate_pathway" = "mmu00030",
                                      "Oxidative_phosphorylation" = "mmu00190",
                                      "Fatty_acid_biosynthesis" = "mmu00061",
                                      "Fatty_acid_degradation" = "mmu00071",
                                      "Fatty_acid_elongation" = "mmu00062",
                                      "Pyruvate_metabolism" = "mmu00620",
                                      "Biosynthesis_of_amino_acids" = "mmu01230",
                                      "Nucleotide_metabolism" = "mmu01232",
                                      "Biosynthesis_of_cofactors" = "mmu01240",
                                      "Apoptosis" = "mmu04210",
                                      "Cell_cycle" = "mmu04110",
                                      "NF_kappa_B_signaling_pathway" = "mmu04064",
                                      "TNF_signaling_pathway" = "mmu04668",
                                      "Cytokine_cytokine_receptor_interaction" = "mmu04060",
                                      "JAK_STAT_signaling_pathway" = "mmu04630"
                                    )

## Genes monitored individually (tsv file, one column for each group of genes)
#MONITORED_GENES = list()
MONITORED_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, 
                                                  "02_MonitoredGenesLists",
                                                  "MonitoredGenes.csv"),
                                       sep = ",",
                                       header = TRUE,
                                       stringsAsFactors = FALSE,
                                       row.names = NULL, fill = TRUE));
MONITORED_GENES = Map('[', MONITORED_GENES, lapply(MONITORED_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

# MONITORED_GENES = c(list(Contamination = CONTAMINATION_GENES), MONITORED_GENES) # Add contamination genes as a group of genes to be monitored

## Genes groups to analyze (csv file, one column for each group of genes)
GROUP_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, 
                                                  "02_MonitoredGenesLists",
                                                  "GeneGroups.csv"),
                                       sep = ",", quote = '"',
                                       header = TRUE,
                                       stringsAsFactors = FALSE,
                                       row.names = NULL, fill = TRUE));
GROUP_GENES = Map('[', GROUP_GENES, lapply(GROUP_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

# ## Genes monitored individually (tsv file, one column for each group of genes)
# CHOSEN_MARKER_GENES = read.table( file.path( PATH_EXPERIMENT_REFERENCE, 
#                                              "02_MonitoredGenesLists",
#                                              "ChosenMarkerGenes.csv"),
#                                   sep = ",",
#                                   header = TRUE,
#                                   stringsAsFactors = FALSE,
#                                   row.names = NULL, fill = TRUE);
# CHOSEN_MARKER_GENES = CHOSEN_MARKER_GENES$Gene

## Genes monitored as modules (tsv file, one column for each group of genes)
MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, "05_Module", "Modules.csv"),
                                        sep = ",",
                                        header = TRUE,
                                        stringsAsFactors = FALSE,
                                        row.names = NULL, fill = TRUE));
MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings
MODULES_GENES[[ "S_PHASE"]] = CELL_CYCLE_SPHASE_GENELIST
MODULES_GENES[[ "G2M_PHASE"]] = CELL_CYCLE_G2MPHASE_GENELIST


