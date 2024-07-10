###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "03d_GlobalHeterogeneity_NoProliferation_NoContamination"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Analysis of Not proliferating Cells with removal of contamining CD4SP cells."

# Seurat Objects from previous analysis step
PATH_GLOBALHETEROGENEITY = list()
for( sample in SAMPLE_SET){
  PATH_GLOBALHETEROGENEITY[[ sample]] = file.path( PATH_EXPERIMENT_OUTPUT, "03a_GlobalHeterogeneity", sample)
}

# Filename of the Seurat Object from previous analysis
FILENAME_SEURAT_OBJECT = "LTaTReg_220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy_seuratObject_final.RDS"

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



#### Filtering / Normalization

# Filters for loading seurat object
LOAD_MIN_CELLS     = 3;    # Retain cells with at least this many features (annotations)
LOAD_MIN_FEATURES  = 200;  # Retain annotations appearing in at least this many cells

# Cells with number of UMIs outside the range will be excluded
FILTER_UMI_MIN     = 851;
FILTER_UMI_MAX     = 21686;

# Cells with number of genes outside the range will be excluded
FILTER_FEATURE_MIN = 603;
FILTER_FEATURE_MAX = 4128;

# Cells with percentage of mitocohondrial genes above threshold will be excluded
FILTER_MITOPCT_MAX = 9.94;

# Cells with percentage of ribosomal genes below threshold will be excluded
FILTER_RIBOPCT_MIN = 11.5;

# Normalization parameters (see Seurat::NormalizeData())
DATA_NORM_METHOD      = "LogNormalize";
DATA_NORM_SCALEFACTOR = 10000;

# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE;
DATA_SCALE        = FALSE;
DATA_VARS_REGRESS = NULL;  # c("nCount_RNA") for UMIs (NULL to ignore)




#### Analysis parameters

# Maximum number of variable features to keep
VARIABLE_FEATURES_MAXNB   = 2000;  # For analysis (PCA)
VARIABLE_FEATURES_SHOWTOP = 200;   # For table in report

# Nearest-neighbor graph construction
FINDNEIGHBORS_K = 30

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION     = 0.8;
FINDCLUSTERS_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results
FINDCLUSTERS_ALGORITHM      = 1;   # 1 = Louvain; 2 = Louvain with multilevel refinement; 3 = SLM; 4 = Leiden

RENAME_CLUSTERS = list()
RENAME_CLUSTERS[[ SAMPLE_ALL]] = list()
RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "0"]] = 4
RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "1"]] = 3
RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "2"]] = 1
RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "3"]] = 0
RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "4"]] = 5
RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "5"]] = 2

# PCA parameters
PCA_NPC              = 50;  # Default number of dimensions to use for PCA (see Seurat::RunPCA())
PCA_PLOTS_NBDIMS     = 3;   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15;  # Number of'top' features to show when plotting PCA loadings

# Dimensionality reduction parameters (TSNE/UMAP)
DIMREDUC_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 5;       # Number of marker genes to show in report and tables (NULL for all)

# Parameter for enrichment analysis in GO Terms
ENRICHMENT_GO_PVALUECUTOFF = 0.05
ENRICHMENT_GO_QVALUECUTOFF = 0.05



#### List of group of clusters of interest

CLUSTER_GROUP_LIST = list( ClusterGroup_NotProliferating_NoContamination = c( 2, 4, 6, 7, 8, 9))


#### Lists of genes of interest

# Number of genes under which genes of a module (MODULES_GENES) are transfered to be analyzed individually (MONITORED_GENES)
MONITORED_GENES_SMALL_MODULE = 5;
MODULES_CONTROL_SIZE = 100;


## Contamination-related genes (txt file, one gene name by line)
#CONTAMINATION_GENES = readLines( file.path( PATH_PROJECT_EXTERNALDATA, "ContaminationGenes.txt"))
CONTAMINATION_GENES=NULL

## Genes monitored individually (tsv file, one column for each group of genes)
MONITORED_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE,
                                                  "02_MonitoredGenesLists",
                                                  "MonitoredGenes.csv"),
                                       sep = ",",
                                       header = TRUE,
                                       stringsAsFactors = FALSE,
                                       row.names = NULL, fill = TRUE));
MONITORED_GENES = Map('[', MONITORED_GENES, lapply(MONITORED_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

# ## Genes monitored individually (tsv file, one column for each group of genes)
# CHOSEN_MARKER_GENES = read.table( file.path( PATH_EXPERIMENT_REFERENCE, 
#                                                   "02_MonitoredGenesLists",
#                                                   "ChosenMarkerGenes.csv"),
#                                        sep = ",",
#                                        header = TRUE,
#                                        stringsAsFactors = FALSE,
#                                        row.names = NULL, fill = TRUE);
# CHOSEN_MARKER_GENES = CHOSEN_MARKER_GENES$Gene

## Genes groups to analyze (csv file, one column for each group of genes)
GROUP_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, 
                                              "02_MonitoredGenesLists",
                                              "GeneGroups.csv"),
                                   sep = ",", quote = '"',
                                   header = TRUE,
                                   stringsAsFactors = FALSE,
                                   row.names = NULL, fill = TRUE));
GROUP_GENES = Map('[', GROUP_GENES, lapply(GROUP_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings


## Genes monitored as modules (tsv file, one column for each group of genes)
MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, "05_Module", "Modules.csv"),
                                        sep = ",",
                                        header = TRUE,
                                        stringsAsFactors = FALSE,
                                        row.names = NULL, fill = TRUE));
MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings
MODULES_GENES[[ "S_PHASE"]] = CELL_CYCLE_SPHASE_GENELIST
MODULES_GENES[[ "G2M_PHASE"]] = CELL_CYCLE_G2MPHASE_GENELIST


