###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "07_Integration_OwenEtAl"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Compare project data with data produced by Owen et al. (PMID:36038290)"

# Retrieve result of CellRanger from previous analysis step
PATH_SEURAT_OBJECT_PROJECT = c( "RUN1_THY" = file.path( PATH_EXPERIMENT_OUTPUT, 
                                                        "03a_GlobalHeterogeneity",
                                                        "All-TregThy",
                                                        "LTaTReg_220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy_seuratObject_final.RDS"))

# Retrieve the result of the Owen et al. 2022 study (PMID 36038290)
PATH_SEURAT_OBJECT_OWENETAL = c( "OwenEtAl_PMID_36038290" = file.path( PATH_EXPERIMENT_REFERENCE,
                                                          "07_PublicDataset", 
                                                          "2022_OwenEtAl_PMID_36038290", 
                                                          "GEO_GSM5851360", 
                                                          "GSM5851360_final_seurat_analysis_object.rds"))

# Retrieve result of CellRanger from previous analysis step when removing contaminating and proliferating cells
PATH_NOCONTAMINATION_NOPROLIFERATION_SEURAT_OBJECT_RDS = file.path( PATH_EXPERIMENT_OUTPUT, 
                                                                           "03d_GlobalHeterogeneity_NoProliferation_NoContamination",
                                                                           "All-TregThy",
                                                                           "ClusterGroup_NotProliferating_NoContamination",
                                                                           "0.8",
                                                                           "LTaTReg_220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy_seuratObject_final.RDS")

#### General

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (for Seurat3 using 'future')
NBCORES = 4;

#### Analysis parameters

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
