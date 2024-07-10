###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples. 
#


#### General

GLOBAL_DESCRIPTION = "LTaTReg - TregThy"

SCIENTIFIC_GROUP = "MILAB"
SCIENTIFIC_PROJECT_NAME = "LTaTReg"
EXPERIMENT_NAME = "220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy"



#### Additional custom variables (for a specific analysis step, tools, ...)

#SAMPLE_ID = ""



#### Input / Output

# Output folder name in data folder (for R session object, lists of cells/genes) 
PATH_PROJECT = file.path( "/mnt/DOSI", 
                                        SCIENTIFIC_GROUP,
                                        "BIOINFO", 
                                        "Projet",
                                        SCIENTIFIC_PROJECT_NAME)

PATH_EXPERIMENT = file.path( PATH_PROJECT, EXPERIMENT_NAME)

PATH_EXPERIMENT_RAWDATA       = file.path( PATH_EXPERIMENT, "00_RawData")
PATH_EXPERIMENT_REFERENCE = file.path( PATH_EXPERIMENT, "01_Reference")
PATH_EXPERIMENT_OUTPUT        = file.path( PATH_EXPERIMENT, "05_Output")



# Create a 'safe' unique prefix for output files
outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME, "_",     
                            EXPERIMENT_NAME, "_"
                            #startTimeFileName, "_",
                            #paramsHash, "_"
                          )



#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



