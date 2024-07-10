###############################################################################
# This file defines SAMPLE parameters as global variables that will be loaded
# before analysis starts. 
#


SAMPLE_ALL = "All-TregThy"
SAMPLE_WT = "WT-TregThy"
SAMPLE_KO = "Lta-TregThy"

SAMPLE_SET = c( SAMPLE_ALL, SAMPLE_WT, SAMPLE_KO)

SAMPLE_COLOR = c( "#00BA38", "#619CFF", "#F8766D")
names( SAMPLE_COLOR) = SAMPLE_SET
