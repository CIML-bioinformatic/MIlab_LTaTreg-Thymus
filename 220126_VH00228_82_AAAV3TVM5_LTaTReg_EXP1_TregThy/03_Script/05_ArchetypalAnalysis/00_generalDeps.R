# ##################################################
# Global declarations and libraries for the analysis
# ##################################################

# Packages for reporting
library( funr)
library( digest)
library( fs)
library( pander)
library( dplyr)

# Package to plot data
library( ggplot2)
library( kableExtra)

# Technology package
library( Seurat)
library( limma)
library( stats)

# Python/R interaction
library( reticulate)
use_python("/usr/bin/python3")
