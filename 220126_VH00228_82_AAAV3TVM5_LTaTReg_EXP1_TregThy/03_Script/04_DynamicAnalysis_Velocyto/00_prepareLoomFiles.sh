#!/bin/bash

# This script aims to execute the velocyto program.
# 

# Provide the folder where the 10x data are stored
export RAW_DATA_PATH=/mnt/DOSI/MILAB/BIOINFO/Projet/LTaTReg/220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy/05_Output/01_CellRanger_FeatureBarcoding/mm10
# Provide the path to the genome annotation file
export GENOME_ANNOTATION_PATH=/mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/01_COMMON_DATA/01_REFERENCE/01_GENOME/CellRanger/mouse/2020-A/refdata-gex-mm10-2020-A/genes/genes.gtf
# Provide the folder where to put the Velocyto results (loom file)
export OUTPUT_PATH=/mnt/DOSI/MILAB/BIOINFO/Projet/LTaTReg/220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy/05_Output/04_DynamicAnalysis_Velocyto
mkdir -p $OUTPUT_PATH/loom

# Ensure Samtools is in the PATH
export PATH=/samtools/bin:$PATH

# Execute Velocyto on 10x data
cd $OUTPUT_PATH
velocyto run10x $RAW_DATA_PATH $GENOME_ANNOTATION_PATH

# Move the Velocyto results (put by the tool in the input folder) in the desired output folder
mv $RAW_DATA_PATH/velocyto/* $OUTPUT_PATH
rmdir $RAW_DATA_PATH/velocyto
