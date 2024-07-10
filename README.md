# Lymphotoxin limits Foxp3+ regulatory T cell development from Foxp3lo precursors via IL-4 signaling

## Article information

**Title:** Lymphotoxin limits Foxp3+ regulatory T cell development from Foxp3lo precursors via IL-4 signaling

**Authors:** 
Alexia Borelli 1, Jérémy C. Santamaria 1, Cloé Zamit 1, Cécile Apert 2,3, Jessica Chevallier 1, Philippe Pierre 1, Rafael J. Argüello 1, Lionel Spinelli 1 and Magali Irla 1*
 
1 Aix-Marseille Univ, CNRS, INSERM, CIML, Centre d'Immunologie de Marseille-Luminy, Turing Centre for Living Systems, Marseille, France
2 Toulouse Institute for Infectious and Inflammatory Diseases (Infinity), INSERM UMR1291 - CNRS UMR5051 - University Toulouse III, Toulouse, France
3 Current affiliation : Microenvironment & Immunity Unit, Institut Pasteur, Paris, France
 
*Correspondence: Magali.Irla@inserm.fr

**Summary:**
Regulatory T cells (Treg), critical players of immune tolerance, develop in the thymus via two distinct developmental pathways involving CD25+Foxp3- and CD25-Foxp3lo precursors. However, the mechanisms regulating the recently identified Foxp3lo precursor pathway remain unclear. Here, we found that the membrane-bound lymphotoxin α1β2 (LTα1β2) heterocomplex is upregulated during Treg development upon TCR/CD28 and IL-2 stimulation. We show that Lta expression limits the maturational development of Treg from Foxp3lo precursors by regulating their proliferation, survival and metabolic profile. Transgenic reporter mice and transcriptomic analyses further revealed that medullary thymic epithelial cells (mTEC) constitute an unexpected source of IL-4. We demonstrate that LTα1β2-lymphotoxin β receptor-mediated interactions with mTEC limit Treg development by down-regulating IL-4 expression in mTEC. Collectively, our findings identify the lymphotoxin axis as the first inhibitory checkpoint of thymic Treg development that fine-tunes the Foxp3lo Treg precursor pathway by limiting IL-4 availability.

**DOI:**

---
---

## Goal of the github
This github project contains the instructions and material to reproduce the analyses reported in the article (and more).
Source code (scripts and dockerfiles) are available in the github repository. Required data and built Docker/Singularity images are available on download from Zenodo. Instructions to reproduce the analyses are provided below.

To reproduce the analysis, you have to first, prepare the environments (see "Prepare the Environments" section below), then execute the analysis step by step (see "Run the analysis" section below).

## Description of the datasets

Single cell suspensions of thymus were obtained by scratching organs through a 70-μm nylon mesh cell strainer with a plastic plunge in PBS 5%BSA 1mM EDTA. Thymic cell suspensions were then incubated with red blood cell lysis buffer for 3 minutes at room temperature and enriched for CD4+ thymocytes by depletion of CD8+ and CD11c+ cells using the AutoMACS® Pro Separator with the Deplete program (Miltenyi). We used cell hashing with hashtag oligonucleotides (HTO) to multiplex the two samples from 2 individual mice. Cell surface staining used for gating cells from the Treg lineage and staining for distinct barcoded anti-mouse CD45 antibody (Biolegend; A0304 and A0305 for Foxp3eGFP and Foxp3eGFPxLta-/- Treg cells, respectively) were performed in PBS 5%BSA 1mM EDTA for 30 min on ice. For each sample, cells from the Treg lineage Live Dead-CD4+CD8-CCR6-CD3e+ expressing either CD25, Foxp3eGFP or both were bulk-sorted with BD FACS Aria II. Sorted cell samples (20,000 Foxp3eGFP Treg cells and 20,000 Foxp3eGFPxLta-/- Treg cells) were pooled with a target of 10,000 captured cells and loaded in a single capture well for subsequent 10x Genomics Single Cell 3’ v3.1 workflow.

Library was performed according to the manufacter’s instructions (single cell 3’ v3.1 protocol, 10x Genomics). Briefly, cells were resuspended in the master mix and loaded together with partitioning oil and gel beads into the chip to generate the gel bead-in-emulsion (GEM). The poly-A RNA from the cell lysate contained in every single GEM was retrotranscripted to cDNA, which contains an Illumina R1 primer sequence, Unique Molecular Identifier (UMI) and the 10x Barcode. The pooled barcoded cDNA was then cleaned up with Silane DynaBeads, amplified by PCR and the appropriated sized fragments were selected with SPRIselect reagent. The pellet and supernatant fractions were separated for subsequent HTO and gene expression library construction. During the library construction, Illumina R2 primer sequence, paired-end constructs with P5 and P7 sequences and a sample index were added. The HTO library was constructed using Truseq D701 and D702 sequences containing i7 indexes. The resulting libraries were pooled and sequenced on an Illumina NextSeq2000 platform with a P2 flow cell (100 cycles). 

The resulting data were saved under the name 220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy.

mRNA and Cell Hashing FASTQ raw files were processed using Cell Ranger v6.0.1 (10X genomics Inc.) software with default parameters to performs alignment, filtering, barcode counting and unique molecular identifier (UMI) counting. Reads were aligned to the mouse mm10 genome. A total number of 4,798 cells were identified with a mean of 44,015 reads per cell and a median of 1,992 genes per cell.


## Prepare the environments

In order to prepare the environment for analysis execution, it is required to:

- Clone the github repository and set the WORKING_DIR environment variable
- Download the pre-processed data
- Install Docker
- Download the Docker images and load them on your system

Below you will find detailed instruction for each of these steps.

---

### Clone the github repository

Use you favorite method to clone this repository in a chosen folder. This will create a folder **LTaTReg** with all the source code. 

---

### Set the WORKING_DIR variable

Then, you must set an environment variable called **WORKING_DIR** with a value set to the path to this folder.

For instance, if you have chosen to clone the Git repository in __"/home/spinellil/workspace"__, then the **WORKING_DIR** variable will be set to __"/home/spinellil/workspace/LTaTReg"__

**On linux:**

```
    export WORKING_DIR=/home/spinellil/workspace/LTaTReg
```

---

### Add you working dir in the code

The code uses variables that are stored in different "parameters" file. One important variable is the PATH_PROJECT which indicate to the code where your project is stored.
You have to modify this variable in the code to reflect your project setup. The dataset folder has a file called **globalParams.R** in the subfolder **03_Script**

```
    LTaTReg
    │
    ├── 220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy
    │   │
    │   └── 03_Script/globalParams.R

```

Edit those files and in each of them locate the line defining the **PATH_PROJECT** variable and change its value to the same value as the **WORKING_DIR** variable you defined before. Then save the files.

```
PATH_PROJECT = "/home/spinellil/workspace/LTaTReg"
```

---


### Download the data

The sample needs its own sub-folder containing the initial data used by the analysis. These data can be downloaded from Zenodo and uncompressed. The Zenodo dataset DOI is @TODO : DOI ZENODO. The initial data from the analysis are the pre-processed data located in the following folder:

* 01_CellRanger_FeatureBarcoding : contains the result of Cell Ranger count analysis from the mRNA fastq

Scripts of the pre-processing steps are provided for Cell Ranger count analysis and corresponding raw data (fastq files) can be downloaded from GEO (see article).

To download and uncompress the data, use the following code:

**On linux:**

```
    cd $WORKING_DIR
    wget @TODO PATH ZENODO -O 220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy_processedData.tar.gz
    tar zxvf 10x_190712_m_moF220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy_processedDataluMemB_processedData.tar.gz
```
 
Once done, you may obtain the following subfolder structure, each of them containing several files.

```
    LTaTReg
    │
    ├── 220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregThy_processedData
    │   │
    │   └── 05_Output/01_CellRanger_FeatureBarcoding
    
```

---

### Install Docker

To install Docker on your system to take advantage of interactive analysis environment with Rstudio, follow the instructions here : https://docs.docker.com/engine/install/

---


### Download the Docker images and load them on your system

Docker image tar files are stored on Zenodo @TODO DOI TO DOCKER TAR. Open a shell command and change dir to the root of the cloned Git repository (WORKING_DIR). Then execute the following commands to download the tarball file, untar it and load the docker images on your system: 

```
    cd $WORKING_DIR
    wget @TODO PATH ZENODO -O LTaTreg_DockerImages.tar.gz
    tar zxvf LTaTreg_DockerImages.tar.gz
    docker load -i milab_ltatreg_r411_seurat4.tar
    docker load -i milab_ltatreg_r42_monocle3.tar
```

---
---

## Run the analysis

### Run the analysis using Docker

The analysis uses the Docker images that contains both the R language, the required R libraries and a Rstudio server.

The study contains 1 sample of single-cell RNA-seq data. The analysis is performed in several steps for which you will find the R script files in the subfolder **03_Script**. 

Each step of analysis generates its own HTML report file and several output files. Some output files of some steps are used by other steps, making a complete workflow of analysis. The output files are stored in the datatset folder in a sub-folder named "05_Output".

In order to start to reproduce the analysis, you have to run the Docker images. We provide two images :

* milab_ltatreg_r411_seurat4 : allows to perform almost all steps of analysis using the Seurat package.
* milab_ltatreg_r42_monocle3 : allows to perform the Monocle 3 analysis step.

To start a docker container from an image, use the following command (on Linux):

```
docker run -d -p 8787:8787 -v /$WORKING_DIR:/$WORKING_DIR -e PASSWORD=<PASSWORD> -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g)  milab_ltatreg_r411_seurat4

docker run -d -p 8888:8787 -v /$WORKING_DIR:/$WORKING_DIR -e PASSWORD=<PASSWORD> -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g)  milab_ltatreg_r42_monocle3

```

where:

* <PASSWORD> is a simple string you will use as password to login into Rstudio

One started, you can open an internet browser and use the URL https://127.0.0.1:8787 for almost all step of analysis and https://127.0.0.1:8788 for the analysis using Monocle 3.

At the login prompt, enter the name of the user session you are connected with and the password you type in place of <PASSWORD>. You are now in a Rstudio environment and the container is able to connect to the **WORKING_DIR** of your system and the folder defined as WORKING_DIR.
Inside you will find the project files. To run the analysis, go into one of the script folder (one of the analysis step), define this folder of the script as the R Working Dir in Rstudio and run the "launch_report_compilation.R" script that will execute the report HTML report.




