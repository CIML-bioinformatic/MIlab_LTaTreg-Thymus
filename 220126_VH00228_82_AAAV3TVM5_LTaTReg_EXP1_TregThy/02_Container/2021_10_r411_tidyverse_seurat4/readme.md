
# r411_tidyverse_seurat4

Based on rocker/tidyverse distribution (includes Rstudio): 
Added:
 - Misc packages for data manipulation, figures, reports
 - Seurat
 - RFutils (ppackge with custom helper functions: VENN, ...)



## Build

docker build . -t milab_ltatreg_r411_seurat4



## Save

docker save milab_ltatreg_r411_seurat4 | gzip > /mnt/DOSI/MILAB/BIOINFO/Projet/LTaTReg/220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregSpl/02_Container/milab_ltatreg_r411_seurat4.tar.gz



## Run Rstudio

Give user details to internal script which sets user and permissions:

```
docker run -d --name milab_ltatreg_r411_seurat4 -p 9090:8787 -e PASSWORD=yourPass -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -v /mnt:/mnt milab_ltatreg_r411_seurat4
```

Then connect to the machine running docker (localhost) on mapped port (8787):
http://127.0.0.1:9090



## Run a command as specified user (not starting Rstudio)

docker run --rm \
           -e PASSWORD=yourPass \
           -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) \
           -v /mnt:/mnt \
           milab_ltatreg_r411_seurat4 \
           /init s6-setuidgid $(whoami) command

