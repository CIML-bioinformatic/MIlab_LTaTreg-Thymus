
# milab_ltatreg_r42_monocle3

Based on rocker/tidyverse distribution (includes Rstudio): 
Added:
 - Misc packages for data manipulation, figures, reports
 - Monocle3



## Build

docker build . -t milab_ltatreg_r42_monocle3



## Save

docker save milab_ltatreg_r42_monocle3 | gzip > /mnt/DOSI/MILAB/BIOINFO/Projet/LTaTReg/220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregSpl/02_Container/milab_ltatreg_r42_monocle3.tar.gz



## Run Rstudio

Give user details to internal script which sets user and permissions:

```
docker run -d --name milab_ltatreg_r42_monocle3 -p 9393:8787 -e PASSWORD=yourPass -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -v /mnt:/mnt milab_ltatreg_r42_monocle3
```

Then connect to the machine running docker (localhost) on mapped port (8787):
http://127.0.0.1:9292



## Run a command as specified user (not starting Rstudio)

docker run --rm \
           -e PASSWORD=yourPass \
           -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) \
           -v /mnt:/mnt \
           milab_ltatreg_r42_monocle3 \
           /init s6-setuidgid $(whoami) command

