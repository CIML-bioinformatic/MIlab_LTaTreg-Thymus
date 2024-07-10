
# Image to run velocyto

The image contains:

 - Python 3.6
 - velocyto
 
 
## BUILD
 
docker build . -t milab_ltatreg_velocyto
 
## SAVE

docker save milab_ltatreg_velocyto | gzip > /mnt/DOSI/MILAB/BIOINFO/Projet/LTaTReg/220126_VH00228_82_AAAV3TVM5_LTaTReg_EXP1_TregSpl/02_Container/milab_ltatreg_velocyto.tar.gz

## RUN

Get an interactive session to run the velocyto scripts:

```
docker run -it --name milab_ltatreg_velocyto -v /mnt:/mnt milab_ltatreg_velocyto /bin/bash
```

