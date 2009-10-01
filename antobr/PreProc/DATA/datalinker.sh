#simple script to create symbolic links from RegCM 
# DATA_LOCATION: specify the dir

DATA_LOCATION=/scratch/DATA/RegCM

mkdir SURFACE
/bin/ln -sf $DATA_LOCATION/SURFACE/* SURFACE/.
mkdir AERGLOB
/bin/ln -sf $DATA_LOCATION/AERGLOB/* AERGLOB/.
mkdir SST
/bin/ln -sf $DATA_LOCATION/SST/* SST/.
mkdir ECWCRP
/bin/ln -sf $DATA_LOCATION/ECWCRP/* ECWCRP/.
mkdir ERA40
/bin/ln -sf $DATA_LOCATION/ERA40/* ERA40/.
mkdir NNRP1
/bin/ln -sf $DATA_LOCATION/NNRP1/* NNRP1/.
mkdir NNRP2
/bin/ln -sf $DATA_LOCATION/NNRP2/* NNRP2/.

#this last added   by hand 
mkdir EIN15 
/bin/ln -sf $DATA_LOCATION/ERAIN150/* EIN15/.

#mkdir FVGCM
#/bin/ln -sf /home/RAID2-D10/RCM3DATA/FVGCM/* FVGCM/.
