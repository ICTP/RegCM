#csh
ln -sf ../../Input/DOMAIN.INFO fort.10
ln -sf ../../PreProc/Terrain/domain.param
ifort -cm -w -w90 -w95 CMAP2RCM.f -L../../Commons/env/liblinux -lnetcdf -o cmap -convert big_endian
./cmap
rm cmap fort.10 domain.param
