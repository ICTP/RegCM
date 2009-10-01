#csh
ln -sf ../../Input/DOMAIN.INFO fort.10
ln -sf ../../PreProc/Terrain/domain.param
ifc -cm -w -w90 -w95 CMAP2RCM.f -L../../Commons/env/liblinux -lnetcdf -o cmap
setenv F_UFMTENDIAN big
./cmap
rm cmap fort.10 domain.param
