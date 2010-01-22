#csh
ln -sf ../../Input/DOMAIN.INFO fort.10
ln -sf ../../PreProc/Terrain/domain.param
ifc -tpp7 -O3 -cm -w -w90 -w95 WM2RCM.f -L../../Commons/env/liblinux -lnetcdf -o wm
setenv F_UFMTENDIAN big
./wm
rm wm fort.10 domain.param
