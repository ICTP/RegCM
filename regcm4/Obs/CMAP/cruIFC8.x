ln -sf ../../Input/DOMAIN.INFO fort.10
ln -sf ../../PreProc/Terrain/domain.param
ifort -cm -w -w90 -w95 CRU2RCM.f -L../../Commons/env/liblinux -lnetcdf -o cru -convert big_endian
./cru
rm cru fort.10 domain.param
