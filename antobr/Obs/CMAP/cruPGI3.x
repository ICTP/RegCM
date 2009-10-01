ln -sf ../../Input/DOMAIN.INFO fort.10
ln -sf ../../PreProc/Terrain/domain.param
pgf77 -byteswapio CRU2RCM.f -L../../Commons/env/liblinuxold -lnetcdf -o cru
./cru
rm cru fort.10 domain.param
