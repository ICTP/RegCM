ln -sf ../../Input/DOMAIN.INFO fort.10
ln -sf ../../PreProc/Terrain/domain.param
pgf77 -byteswapio WM2RCM.f -L../../Commons/env/liblinux -lnetcdf -o wm
./wm
rm wm fort.10 domain.param
