dset ^PLEV_VAR
title model-output on P-plane
options big_endian
undef -1.e34
pdef   49   32 lcc   45.39   13.48   24.50   16.00   30.00   60.00   13.48  60000.  60000.
xdef  175 linear   -9.60  0.2703
ydef   78 linear   34.70  0.2703
zdef 11 levels 1000 925 850 700 500 400 300 250 200 150 100
tdef 60 linear 06z01Jul1994 6hr
vars 13
h 11 0 p-level height
t 11 0 air temperature
u 11 0 zonal wind
v 11 0 meridional wind
omega 11 0 vertical velocity of p-level
q 11 0 specific moisture
qc 11 0 cloud water
ps 0 99 surface pressure
slp 0 99 sea level pressure
rain 0 99 precipitation
tg 0 99 temperature of lower soil layer
swt 0 99 total soil water in mm H2O
rno 0 99 accumulated infiltration
endvars
