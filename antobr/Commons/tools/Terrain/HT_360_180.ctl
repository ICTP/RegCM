dset ^HT_360_180.dat
title RegCM domain information
options big_endian
undef -32767.
xdef 360 linear  -179.5  1.0000
ydef 180 linear  -89.50  1.0000
zdef 1 levels 1000.00
tdef 1 linear 00z01Jan2001 1mo
vars 2
ht      0 99 surface elevation         
htsd    0 99 surface elevation std. dev
endvars
