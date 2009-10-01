#!/bin/csh -f
make clean
make SST_ERSST
./SST_ERSST
/bin/rm -f SST_ERSST*.o SST_ERSST
make ICBC
./ICBC
/bin/rm -f ICBC.o ICBC SST.RCM
