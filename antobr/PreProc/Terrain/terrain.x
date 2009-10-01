#!/bin/csh -f
make clean
make
./terrain
/bin/rm -f terrain*.o terrain
