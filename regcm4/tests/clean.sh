#!/bin/sh
# script to be execute before regcm execution ! 
# it accepts one parameter: month or day 
#
# it should work for all the regcm benchmarks
if test $# -ne 1
then
   # Not enough arguments
   echo -e "Usage:\n$0  month/day  " 1>&2
   exit 1
fi

# clean output directory
rm -rf output
rm *SAV*

# clean all possible useless file for such kind of run

# copy the right input
cp regcm.in-one-day regcm.in 

# 
mkdir output
#

echo -e " ready to execute " 
