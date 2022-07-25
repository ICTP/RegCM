#!/bin/bash

if [ $# -lt 2 ]
then
  echo "Usage : $0 domain.clm.file_2d.nc domain_SRF.file.nc.ctl [zlev]"
  echo
  echo "zlev defaults to levgrnd"
  exit -1
fi

file=$1
srfctl=$2
zlev=levgrnd
if [ $# -eq 3 ]
then
  zlev=$3
  echo "Selecting levels $zlev"
fi
xdefstr=`cat $srfctl | grep "^xdef"`
ydefstr=`cat $srfctl | grep "^ydef"`
pdefstr=`cat $srfctl | grep "^pdef"`
ctlfile=${file}.ctl
domain=`echo $file | cut -d '.' -f 1`
ncdump -h $file > hdump
level=(`ncks -H -C -s "%g " -v $zlev ${file}`)
months=( jan feb mar apr may jun jul aug sep oct nov dec )
yearmon=`echo $file | cut -d '.' -f 5 | cut -d '_' -f 1`
year=`echo $yearmon | cut -b 1-4`
mon=`echo $yearmon | cut -b 5-6`
vars=(`cat hdump | grep "(iy, jx)" | \
       awk '{print $2}' | sed -e 's/(iy,//'`)
timevars=(`cat hdump | grep "time, iy, jx)" | grep time | \
           awk '{print $2}' | sed -e 's/(time,//'`)
groundvars=(`cat hdump | grep "time, $zlev, iy, jx)" | \
             grep time | awk '{print $2}' | sed -e 's/(time,//'`)
nvars=$(( ${#vars[@]} + ${#timevars[@]} + ${#groundvars[@]} ))
nlevs=${#level[@]}

echo dset ^$file > $ctlfile
echo dtype netcdf >> $ctlfile
echo undef 1e+20_FillValue >> $ctlfile
echo title ICTP Regional Climatic model V5 >> $ctlfile
echo $pdefstr >> $ctlfile
echo $xdefstr >> $ctlfile
echo $ydefstr >> $ctlfile
echo zdef  $nlevs levels ${level[@]} >> $ctlfile 
echo tdef  1 linear 00Z15${months[$(( $mon - 1 ))]}${year} 1mon >> $ctlfile
echo vars $nvars >> $ctlfile
for var in ${vars[@]}
do
  longname=`cat hdump | grep -P "\t$var:long_name" | cut -d '"' -f 2```
  units=`cat hdump | grep -P "\t$var:units" | cut -d '"' -f 2```
  echo "$var=>$var 0 y,x $longname ($units)" >> $ctlfile
done
for var in ${timevars[@]}
do
  longname=`cat hdump | grep -P "\t$var:long_name" | cut -d '"' -f 2```
  units=`cat hdump | grep -P "\t$var:units" | cut -d '"' -f 2```
  echo "$var=>$var 0 t,y,x $longname ($units)" >> $ctlfile
done
for var in ${groundvars[@]}
do
  longname=`cat hdump | grep -P "\t$var:long_name" | cut -d '"' -f 2```
  units=`cat hdump | grep -P "\t$var:units" | cut -d '"' -f 2```
  echo "$var=>$var $nlevs t,z,y,x $longname ($units)" >> $ctlfile
done
echo endvars >> $ctlfile
rm -f hdump
