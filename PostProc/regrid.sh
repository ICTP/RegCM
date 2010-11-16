#!/bin/bash

CDO=`which cdo`
BC=`which bc`

function calc
(
  scale=2; echo "$@" | $BC
)

if [ -z "$CDO" ]
then
  echo "Cannot find a cdo executable in Your path."
  echo
  echo 'Please install cdo (https://code.zmaw.de/projects/cdo)'
  echo 'If on CINECA SP6, You may want to:'
  echo 
  echo '            module load autoload cdo'
  echo
  exit 1
fi

if [ -z "$BC" ]
then
  echo "Cannot find a bc executable in Your path."
  echo
  echo 'Please install bc (http://www.gnu.org/software/bc)'
  echo
  exit 1
fi

if [ $# -ne 4 ]
then
  echo "Not enough arguments."
  echo ""
  echo "RegCM regridder example usage:"
  echo ""
  echo "                `basename $0` inputfile.nc griddes method"
  echo ""
  echo "where griddes is the regular latlon grid description in the form:"
  echo
  echo "                minlat,maxlat,deltalat minlon,maxlon,deltalon"
  echo
  echo "and method is one in:"
  echo
  echo "                bil bic dis nn"
  echo
  echo "aka bilinear, bicubic, distance weighted and nearest neighbor"
  echo
  echo "Example:"
  echo "         `basename $0` EUROPE_SRF.1990060100.nc 30.0,48.0,0.1 8.0,21.0,0.1 bil"
  echo ""
  exit 1
fi

infile=$1
latdes=$2
londes=$3
method=$4

declare -a lat
declare -a lon

lat=(`echo "${latdes//,/ }"`)
lon=(`echo "${londes//,/ }"`)

x1=`echo ${lon[0]}`
x2=`echo ${lon[1]}`
x3=`echo ${lon[2]}`
y1=`echo ${lat[0]}`
y2=`echo ${lat[1]}`
y3=`echo ${lat[2]}`

xx1=`echo "($x1+0.5)/1" | $BC`
xx2=`echo "($x2+0.5)/1" | $BC`
[ $xx1 -lt -180 ] && x1=`calc "$x1+360.0"`
[ $xx2 -lt -180 ] && x2=`calc "$x2+360.0"`
[ $xx1 -gt  180 ] && x1=`calc "$x1-360.0"`
[ $xx2 -gt  180 ] && x2=`calc "$x2-360.0"`

if [ $xx2 -lt 0 -a $xx1 -gt 0 ]
then
  xsize=`calc "((180.0-($x1)+($x2))/($x3))+1"`
else
  xsize=`calc "((($x2)-($x1))/($x3))+1"`
fi
ysize=`calc "((($y2)-($y1))/($y3))+1"`

echo "gridtype = lonlat" > griddes
echo "xsize = $xsize" >> griddes
echo "ysize = $ysize" >> griddes
echo "xfirst = $x1" >> griddes
echo "xinc = $x3" >> griddes
echo "yfirst = $y1" >> griddes
echo "yinc = $y3" >> griddes

$CDO remap$method,griddes $infile `basename $infile .nc`_lonlat.nc

rm -f griddes
exit 0
