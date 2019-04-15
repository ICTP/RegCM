#!/bin/bash
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#
#    ICTP RegCM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# This scripts prepares a dataset downloaded from dss.ucar.edu
#
#  ds083.2 NCEP FNL Operational Model Global Tropospheric Analyses,
#          continuing from July 1999
#
# transforming from the original GRIB1 or GRIB2 format to NetCDF format
# for the RegCM ICBC to prepare input for the RegCM model.
#

G1LIST=(PRES_3_SFC_10 HGT_3_ISBL_10 R_H_3_ISBL_10 \
        TMP_3_ISBL_10 U_GRD_3_ISBL_10 V_GRD_3_ISBL_10)
T1VAR="TMP_3_ISBL_10"
RH1VAR="R_H_3_ISBL_10"
G1VAR="HGT_3_SFC_10"
G1LAT=lat_3
G1LON=lon_3
G1LLV="$G1LAT $G1LON"
G2LIST=(PRES_P0_L1_GLL0 HGT_P0_L100_GLL0 RH_P0_L100_GLL0 \
        TMP_P0_L100_GLL0 UGRD_P0_L100_GLL0 VGRD_P0_L100_GLL0)
T2VAR="TMP_P0_L100_GLL0"
RH2VAR="RH_P0_L100_GLL0"
G2VAR="HGT_P0_L1_GLL0"
G2LAT=lat_0
G2LON=lon_0
G2LLV="$G2LAT $G2LON"

transform ()
{
  local GLIST=()
  local GLLV=""
  local GLAT=""
  local GLON=""
  local TVAR=""
  local RHVAR=""
  if [ $2 -eq 1 ]
  then
    GLIST=(${G1LIST[@]})
    GLLV=$G1LLV
    GLAT=$G1LAT
    GLON=$G1LON
    TVAR=$T1VAR
    RHVAR=$RH1VAR
  else
    GLIST=(${G2LIST[@]})
    GLLV=$G2LLV
    GLAT=$G2LAT
    GLON=$G2LON
    TVAR=$T2VAR
    RHVAR=$RH2VAR
  fi
  echo Using $GLAT $GLON $TVAR $RHVAR
  local NLIST=(ps hga rha ta ua va)
  local XVAR=`echo ${GLIST[@]} | sed -e 's/\([ ]\+\)/,/g'`
  $G2NC $1 -e grb > /dev/null 2>&1
  local VVAR1=`ncdump -h $1.nc | grep $TVAR | \
               head -1 | sed -e 's/.\+(//' -e 's/) ;//' \
	                     -e 's/,//g' -e "s/ $GLLV//"`
  local VVAR2=`ncdump -h $1.nc | grep $RHVAR | \
               head -1 | sed -e 's/.\+(//' -e 's/) ;//' \
	                     -e 's/,//g' -e "s/ $GLLV//"`
  local CVAR=`echo $VVAR1 $VVAR2 $GLLV | sed -e 's/\([ ]\+\)/,/g'`
  ncks -v $CVAR,$XVAR $1.nc tmp.nc
  mv tmp.nc $1.nc
  $GRNM -d $GLAT,lat -d $GLON,lon -d $VVAR1,lev -d $VVAR2,rhlev  \
        -v $GLAT,lat -v $GLON,lon -v $VVAR1,lev -v $VVAR2,rhlev  ${1}.nc
  for (( i = 0; i < ${#GLIST[@]}; i ++ ))
  do
    $GRNM -v ${GLIST[$i]},${NLIST[$i]} ${1}.nc
  done
  year=`echo $1 | cut -b 5-8` # fnl_YYYYMMDD_HH_MM
  mkdir -p $year
  mv ${1}.nc $year
}

orog ()
{
  local GVAR=""
  local GLLV=""
  local GLAT=""
  local GLON=""
  local TVAR=""
  local RHVAR=""
  if [ $2 -eq 1 ]
  then
    GVAR=$G1VAR
    GLLV=$G1LLV
    GLAT=$G1LAT
    GLON=$G1LON
    TVAR=$T1VAR
    RHVAR=$RH1VAR
  else
    GVAR=$G2VAR
    GLLV=$G2LLV
    GLAT=$G2LAT
    GLON=$G2LON
    TVAR=$T2VAR
    RHVAR=$RH2VAR
  fi
  [ -d fixed ] || mkdir fixed
  [ -f fixed/fixed_orography.nc ] && return
  local NVAR=orog
  $G2NC $1 -e grb > /dev/null 2>&1
  local VVAR1=`ncdump -h $1.nc | grep $TVAR | \
               head -1 | sed -e 's/.\+(//' -e 's/) ;//' \
	                     -e 's/,//g' -e 's/ lat_3 lon_3//'`
  local VVAR2=`ncdump -h $1.nc | grep $RHVAR | \
               head -1 | sed -e 's/.\+(//' -e 's/) ;//' \
	                     -e 's/,//g' -e 's/ lat_3 lon_3//'`
  local CVAR=`echo $VVAR1 $VVAR2 lat_3 lon_3 | sed -e 's/\([ ]\+\)/,/g'`
  ncks -v $CVAR,$GVAR $1.nc tmp.nc
  mv tmp.nc $1.nc
  $GRNM -d $GLAT,lat -d $GLON,lon -d $VVAR1,lev -d $VVAR2,rhlev \
        -v $GLAT,lat -v $GLON,lon -v $VVAR1,lev -v $VVAR2,rhlev \
	-v $GVAR,$NVAR ${1}.nc
  mv ${1}.nc fixed/fixed_orography.nc
}

echo
echo "This scripts prepares a dataset downloaded from dss.ucar.edu"
echo
echo  "ds083.2 NCEP FNL Operational Model Global Tropospheric Analyses,"
echo  "        continuing from July 1999"
echo
echo "transforming from the original GRIB1 or GRIB2 format to NetCDF format"
echo "for the RegCM ICBC to prepare input for the RegCM model."
echo

G2NC=`which ncl_convert2nc`
if [ "X$G2NC" = "X" ]
then
  echo NCL is not currently installed on this system.
  echo This script rely on it. Bailing out.
  exit 2
fi

GRNM=`which ncrename`
if [ "X$GRNM" = "X" ]
then
  echo NCO is not currently installed on this system.
  echo This script rely on it. Bailing out.
  exit 3
fi

WGRB=`which wgrib`
if [ "X$WGRB" = "X" ]
then
  echo wgrib is not currently installed on this system.
  echo This script rely on it. Bailing out.
  exit 4
fi

fnls=`ls -1 fnl_* 2> /dev/null | head -1`
if [ "X$fnls" = "X" ]
then
  echo "Nothing to do."
  exit 0
fi

ff=`ls -1 fnl_* | head -1`
INPVER=1
XTEST=`$WGRB -s $ff`
[ -z "$XTEST" ] && INPVER=2

echo "Extracting fixed orography."
orog $ff $INPVER

echo "Converting to netCDF."
fnlsg=`ls -1 fnl_*`
for file in $fnlsg
do
  INPVER=1
  XTEST=`$WGRB -s $file`
  [ -z "$XTEST" ] && INPVER=2
  transform $file $INPVER
  rm -f $file
done

echo "Done"
exit 0
