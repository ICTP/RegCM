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
# This scripts prepares a dataset downloaded from dss.ucar.edu, either
#
#  ds082.0 NCEP FNL Operational Model Global Tropospheric Analyses,
#          July 1976 to April 1997
#  ds083.0 NCEP FNL Operational Model Global Tropospheric Analyses,
#          April 1997 to June 2007 
#  ds083.2 NCEP FNL Operational Model Global Tropospheric Analyses,
#          continuing from July 1999
#
# transforming from the original GRIB1 or GRIB2 format to NetCDF format
# for the RegCM ICBC to prepare input for the RegCM model.
#

transform ()
{
  local GLIST=(TMP_3_SFC_10 PRES_3_SFC_10 HGT_3_ISBL_10 R_H_3_ISBL_10 \
         TMP_3_ISBL_10 U_GRD_3_ISBL_10 V_GRD_3_ISBL_10)
  local NLIST=(ts ps hga rha ta ua va)
  local XVAR=`echo ${GLIST[@]} | sed -e 's/\([ ]\+\)/,/g'`
  local CVAR="lon_3,lat_3,lv_ISBL2,lv_ISBL6"
  $G2NC $1 -u time -e grb -v ${CVAR},${XVAR} > /dev/null 2>&1
  $GRNM -d lat_3,lat -d lon_3,lon -d lv_ISBL2,lev -d lv_ISBL6,rhlev  \
        -v lat_3,lat -v lon_3,lon -v lv_ISBL2,lev -v lv_ISBL6,rhlev  ${1}.nc
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
  [ -d fixed ] || mkdir fixed
  [ -f fixed/fixed_orography.nc ] && return
  local GVAR=HGT_3_SFC_10
  local NVAR=orog
  local CVAR="lon_3,lat_3,lv_ISBL2,lv_ISBL6"
  $G2NC $1 -u time -e grb -v ${CVAR},${GVAR} > /dev/null 2>&1
  $GRNM -d lat_3,lat -d lon_3,lon -d lv_ISBL2,lev -d lv_ISBL6,rhlev \
        -v lat_3,lat -v lon_3,lon -v lv_ISBL2,lev -v lv_ISBL6,rhlev ${1}.nc
  $GRNM -v $GVAR,$NVAR ${1}.nc
  mv ${1}.nc fixed/fixed_orography.nc
}

echo
echo "This scripts prepares a dataset downloaded from dss.ucar.edu, either"
echo
echo  "ds082.0 NCEP FNL Operational Model Global Tropospheric Analyses,"
echo  "        July 1976 to April 1997"
echo  "ds083.0 NCEP FNL Operational Model Global Tropospheric Analyses,"
echo  "        April 1997 to June 2007"
echo  "ds083.2 NCEP FNL Operational Model Global Tropospheric Analyses,"
echo  "        continuing from July 1999"
echo
echo "transforming from the original GRIB1 or GRIB2 format to NetCDF format"
echo "for the RegCM ICBC to prepare input for the RegCM model."
echo
if [ ! -f datafiles.tar ]
then
  echo Cannot find the downloaded archive file \'datafiles.tar\'
  echo Please execute this script in the directory where this file is.
  exit 1
fi

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
  exit 2
fi

echo "Extracting..."
tar xvf datafiles.tar > /dev/null

echo "Extracting fixed orography."
orog `ls -1 fnl_* | head -1`

echo "Converting to netCDF."
fnlsg=`ls -1 fnl_*`
for file in $fnlsg
do
  transform $file
  rm -f $file
done

echo "Done"
exit 0
