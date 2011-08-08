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

listfirst()
{
  fnlsg=`ls -1 fnl_* | head -1`
  fnlsnc=$fnlsg.nc
}

listall ()
{
  fnlsg=`ls -1 fnl_*`
  fnlsnc=""
  for file in $fnlsg
  do
    fnlsnc="$fnlsnc $file.nc"
  done
}

transform ()
{
  echo "Doing ${2}"
  $G2NC $fnlsg -u time -e grb -v lon_3,lat_3,lv_ISBL2,${1} > /dev/null 2>&1
  for file in $fnlsnc
  do
    mv $file ${2}_$file
    $GRNM -v ${1},${2} ${2}_$file
    $GRNM -d lat_3,lat ${2}_$file
    $GRNM -d lon_3,lon ${2}_$file
    $GRNM -d lv_ISBL2,plev ${2}_$file
    $GRNM -v lat_3,lat ${2}_$file
    $GRNM -v lon_3,lon ${2}_$file
    $GRNM -v lv_ISBL2,plev ${2}_$file
  done
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

listfirst
transform HGT_3_SFC_10 orog
mv *.nc fixed_orograpy.nc

listall
transform TMP_3_SFC_10 ts
transform PRES_3_SFC_10 ps
transform HGT_3_ISBL_10 hga
transform R_H_3_ISBL_10 rha
transform TMP_3_ISBL_10 ta
transform U_GRD_3_ISBL_10 ua
transform V_GRD_3_ISBL_10 va

rm -f $fnlsg
exit 0
