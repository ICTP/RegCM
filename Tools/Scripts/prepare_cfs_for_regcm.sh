#!/bin/bash

WORKDIR=/scratch/ggiulian/test/cfsdata/GRIB
CFSDIR=/scratch/ggiulian/test/cfsdata/CFS

run=00     # The run (one in 00-06-12-18)
emem=01    # The ensemble member (one in 01 02 03 04)

download=1 # To also download the data
do_pack=1  # Keep this to lower disk usage

DWN=`which wget`
if [ "X$DWN" == "X" ]
then
  echo "Need wget to download data"
  exit 1
fi
WG2=`which wgrib2`
if [ "X$WG2" == "X" ]
then
  echo "Need wgrib2 to convert GRIB2 to NetCDF"
  exit 1
fi

mkdir -p $WORKDIR
cd $WORKDIR

today=`date -u +%Y%m%d`
year=`date -u +%Y`
month=`date -u +%m`
day=`date -u +%d`

# TO TEST FOR A PARTICULAR DATE

#today=20170522
#year=2017
#month=05
#day=22

cmdline="-nc_table prc.table -nc4 -netcdf"
pgbf=pgbf
pgbanl=pgbanl
flxf=flxf

# Download data
if [ $download -eq 1 ]
then
  PROTO=http
  SITE=nomads.ncep.noaa.gov
  SOURCE=pub/data/nccf/com/cfs/prod/cfs
  DSITE=${PROTO}://${SITE}/${SOURCE}/cfs.${today}/${run}/6hrly_grib_${emem}
  echo $DSITE
  for (( d = 0; d < 32; d ++ ))
  do
    for (( h = 0; h < 24; h += 6 ))
    do
      dd=$(( $d * 24 + $h ))
      d1=`date -u +%Y%m%d%H \
	      --date "$year-$month-$day $run:00:00 UTC + $dd hours"`
      if [ $d -eq 0 -a $h -eq 0 ]
      then
        wget -c $DSITE/${pgbanl}.${emem}.$year$month$day$run.grb2
        wget -c $DSITE/${flxf}${d1}.${emem}.${today}${run}.grb2
        continue
      fi
      for file in ${flxf} ${pgbf}
      do
        wget -c $DSITE/${file}${d1}.${emem}.${today}${run}.grb2
      done
    done
  done
fi

# Change data format

echo '# PROCTABLE'                                     >  prc.table
echo '$lev_type  100:level:pressure level:mb:0.01'     >> prc.table
echo '$nlev      20'                                   >> prc.table
echo '$levs  1000 975 950 925 900 850 800 700 600 500' >> prc.table
echo '       400 300 250 200 150 100 70 30 20 10'      >> prc.table
if [ $do_pack -eq 1 ]
then
  echo 'UGRD:*:uwnd:short:-140:175'                    >> prc.table
  echo 'VGRD:*:vwnd:short:-140:175'                    >> prc.table
  echo 'TMP:*:air:short:137:362.5'                     >> prc.table
  echo 'HGT:*:hgt:short:-1500:90000'                   >> prc.table
  echo 'RH:*:rhum:short:-25:125'                       >> prc.table
  echo 'TMP:surface:sst:short:250:350'                 >> prc.table
else
  echo 'UGRD:*:uwnd'                                   >> prc.table
  echo 'VGRD:*:vwnd'                                   >> prc.table
  echo 'TMP:*:air'                                     >> prc.table
  echo 'HGT:*:hgt'                                     >> prc.table
  echo 'RH:*:rhum'                                     >> prc.table
  echo 'TMP:surface:sst'                               >> prc.table
fi
echo '#EOF'                                            >> prc.table

# Prepare
cat ${flxf}*.grb2 > flx.grb2
cat ${pgbanl}.${emem}.$year$month$day$run.grb2 ${pgbf}*.grb2 > plv.grb2
cat ${flxf}${today}${run}.${emem}.${today}${run}.grb2 >> \
	${pgbanl}.${emem}.$year$month$day$run.grb2

START_TIME=${today}${run}
PLVDEST=${CFSDIR}/${START_TIME}/${emem}/PLEV
SSTDEST=${CFSDIR}/${START_TIME}/${emem}/SST
mkdir -p $PLVDEST
mkdir -p $SSTDEST

# Process
$WG2 -match 'TMP:surface'        \
  $cmdline ${SSTDEST}/sst.${START_TIME}.nc  flx.grb2
$WG2 -match ':TMP' -match 'mb:'  \
  $cmdline ${PLVDEST}/air.${START_TIME}.nc  plv.grb2
$WG2 -match ':RH' -match 'mb:'   \
  $cmdline ${PLVDEST}/rhum.${START_TIME}.nc plv.grb2
$WG2 -match ':UGRD' -match 'mb:' \
  $cmdline ${PLVDEST}/uwnd.${START_TIME}.nc plv.grb2
$WG2 -match ':VGRD' -match 'mb:' \
  $cmdline ${PLVDEST}/vwnd.${START_TIME}.nc plv.grb2
$WG2 -match ':HGT' -match 'mb:'  \
  $cmdline ${PLVDEST}/hgt.${START_TIME}.nc  plv.grb2

# Cleanup
rm -f flx.grb2 plv.grb2 prc.table
rm -f ${flxf}*.grb2 ${pgbf}*.grb2
