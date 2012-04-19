#!/bin/bash

CDO=`which cdo`
NCRCAT=`which ncrcat`
NCKS=`which ncks`
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

if [ -z "$NCRCAT" ]
then
  echo "Cannot find a ncrcat executable in Your path."
  echo
  echo 'Please install nco (http://nco.sourceforge.net)'
  echo 'If on CINECA SP6, You may want to:'
  echo 
  echo '            module load autoload nco'
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

if [ $# -lt 2 ]
then
  echo "Not enough arguments."
  echo ""
  echo "RegCM averager example usage:"
  echo ""
  echo "       `basename $0` timespan inputfile1.nc [inputfile2.nc ...]"
  echo ""
  echo "where timespan defines a time window i.e. one in:"
  echo
  echo "         day month year seasD seasJ"
  echo
  echo "where seasD and seasJ identifies seasons starting in Dec or Jan"
  echo
  echo "Example:"
  echo "       `basename $0` month EUROPE_SRF.1990060100.nc EUROPE_SRF.1990060200.nc"
  echo ""
  exit 1
fi

tspan=$1
shift

declare -a infiles
infiles=($@)
nfiles=$#

internal=0
if [ "$tspan" == "day" ]
then
  operator="daymean"
  internal=1
elif [ "$tspan" == "month" ]
then
  operator="timmean"
  internal=1
elif [ "$tspan" == "year" ]
then
  operator="yearmean"
elif [ "$tspan" == "seasD" ]
then
  seastart="December"
  operator="seasmean"
elif [ "$tspan" == "seasJ" ]
then
  seastart="January"
  operator="seasmean"
else
  echo "Unknown timespan. I understand just the following:"
  echo
  echo "         day month year seasD seasJ"
  echo
  exit 1
fi

# Assumptions: Each files last timestep "belongs" to the next timeframe
# We just extract it and "use" in the next mean operation.
# The first timeframe in a run (not a restarted one) has the timestep
# at distance 0 from reference.

if [ $internal -eq 1 ]
then
  echo "Averaging..."
  for (( i = 0; i < $nfiles; i ++ ))
  do
    file=${infiles[$i]}
    xtimes=`$NCKS -v time -H -s "%g\n" $file`
    set -- $xtimes
    NT=$(( $# - 1 ))
    first=$1
    if [ $first -eq 0 ]
    then
      $NCKS -d time,0,$(( $NT - 1 )) $file tmp.nc
      $CDO $operator tmp.nc `basename ${file} .nc`_${operator}.nc
      rm -f tmp.nc
      $NCKS -d time,$NT $file tmp.nc
    else
      $NCKS -d time,0,$(( $NT - 1 )) $file tmp1.nc
      $NCRCAT tmp.nc tmp1.nc tmp2.nc && rm -f tmp.nc tmp1.nc
      $CDO $operator tmp2.nc `basename ${file} .nc`_${operator}.nc
      rm -f tmp2.nc
      $NCKS -d time,$NT $file tmp.nc
    fi
  done
  rm -f tmp.nc
  echo "Done"
else
# Annual or seasonal means.
# We neeed first to merge files, then process.
  echo "Merging..."
# Strip last timestp from last file
  lastf=$(( $nfiles - 2 ))
  $NCKS -d time,0,$(( $NT - 1 )) ${infiles[$lastf]} tmp.nc
  $NCRCAT ${infiles[0:$lastf-1]} tmp.nc tmp1.nc && rm tmp.nc
  echo "Averaging..."
  CDO_SEASON_START=$seastart $CDO $operator tmp1.nc ${base}_${operator}.nc
  rm -f tmp1.nc
fi

exit 0
