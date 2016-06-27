#!/bin/bash

CDO=`which cdo`

if [ -z "$CDO" ]
then
  echo "Cannot find a cdo executable in Your path."
  echo
  echo 'Please install cdo (https://code.zmaw.de/projects/cdo)'
  echo
  exit 1
fi

###############################################################################
# This is the master script for pre-processing the CMIP5 emissions
###############################################################################

echo "EMISSION INTERPOLATION SCRIPT READY."

if [ $# -lt 1 ]
then
  echo "Not enough arguments."
  echo "Need input namelist file to work."
  echo "Usage: "
  echo "         "`basename $0` " regcm.in"
  echo
  exit 1
fi

INPFILE=$1
DOMNAME=`cat $INPFILE | grep domname | grep -v coarse_domname | \
	 cut -d "=" -f 2 | tr "'" " " | sed -e 's/ //g' -e 's/,//'`
RCMINPDIR=`cat $INPFILE | grep dirter | cut -d "=" -f 2 | tr "'" " " | \
           sed -e 's/ //g' -e 's/,//'`
EMISSDIR=`cat $INPFILE | grep inpglob | cut -d "=" -f 2 | tr "'" " " | \
           sed -e 's/ //g' -e 's/,//'`

## global emissions file locations / output of creation ##
data_dir="$EMISSDIR/BIONOX"
## grid of RCPs location
CMIP5_dir="$EMISSDIR/RCP_EMGLOB_PROCESSED/grids"
## grid of REGMC location
REGCM_dir="$RCMINPDIR"
#output_directory
out_dir="$RCMINPDIR"

###############################################################################
#   INTERPOLATION MAP GENERAL
###############################################################################

#GENERATE WEIGHTS
file_list=(`ls ${data_dir}/*.nc`)

$CDO gencon,$REGCM_dir/${DOMNAME}_grid.nc -setgrid,$CMIP5_dir/bionit_grid.nc \
      ${file_list[0]} remapweights.nc

chfiles=""
for file in ${file_list[*]}
do
  ofile=`basename $file`
  echo "Producing $ofile..."
  $CDO remap,$REGCM_dir/${DOMNAME}_grid.nc,remapweights.nc $file $out_dir/$ofile
  chfiles="$chfiles $out_dir/$ofile"
done

# here the final naming has to be interactive with regcm.in
rm -f $out_dir/${DOMNAME}_BIONEM.nc
$CDO merge $chfiles $out_dir/${DOMNAME}_BIONEM.nc 

echo 'Cleanup...'
for file in ${file_list[*]}
do
  ofile=`basename $file`
  rm -f $out_dir/$ofile
done
rm -f remapweights.nc
echo 'Done'
exit 0
