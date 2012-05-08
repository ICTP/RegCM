#! /bin/tcsh

#set echo verbose
###############################################################################
# This is the master script for pre-processing the CMIP5 emissions
###############################################################################

## global emissions file locations / output of creation ##
set data_dir=/home/netapp-clima/users/fsolmon/ARGO/RCP_EMGLOB_PROCESSED/creation/output
## grid of REGMC location
set REGCM_dir=/home/netapp-clima/users/fsolmon/ARGO/RCP_EMGLOB_PROCESSED/interpolation/grids
## grid of RCPs location
set CMIP5_dir=/home/netapp-clima/users/fsolmon/ARGO/RCP_EMGLOB_PROCESSED/interpolation/grids
#output_directory
set out_dir=/home/netapp-clima/users/fsolmon/ARGO/RCP_EMGLOB_PROCESSED/interpolation/output


###############################################################################
#   INTERPOLATION MAP GENERAL
###############################################################################


if ( ! -d $out_dir ) mkdir $out_dir
rm -f $out_dir/*.nc

#GENERATE WEIGHTS
set file_list=`ls ${data_dir}/*.nc`
echo $file_list
cdo gencon,$REGCM_dir/REGCM_grid.nc  -setgrid,$CMIP5_dir/CMIP5_grid.nc $file_list[1] remapweights.nc
foreach ifile ($data_dir/*.nc)
  set ofile=`basename $ifile`
  cdo remap,$REGCM_dir/REGCM_grid.nc,remapweights.nc $ifile $out_dir/$ofile
end

# here the final naming has to be interactive with regcm.in
cdo merge $out_dir/*.nc $out_dir/CHEMISS.nc 
exit

