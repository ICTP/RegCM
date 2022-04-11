#!/bin/bash



######################################################################

function print_header
{

    cat <<EOF


  This is `basename $0` part of the RegCM version 4
   GIT Revision: `git describe --abbrev=4 --dirty --always --tags`

 : this run start at    : `date '+%Y-%m-%d %H:%M:%S%z'` 
 : it is submitted by   : `whoami`
 : it is running on     : $HOST
 : in directory         : $PWD 



EOF
}


######################################################################

function help
{
    more <<EOF

    NAME: `basename $0` 

    SYNOPSIS: `basename $0` [-i] namelist [options]

    MANDATORY ARGUMENTS:
    namelist: set namelist

    OTHER OPTIONS
    -i (, --input=) namelist: set namelist [for compatibility reasons]
    -e (, --emissions=) name: set path to emissions inventories [default: RCP database]
    -v (, --verbose=) level: set verbose level (values are 0; 1 [default]; 2, 3; 4)
    -c (, --cdo=) path: set cdo path name if unusual
    -o (, --output=) outout: set path to interpolated emissions file [default: "dirter" directory]
    -h (,--help): print this help 

    
EOF
exit 0
}
######################################################################



######################################################################

function get_args
{

    while [ "$#" -ge 1 ] 
      do
      case $1 in 
	  "-h"|"--help")   
	  help;
          exit;  
	  ;;
	  "-v"|"--verbose=")
	  if [ -n "$2" ]; then 
	      shift 1 ; VERBOSE=$1 ; 
          fi
	  ;;
	  "-i"|"--input=")
	  if [ -n "$2" ]; then 
	      shift 1 ; NAMELIST=$1 ;  
	  fi
	  ;;
	  "-e"|"--emissions=")
	  if [ -n "$2" ]; then 
	      shift 1 ; data_dir=$1 ;  
	  fi
	  ;;
	  "-o"|"--output=")
	  if [ -n "$2" ]; then 
	      shift 1 ; out_dir=$1 ;  
	  fi
	  ;;
	  "-c"|"--cdo=")
	  if [ -n "$2" ]; then 
	      shift 1 ; CDO=$1 ;  
          fi
	  ;;
  	  -*)
	  echo -e "ERROR: unknown option. Use -h for help\n"; exit;
	  ;;
	  *)
          NAMELIST=$1 ;  
	  ;;
       esac
       shift
       done


# set default values
VERBOSE=${VERBOSE:=1}
CDO=${CDO:=`which cdo 2> /dev/null`}


}
######################################################################



######################################################################

function test_settings
{

# set verbose options (default is 1)
case ${VERBOSE} in
    "0")
    CDOOPTIONS="-s"
    ;;
    "2")
    set -v
    ;;
    "3")
    set -x;
    ;;
    "4")
    set -x;
    CDOOPTIONS="-v"
    ;;
esac


# test namelist
if [ ! -f "${NAMELIST}" ] ; then
  echo "Cannot find namelist ${NAMELIST}: no such file or directory."
  echo "Please set it with -i option. Use -h for help"
  exit 1
else
  # set domain name
  DOMNAME=`cat $NAMELIST | grep domname | grep -v coarse_domname | \
           cut -d "=" -f 2 | tr "'" " " | \
           cut -d "," -f 1 | sed -e 's/ //g' `

  # set model simulation input directory
  REGCM_dir=`cat $NAMELIST | grep dirter | cut -d "=" -f 2 | tr "'" " " | \
           cut -d "," -f 1 | sed -e 's/ //g' `
  out_dir=${out_dir:="${REGCM_dir}"}

  # set emission inventories path (default directory is global RCP directory)
  EMISSDIR=`cat $NAMELIST | grep inpglob | cut -d "=" -f 2 | tr "'" " " | \
           cut -d "," -f 1 | sed -e 's/ //g' `
  data_dir=${data_dir:="${EMISSDIR}/RCP_EMGLOB_PROCESSED/iiasa/"}

  # set model chemistry type
  CHEMTYPE=`cat $NAMELIST | grep chemsimtype | cut -d "=" -f 2 | tr "'" " " | \
           cut -d "," -f 1 | sed -e 's/ //g' `
fi

# test CDO availability
if [  x${CDO}x == "xx" ] ; then
  echo "Cannot find a cdo executable in your path."
  echo
  echo 'Please install cdo (https://code.zmaw.de/projects/cdo)'
  echo
  exit 1
fi

# test model domain name
if [ x${DOMNAME}x == "xx" ] ; then 
  echo "model input path is not set." 
  echo "Please check out your namelist."
  exit 1
fi

# test model simulation path
if [ ! -d "${REGCM_dir}" ] ; then
  echo "Cannot find model input path ${REGCM_dir}: no such file or directory."
  echo "Please check out your namelist."
  exit 1
fi

# test emission data path
if [ ! -d "${data_dir}" ] ; then
  echo "Cannot find emission data path ${data_dir}: no such file or directory."
  echo "Please check out your namelist or rerun -d option. Use -h for help"
  exit 1
fi

# test interpolated data path
if [ ! -d "${out_dir}" ] ; then
  echo "Cannot find interpolated emission data path ${out_dir}: no such file or directory."
  echo "Please  check out your namelist or rerun -o option. Use -h for help"
  exit 1
fi

if [ ${VERBOSE} -ge 1 ] ; then
    echo "EMISSION INTERPOLATION SCRIPT READY."
fi


}
######################################################################


######################################################################

function set_filelist
{

# define species for chemistry mode

# set emissions file list
file_list=()
file_list=(`ls ${data_dir}/soil_erodibility_factor.nc 2> /dev/null`) 


}
######################################################################



######################################################################

######################################################################




######################################################################

function emissions_interpolate
{

    # set file list
    set_filelist

    # create weights from first specie file
    $CDO gencon,$REGCM_dir/${DOMNAME}_grid.nc -setgrid,soil_erodibility_factor.nc \
	soil_erodibility_factor.nc remapweights.nc
    #$CDO -O $CDOOPTIONS remap,$REGCM_dir/${DOMNAME}_grid.nc,remapweights.nc ${file_list[0]} $REGCM_dir/${DOMNAME}_DUST.nc
    $CDO -O $CDOOPTIONS remap,$REGCM_dir/${DOMNAME}_grid.nc,remapweights.nc soil_erodibility_factor.nc $REGCM_dir/${DOMNAME}_DUST_TMP1.nc


    $CDO gencon,$REGCM_dir/${DOMNAME}_grid.nc -setgrid,aeolian_roughness_lenght.nc \
        aeolian_roughness_lenght.nc remapweights.nc

    $CDO -O $CDOOPTIONS remap,$REGCM_dir/${DOMNAME}_grid.nc,remapweights.nc aeolian_roughness_lenght.nc $REGCM_dir/${DOMNAME}_DUST_TMP2.nc

    # merge averything
    $CDO $CDOOPTIONS  merge $REGCM_dir/${DOMNAME}_DUST_TMP1.nc $REGCM_dir/${DOMNAME}_DUST_TMP2.nc $REGCM_dir/${DOMNAME}_DUSTPARAM.nc

}
######################################################################


######################################################################

function cleanup
{

    if [ ${VERBOSE} -ge 1 ] ; then
	echo 'Cleanup...'
    fi
    if [ ${VERBOSE} -ge 1 ] ; then
	echo 'Done'
    fi

}
######################################################################




###############################################################################

# This is the master script for pre-processing the emissions

###############################################################################


get_args "$@" 

print_header

test_settings

emissions_interpolate

cleanup

exit 0

