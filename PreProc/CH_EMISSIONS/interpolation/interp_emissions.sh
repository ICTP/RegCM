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
case $CHEMTYPE in
  "DUST")
  echo "WARNING: no need to interpolate emissions for $CHEMTYPE simulation, abort computation"
  exit 0
  ;;    
  "DU12")
  echo "WARNING: no need to interpolate emissions for $CHEMTYPE simulation, abort computation"
  exit 0
  ;;    
  "SSLT")
  echo "WARNING: no need to interpolate emissions for $CHEMTYPE simulation, abort computation"
  exit 0
  ;;
  "SULF")
  SPECIELIST=(SO2)
  ;;
  "CARB")
  SPECIELIST=(BC OC)
  ;;
  "SUCA")
  SPECIELIST=(BC OC SO2)
  ;;
  "SUCE")
  SPECIELIST=(BC OC SO2)
  ;;
  "AERO")
  SPECIELIST=(BC OC SO2)
  ;;
  "CBMZ")  
  SPECIELIST=(ALD2 AONE BC C2H6 CH3OH CH4 CO ETH HCHO NH3 NOx OC OLEI OLET PAR RCOOH SO2 TOL XYL ISOP_BIO) 
  ;;
  "DCCB")  
  SPECIELIST=(ALD2 AONE BC C2H6 CH3OH CH4 CO ETH HCHO NH3 NOx OC OLEI OLET PAR RCOOH SO2 TOL XYL ISOP_BIO) 
  ;;
  *) 
  echo "Cannot find chemistry simulation type: \"${CHEMTYPE}\""
  echo "Please  check out your namelist"
  exit 1
  ;;
esac  

# set emissions file list
n=${#SPECIELIST[@]}
file_list=()
for (( i = 0; i < $n; i++ )) ; do
    file_list+=(`ls ${data_dir}/*_"${SPECIELIST[i]}"[._]*nc 2> /dev/null`) 
done

}
######################################################################



######################################################################

function init_timeframe
{

  # set model starting date for emission processing
  START=`cat $NAMELIST | grep gdate1 | cut -d "=" -f 2 | tr "'" " " | \
           cut -d "," -f 1 | sed -e 's/ //g' `
  START=${START:0:4}-${START:4:2}-${START:6:2}T00:00:00

  # set model ending date
  END=`cat $NAMELIST | grep gdate2 | cut -d "=" -f 2 | tr "'" " " | \
           cut -d "," -f 1 | sed -e 's/ //g' `
  END=${END:0:4}-${END:4:2}-${END:6:2}T00:00:00


}
######################################################################




######################################################################

function emissions_interpolate
{

    # set file list
    set_filelist

    # create weights from first specie file
    $CDO gencon,$REGCM_dir/${DOMNAME}_grid.nc -setgrid,${file_list[0]} \
	${file_list[0]} remapweights.nc

    # interpolate all species onto model grid 
    init_timeframe
    chfiles=""
    for file in ${file_list[*]}; do
	ofile=`basename $file`
	if [ ${VERBOSE} -ge 1 ] ; then
	    echo "Producing $ofile..."
	fi
	
	# adjust data to simulation time period
	$CDO $CDOOPTIONS seldate,$START,$END $file $out_dir/tfile
	
	# grid interpolation
	$CDO -O $CDOOPTIONS remap,$REGCM_dir/${DOMNAME}_grid.nc,remapweights.nc $out_dir/tfile $out_dir/$ofile

	chfiles+="$out_dir/$ofile "
    done    
    
    rm -f $out_dir/${DOMNAME}_CHEMISS.nc
#REMAP_EXTRAPOLATE="on"                                          
    $CDO $CDOOPTIONS -O merge $chfiles $out_dir/${DOMNAME}_CHEMISS.nc     

}
######################################################################


######################################################################

function cleanup
{

    if [ ${VERBOSE} -ge 1 ] ; then
	echo 'Cleanup...'
    fi
    for file in ${file_list[*]}
      do
      ofile=`basename $file`
      rm -f $out_dir/$ofile 
    done
    rm -f remapweights.nc $out_dir/tfile

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

