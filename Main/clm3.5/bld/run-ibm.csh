#! /usr/bin/csh -f

#-----------------------------------------------------------------------
## IBM
##------------
##
## This is an example script to build and run the default CLM configuration
## on an IBM SP.  This is setup to run on NCAR's machine bluevista.
##
## To submit this on bluevista:
##
## bsub < run-ibm.csh
##
##

## Setting LSF options for batch queue submission.
#BSUB -a poe                    # use poe for multiprocessing
#BSUB -x                        # exclusive use of node (not_shared)
## Number of tasks and tasks per node (CHANGE THIS IF YOU TURN smp on)
#BSUB -n 16                     # total number of MPI-tasks (processors) needed
#BSUB -R "span[ptile=16]"       # max number of tasks (MPI) per node
#BSUB -o out.%J                 # output filename
#BSUB -e out.%J                 # error filename
#BSUB -q regular                # queue
#BSUB -W 0:10                   # wall clock limit
#BSUB -P xxxxxxxx               # Project number to charge to (MAKE SURE YOU CHANGE THIS!!!)

## POE Environment.  Set these for interactive jobs.  They're ignored in batch submission.
## MP_NODES is the number of nodes.  
setenv MP_NODES 1
setenv MP_TASKS_PER_NODE 16
setenv MP_EUILIB us
setenv MP_RMPOOL 1

unsetenv MP_PROCS

setenv MP_STDINMODE 0

# should be set equal to (CPUs-per-node / tasks_per_node)
# Only activated if smp=on below
setenv OMP_NUM_THREADS 4

## suggestion from Jim Edwards to reintroduce XLSMPOPTS on 11/13/03
setenv XLSMPOPTS "stack=256000000"
setenv AIXTHREAD_SCOPE S
setenv MALLOCMULTIHEAP true
setenv OMP_DYNAMIC false
## Do our best to get sufficient stack memory
limit stacksize unlimited

## netCDF stuff
setenv INC_NETCDF /usr/local/include
setenv LIB_NETCDF /usr/local/lib64/r4i4

## ROOT OF CLM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CLM distribution.
## (the root directory contains the subdirectory "src")
set clmroot   = /fis/cgd/.......       # (MAKE SURE YOU CHANGE THIS!!!)

## ROOT OF CLM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CLM distribution.
setenv CSMDATA /fs/cgd/csm/inputdata/lnd/clm2       # (MAKE SURE YOU CHANGE THIS!!!)

## Configuration settings:
set spmd     = on       # settings are [on   | off       ] (default is off)
set smp      = off      # settings are [on   | off       ] (default is off)
set maxpft   = 4        # settings are 4->17               (default is 4)
set rtm      = off      # settings are [on   | off       ] (default is off)   
## IF YOU CHANGE ANY OF THE CONFIGURATION SETTINGS -- DELETE THE $blddir/config_cache.xml AND RESUBMIT
## (see below)
#--------------------------------------------------------------------------------------------

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CLM configuration scripts.
set case    = clmrun
set wrkdir  = /ptmp/$LOGNAME
set blddir  = $wrkdir/$case/bld
set rundir  = $wrkdir/$case
set cfgdir  = $clmroot/bld
set usr_src = $clmroot/bld/usr.src

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## Build (or re-build) executable
set flags = "-maxpft $maxpft -rtm $rtm -usr_src $usr_src"
if ($spmd == on ) set flags = "$flags -spmd"
if ($spmd == off) set flags = "$flags -nospmd"
if ($smp  == on ) set flags = "$flags -smp"
if ($smp  == off) set flags = "$flags -nosmp"

set fsurdat="surfdata_64x128_c070501.nc"

echo "cd $blddir"
cd $blddir                  || echo "cd $blddir failed" && exit 1
## Check if config_cache.xml file exists -- if so just run make -- if NOT then run configure.
## IF YOU CHANGE ANY OF THE CONFIGURATION SETTINGS -- DELETE THE $blddir/config_cache.xml AND RESUBMIT
#--------------------------------------------------------------------------------------------
if ( ! -f $blddir/config_cache.xml ) then
    echo "flags to configure are $flags"
    $cfgdir/configure $flags    || echo "configure failed" && exit 1
    echo "Building CLM in $blddir ..."
    gmake -j8 >&! MAKE.out      || echo "CLM build failed: see $blddir/MAKE.out" && exit 1
else
    echo "Re-building CLM in $blddir ..."
    rm -f Depends
    gmake -j8 >&! REMAKE.out      || echo "CLM build failed: see $blddir/REMAKE.out" && exit 1
endif

## Create the namelist
cd $rundir                      || echo "cd $blddir failed" && exit 1

cat >! lnd.stdin << EOF
 &clm_inparm
 caseid         = '$case'
 ctitle         = '$case'
 finidat        = '$CSMDATA/inidata_3.1/offline/clmi_0000-01-01_064x128_c070403.nc'
 fsurdat        = "$CSMDATA/surfdata/$fsurdat"
 fatmgrid       = "$CSMDATA/griddata/griddata_64x128_060829.nc"
 fatmlndfrc     = "$CSMDATA/griddata/fracdata_64x128_USGS_070110.nc"
 fpftcon        = '$CSMDATA/pftdata/pft-physiology.c070207'
 offline_atmdir = "$CSMDATA/NCEPDATA.Qian.T62.c051024"
 frivinp_rtm    = '$CSMDATA/rtmdata/rdirc.05.061026'
 nsrest         =  0
 nelapse        =  48
 dtime          =  1800
 start_ymd      =  20030101
 start_tod      =  0
 irad           = -1
 wrtdia         = .true.
 mss_irt        =  0
 hist_dov2xy    = .true.
 hist_nhtfrq    =  -24
 hist_mfilt     =  1
 hist_crtinic   = 'MONTHLY'
 /
 &prof_inparm
 /
EOF

## Run CLM 

cd $rundir                    || echo "cd $rundir failed" && exit 1
echo "running CLM in $rundir"

setenv LID "`date +%y%m%d-%H%M%S`"

if ($spmd == on) then

  mpirun.lsf $blddir/clm >&! clm.log.$LID       || echo "CLM run failed" && exit 1
else 
  $blddir/clm  >&! clm.log.$LID                 || echo "CLM run failed" && exit 1
endif

exit 0
