#!/bin/bash

#SBATCH --account             CMPNS_ictpclim
#SBATCH --job-name            MED-12_POST
#SBATCH --mail-type           FAIL
#SBATCH --mail-user           ggiulian@ictp.it
#SBATCH --nodes               1
#SBATCH --ntasks-per-node     112
#SBATCH --partition           dcgp_usr_prod
#SBATCH --qos=qos_lowprio
#SBATCH --time                01:30:00

#
# CORDEX-CMIP6 post-processing of one month of RegCM output.
#
#   regcmcordex.sh DATADIR OUTPUT YYYYMMDDHH   process the RegCM files of a
#                                               month (e.g. 1980020100), then
#                                               remove them
#   regcmcordex.sh DATADIR OUTPUT finalize     assemble the archive files
#                                               left in segments, once the
#                                               whole experiment is processed
#
# Each month is written as segments (OUTPUT/.../<var>/.segments), and the
# archive-size files (1 year sub-daily, 5 years daily, 10 years monthly) are
# assembled as soon as their segments are all there. The months can be
# submitted in any order and run at the same time. With first_year and
# last_year set to the years of the experiment, its first and last files are
# also assembled without the finalize step.
#

source /leonardo/home/userexternal/ggiulian/modules_new

[ $# -lt 3 ] && exit 1

datadir=$1
output=$2
idate=$3

echo 'Input directory : '$datadir
echo 'Output directory : '$output
echo 'Processing date(s): '$idate

pycordex=/leonardo/home/userexternal/ggiulian/RegCM/Tools/Scripts/pycordexer
mail=ggiulian@ictp.it
domain=MED-12
global=ERA5
experiment=evaluation
ensemble=r1i1p1f1
notes="None"
proc=20
regcm_model=RegCM-ES
regcm_release=1.1
regcm_version_id=v1-r1
# Years of the experiment (leave empty if not known: run finalize at the end)
first_year=1980
last_year=2024

chunkargs=""
[ -n "$first_year" ] && chunkargs="$chunkargs --first-year $first_year"
[ -n "$last_year" ] && chunkargs="$chunkargs --last-year $last_year"

if [ "$idate" == "finalize" ]
then
  # Assemble the files still in segments (the edges of the experiment if the
  # years are not given, or chunks with months missing: these are reported)
  $pycordex/cordex_finalize.py $chunkargs $output
  exit $?
fi

allargs="-m $mail -d $domain -g $global -e $experiment -b $ensemble \
         -n "$notes" -o $output -p $proc --regcm-model-name $regcm_model \
         -r $regcm_release --regcm-version-id $regcm_version_id $chunkargs"

srffile=$datadir/*_SRF.${idate}*.nc
stsfile=$datadir/*_STS.${idate}*.nc
radfile=$datadir/*_RAD.${idate}*.nc
atmfile=$datadir/*_ATM.${idate}*.nc

srfvars=tas,pr,evspsbl,huss,hurs,ps,psl,sfcWind,uas,vas,clt,rsds,rlds
srfvars=$srfvars,ts,prc,prhmax,prsn,tauu,tauv,zmla,prw,rsus,rlus,hfss,hfls
srfvars=$srfvars,ua50m,ua100m,ua150m,va50m,va100m,va150m,ta50m,hus50m
srfvars=$srfvars,mrros,mrro,cape,cin,li,evspsblpot,z0,hfso,snm,mrsos
stsvars=prmean,psmean,tasmean,tasmax,tasmin,sfcWindmax,sundmean,wsgsmax
radvars=clwvi,clivi,rlut,rsut,rsdt,clh,clm,cll,cld
atmvars=ua,va,ta,hus,zg,wa,mrsol,mrso,snw,tsl,cli,clw

pids=""
err=""
$pycordex/pycordexer.py $allargs $srffile $srfvars & pids+="$! "
$pycordex/pycordexer.py $allargs $radfile $radvars & pids+="$! "
$pycordex/pycordexer.py $allargs $stsfile $stsvars & pids+="$! "
$pycordex/pycordexer.py $allargs $atmfile $atmvars & pids+="$! "

for p in $pids; do wait $p || err=$?; done
[[ $err ]] && exit -1

# The month is now in the segments (or in the assembled files): the RegCM
# files can be removed
rm $srffile $radfile $stsfile $atmfile
