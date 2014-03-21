#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Compute daily averages out of a RegCM output file"
  echo "Usage:"
  echo "        $0 regcmfile.nc [outfile.nc]"
  echo
  exit 1
fi

if [ ! -f $1 ]
then
  echo "RegCM file $1 not found."
  exit 1
else
  rcmfile=$1
fi

if [ $# -eq 3 ]
then
  outfile=$2
else
  outfile=`basename $rcmfile .nc`_day.nc
fi

times=(`ncks -C -v time -H -s "%f " $rcmfile`)
nt=${#times[@]}
freq=`echo "${times[1]}-${times[0]}" | bc | cut -d "." -f 1`
tpd=$(( 24 / $freq ))
ndays=$(( $nt/$tpd ))
ncdump -v time_bnds $rcmfile 2> /dev/null > /dev/null
not_time_bnds=$?

for (( iday = 0; iday < $ndays; iday ++ ))
do
  tstart=$(( $iday*$tpd ))
  tstop=$(( $tstart+$tpd-1 ))
  cday=`printf %03d $iday`
  tmp1="$$.xtmp$cday.nc"
  tmp2="$$.ytmp$cday.nc"
  tmp3="tmp$cday.nc"
  ncks -h -O -d time,$tstart,$tstop $rcmfile $tmp1
  ncra -h -O $tmp1 $tmp2
  if [ $not_time_bnds -eq 1 ]
  then
    tb=(`ncks -H -s "%f " -v time -C $tmp1`) && rm $tmp1
    bb1=`echo ${tb[0]}-$freq | bc`f
    bb2=`echo ${tb[-1]}`f
    ncap2 -h -O -s "defdim(\"time_bounds\",2); \
                 time_bnds[\$time,\$time_bounds] = 0.0f; \
                 time_bnds(:,:)={$bb1,$bb2}; \
                 time@bounds=\"time_bnds\";
                 time_bnds@units=time@units;
                 time_bnds@calendar=time@calendar;" $tmp2 $tmp3 && rm $tmp2
  else
    tb=(`ncks -H -s "%f " -v time_bnds -C $tmp1`) && rm $tmp1
    ncap2 -h -O -s "time_bnds(:,:)={${tb[0]},${tb[-1]}};" \
              $tmp2 $tmp3 && rm $tmp2
  fi
done
ncrcat -h -O tmp???.nc $outfile
rm -f tmp???.nc

exit 0
