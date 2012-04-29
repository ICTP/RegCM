#!/bin/bash

#############################################################
# Load modules                                              #
#############################################################

. /etc/profile.d/modules.sh

# --- remove all defined modules ---
module --long list >& .log
lstlen=`cat .log | wc -l`
lstlen=$((lstlen-2))
lstmod=`cat .log | tail -n $lstlen | awk '{print $1}'`
for mod in $lstmod
do
  module rm $mod
done

# --- define required modules ---
module load cdo
module load nco

#############################################################
# Parameters                                                #
#############################################################

regcmin="regcm.in_CAS50km"
regcmout="regcmout.txt"
romsin="cas.in"

#############################################################
# RegCM                                                     #
#############################################################

# get date stamp
dstamp=`date +"%d-%m-%y_%H:%M"`

# backup original files
cp $romsin ${romsin/.in/}_$dstamp.in
mv $regcmout ${regcmout/.txt/}_$dstamp.txt
cp $regcmin ${regcmin}_$dstamp

# get latest SAV file name and query new restart date
fsav=`ls -al output/*_SAV.* | tail -n 1 | awk '{print $9}'`
#fsav=`ls -al output/*_TMPSAV.* | tail -n 1 | awk '{print $9}'`
mdate1=`echo "$fsav" | awk -F. '{print $2}'` 
prefix=`cat $regcmin | grep domname | awk -F\' '{print $2}'`
echo "[debug] -- model prefix is $prefix"
echo "[debug] -- new restart date is $mdate1"

# link TMPSAV as SAV
#if [ -f "output/${prefix}_SAV.$mdate1" ]; then
#  echo "[debug] -- TMPSAV file is already linked as ${prefix}_SAV.$mdate1"
#else
#  cd output
#  rm -f ${prefix}_SAV.$mdate1
#  ln -s ${prefix}_TMPSAV.$mdate1 ${prefix}_SAV.$mdate1
#  cd - >& /dev/null
#fi

# mdate1
# check restart time is same or not?
str1=`cat $regcmin | grep "mdate1" | tail -n 1`
val1=`echo "$str1" | awk -F= '{print $2}' | tr -d ' '`
val2=$mdate1
if [ "$val1" == "$val2," ]; then
  echo "[debug] -- the restart time is already changed. do not change it again!"
else
  str2=${str1/$val1/$val2,}
  cat $regcmin | sed "s/$str1/$str2/g" > .tmp
  mv .tmp $regcmin
  echo "[debug] -- mdate1 is changed to '$val2' in '$regcmin' file."
fi

# ifrest 
# check that the run is restarted before or not?
str1=`cat $regcmin | grep "ifrest"`
if [ -n "`echo $str1 | grep "true"`" ]; then
  echo "[debug] -- the simulation is restarted before. do not change it again!"
else
  val1=`echo "$str1" | awk -F= '{print $2}'`
  str2=${str1/$val1/ .true. ,}
  cat $regcmin | sed "s/$str1/$str2/g" > .tmp
  mv .tmp $regcmin
  echo "[debug] -- ifrest is changed to '.true.' in '$regcmin' file."
fi

#############################################################
# ROMS                                                      #
#############################################################

# get list of dates from the ROMS restart file
cdo -s showdate output/ocean_rst.nc | tr " " "\n" | grep - >&.tmp

# find time indices to split data
yy=${mdate1:0:4}
mm=${mdate1:4:2}
dd=${mdate1:6:2}
dstr="$yy-$mm-$dd"
lno=`awk -v dstr=$dstr '{if($1==dstr) print NR}' .tmp`

# create restart file
if [ $((lno-1)) -eq "0" ]; then
  echo "[debug] -- the restart file has just created. do not create it again!"
else
  mv output/ocean_rst.nc output/ocean_rst_$dstamp.nc
  ncks -d ocean_time,$((lno-1)) output/ocean_rst_$dstamp.nc output/ocean_rst.nc
  echo "[debug] -- restart file is created for time step '$lno'."
  echo "[debug] -- old one is saved as 'output/ocean_rst_$dstamp.nc'."
fi

# modify parameters
# NRREC
str1=`cat $romsin | grep "NRREC" | head -n 1`
str2=`echo $str1 | awk -F! '{print $2}'`
str1=${str1/\!$str2/""}
val1=`echo "$str1" | awk -F"==" '{print $2}' | tr -d ' '`
val2="-1"
if [ "$val1" == "$val2" ]; then
  echo "[debug] -- the restart record (NRREC) is already changed. do not change it again!"
else
  str2=${str1/$val1/$val2}
  cat $romsin | sed "s/$str1/$str2 \!$val1/g" > .tmp
  mv .tmp $romsin
  echo "[debug] -- NRREC is changed to '$val2' in '$romsin' file."
fi

# LDEFOUT
str1=`cat $romsin | grep "LDEFOUT" | head -n 1`
str2=`echo $str1 | awk -F! '{print $2}'`
str1=${str1/\!$str2/""}
val1=`echo "$str1" | awk -F"==" '{print $2}' | tr -d ' '`
val2="F"
if [ "$val1" == "$val2" ]; then
  echo "[debug] -- the LDEFOUT option is already changed. do not change it again!"
else
  str3="     `echo "$str1" | awk -F"==" '{print $1}' | tr -d ' '` == $val2 "
  cat $romsin | sed "s/$str1/$str3/g" > .tmp
  mv .tmp $romsin
  echo "[debug] -- LDEFOUT is changed to '$val2' in '$romsin' file."
fi

# ININAME
str1=`cat $romsin | grep "ININAME" | head -n 1`
str2=`echo $str1 | awk -F! '{print $2}'`
str1=${str1/\!$str2/""}
val1=`echo "$str1" | awk -F"==" '{print $2}' | tr -d ' '`
val2="output/ocean_rst.nc"
if [ "$val1" == "$val2" ]; then
  echo "[debug] -- the restart file (ININAME) is already changed. do not change it again!"
else
  str2=${str1/$val1/$val2}
  cat $romsin | sed "s:$str1:$str2:g" > .tmp
  mv .tmp $romsin
  echo "[debug] -- ININAME is changed to '$val2' in '$romsin' file."
fi

