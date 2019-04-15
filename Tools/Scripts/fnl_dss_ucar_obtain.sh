#!/bin/bash
#################################################################
# Bash Script to retrieve online Data files of 'ds083.2',
# This script uses 'wget' to download data.
#
# Highlight this script by Select All, Copy and Paste it into a file;
# make the file executable and run it on command line.
#
# You need pass in your password as a parameter to execute
# this script; or you can set an evnironment variable RDAPSWD
# if your Operating System supports it.
#
# Contact baseball@ucar.edu (Gregg Walters) for further assistance.
#################################################################

LOGIN="SET_YOUR_LOGIN"
PASSWORD="SET_YOUR_PASS"

if [ "X$PASSWORD" = "XSET_YOUR_PASS" -o "X$LOGIN" = "XSET_YOUR_LOGIN" ]
then
  echo "Please edit this script and set your login and password for " \
       "dss.ucar.edu in it"
  exit 1
fi

if [ $# -ne 3 ]
then
  echo
  echo "Usage : `basename $0` year month day"
  echo
  echo with year in the format YYYY, month in the format MM and day in DD
  echo
  echo Example: `basename $0` 2001 01 02
  exit 1
fi

year=$1
month=$2
day=$3

v=`wget -V |grep 'GNU Wget ' | cut -d ' ' -f 3`
a=`echo $v | cut -d '.' -f 1`
b=`echo $v | cut -d '.' -f 2`
xc=$((100 * $a + $b))
if [ $xc -gt 109 ]
then
  opt='wget --no-check-certificate'
else
  opt='wget'
fi

opt1='-O /dev/null --save-cookies auth.dss_ucar_edu --post-data'
opt2="email=$LOGIN&passwd=$PASSWORD&action=login"
$opt $opt1="$opt2" https://dss.ucar.edu/cgi-bin/login
opt1="-N --load-cookies auth.dss_ucar_edu"
opt2="$opt $opt1 http://dss.ucar.edu/dsszone/ds083.2/"

gribf=grib1
[ $year -gt 2007 ] && gribf=grib2
[ $year -eq 2007 -a $month -eq 12 -a $day -gt 6 ] && gribf=grib2

if [ $year -eq 2007 -a $month -eq 12 -a $day -eq 6 ]
then
  filelist=" grib1/2007/2007.12/fnl_20071206_00_00 \
             grib1/2007/2007.12/fnl_20071206_06_00 \
	     grib2/2007/2007.12/fnl_20071206_12_00 \
	     grib2/2007/2007.12/fnl_20071206_18_00 "
else
  filelist=" $gribf/$year/$year.$month/fnl_$year$month${day}_00_00 \
             $gribf/$year/$year.$month/fnl_$year$month${day}_06_00 \
	     $gribf/$year/$year.$month/fnl_$year$month${day}_12_00 \
	     $gribf/$year/$year.$month/fnl_$year$month${day}_18_00 "
fi

for file in $filelist
do
 syscmd="$opt2$file"
 echo "$syscmd ..."
 $syscmd
done

rm -f auth.dss_ucar_edu
exit 0
