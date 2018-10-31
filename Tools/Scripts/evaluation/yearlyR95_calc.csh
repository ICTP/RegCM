#! /bin/csh

cdo timmin pre2_1975-2004OK.nc minfile_1975-2004.nc
cdo timmax pre2_1975-2004OK.nc maxfile_1975-2004.nc
cdo timpctl,95 pre2_1975-2004OK.nc minfile_1975-2004.nc maxfile_1975-2004.nc 95pctlref_1975-2004.nc

set year = 1975
set ymax = 2004
#1.here selecting each year only because doing this it needs less computational memory 
#2.we want to consider only precipitation > 1 mm for the r95 calculation (pre_1975-2004OK.nc is the daily precip in mm)
while ( $year <= $ymax )

cdo -selyear,$year pre_1975-2004OK.nc pre_$year\OK.nc
ncap2 -v -s 'pre2[time,lat,lon]=0.0; where (pre>1) pre2=pre' pre_$year\OK.nc PRE2_$year.nc
cdo sub PRE2_$year.nc PRE2_$year.nc zerofile_$year.nc
cdo add zerofile_$year.nc 95pctlref_1975-2004.nc 95pctl_refx$year.nc
cdo eca_r95ptot PRE2_$year.nc 95pctl_refx$year.nc r95ptot1_$year.nc
@ year++
end

ncrcat r95ptot1_* r95ptot1_1975-2004.nc
#
#to calculate the trend of the R95 time series:

cdo trend r95ptot1_1975-2004.nc atrend.nc r95trend.nc
