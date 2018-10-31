#! /bin/csh

set year = 1975
set ymax = 2004
#1.here selecting each year only because doing this it needs less computational memory 
#2.we want to consider only precipitation > 1 mm for the r95 calculation (pre_1975-2004OK.nc is the daily precip in mm)
while ( $year <= $ymax )

cdo -selyear,$year pre_1975-2004OK.nc pre_$year\OK.nc
ncap2 -v -s 'pre2[time,lat,lon]=0.0; where (pre>1) pre2=pre' pre_$year\OK.nc pre2_$year\OK.nc
@ year++
end

ncrcat pre2_*\OK.nc pre2_1975-2004OK.nc
cdo timmin pre2_1975-2004OK.nc minfile_1975-2004.nc
cdo timmax pre2_1975-2004OK.nc maxfile_1975-2004.nc
cdo timpctl,95 pre2_1975-2004OK.nc minfile_1975-2004.nc maxfile_1975-2004.nc 95pctlref_1975-2004.nc
#this because "eca_r95ptot" needs as input files, files with the same number of timesteps:
cdo duplicate,10800 95pctlref_1975-2004.nc 95pctl_1975-2004.nc
cdo eca_r95ptot pre2_1975-2004OK.nc 95pctl_1975-2004.nc r95ptot1_1975-2004.nc 

#for the future: do the same, except that you don't have to calculate the "95pctl.nc", since you should use the reference one!
