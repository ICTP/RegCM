#! /bin/csh

#annual
cdo -L -z zip -f nc4 -ydrunmean,5  tmax_1979-2017.nc tmax-normy.nc
cdo -L -z zip -f nc4 -duplicate,38 tmax-normy.nc tmax-normy-dup.nc
cdo -L -z zip -f nc4 -eca_hwdi,5,5 tmax_1979-2017.nc tmax-normy-dup.nc heatw55.nc 
cdo -L -z zip -f nc4 -eca_hwdi,6,6 tmax_1979-2017.nc tmax-normy-dup.nc heatw66.nc 
cdo -L -z zip -f nc4 -eca_hwdi,7,7 tmax_1979-2017.nc tmax-normy-dup.nc heatw77.nc
