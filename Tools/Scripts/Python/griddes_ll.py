#!/usr/bin/env python3

import sys

if len(sys.argv) < 2:
    print('Need netcdf filename.')
    sys.exit(-1)

from netCDF4 import Dataset
from math import floor, ceil, fmod

ds = Dataset(sys.argv[1],'r')

try:
    res = float(sys.argv[2])
except:
    res = 0.5

try:
    lat = ds.variables['xlat']
    lon = ds.variables['xlon']
except:
    try:
        lat = ds.variables['lat']
        lon = ds.variables['lon']
    except:
        print('No recognizable coordinates.')
        sys.exit(-1)


la1 = res*floor(min(lat[1,:])/res)
la2 = res*ceil(max(lat[-1,:])/res)
m1 = lon[1,1]
m2 = lon[-1,1]
if m1*m2 < 0:
    m1 = fmod(m1 + 360.0,360.0)
    m2 = fmod(m2 + 360.0,360.0)
    lo1 = res*floor(min(m1,m2)/res)
else:
    lo1 = res*floor(min(m1,m2)/res)
m1 = lon[1,-1]
m2 = lon[-1,-1]
if m1*m2 < 0:
    m1 = fmod(m1 + 360.0,360.0)
    m2 = fmod(m2 + 360.0,360.0)
    lo2 = res*floor(max(m1,m2)/res)
else:
    lo2 = res*floor(max(m1,m2)/res)

ds.close( )

nlat = int((la2-la1)/res) + 1
if lo2 > lo1:
    nlon = int((lo2-lo1)/res) + 1
else:
    lo2 = lo2 + 360
    nlon = int((lo2-lo1)/res) + 1

print('gridtype  = lonlat')
print('gridsize  = ',nlon*nlat)
print('datatype  = float')
print('xsize     = ',nlon)
print('ysize     = ',nlat)
print('xname     = lon')
print('xlongname = "longitude"')
print('xunits    = "degrees_east"')
print('yname     = lat')
print('ylongname = "latitude"')
print('yunits    = "degrees_north"')
print('xfirst    = ',lo1)
print('xinc      = ',res)
print('yfirst    = ',la1)
print('yinc      = ',res)

sys.exit(0)
