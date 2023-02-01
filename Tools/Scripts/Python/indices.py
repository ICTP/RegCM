#!/usr/bin/env python

from netCDF4 import Dataset, num2date
import numpy as np
from bisect import bisect
import sys

earth_radius = 3958.75
def get_distances(alat,alon,xlat,xlon):
    lat_dif = np.radians(alat - xlat)
    long_dif = np.radians(alon - xlon)
    sin_d_lat = np.sin(lat_dif / 2.)
    sin_d_long = np.sin(long_dif / 2.)
    step_1 = ((sin_d_lat ** 2) + (sin_d_long ** 2) * 
              np.cos(np.radians(alat)) * np.cos(np.radians(xlat)))
    step_2 = 2 * np.arctan2(np.sqrt(step_1), np.sqrt(1-step_1))
    dist = step_2 * earth_radius
    return dist

try:
    ncfile = Dataset(sys.argv[1], 'r')
except:
    if len(sys.argv) > 1:
        print('Cannot open data file '+sys.argv[1])
    else:
        print('Need at least a netcdf file name!')
    sys.exit(-1)

try:
    xlat = np.array((sys.argv[2],),np.float64)
    xlon = np.array((sys.argv[3],),np.float64)
except:
    if len(sys.argv) == 4:
        print('Cannot convert args to lat,lon values: '+sys.argv[2:])
    else:
        print('Need lat lon coordinates!')
    sys.exit(-2)

lats = ncfile.variables['xlat'][:]
lons = ncfile.variables['xlon'][:]
vtime = ncfile.variables['time']
xtime = vtime[:]
times = num2date(xtime, units=vtime.units, calendar=vtime.calendar)

alat = None
alon = None

if lats.ndim == 1:
    alon,alat = np.meshgrid(lons,lats)
else:
    alat = lats
    alon = lons

dd = get_distances(alat,alon,xlat,xlon)

# Return nearest point.
[jj,ii] = np.unravel_index(dd.argmin(), dd.shape)

print(jj,ii)

ncfile.close()
