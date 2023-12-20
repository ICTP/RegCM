#!/usr/bin/env python3

import sys

if len(sys.argv) < 2:
    print('Need netcdf filename.')
    sys.exit(-1)

from netCDF4 import Dataset
import numpy as np

ds = Dataset(sys.argv[1],'r')

try:
    res = float(sys.argv[2])
except:
    res = 0.5

lat = ds.variables['xlat'][:]
lon = ds.variables['xlon'][:]

nx,ny = np.shape(lat)

print('Nx             = ', nx)
print('Ny             = ', ny)
print('Lon 1,1        = ',lon[0,0])
print('Lon 1,ny       = ',lon[-1,0])
print('Lon nx,ny      = ',lon[-1,-1])
print('Lon nx,1       = ',lon[0,-1])
print('Lat 1,1        = ',lat[0,0])
print('Lat 1,ny       = ',lat[-1,0])
print('Lat nx,ny      = ',lat[-1,-1])
print('Lat nx,1       = ',lat[0,-1])
print('Pole Lat,Lon   = (',
        ds.grid_north_pole_latitude,', ',
        ds.grid_north_pole_longitude,')')
print('Center Lat,Lon = (',
        ds.latitude_of_projection_origin,', ',
        ds.longitude_of_projection_origin,')')

sys.exit(0)
