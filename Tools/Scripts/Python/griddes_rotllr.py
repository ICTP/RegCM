#!/usr/bin/env python3

import sys

def scoords(lon,lat):
    x = lon
    y = lat
    if ( lon < 0.0 ):
        x = 360+x
    return('('+"{0:0.2f}".format(x)+', '+"{0:0.2f}".format(y)+')')

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

try:
    lat = ds.variables['xlat'][:]
    lon = ds.variables['xlon'][:]
except:
    try:
        lat = ds.variables['lat'][:]
        lon = ds.variables['lon'][:]
    except:
        print('No coordinates...')
        sys.exit(-1)

ny, nx = np.shape(lat)
hny, hnx = ny//2 , nx//2

print('Nx             = ', nx)
print('Ny             = ', ny)
print('------------------------------------------------------------')
print('Coordinates in actual coordinates of corners:')
print('------------------------------------------------------------')
print('TLC            = ',scoords(lon[-1,0],lat[-1,0]))
print('CNB            = ',scoords(lon[-1,hnx],lat[-1,hnx]))
print('TRC            = ',scoords(lon[-1,-1],lat[-1,-1]))
print('CWD            = ',scoords(lon[hny,0],lat[hny,0]))
print('CPD            = ',scoords(lon[hny,hnx],lat[hny,hnx]))
print('CED            = ',scoords(lon[hny,-1],lat[hny,-1]))
print('BLC            = ',scoords(lon[0,0],lat[0,0]))
print('CSB            = ',scoords(lon[0,hnx],lat[0,hnx]))
print('BRC            = ',scoords(lon[0,-1],lat[0,-1]))
print('------------------------------------------------------------')
print('Pole Lat,Lon   = (',
        ds.grid_north_pole_latitude,', ',
        ds.grid_north_pole_longitude,')')
print('Center Lat,Lon = (',
        ds.latitude_of_projection_origin,', ',
        ds.longitude_of_projection_origin,')')

sys.exit(0)
