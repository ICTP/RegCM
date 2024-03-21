#!/usr/bin/env python3

import sys
import numpy as np
from netCDF4 import Dataset

try:
    ds = Dataset(sys.argv[1],"r")
except:
    print('I do need the input file path to work on.')
    print(sys.argv[0]+' regcm_[ATM or SRF].nc file(s)')
    sys.exit(1)

if ds.projection != 'ROTLLR':
    print('Written for Rotated Lat-Lon projection (ROTLLR)')
    sys.exit(1)

plat = ds.grid_north_pole_latitude
plon = ds.grid_north_pole_longitude
xlat = ds.variables['xlat'][:]
xlon = ds.variables['xlon'][:]
ds.close( )

def rotate(u,v):
    if np.abs(plat - 90.0) < 0.001:
        return u, v
    plam = np.radians(plon)
    pphi = np.radians(plat)
    phi = np.radians(xlat)
    lam = np.where(abs(xlat)>89.99999, 0.0, np.radians(xlon))
    dlam = plam - lam
    f1 = np.cos(pphi)*np.sin(dlam)
    f2 = np.cos(phi)*np.sin(pphi) - np.sin(phi)*np.cos(pphi)*np.cos(dlam)
    delta = np.arctan(f1/f2)
    return np.cos(delta)*u+v*np.sin(delta), -u*np.sin(delta)+v*np.cos(delta)

for fname in sys.argv[1:]:
    print('Rotating vectors in '+fname+'.....')
    ds = Dataset(fname,"r+")
    if "uas" in ds.variables: # SRF file
        u = ds.variables['uas']
        v = ds.variables['vas']
        tu = ds.variables['tauu']
        tv = ds.variables['tauv']
        u100 = ds.variables['ua100m']
        v100 = ds.variables['va100m']

        u.standard_name = 'eastward_wind'
        u.long_name = 'Eastward Near-Surface Wind'
        v.standard_name = 'northward_wind'
        v.long_name = 'Northward Near-Surface Wind'
        tu.long_name = "Surface Downward Eastward Wind Stress"
        tu.standard_name = "surface_downward_eastward_stress"
        tv.long_name = "Surface Downward Northward Wind Stress"
        tv.standard_name = "surface_downward_northward_stress"
        u100.long_name = "Eastward Wind at 100m" ;
        u100.standard_name = "eastward_wind" ;
        v100.long_name = "Northward Wind at 100m" ;
        v100.standard_name = "northward_wind" ;

        dims = np.shape(u)
        for t in range(dims[0]):
            ux = u[t,Ellipsis]
            vx = v[t,Ellipsis]
            tux = tu[t,Ellipsis]
            tvx = tv[t,Ellipsis]
            u0x = u100[t,Ellipsis]
            v0x = v100[t,Ellipsis]
            u[t,Ellipsis], v[t,Ellipsis] = rotate(ux,vx)
            tu[t,Ellipsis], tv[t,Ellipsis] = rotate(tux,tvx)
            u100[t,Ellipsis], v100[t,Ellipsis] = rotate(u0x,v0x)
    elif "ua" in ds.variables: # ATM file
        u = ds.variables['ua']
        v = ds.variables['va']
        u.standard_name = 'eastward_wind'
        u.long_name = 'Eastward Wind'
        v.standard_name = 'northward_wind'
        v.long_name = 'Northward Wind'
        dims = np.shape(u)
        for t in range(dims[0]):
            for z in range(dims[1]):
                ux = u[t,z,Ellipsis]
                vx = v[t,z,Ellipsis]
                u[t,z,Ellipsis], v[t,z,Ellipsis] = rotate(ux,vx)
    else:
        print('Skip: Nothing to rotate here.')
    ds.close( )
