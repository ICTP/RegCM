#!/usr/bin/env python3

import os
import sys
import cf_xarray as cfxr
import xarray as xr
from cartopy import crs as ccrs
from matplotlib import pyplot as plt

infile = sys.argv[1]
vname = os.path.basename(infile).split('_')[0]
ds = xr.open_dataset(infile)
grid_mapping = ds.cf["grid_mapping"]

if grid_mapping.grid_mapping_name == 'rotated_latitude_longitude':
    pole_latitude = grid_mapping.grid_north_pole_latitude
    pole_longitude = grid_mapping.grid_north_pole_longitude
    projection = ccrs.RotatedPole(pole_longitude, pole_latitude)
else:
    print('Implemented only for Rotated Latitude/Longitude')
    sys.exit(-1)

p = ds[vname][0,:,:].plot(subplot_kws=dict(projection=projection,
                      transform=ccrs.PlateCarree()))
p.axes.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
p.axes.coastlines(resolution='50m')
plt.show( )
