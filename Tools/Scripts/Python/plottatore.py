#!/usr/bin/env python3

import sys
import os
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('cairo')

def usage(msg):
    print(msg)
    print(sys.argv[0]+': plot a variable from NetCDF 2D file')
    print('Usage: '+sys.argv[0]+' filename.nc variable')
    sys.exit(1)

def main(filename,variable):
    lat = None
    lon = None
    var = None
    try:
        ds = xr.open_dataset(filename)
        lon = ds.xlon[:]
        lat = ds.xlat[:]
        var = -ds[variable]
    except Exception as e:
        usage('File or variable error : '+str(e))
    if len(var.shape) == 4:
        pvar = var[0,0,Ellipsis]
    elif len(var.shape) == 3:
        pvar = var[0,Ellipsis]
    else:
        pvar = var
    ax = plt.axes(projection=ccrs.PlateCarree())
    if variable == 'pr':
        pvar.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(),
                             x='xlon', y='xlat', vmin=-5e-5, vmax=5e-5)
    else:
        pvar.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(),
                             x='xlon', y='xlat')
    ax.coastlines( )
    ax.gridlines(draw_labels=True)
    plt.savefig(os.path.join('plots',variable+'.pdf'), format='pdf',
		papertype='a4',orientation='landscape',transparent=True)
    #plt.show()

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2])
