#!/usr/bin/python

"""
 Purpose: Draw a base map depicting the domain
          Selected projection: Miller
 Date:    April 12, 2017
 Author:  Susanna Strada

 REFERENCES:
 1. Basemap Tool
       http://basemaptutorial.readthedocs.org/en/latest/index.html
       http://matplotlixb.org/basemap/api/basemap_api.html#module-mpl_toolkits.basemap
 2. For Lambert Proj
       Conformal:  http://matplotlib.org/basemap/users/lcc.html
       Equal-Area: http://matplotlib.org/basemap/users/laea.html
"""

######################################################
# Import modules you need
#-----------------------------------------------------
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
from netCDF4 import Dataset as nc
import numpy as np
import os


######################################################
# Specify directories & job's name-period
#-----------------------------------------------------
dirnc = 'input'
domname = 'mm5'

######################################################
### Open & read NCDF files & vars 
#-----------------------------------------------------
# RegCM file
RCMf    = nc(os.path.join(dirnc,domname+'_DOMAIN000.nc'), mode='r')
lat     = RCMf.variables['xlat'][:,:]
lon     = RCMf.variables['xlon'][:,:]
topo    = RCMf.variables['topo'][:,:]
mask    = RCMf.variables['mask'][:,:]
RCMf.close()

lat_start  = np.min(lat[:,:])
lat_end    = np.max(lat[:,:])
lon_start  = np.min(lon[:,:])
lon_end    = np.max(lon[:,:])

######################################################
### Build a map using a specific projection 
### e.g. mill = Miller --> Ref.: http://matplotlib.org/basemap/users/mill.html
#-----------------------------------------------------
### MAp over the MedCORDEX domain  
def map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end, fontsize=10):
    """
    use: map = map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end) # to create a basemap object
    """
    m = Basemap(ax=ax, llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end,urcrnrlat=lat_end,
                resolution='i', area_thresh=10000., projection='mill', lon_0=46.47, lat_0=11.39, lat_ts=0)

    m.drawcoastlines(color='k',linewidth=1, zorder=10)
    m.drawcountries(color='k',linewidth=0.5, zorder=11)
    m.drawparallels(range(10, 80, 5), labels=[1,0,0,0], fontsize=fontsize, dashes=[1, 2],linewidth=1, color='k', zorder=12)
    m.drawmeridians(range(-30, 80, 5),labels=[0,0,0,1], fontsize=fontsize, dashes=[1, 2],linewidth=1, color='k', zorder=12)
    m.fillcontinents(color='coral',lake_color='aqua', zorder=1)

    return m

### Test if everything works fine
fig_MED = plt.figure(1, figsize=(12,8))
ax      = fig_MED.add_subplot(1,1,1)
m       = map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end)
## Print out 
plt.show()

# End
