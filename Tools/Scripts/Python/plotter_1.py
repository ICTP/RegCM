#!/usr/bin/env python

from netCDF4 import Dataset, num2date
import numpy as np
import sys
import math
from mpl_toolkits.basemap import Basemap
from scipy import ndimage
import matplotlib.pyplot as plt
import colormaps as cmaps

def block_mean(ar, fact):
   assert isinstance(fact, int), type(fact)
   sx, sy = ar.shape
   X, Y = np.ogrid[0:sx, 0:sy]
   regions = sy/fact * (X/fact) + Y/fact
   res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
   res.shape = (sx/fact, sy/fact)
   return res

try:
    ncfile = Dataset(sys.argv[1], 'r')
except:
    if len(sys.argv) > 1:
        print('Cannot open data file '+sys.argv[1])
    else:
        print('Need at least a netcdf file name!')
    sys.exit(-1)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

lats = ncfile.variables['xlat'][:]
lons = ncfile.variables['xlon'][:]
llcrnrlon = np.amin(lons)
urcrnrlon = np.amax(lons)
llcrnrlat = np.amin(lats)
urcrnrlat = np.amax(lats)
m =  Basemap(llcrnrlon = llcrnrlon, llcrnrlat = llcrnrlat,
             urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat,
             resolution = 'h', projection = 'cyl')

m.drawlsmask(land_color='#00441b', ocean_color='#8be5e5', lakes=True)
#m.drawlsmask(land_color='#aaaaaa', ocean_color='#ffffff', lakes=True)
m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(0., 90, 10.)
m.drawparallels(parallels, labels = [1, 0, 0, 0], fontsize = 10)
meridians = np.arange(-180., 180., 10.)
m.drawmeridians(meridians, labels = [0, 0, 0, 1], fontsize = 10)

vtime = ncfile.variables['time']
xtime = vtime[:]
times = num2date(xtime, units=vtime.units, calendar=vtime.calendar)

temp = ncfile.variables['tas'][-1,0,:,:]
minv = int(math.floor(np.amin(temp)) - 1.0)
maxv = int(math.ceil(np.amax(temp)) + 1.0)
clevs = range(minv,maxv)

#u = ncfile.variables['uas'][-1,0,:,:]
#v = ncfile.variables['vas'][-1,0,:,:]

ncfile.close()

cs = m.contourf(lons, lats, temp, clevs)
cbar = m.colorbar(cs, location = 'bottom', pad = "8%")
cbar.set_label('Temperature')
#cs = m.barbs(block_mean(lons,10), block_mean(lats,10),
#             block_mean(u,10), block_mean(v,10))
#plt.title('Wind for ' + str(times[-1]))
plt.title('Temperature for ' + str(times[-1]))
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='magma', cmap=cmaps.magma)
plt.set_cmap(cmaps.magma)

plt.show()
