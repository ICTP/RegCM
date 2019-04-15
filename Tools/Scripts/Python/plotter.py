#!/usr/bin/env python

from netCDF4 import Dataset, num2date
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

try:
    ncfile = Dataset(sys.argv[1], 'r')
except:
    if len(sys.argv) > 1:
        print('Cannot open data file '+sys.argv[1])
    else:
        print('Need at least a netcdf file name!')
    sys.exit(-1)

months = [ 6, 7, 8 ]
cmonth = 'JJA'

fig = plt.figure(figsize=(8, 8))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]
llcrnrlon = np.amin(lons)
urcrnrlon = np.amax(lons)
llcrnrlat = np.amin(lats)
urcrnrlat = np.amax(lats)
m =  Basemap(llcrnrlon = llcrnrlon, llcrnrlat = llcrnrlat,
             urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat,
             resolution = 'l', projection = 'cyl')

m.drawlsmask(land_color='#00441b', ocean_color='#8be5e5', lakes=True)
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
sm = np.array([x.month for x in times])
mm = np.array(months)
ii = np.array([i for i, x in enumerate(sm) if x in mm])
y1 = times[0].year
y2 = times[-1].year

temp = np.mean(ncfile.variables['tas'][ii,:,:], axis=0)
ncfile.close()

clevs = range(np.amin(temp), np.amax(temp))
cs = m.contourf(lons, lats, temp, clevs, cmap = plt.cm.YlOrRd)
cbar = m.colorbar(cs, location = 'bottom', pad = "8%")
cbar.set_label('Temperature')
plt.title('Average '+cmonth+' Temperature from ' + repr(y1) + ' to ' + repr(y2))
plt.show()
