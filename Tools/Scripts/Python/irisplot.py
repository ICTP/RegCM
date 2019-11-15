#!/usr/bin/env python3.7

import numpy as np
import cartopy.feature as cfeature
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import iris
import iris.plot as iplt
import sys
import os

interactive = True

SMALL_SIZE = 5
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

if not interactive:
    matplotlib.use('PDF')

seasons = [ 'DJF', 'MAM', 'JJA', 'SON' ]

def plotseason(data,seas,clev):
    plt.subplot(2, 2, seas)
    cf = iplt.contourf(data[seas-1],clev,cmap='seismic',extend="both")
    plt.title(seasons[seas-1])
    plt.gca().add_feature(cfeature.BORDERS)
    plt.gca().coastlines()
    #plt.gca().gridlines()
    return cf

fname = sys.argv[1]
biases = iris.load_cube(fname)

up = np.max(biases.data)
down = abs(np.min(biases.data))

best = min(up,down)/2.0

contour_levels = np.linspace(-best,best,100)

fig = plt.figure(dpi=100)

for i in range(1,5):
    cf = plotseason(biases,i,contour_levels)

fig.subplots_adjust(bottom=0.15)

# make an axes to put the shared colorbar in
colorbar_axes = plt.gcf().add_axes([0.35, 0.05, 0.3, 0.05])
colorbar = plt.colorbar(cf, colorbar_axes, orientation='horizontal')
colorbar.set_label('%s' % biases.units)

# limit the colorbar to 7 tick marks
colorbar.locator = matplotlib.ticker.MaxNLocator(7)
colorbar.formatter = matplotlib.ticker.FuncFormatter(lambda x, pos: "%.1e" % x)
colorbar.update_ticks()

title = 'Seasonal bias for %s' % biases.long_name
plt.suptitle('Seasonal bias for %s' % biases.long_name)

if interactive:
    iplt.plt.show( )
else:
    imgspec = 'pdf'
    figname = os.path.basename(fname)
    figname = os.path.splitext(figname)[0]+"."+imgspec
    iplt.plt.savefig(figname, format=imgspec,
        papertype='a4',orientation='landscape',transparent=True)
    iplt.plt.draw( )
