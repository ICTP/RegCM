#!/usr/bin/python25.6

"""
 Purpose: Plot Monthly AVG surface air temperature (RegCM var: tas) and precipitation (RegCM var: pr) 
          To run the script: 
		python prec_tas_RegCM_SouthAmerica.py 
 Date:    Sept 26, 2018
 Author:  S. STRADA
"""

# REFERENCES:
# Read NetCDF file
#       http://www.hydro.washington.edu/~jhamman/hydro-logic/blog/2013/10/12/plot-netcdf-data/
# Python colormaps
#	https://matplotlib.org/users/colormaps.html


######################################################
# Import modules you need
#-----------------------------------------------------
from netCDF4 import Dataset as nc
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab
import os
from matplotlib.colors import Normalize

# Import functions
from module_midpoint_normalize import MidpointNormalize
from module_Basemap_RegCM_domain import map_RegCMdomain
# Call the map function 
# m = map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end, lon0, lat0, fontsize, dparall, dmerid)


######################################################
# Specify directories, file or other info
#-----------------------------------------------------
dirRegCM = '/home/clima-archive4/CORDEX2/monthly_original/SOUTH_AMERICA/South_America_ERAINT_evaluation/'
year     = '1985'
month    = '03'

# Needed info to draw a domain map
# Latitude and longitude of projection origin (from RegCM output file)
lat0     = -22.
lon0     = -59.
# Step to draw parallels and meridians
dparall  = 20
dmerid   = 20


######################################################
### Open & read NCDF files & vars 
#-----------------------------------------------------
inf_mod   = nc(os.path.join(dirRegCM,'SACOR_SRF.'+year+month+'.nc'), mode='r') 
lat       = inf_mod.variables['xlat'][:,:]
lon       = inf_mod.variables['xlon'][:,:]
mask      = inf_mod.variables['mask'][:,:]
PREC      = inf_mod.variables['pr'][0,:,:] # Important: you have a single element on the time dimension, but you should put that axis equal to 0
TAS       = inf_mod.variables['tas'][0,0,:,:] 
inf_mod.close()

# Print MIN-MAX discarding all NaN values, if any
print("Temp MIN", np.nanmin(TAS), "Temp MAX", np.nanmax(TAS))
print("Prec MIN", np.nanmin(PREC), "Prec MAX", np.nanmax(PREC))

# Convert precipitation from kg/m2/s to mm/day
PREC      = PREC * 86400.
print("Prec MIN", np.nanmin(PREC), "Prec MAX", np.nanmax(PREC))
PREC_ma   = ma.masked_where(mask == 0, PREC)

# Convert temperature from Kelvin to Celsius
TAS       = TAS - 273.15
print("Temp MIN", np.nanmin(TAS), "Temp MAX", np.nanmax(TAS))

# Compute a domain average discarding NaN values, if any, and print the result 
PREC_avg  = np.nanmean(PREC[:,:])
print("Prec. Domain-AVG", PREC_avg)
TAS_avg   = np.nanmean(TAS[:,:])
print("Temp. Domain-AVG", TAS_avg)

# Identify MIN-MAX lat-lon (to be used to build the domain map)
lat_start = np.min(lat[:])
lat_end   = np.max(lat[:])
lon_start = np.min(lon[:])
lon_end   = np.max(lon[:])
print("Lat start and end; Lon start and end ", lat_start, lat_end, lon_start, lon_end)


######################################################
### Create a figure environment (you can adjust the figsize)
#-----------------------------------------------------
fig_MAP = plt.figure(1, figsize=(12,8))


######################################################
### Make plots
#-----------------------------------------------------
# TEMPERATURE
# Define contour levels
clev_TAS     = [-16., -12., -8., -4., -2., 2., 4., 8., 12., 16.]
norm_TAS     = colors.BoundaryNorm(boundaries=clev_TAS, ncolors=256)
cmap_TAS     = plt.cm.bwr

# Figure title
plt.suptitle('Monthly AVG ('+year+'/'+month+'), RegCM (25 km)', fontsize=15)

# TAS: Raster plot 
ax = fig_MAP.add_subplot(2,2,1)                   # Numbers specify: Nrows, Ncols, FigNum
ax.set_title('Surf. Air Temp. ($^{\circ}$C), Raster plot', fontsize=14)
plt.ylabel('Latitude (degrees)', fontsize=12, labelpad=29)
plt.xlabel('Longitude (degrees)', fontsize=12, labelpad=18)
m            = map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end, lon0, lat0, 10, dparall, dmerid)
lon_S, lat_S = m(lon, lat) # Compute map projection coordinates
#plotTAS1     = m.pcolormesh(lon_S, lat_S, TAS, norm=MidpointNormalize(midpoint=0.), ax=ax, zorder=2, cmap=cmap_TAS)
plotTAS1     = m.pcolormesh(lon_S, lat_S, TAS, norm=norm_TAS, ax=ax, zorder=2, cmap=cmap_TAS)
# Add colorbar
cbarD = m.colorbar(plotTAS1, spacing='proportional', ticks=clev_TAS, boundaries=clev_TAS, pad="2%")
cbarD.set_label('$^{\circ}$C')


# TAS: Filled contour plot 
ax = fig_MAP.add_subplot(2,2,2)                   # Numbers specify: Nrows, Ncols, FigNum
ax.set_title('Surf. Air Temp. ($^{\circ}$C), Filled contour plot', fontsize=14)
plt.ylabel('Latitude (degrees)', fontsize=12, labelpad=29)
plt.xlabel('Longitude (degrees)', fontsize=12, labelpad=18)
m            = map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end, lon0, lat0, 10, dparall, dmerid)
lon_S, lat_S = m(lon, lat) # compute map proj coordinates
plotTAS2     = m.contourf(lon_S, lat_S, TAS, norm=norm_TAS, levels=clev_TAS, ax=ax, zorder=2, cmap=cmap_TAS, extend='both')
# Add colorbar
cbarD = m.colorbar(plotTAS2, spacing='proportional', ticks=clev_TAS, boundaries=clev_TAS, pad="2%")
cbarD.set_label('$^{\circ}$C')


#-----------------------------------------------------
# PRECIPITATION
# Define contour levels
clev_PR      = [0.0, 2.0, 4.0, 6., 8.0, 10.0, 15., 20., 30.0]
norm_PR      = colors.BoundaryNorm(boundaries=clev_PR, ncolors=256)
cmap_PR      = plt.cm.Blues

# PR: Raster plot 
ax = fig_MAP.add_subplot(2,2,3)                   # Numbers specify: Nrows, Ncols, FigNum
ax.set_title('Prec. (mm day$^{-1}$), Raster plot', fontsize=14)
plt.ylabel('Latitude (degrees)', fontsize=12, labelpad=29)
plt.xlabel('Longitude (degrees)', fontsize=12, labelpad=18)
m            = map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end, lon0, lat0, 10, dparall, dmerid)
lon_S, lat_S = m(lon, lat) # Compute map projection coordinates
plotPR1      = m.pcolormesh(lon_S, lat_S, PREC_ma, norm=norm_PR, vmin=0.0, vmax=30.0, ax=ax, zorder=2, cmap=cmap_PR)
# Add colorbar
cbarD = m.colorbar(plotPR1, spacing='proportional', ticks=clev_PR, boundaries=clev_PR, pad="2%")
cbarD.set_label('mm day$^{-1}$')


# PR: Filled contour plot 
ax = fig_MAP.add_subplot(2,2,4)                   # Numbers specify: Nrows, Ncols, FigNum
ax.set_title('Prec. (mm day$^{-1}$), Filled contour plot', fontsize=14)
plt.ylabel('Latitude (degrees)', fontsize=12, labelpad=29)
plt.xlabel('Longitude (degrees)', fontsize=12, labelpad=18)
m            = map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end, lon0, lat0, 10, dparall, dmerid)
lon_S, lat_S = m(lon, lat) # compute map proj coordinates
plotPR2      = m.contourf(lon_S, lat_S, PREC, norm=norm_PR, levels=clev_PR, ax=ax, zorder=2, cmap=cmap_PR, extend='both')
# Add colorbar
cbarD = m.colorbar(plotPR2, spacing='proportional', ticks=clev_PR, boundaries=clev_PR, pad="2%")
cbarD.set_label('mm day$^{-1}$')

#-----------------------------------------------------
plt.tight_layout(pad=3.5, h_pad=2.0, w_pad=2.0)        # To adjust spacing between subplots and inside the figure
plt.show() # To pop up the plot

