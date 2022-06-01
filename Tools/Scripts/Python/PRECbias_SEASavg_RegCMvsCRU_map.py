#!/usr/bin/python25.6

"""
 Purpose: Plot seasonal bias in precipitation rate using CRU observations and RegCM output
          Period: 1980-1989
          Used variables: pre (CRU; mm/month) and pr (RegCM; kg/m2/s) 
 Author:  S. STRADA
 Date:    Sept. 27, 2018
"""

# REFERENCES:
# Read NetCDF file
#       http://www.hydro.washington.edu/~jhamman/hydro-logic/blog/2013/10/12/plot-netcdf-data/
# Python colormaps
#       https://matplotlib.org/users/colormaps.html


######################################################
# Import modules you need
#-----------------------------------------------------
from netCDF4 import Dataset as nc
import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab
from matplotlib.colors import Normalize
from cdo import *   # python version
cdo = Cdo()

# Import functions
from module_Basemap_RegCM_domain import map_RegCMdomain
# Call the map function 
# m = map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end, lon0, lat0, fontsize, dparall, dmerid)


######################################################
### Files and directories 
#-----------------------------------------------------
dirIN  = 'files'
dirOUT = 'files'
infOBS = os.path.join(dirIN,'cru_ts4_1980_1989_pre_SEASavg.nc')
infMOD = os.path.join(dirIN,'SACOR_SRF_1980_1989_SEASavg.nc')

# Needed info to draw a domain map
# Latitude and longitude of projection origin (from RegCM output file)
lat0      = -22.
lon0      = -59.
# Boundaries of the domain map
lat_start = -60. 
lat_end   = 20. 
lon_start = -112. 
lon_end   = -7.0
# Step to draw parallels and meridians
dparall   = 20
dmerid    = 20


######################################################
### Open & read NCDF files & vars 
#-----------------------------------------------------
### OBS 
OBSf     = nc(infOBS, mode='r') 
lat      = OBSf.variables['lat'][:]
lon      = OBSf.variables['lon'][:]
PR_OBS   = OBSf.variables['pre'][:,:,:]
OBSf.close()

# Print MIN-MAX discarding all NaN values, if any
print("Prec MIN", np.nanmin(PR_OBS), "Prec MAX", np.nanmax(PR_OBS))
print("Prec JJA MIN", np.nanmin(PR_OBS[2,:,:]), "Prec MAX", np.nanmax(PR_OBS[2,:,:]))
print("Prec DJF MIN", np.nanmin(PR_OBS[0,:,:]), "Prec MAX", np.nanmax(PR_OBS[0,:,:]))


### RegCM
# Use CDO in your Python script to regrid RegCM output on the CRU grid: from 25km to 50 km 
outMOD = os.path.join(dirOUT,'SACOR_SRF_1980_1989_SEASavg_regridCRU.nc')
outMOD = cdo.remapdis(infOBS, input = " " + infMOD, output = outMOD)

MODf   = nc(outMOD, mode='r')
PR_MOD = MODf.variables['pr'][:,:,:]
MODf.close()

# Convert precipitation from kg/m2/s to mm/month 
# We assume 30 days per month (approximation)
PR_MOD      = PR_MOD*30.*86400. 
print("Prec MIN", np.nanmin(PR_MOD), "Prec MAX", np.nanmax(PR_MOD))

# Print MIN-MAX discarding all NaN values, if any
print("Prec JJA MIN", np.nanmin(PR_MOD[2,:,:]), "Prec MAX", np.nanmax(PR_MOD[2,:,:]))
print("Prec DJF MIN", np.nanmin(PR_MOD[0,:,:]), "Prec MAX", np.nanmax(PR_MOD[0,:,:]))


######################################################
### Compute the seasonal bias at each grid cell
#-----------------------------------------------------
Prec_DIFF = PR_MOD[:,:,:] - PR_OBS[:,:,:]
print("Prec Bias: MIN = ", np.nanmin(Prec_DIFF), ", MAX = ", np.nanmax(Prec_DIFF))


######################################################
### Create figure environment
#-----------------------------------------------------
figMAP = plt.figure(1, figsize=(12,12))


######################################################
### Make plots
#-----------------------------------------------------
clevs    = [-50., -30., -10., -5., 5., 10., 30., 50.]
norm     = colors.BoundaryNorm(boundaries=clevs, ncolors=256)
cmap     = plt.cm.BrBG                
cmap.set_bad(color='gray')

seas={1:"DJF",2:"MAM",3:"JJA",4:"SON"} # set dictionary

# Set up figure title
plt.suptitle('Seasonal Bias: RegCM - CRU4.0 [1980-1989]', fontsize=15)

ii=1
for jj in range(0, 4): #loop over the 4 seasons

     ax = plt.subplot(2,2,ii) # set up panel plot
     m            = map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end, lon0, lat0, 10, dparall, dmerid)
     lon_S, lat_S = np.meshgrid(lon, lat)
     x, y         = m(lon_S, lat_S) # compute map proj coordinates.
     plt.ylabel('Latitude (degrees)', fontsize=10, labelpad=29)
     plt.xlabel('Longitude (degrees)', fontsize=10, labelpad=18)
     plot         = m.pcolormesh(x, y, Prec_DIFF[jj,:,:], norm=norm, vmin=-50., vmax=50., cmap=cmap)

     ax.text(0.020,0.9, " %s" % (seas[jj+1]), transform=ax.transAxes, fontsize=12, bbox=dict(boxstyle='square',
            fc='w', alpha=1.0), zorder=100) # draw label

     ii=ii+1

# Set up label
cax = figMAP.add_axes([0.2, 0.035, 0.6, 0.015])
art = plt.colorbar(plot, cax, ticks=clevs, spacing='uniform', orientation='horizontal')
art.set_label('Precipitation (mm month$^{-1}$)', fontsize=12)
art.ax.tick_params(labelsize=10)

# Adjust layout
plt.tight_layout(pad=3.5, w_pad=1.5, h_pad=1.5)

# Show
plt.show()
