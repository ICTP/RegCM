#!/usr/bin/python25.6

"""
 Purpose: Plot the seasonal cycle of the precipitation in RegCM output and CRU observations
 Author:  S. STRADA
 Date:    Sept. 27, 2018
"""

# REFERENCES:
# Read NetCDF file
#       http://www.hydro.washington.edu/~jhamman/hydro-logic/blog/2013/10/12/plot-netcdf-data/

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
import os
import pylab
from matplotlib.colors import Normalize
from cdo import *   # python version
cdo = Cdo()


######################################################
### Files and directories 
#-----------------------------------------------------
dirIN  = 'files'
infOBS = os.path.join(dirIN,'cru_ts4_1980_1989_pre_MONavg.nc')
infOBS_SA = os.path.join(dirIN,'cru_ts4_1980_1989_pre_MONavg_SA.nc')
infMOD = os.path.join(dirIN,'SACOR_SRF_1980_1989_MONavg.nc')
infMOD_SA = os.path.join(dirIN,'SACOR_SRF_1980_1989_MONavg_SA.nc')


######################################################
### Open & read NCDF files & vars 
#-----------------------------------------------------
# Use CDO in Python to select a lat-lon box
infOBS_SA = cdo.sellonlatbox("-112,-7,-60,20 ", input = infOBS, output = infOBS_SA)
### OBS 
OBSf      = nc(infOBS_SA, mode='r')
PR_OBS    = OBSf.variables['pre'][:,:,:]
OBSf.close()
# Print MIN-MAX discarding all NaN values, if any
print("Prec MIN", np.nanmin(PR_OBS), "Prec MAX", np.nanmax(PR_OBS))


### RegCM
infMOD_SA = cdo.sellonlatbox("-112,-7,-60,20 ", input = infMOD, output = infMOD_SA)
MODf      = nc(infMOD_SA, mode='r')
PR_MOD    = MODf.variables['pr'][:,:,:]
mask      = MODf.variables['mask'][:,:]
MODf.close()
print("Prec South-America MIN", np.nanmin(PR_MOD), "Prec South-America MAX", np.nanmax(PR_MOD))

# Convert precipitation from kg/m2/s to mm/month 
# We assume 30 days per month (approximation)
PR_MOD = PR_MOD*30.*86400.
print("Prec South-America MIN", np.nanmin(PR_MOD), "Prec South-America MAX", np.nanmax(PR_MOD))
PR_MOD_ma = np.empty_like(PR_MOD)
for t in range(0, 12):
    PR_MOD_ma[t,:,:] = ma.masked_where(mask == 0, PR_MOD[t,:,:])


######################################################
### Compute GLOBAL ANNUAL mean at each grid cell
#-----------------------------------------------------
PRobs_ANNcycle = np.ma.empty((12), dtype='float')
PRobs_ANNcycle.fill(np.nan)
PRmod_ANNcycle = np.ma.empty((12), dtype='float')
PRmod_ANNcycle.fill(np.nan)
PR_DIFF        = np.ma.empty((12), dtype='float')
PR_DIFF.fill(np.nan)

for t in range(0, 12):
   PRobs_ANNcycle[t]  = np.nanmean(PR_OBS[t,:,:])
   PRmod_ANNcycle[t]  = np.nanmean(PR_MOD_ma[t,:,:])
   PR_DIFF[t]         = PRmod_ANNcycle[t] - PRobs_ANNcycle[t]


######################################################
### Create figure to make a plot
#-----------------------------------------------------
fig_XY = plt.figure(1, figsize=(10,8))


######################################################
### Make plots
#-----------------------------------------------------
## Precipitation
ax = fig_XY.add_subplot(1,1,1)                   # Numbers specify: Nrows, Ncols, FigNum
x_axis = np.arange(1,13,1)
#for t in range(0, 12):
#   plt.scatter([t+1], PRobs[t], color="green")
#   plt.scatter([t+1], PRrcm[t], color="blue")
plt.plot(x_axis, PRmod_ANNcycle, color="blue", ls=':', lw=1.5, label="RegCM")
plt.plot(x_axis, PRobs_ANNcycle, color="green", ls=':', lw=1.8, label="CRU")
plt.plot(x_axis, PR_DIFF, color="red", ls='-', lw=2.0, label="RegCM - CRU")
plt.title('Seasonal cycle of Precipitation', fontsize=15)
plt.xlabel('Month', fontsize=12)
#plt.ylabel('$\lambda$E (MJoule m$^{-2}$ day$^{-1}$)', fontsize=20)
plt.ylabel('Prec. (mm month$^{-1}$)', fontsize=12)
plt.axis([0., 13., -100.0, 250.])
plt.xticks(np.arange(0., 13., 1.0), (['','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','']), fontsize=12)
plt.yticks(np.arange(-100.0, 250., 20.0), fontsize=12)
plt.axhline(y=0.0, color='k', linestyle='-')
plt.grid(True)
plt.legend(loc='best', fontsize=12)


#-----------------------------------------------------
plt.show()

# End
