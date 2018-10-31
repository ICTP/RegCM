#/usr/bin/python25.6

"""
 Purpose: Plot cross section of relative humidity
          Used variable: rh (RegCM) 
 Author:  Susanna STRADA
 Date:    Sept. 27, 2018
"""

# REFERENCES:
# Read NetCDF file
#       http://www.hydro.washington.edu/~jhamman/hydro-logic/blog/2013/10/12/plot-netcdf-data/
# Basemap and hatching
# 	https://stackoverflow.com/questions/41664850/hatch-area-using-pcolormesh-in-basemap
# WRF Python
# 	http://wrf-python.readthedocs.io/en/latest/plot.html#vertical-cross-section
#	http://wrf-python.readthedocs.io/en/latest/user_api/index.html#interpolation-routines
# 	https://groups.google.com/a/ucar.edu/forum/#!msg/wrfpython-talk/Vu6ydwl2U20/BV4TWf7rCQAJ


######################################################
# Import modules you need
#-----------------------------------------------------
from netCDF4 import Dataset as nc
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
import os

# Import function to create the domain
from module_Basemap_RegCM_domain import map_RegCMdomain
from wrf import to_np, CoordPair, vertcross, latlon_coords, interpline, interp2dxy, xy
from module_find_index_closest_point import find_ind_latlon2d


######################################################
# Specify directories & job's name-period
#-----------------------------------------------------
dirIN  = 'files'
infMOD = os.path.join(dirIN,'SACOR_ATM.198503_pressure.nc')

## Needed info to draw a domain map
# Latitude and longitude of projection origin (from RegCM output file)
lat0     = -22.
lon0     = -59.
# Boundaries of the domain map
lat_start = -60.
lat_end   = 20.
lon_start = -112.
lon_end   = -7.0
# Step to draw parallels and meridians
dparall  = 20
dmerid   = 20

## Start and end coordinates for the cross section
start_lat = -40.0; start_lon = -80.0
end_lat = -20.0; end_lon = -60.0   


######################################################
### Open & read NCDF files & vars 
#-----------------------------------------------------
inf_mod   = nc(infMOD, mode='r')
lat       = inf_mod.variables['xlat'][:,:]
lon       = inf_mod.variables['xlon'][:,:]
PS  	  = inf_mod.variables['ps'][0,:,:]
plev      = inf_mod.variables['plev'][:]
RH        = inf_mod.variables['rh'][0,:,:,:] # Important: you have a single element on the time dimension, but you should put that axis equal to 0
inf_mod.close()

# Print MIN-MAX discarding all NaN values, if any
print("RH MIN", np.nanmin(RH), "RH MAX", np.nanmax(RH))

# Get number of pressure level
nlev = len(plev)


######################################################
### Compute the cross section 
#-----------------------------------------------------
## 1. Find the starting and ending point of the transect
start_i, start_j = find_ind_latlon2d(lat, start_lat, lon, start_lon)
end_i, end_j     = find_ind_latlon2d(lat, end_lat, lon, end_lon)
print(start_i, start_j)
print(lat[start_i, start_j], lon[start_i, start_j])
print(end_i, end_j)
print(lat[end_i, end_j], lon[end_i, end_j])
start = (start_j, start_i)
end   = (end_j, end_i)

## 2. Compute the vertical cross-section interpolation
xy   	   = xy(PS, start_point=start, end_point=end)
rh_cross   = interp2dxy(RH, xy, meta=False)
print(np.nanmin(rh_cross), np.nanmax(rh_cross))

## 3. Interpolate LAT-LON along the cross section and build LAT-LON pairs (to put them on the X-axis in the final plot)
lat_interp   = interpline(lat, start_point=CoordPair(x=start_j, y=start_i), end_point=CoordPair(x=end_j, y=end_i), latlon=True, meta=False)
lon_interp   = interpline(lon, start_point=CoordPair(x=start_j, y=start_i), end_point=CoordPair(x=end_j, y=end_i), latlon=False, meta=False)
latlon_pairs = np.dstack((lat_interp, lon_interp))
ps_interp    = interpline(PS, start_point=CoordPair(x=start_j, y=start_i), end_point=CoordPair(x=end_j, y=end_i), latlon=False, meta=False)


######################################################
### Create figure to make a plot
#-----------------------------------------------------
fig_MAP = plt.figure(1, figsize=(15,10))


######################################################
### Make plots
#-----------------------------------------------------
clevs = [5.0, 10.0, 20., 30., 40., 50., 60., 70., 80., 90.]
norm  = colors.BoundaryNorm(boundaries=clevs, ncolors=256)
cmap  = "BrBG" 

#plt.title("Vertical Cross Section of Relative Humidity (%)", fontsize=15)

# Plot an horizontal map
ax = fig_MAP.add_subplot(1,2,1)                   # Numbers specify: Nrows, Ncols, FigNum
ax.set_title('RH at '+str(plev[2])+' (%) ', fontsize=12)
plt.ylabel('Latitude (degrees)', fontsize=10, labelpad=33)
plt.xlabel('Longitude (degrees)', fontsize=10, labelpad=24)
m            = map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end, lon0, lat0, 10, dparall, dmerid)
lon_S, lat_S = m(lon, lat) 
plot_RH      = m.contourf(lon_S, lat_S, RH[2,:,:], norm=norm, levels=clevs, ax=ax, zorder=1, cmap=cmap, extend='both')
# Draw the transect line
point_x, point_y = m([start_lon, end_lon], [start_lat, end_lat])
m.plot([point_x[0], point_x[1]], [point_y[0], point_y[1]], color="blue",
        marker="o", zorder=3, ax=ax)
# Add colorbar
cbarD = m.colorbar(plot_RH, pad="2%")
cbarD.set_label('%')

#-----------------------------------------------------
# Plot the vertical cross section

nxy    = np.shape(rh_cross)[1]
print(nxy)
x      = np.arange(0,nxy,1)
X,Z    = np.meshgrid(x, plev)
x_axis = np.arange(0,nxy,1)

y1     = np.empty(nxy,'float')
y1.fill(1000.)
y2     = np.minimum(ps_interp[:]/100.,y1)

ax = fig_MAP.add_subplot(1,2,2)                   # Numbers specify: Nrows, Ncols, FigNum
ax.set_title('RH (%)', fontsize=12)
plt.ylabel('Vertical levels (hPa)', fontsize=12, labelpad=33)
plt.xlabel('Latitude, Longitude (degrees)', fontsize=10, labelpad=24)
cross_RH = ax.contourf(x_axis, plev, to_np(rh_cross[:,:]), levels=clevs, cmap=cmap, extend='both')

# Add the orography
YY = np.arange(0, plev.shape[0], 2)
XX = np.arange(0, x_axis.shape[0], 2)
points = np.meshgrid(YY, XX)
xx, yy = np.meshgrid(YY, XX)
x, y   = m(xx, yy) 
ax.fill_between(x_axis, y2, 1000.0, facecolor='black')
plt.gca().invert_yaxis()

# Add the color bar
plt.colorbar(cross_RH, ax=ax)

# Set the x-ticks to use latitude and longitude labels.
coord_pairs = latlon_pairs[0,:,:]
x_ticks     = np.arange(coord_pairs.shape[0])
x_labels    = coord_pairs
ax.set_xticks(x_ticks[::5])
ax.set_xticklabels(x_labels[::5], rotation=70, fontsize=8)

# Set the y-ticks 
vert_vals = plev[::-1]
ax.set_yticks(vert_vals)
ax.set_yticklabels(vert_vals[:], fontsize=8)


#-----------------------------------------------------
plt.tight_layout(pad=3.5, h_pad=1.5, w_pad=1.2)        # To adjust spacing between subplots and inside the figure
plt.show()
