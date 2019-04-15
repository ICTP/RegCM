#!/usr/bin/python2.6

""" Here a comment starts, with 3 quotation marks. In the same way, the comment ends ... 

 Purpose: Draw a base map of the CORDEX domain 
          Selected projection: Lambert Conformal Projection
 Date:    Sept. 26, 2018
 Author:  S. STRADA

 REFERENCES:
 Basemap Tool
	http://basemaptutorial.readthedocs.org/en/latest/index.html
	https://matplotlib.org/basemap/

""" 


######################################################
# Import modules you need
#-----------------------------------------------------
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
import numpy as np


######################################################
### Python fuction to build a map using a specific projection 
#-----------------------------------------------------
def map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end, lon0, lat0, fontsize, dparall, dmerid):
    """
    How to call the function in a script to create a basemap object : 
  	1. Import function to create the domain
 	from module_RegCM_domain import basemap_RegCMdomain

	2. Call the function and pass to it all needed variables
	map = map_RegCMdomain(ax, lat_start, lat_end, lon_start, lon_end, lon_end, lon_0, lat0, fontsize)) 

	Setup Miller Cyclindrical Projection
	 --> llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon are the lat/lon values of the lower left and upper right corners of the map
	 --> resolution = 'i' means intermediate coastline resolution
	 --> area_thresh=1000 means don't plot coastline features less than 1000 km^2 in area (pay attention to this if you need to plot small islands!)

    """
    m = Basemap(ax=ax, llcrnrlon=lon_start, llcrnrlat=lat_start, urcrnrlon=lon_end, urcrnrlat=lat_end,
                resolution='i', area_thresh=1000., projection='mill', lon_0=lon0, lat_0=lat0, lat_ts=0)

    m.drawcoastlines(color='k',linewidth=1, zorder=10)
    m.drawcountries(color='k',linewidth=0.5, zorder=11)
    m.drawparallels(range(-90, 90, dparall), labels=[1,0,0,0], fontsize=fontsize, dashes=[1, 2],linewidth=1, color='k', zorder=12)
    m.drawmeridians(range(-180, 180, dmerid),labels=[0,0,0,1], fontsize=fontsize, dashes=[1, 2],linewidth=1, color='k', zorder=12)

    return m
