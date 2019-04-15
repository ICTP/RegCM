#!/usr/bin/python2.6

""" 
 Purpose: Find indices of the closest grid-cell once given lat-lon coordinates
 Date:    Sept. 26, 2018
 Author:  S. STRADA
""" 


######################################################
# Import modules you need
#-----------------------------------------------------
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
import numpy as np


######################################################
### Python function to find indices of a specific lat-lon points based on a 2D lat-lon grid  
#-----------------------------------------------------
def find_ind_latlon2d(lats, lat_pt, lons, lon_pt):
  a = abs(lats - lat_pt) + abs(lons - lon_pt)
  i,j = np.unravel_index(a.argmin(),a.shape)
  return i, j
