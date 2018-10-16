#!/usr/bin/python25.6

"""
 Purpose: Define a way to explicitly set the mid point of a colorbar
          Use in a plot, e.g.: 
		norm=MidpointNormalize(midpoint=0.)
 Date:    Sept 26, 2018
 Author:  S. STRADA
"""


######################################################
# Import modules you need
#-----------------------------------------------------
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.colors import Normalize

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

