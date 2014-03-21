#!/usr/bin/env python

import os, sys
from read import RegCMReader, CRUReader
from plot import Plotter

usage = """usage: ./monthly_mean.py regcm_model_dir beg_yr end_yr model_nc_variable 'observ_file' observ_nc_variable

Parameters:
regcm_model_dir: a directory where RegCM model data are located
beg_yr: beginning of a year range to calculate monthly mean
end_yr: end of the year range for monthly mean
model_nc_variable: a variable from the RegCM netCDF files 
'observ_file': a CRU netCDF file
model_nc_variable: a variable from the CRU netCDF file

Example: 
./app.py data/ 1990 1995 t2m 'obs/CRUTMP.CDF' TMP
"""

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print usage
        sys.exit(1)
    nc_dir = sys.argv[1]
    beg_yr = int(sys.argv[2])
    end_yr = int(sys.argv[3])
    nc_var = sys.argv[4]
    obs_pattern = sys.argv[5]
    obs_nc_var = sys.argv[6]

    for yr in xrange(beg_yr, end_yr):
        pattern = os.path.join(nc_dir, '*' + str(yr) + '*.nc')
#        for 
        r = RegCMReader(pattern)
    value = r.get_value(nc_var).mean()
    time_limits = value.get_limits('time')
    crd_limits = value.get_latlonlimits()

    obs_r = CRUReader(obs_pattern)
    obs_value = obs_r.get_value(obs_nc_var, imposed_limits={'time': time_limits}, latlon_limits=crd_limits).mean()
    if obs_nc_var == "TMP":
        obs_value.to_K()

    value.regrid(obs_value.latlon)
    diff = obs_value - value
    plt = Plotter(diff)
    plt.plot(levels = (-5, 5))
    plt.show()
    plt.save('image', format='png')
    plt.close()
