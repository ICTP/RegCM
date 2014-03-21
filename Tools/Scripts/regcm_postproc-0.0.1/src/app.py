#!/usr/bin/env python
#    RegCM postprocessing tool
#    Copyright (C) 2014 Aliou, Addisu, Kanhu, Andrey

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
from read import RegCMReader, CRUReader
from plot import Plotter

usage = """usage: ./app.py 'model_file_pattern' model_nc_variable 'observ_glob_pattern' observ_nc_variable

Parameters:
'model file pattern': a glob pattern for one or more netCDF files made from RegCM program (in apostrophes)
model_nc_variable: a variable from the RegCM netCDF files 
'observ_glob_pattern': a glob pattern for one or more CRU netCDF files (also in apostrophes)
model_nc_variable: a variable from the CRU netCDF files 

Example: 
./app.py 'data/Africa_SRF.1970*.nc' t2m 'obs/CRUTMP.CDF' TMP
"""

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print usage
        sys.exit(1)
    pattern = sys.argv[1]
    nc_var = sys.argv[2]
    obs_pattern = sys.argv[3]
    obs_nc_var = sys.argv[4]

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
