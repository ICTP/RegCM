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

import glob
import numpy as np
from netCDF4 import Dataset, num2date, date2num
from value import Value
from plot import Plotter

class Reader(object):

    """The generic class that reads the netCDF file. Can be subclassed with different netCDF layouts.
    """

    def __init__(self, pattern):
        self._files = []
        self._ds = []
        for f in sorted(glob.glob(pattern)):
            self._files.append(f)
            self._ds.append(Dataset(f))

    def get_value(self, var_name, imposed_limits=None, latlon_limits=None):
        # initialize limits dictionary if its not there
        if imposed_limits is None:
            imposed_limits = {}
        # if var_name's dimensions do not include lat and lon, we need to translate latlon_limits to the limits of actual var_name's dimensions
        if latlon_limits is not None:
            imposed_limits.update(self._translate_latlon_limits(latlon_limits))
        else:
            latlon_limits = self._get_latlon_limits()
        value = Value()
        # get data from all the datasets
        for f, ds in zip(self._files, self._ds):
            # data variable (we might need to impose the limits on it later on)
            print 'Reading data from file ' + f
            data = ds.variables[var_name]
            limits = {}
            # dimension names and values
            dim_names = data.dimensions
            dims = [ds.variables[name] for name in dim_names]
            title = getattr(data, 'long_name', None)
            # units (must be according to standart)
            units = getattr(data, 'units', None)
            for idim, (name, dim) in enumerate(zip(dim_names, dims)):
                dim_data = dim[:]
                if name in imposed_limits.keys():
                    # treat time differently
                    if name == "time":
                        imposed_limits[name] = date2num(imposed_limits[name], units=dim.units, calendar=getattr(dim, 'calendar', 'gregorian'))
                    # impose a limit on the dimension
                    min_dim, max_dim = imposed_limits[name]
                    dim_idx = np.where((dim_data >= min_dim) & (dim_data <= max_dim))
                    dim_data = dim_data[dim_idx]
                    data = data[dim_idx]
                    limits[name] = imposed_limits[name]
                else:
                    data = data[:]
                    limits[name] = [np.min(dim_data), np.max(dim_data)]
                # treat time differently
                if name == 'time':
                    limits[name] = num2date(limits[name], units=dim.units, calendar=getattr(dim, 'calendar', 'gregorian'))
                    dims[idim] = num2date(dim_data, units=dim.units, calendar=getattr(dim, 'calendar', 'gregorian'))
                else:
                    dims[idim] = dim_data
                axes_tuple = tuple(range(1, len(dims)) + [0,])
                data = np.transpose(data, axes_tuple)
            latlon = self._get_latlon_within_limits(latlon_limits)
            value_i = Value(data, title, units, dims, dim_names, latlon, limits, latlon_limits)
            value.update(value_i)
        return value

    def _get_latlon_limits(self):
        latlon_limits = {}
        # in hope that all datasets have the same latitude and longitude points
        ds = self._ds[0]
        for crd, crd_name in self.crd_names.iteritems():
            crd_value = ds.variables[crd_name][:]
            latlon_limits[crd] = [np.min(crd_value), np.max(crd_value)]
        return latlon_limits


class RegCMReader(Reader):
    crd_names = {'lat': 'xlat', 'lon': 'xlon'}

    def _get_latlon_within_limits(self, latlon_limits):
        # in hope that all datasets have the same latitude and longitude points
        latlon = []
        ds = self._ds[0]
        for crd in ['lat', 'lon']:
            crd_name = self.crd_names[crd]
            crd_value = ds.variables[crd_name][:]
            crd_shape = crd_value.shape
            min_crd, max_crd = latlon_limits[crd]
            crd_idx = np.where((crd_value >= min_crd) & (crd_value <= max_crd))
            latlon.append(crd_value[crd_idx].reshape(crd_shape))
        return latlon

class CRUReader(Reader):
    crd_names = {'lat': 'lat', 'lon': 'lon'}

    def _translate_latlon_limits(self, latlon_limits):
        # in CRU file a variable depends on actual lat and lon
        for limit in latlon_limits.keys():
            assert limit in self.crd_names.keys()
        return latlon_limits

    def _get_latlon_within_limits(self, latlon_limits):
        # in hope that all datasets have the same latitude and longitude points
        latlon = []
        ds = self._ds[0]
        for crd in ['lat', 'lon']:
            crd_name = self.crd_names[crd]
            crd_value = ds.variables[crd_name][:]
            min_crd, max_crd = latlon_limits[crd]
            crd_idx = np.where((crd_value >= min_crd) & (crd_value <= max_crd))
            latlon.append(crd_value[crd_idx])
        lat, lon = np.meshgrid(latlon[0], latlon[1])
        return [lat.T, lon.T]

