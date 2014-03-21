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
from datetime import datetime
import numpy as np
from netCDF4 import Dataset, num2date, date2num
from scipy.interpolate import RectBivariateSpline

class Value(object):

    def __init__(self, data=None, title=None, units=None, dims=None, 
                       dim_names=None, latlon=None, limits=None, latlon_limits=None):
        self.data = data
        self.title = title
        self.units = units
        self.dims = dims
        self.dim_names = dim_names
        self.latlon = latlon
        self.limits = limits
        self.latlon_limits = latlon_limits

    def update(self, value):
        # title update only if value is empty
        if self.title is None:
            self.title = value.title
        else:
            assert self.title == value.title
        # the same goes with units
        if self.units is None:
            self.units = value.units
        else:
            assert self.units == value.units
        # and with dim_names
        if self.dim_names is None:
            self.dim_names = value.dim_names
        else:
            assert self.dim_names == value.dim_names
        # and with latlon
        if self.latlon is None:
            self.latlon = value.latlon
        else:
            assert np.all(self.latlon[0] == value.latlon[0])
            assert np.all(self.latlon[1] == value.latlon[1])
        # and with latlon_limits
        if self.latlon_limits is None:
            self.latlon_limits = value.latlon_limits
        else:
            assert self.latlon_limits == value.latlon_limits
        # data(over time axis), dims(time) and limits(time) should be expanded
        if self.data is None:
            self.data = value.data
        else:
            # hope that time axis is 0 (can check for it in dims)
            self.data = np.vstack((self.data, value.data))
        # dims
        if self.dims is None:
            self.dims = value.dims
        else:
            # also hope that time axis is 0
            self.dims[0] = np.hstack((self.dims[0], value.dims[0]))
        # limits
        if self.limits is None:
            self.limits = value.limits
        else:
            # hope that time goes increasing
            self.limits['time'] = [self.limits['time'][0], value.limits['time'][1]]

    def regrid(self, latlon):
        # ugly hack
        data = []
        if len(self.data.shape) == 2:
            data_chunks = [self.data, ]
        else:
            data_chunks = self.data
        for data_chunk in data_chunks:
            lut = RectBivariateSpline(self.latlon[0][:,0], self.latlon[1][0], data_chunk)
            new_shape = latlon[0].shape
            data.append(lut.ev(latlon[0].ravel(),latlon[1].ravel()).reshape(new_shape))
        if len(data) == 1:
            self.data = data[0]
        else:
            self.data = np.array(data)
        self.latlon = latlon

    def get_limits(self, name):
        return self.limits[name]

    def get_latlonlimits(self):
        return self.latlon_limits

    def mean(self):
        self.data = np.mean(self.data, axis=0)
        return self

    def __sub__(self, other):
        assert self.data.shape == other.data.shape
        return Value(self.data - other.data, self.title, self.units, self.dims, self.dim_names, self.latlon, self.limits, self.latlon_limits)

    def __abs__(self):
        return Value(np.abs(self.data), self.title, self.units, self.dims, self.dim_names, self.latlon, self.limits, self.latlon_limits)

    def to_C(self):
        self.data -= 273.15

    def to_K(self):
        self.data += 273.15

    def m_average(self, names, data, times):
        fdate = datetime.date(times[0])
        edate = datetime.date(times[1])
        ac = fdate.isocalendar()
        bc = edate.isocalendar()
        alist = []
        blist = []
        for i in ac:
            alist.append(i)

        for j in bc:
            blist.append(j)

        lcount = [0]
        j = 0
        count = 0
        maverage = []
        temp = []
        while alist[0] <= blist[0]:
                if alist[0] == blist[0] and alist[1] == blist[1] and alist[2] == blist[2]:
                   break
                while  alist[2] <= 11:
                       d1 = datetime.date(alist[0], alist[2], alist[1])
                       alist[2] += 1
                       d2 = datetime.date(alist[0], alist[2], alist[1])
                       diff = (d2 -d1).days
                       for i in range(1,diff+1):
                           lcount.append(i)
                           count = count + diff
                           chunk = lcount[0:diff]
                temp = llcount.append(chunk)

                alist[2] = 1
                alist[0] += 1
                m_mean = np.mean(self.data[temp], axis=0)
                maverage.append(m_mean)
        return maverage

