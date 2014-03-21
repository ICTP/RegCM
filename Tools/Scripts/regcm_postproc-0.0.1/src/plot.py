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

import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from value import Value

class Plotter(object):

    def __init__(self, value):
        self._value = value
        self.lat, self.lon = value.latlon

    def plot(self, coastlines=True,
                   countries=True,
                   places=True,
                   title=None,
                   levels = None):
        if levels is not None:
            l_min, l_max = levels
            l = (l_max - l_min) / 10
            levels = range(l_min, l_max + l, l)
        projection = ccrs.PlateCarree()
        self.fig, self.ax = plt.subplots(subplot_kw={'projection': projection})

        if coastlines:
            self.ax.coastlines('10m')
        if countries:
            countries = cfeature.NaturalEarthFeature(
                        scale='110m', category='cultural', name='admin_0_countries')
            self.ax.add_feature(countries, color='r', alpha=0.1)
        if places:
            places = cfeature.NaturalEarthFeature(
                        scale='110m', category='cultural', name='populated_places')
            self.ax.add_feature(places, color='b', hatch='o')

        cx = self.ax.contourf(self.lon, self.lat, self._value.data, transform=ccrs.PlateCarree(),cmap='bwr', levels=levels)

        # To mask out OCEAN or LAND
        #ax.add_feature(cfeature.OCEAN)
        #ax.add_feature(cfeature.LAND)

        self.ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                        linewidth=1, color='blue', alpha=0.5, linestyle='-')

        self.fig.colorbar(cx)

        times = self._value.limits['time']
        plt.title(self._value.title + ' [' + self._value.units + ']\n' +
            'mean between ' + str(times[0]) + ' and ' + str(times[1]) + '\n')

    def show(self):
        plt.show()

    def save(self, filename, format):
        plt.savefig(filename + '.' + format)

    def close(self):
        plt.close(self.fig)


if __name__ == "__main__":
    pass