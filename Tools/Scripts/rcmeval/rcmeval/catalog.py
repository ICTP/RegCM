#
#    This file is part of ICTP RegCM.
#
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
"""Reference datasets (observations and reanalyses) described in a catalog
YAML file, one per site:

  datasets:
    NAME:
      root: /path                    # prepended to the file templates
      frequency: mon | day           # of the files
      surface: all | land | ocean    # where the dataset has data
      variables:
        VAR:                         # CMIP name of the model variable
          file: "{year}/x_{year}_{month:02d}.nc"   # glob allowed
          name: x                    # variable in the file
          factor: 1.0                # value = file*factor + offset,
          offset: 0.0                #   giving the target units
          per_day: false             # monthly totals: divide by days
          levels: false              # true: pressure levels (hPa)
          time: index                # optional: time axis is the month
                                     #   index (1-12) in a yearly file

{year} and {month} in the templates select the files of the period; a
template without them is a single file. Time is taken from the files, or
from the template month when the file has no time axis."""

import calendar
import glob
import os
import numpy as np
import xarray as xr
import yaml

from .grid import regular
from .cache import files_signature

_LAT = ('lat', 'latitude', 'y')
_LON = ('lon', 'longitude', 'x')
_LEV = ('pressure_level', 'level', 'plev', 'lev', 'isobaricInhPa')
_TIME = ('time', 'valid_time', 't')


def load_catalog(path):
    with open(path) as f:
        return yaml.safe_load(f)['datasets']


class Reference:

    def __init__(self, name, spec, var):
        if var not in spec.get('variables', {}):
            raise KeyError(f'{var} not in dataset {name}')
        self.name = name
        self.var = var
        self.root = spec.get('root', '')
        self.frequency = spec.get('frequency', 'mon')
        self.surface = spec.get('surface', 'all')
        self.vs = spec['variables'][var]
        self._grid = None

    def files(self, years, months=None):
        tmpl = os.path.join(self.root, self.vs['file'])
        out = []
        if '{month' in tmpl:
            for y in years:
                for m in (months or range(1, 13)):
                    out += [(y, m, f) for f in sorted(glob.glob(
                        tmpl.format(year=y, month=m)))]
        elif '{year' in tmpl:
            for y in years:
                out += [(y, None, f) for f in sorted(glob.glob(
                    tmpl.format(year=y)))]
        else:
            out = [(None, None, f) for f in sorted(glob.glob(tmpl))]
        return out

    def sources(self, years, months=None):
        return files_signature([f for _, _, f in self.files(years, months)])

    def _norm(self, da, y, m):
        ren = {}
        for d in da.dims:
            if d in _LAT:
                ren[d] = 'lat'
            elif d in _LON:
                ren[d] = 'lon'
            elif d in _LEV:
                ren[d] = 'plev'
            elif d in _TIME:
                ren[d] = 'time'
        da = da.rename(ren)
        da = da.drop_vars([c for c in da.coords if c not in da.dims],
                          errors='ignore')
        drop = [d for d in da.dims if d not in ('time', 'lat', 'lon', 'plev')
                and da.sizes[d] == 1]
        if drop:
            da = da.squeeze(drop, drop=True)
        if 'time' not in da.dims:
            da = da.expand_dims(time=[np.datetime64(f'{y:04d}-{m:02d}-01')])
        if da['lon'].max() > 180.0:
            da = da.assign_coords(lon=((da['lon'] + 180.0) % 360.0) - 180.0)
        da = da.sortby('lon').sortby('lat')
        return da

    def open(self, years, bbox=None, months=None):
        """Data of the period in the target units, cropped to bbox
        (lat0, lat1, lon0, lon1)."""
        parts = []
        for y, m, f in self.files(years, months):
            if self.vs.get('time') == 'index':
                ds = xr.open_dataset(f, chunks={}, decode_times=False)
                tname = next(d for d in ds[self.vs['name']].dims
                             if d in _TIME)
                n = ds.sizes[tname]
                ds = ds.assign_coords({tname: np.array(
                    [np.datetime64(f'{y:04d}-{k + 1:02d}-01')
                     for k in range(n)])})
            else:
                ds = xr.open_dataset(f, chunks={})
            da = self._norm(ds[self.vs['name']], y, m)
            if bbox is not None:
                da = da.sel(lat=slice(bbox[0], bbox[1]),
                            lon=slice(bbox[2], bbox[3]))
            yrs = da['time.year']
            keep = yrs.isin(list(years))
            if months:
                keep = keep & da['time.month'].isin(list(months))
            da = da.sel(time=keep)
            if da.sizes['time']:
                parts.append(da)
        if not parts:
            raise FileNotFoundError(f'No {self.var} data in {self.name} '
                                    f'for {years[0]}-{years[-1]}')
        da = xr.concat(parts, dim='time') if len(parts) > 1 else parts[0]
        if self.vs.get('levels'):
            da = da.assign_coords(plev=da['plev'].astype(float))
        f = float(self.vs.get('factor', 1.0))
        o = float(self.vs.get('offset', 0.0))
        da = da.astype('float64') * f + o
        if self.vs.get('per_day'):
            ndays = xr.DataArray(
                [calendar.monthrange(int(t.dt.year), int(t.dt.month))[1]
                 for t in da['time']], dims='time')
            da = da / ndays
        da.name = self.var
        return da

    def monthly(self, years, bbox=None, months=None):
        da = self.open(years, bbox, months)
        if self.frequency != 'mon':
            da = da.resample(time='MS').mean()
        else:
            da = da.assign_coords(time=da['time'].values.astype(
                'datetime64[M]').astype('datetime64[ns]'))
        return da.load()

    def daily(self, years, bbox=None):
        if self.frequency != 'day':
            raise ValueError(f'{self.name} has no daily data')
        return self.open(years, bbox).load()

    def grid(self, da):
        return regular(da['lat'].values, da['lon'].values, self.name)
