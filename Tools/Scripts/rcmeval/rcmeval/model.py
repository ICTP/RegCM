#
#    This file is part of ICTP RegCM.
#
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
"""Readers of RegCM simulations: the RegCM output files
(DOMAIN_TYPE.YYYYMMDDHH.nc) and the CORDEX archive written by pycordexer.

A run gives its grid and, for a variable and a period, the monthly means,
the daily means and the monthly mean profiles on pressure levels, all in
the target units (units.py)."""

import datetime as dt
import glob
import os
import re
import numpy as np
import netCDF4 as nc
import xarray as xr

from .grid import Grid
from .units import to_target, qsat
from .cache import files_signature

RD = 287.04
CPD = 1004.64
G = 9.80665

_RX = re.compile(r'^(?P<dom>.+)_(?P<typ>[A-Z]{3})\.(?P<date>\d{10})\.nc$')
MEAN_PREF = ('SRF', 'RAD', 'STS', 'ATM')
DAILY_PREF = ('STS', 'SRF')
PROFILE_VARS = ('ta', 'hus', 'hur', 'ua', 'va', 'wap')


def _months(years, months=None):
    out = []
    for y in years:
        for m in (months or range(1, 13)):
            out.append((y, m))
    return out


def _center_time(ds):
    """Time at the middle of the averaging interval when bounds exist."""
    tb = [v for v in ('time_bnds', 'time_bounds') if v in ds.variables]
    if tb:
        b = ds[tb[0]].values
        mid = b[:, 0] + (b[:, 1] - b[:, 0]) / 2
        ds = ds.assign_coords(time=mid).drop_vars(tb[0])
    return ds


def _squeeze(da, keep):
    """Drop singleton dimensions (heights, ensemble members) other than
    those in keep."""
    drop = [d for d in da.dims if d not in keep and da.sizes[d] == 1]
    return da.squeeze(drop, drop=True) if drop else da


def _monthly(da):
    m = da.resample(time='MS').mean()
    return m.dropna('time', how='all')


def _interp_p(pfull, f, plev):
    """Linear in log p interpolation of f (k, n) at the pressures pfull
    (k, n, increasing with k) to plev (Pa); NaN outside."""
    lp = np.log(pfull)
    out = np.full((len(plev), f.shape[1]), np.nan)
    for i, p in enumerate(np.log(plev)):
        k = (lp < p).sum(axis=0)
        ok = (k >= 1) & (k <= lp.shape[0] - 1)
        c = np.flatnonzero(ok)
        kk = k[c]
        w = (p - lp[kk - 1, c]) / (lp[kk, c] - lp[kk - 1, c])
        out[i, c] = f[kk - 1, c] + w * (f[kk, c] - f[kk - 1, c])
    return out


class RegCMRun:
    """Simulation in the RegCM output format."""

    kind = 'regcm'

    def __init__(self, name, path, boundary=10):
        self.name = name
        self.path = path
        self.boundary = boundary
        self.files = {}
        for f in sorted(glob.glob(os.path.join(path, '*_???.??????????.nc'))):
            m = _RX.match(os.path.basename(f))
            if m:
                d = dt.datetime.strptime(m['date'], '%Y%m%d%H')
                self.files.setdefault(m['typ'], []).append((d, f))
                self.domain = m['dom']
        if not self.files:
            raise FileNotFoundError(f'No RegCM output files in {path}')
        self.vars = {}
        for typ, fl in self.files.items():
            with nc.Dataset(fl[0][1]) as d:
                for v in d.variables:
                    if d[v].ndim >= 3 and v not in ('xlat_bnds',
                                                     'xlon_bnds'):
                        self.vars.setdefault(v, []).append(typ)
        self._grid = None

    def has(self, var):
        return var in self.vars or (var in ('hur', 'wap') and
                                    'ta' in self.vars)

    def grid(self):
        if self._grid is None:
            typ = 'SRF' if 'SRF' in self.files else next(iter(self.files))
            with nc.Dataset(self.files[typ][0][1]) as d:
                lat = d['xlat'][:]
                lon = d['xlon'][:]
                if 'areacella' in d.variables:
                    area = d['areacella'][:]
                else:
                    ds = float(d.getncattr('grid_size_in_meters'))
                    area = np.full(lat.shape, ds * ds)
                land = (np.asarray(d['mask'][:]) > 0.5).astype(float) \
                    if 'mask' in d.variables else None
            valid = np.ones(lat.shape, bool)
            b = self.boundary
            if b > 0:
                valid[:b, :] = valid[-b:, :] = False
                valid[:, :b] = valid[:, -b:] = False
            self._grid = Grid(lat, lon, area, valid, land, self.name)
        return self._grid

    def _select(self, typ, years):
        start = dt.datetime(min(years), 1, 1) - dt.timedelta(days=32)
        end = dt.datetime(max(years) + 1, 1, 1)
        return [f for d, f in self.files[typ] if start <= d < end]

    def _type(self, var, pref):
        for t in pref:
            if t in self.vars.get(var, []):
                return t
        if var in self.vars:
            return self.vars[var][0]
        raise KeyError(f'{var} not in {self.name}')

    def sources(self, var, years, kind='mean'):
        if var in ('hur', 'wap'):
            return files_signature(self._select('ATM', years))
        pref = DAILY_PREF if kind == 'daily' else MEAN_PREF
        return files_signature(self._select(self._type(var, pref), years))

    def _open(self, var, typ, years):
        files = self._select(typ, years)

        def keep(ds):
            v = [var] + [x for x in ('time_bnds',) if x in ds.variables]
            return _center_time(ds[v])
        ds = xr.open_mfdataset(files, combine='nested', concat_dim='time',
                               data_vars='minimal', coords='minimal',
                               compat='override', preprocess=keep)
        da = _squeeze(ds[var], ('time', 'kz', 'plev'))
        da = da.sel(time=da['time.year'].isin(list(years)))
        return to_target(da, var)

    def monthly(self, var, years):
        """Monthly means (time, y, x)."""
        da = self._open(var, self._type(var, MEAN_PREF), years)
        return _monthly(da).load()

    def daily(self, var, years):
        """Daily means (time, y, x)."""
        typ = self._type(var, DAILY_PREF)
        da = self._open(var, typ, years)
        if typ != 'STS':
            da = da.resample(time='1D').mean()
        return da.load()

    def profile_monthly(self, var, years, plev):
        """Monthly mean of var interpolated to the pressures plev (hPa) at
        every output time of the ATM files: (time, plev, y, x)."""
        need = {'ta': ['ta'], 'hus': ['hus'], 'ua': ['ua'], 'va': ['va'],
                'hur': ['ta', 'hus'], 'wap': ['wa', 'ta', 'hus']}[var]
        pa = np.asarray(plev, float) * 100.0
        acc = {}
        for f in self._select('ATM', years):
            with xr.open_dataset(f) as d:
                d = d.sel(time=d['time.year'].isin(list(years)))
                if d.sizes['time'] == 0:
                    continue
                if 'pai' in d:
                    p = 1.0e5 * d['pai'].values ** (CPD / RD)
                else:
                    p = (d['ptop'].values * 100.0 + d['sigma'].values[
                        None, :, None, None] * (d['ps'].values[:, None] -
                                                d['ptop'].values * 100.0))
                flds = {v: d[v].values for v in need}
                nt, nk = p.shape[:2]
                ny, nx = p.shape[2:]
                for it in range(nt):
                    pk = p[it].reshape(nk, -1)
                    if var in ('ta', 'hus', 'ua', 'va'):
                        x = flds[var][it].reshape(nk, -1)
                    elif var == 'hur':
                        t = flds['ta'][it].reshape(nk, -1)
                        q = flds['hus'][it].reshape(nk, -1)
                        x = 100.0 * q / qsat(t, pk)
                    else:
                        t = flds['ta'][it].reshape(nk, -1)
                        q = flds['hus'][it].reshape(nk, -1)
                        rho = pk / (RD * t * (1.0 + 0.608 * q))
                        x = -rho * G * flds['wa'][it].reshape(nk, -1)
                    if var == 'hus':
                        x = x * 1000.0
                    y = _interp_p(pk, x, pa).reshape(len(pa), ny, nx)
                    tt = d['time'].values[it]
                    mkey = np.datetime64(str(tt)[:7] + '-01')
                    s, n = acc.get(mkey, (0.0, 0.0))
                    good = np.isfinite(y)
                    acc[mkey] = (s + np.where(good, y, 0.0),
                                 n + good.astype(float))
        if not acc:
            raise KeyError(f'No ATM data for {var} in {self.name}')
        times = sorted(acc)
        data = np.stack([np.where(acc[t][1] > 0,
                                  acc[t][0] / np.maximum(acc[t][1], 1.0),
                                  np.nan) for t in times])
        return xr.DataArray(data, dims=('time', 'plev', 'y', 'x'),
                            coords={'time': np.array(times),
                                    'plev': np.asarray(plev, float)},
                            name=var)


_FREQS = ('mon', 'day', '6hr', '3hr', '1hr')


class CordexRun:
    """Simulation in the CORDEX archive layout (files var_*_freq_*.nc
    anywhere under path)."""

    kind = 'cordex'

    def __init__(self, name, path, boundary=0):
        self.name = name
        self.path = path
        self.boundary = boundary
        self._grid = None
        self._index = {}
        for f in glob.glob(os.path.join(path, '**', '*.nc'), recursive=True):
            parts = os.path.basename(f)[:-3].split('_')
            freq = next((p for p in parts if p in _FREQS + ('fx',)), None)
            if freq:
                self._index.setdefault((parts[0], freq), []).append(f)
        if not self._index:
            raise FileNotFoundError(f'No CORDEX files under {path}')

    def has(self, var):
        return any(v == var for v, _ in self._index) or (
            var in PROFILE_VARS)

    def _files(self, var, years, pref):
        for fr in pref:
            fl = self._index.get((var, fr))
            if fl:
                sel = []
                for f in sorted(fl):
                    rng = os.path.basename(f)[:-3].split('_')[-1]
                    y0, y1 = int(rng[:4]), int(rng.split('-')[-1][:4])
                    if y1 >= min(years) and y0 <= max(years):
                        sel.append(f)
                if sel:
                    return sel, fr
        raise KeyError(f'{var} not in {self.name}')

    def sources(self, var, years, kind='mean'):
        pref = ('day',) + _FREQS if kind == 'daily' else _FREQS
        return files_signature(self._files(var, years, pref)[0])

    def grid(self):
        if self._grid is None:
            f = next(iter(self._index.values()))[0]
            with nc.Dataset(f) as d:
                lat = d['lat'][:]
                lon = d['lon'][:]
            area = None
            land = None
            for v, tgt in (('areacella', 'area'), ('sftlf', 'land')):
                fx = self._index.get((v, 'fx'))
                if fx:
                    with nc.Dataset(fx[0]) as d:
                        a = np.asarray(d[v][:], float)
                    if tgt == 'area':
                        area = a
                    else:
                        land = a / 100.0 if a.max() > 1.5 else a
            if area is None:
                from .grid import EARTH_RADIUS
                dy = np.gradient(np.radians(lat), axis=0) * EARTH_RADIUS
                dx = (np.gradient(np.radians(lon), axis=1) * EARTH_RADIUS *
                      np.cos(np.radians(lat)))
                area = np.abs(dx * dy)
            valid = np.ones(lat.shape, bool)
            b = self.boundary
            if b > 0:
                valid[:b, :] = valid[-b:, :] = False
                valid[:, :b] = valid[:, -b:] = False
            self._grid = Grid(lat, lon, area, valid, land, self.name)
        return self._grid

    def _open(self, var, years, pref):
        files, fr = self._files(var, years, pref)
        ds = xr.open_mfdataset(files, combine='by_coords',
                               data_vars='minimal', coords='minimal',
                               compat='override')
        da = _squeeze(ds[var], ('time',))
        da = da.sel(time=da['time.year'].isin(list(years)))
        return to_target(da, var), fr

    def monthly(self, var, years):
        da, fr = self._open(var, years, _FREQS)
        return (da if fr == 'mon' else _monthly(da)).load()

    def daily(self, var, years):
        da, fr = self._open(var, years, ('day',) + _FREQS[2:])
        if fr != 'day':
            da = da.resample(time='1D').mean()
        return da.load()

    def profile_monthly(self, var, years, plev):
        name = {'wap': 'wa'}.get(var, var)
        out = []
        for p in plev:
            da, fr = self._open(f'{name}{int(p)}', years, _FREQS)
            da = to_target(da, var, da.attrs.get('units'))
            out.append(da if fr == 'mon' else _monthly(da))
        return xr.concat(out, dim=xr.DataArray(np.asarray(plev, float),
                                               dims='plev')).transpose(
            'time', 'plev', ...).load()


def open_run(name, spec):
    fmt = spec.get('format', 'regcm')
    if fmt == 'regcm':
        return RegCMRun(name, spec['path'], spec.get('boundary', 10))
    if fmt == 'cordex':
        return CordexRun(name, spec['path'], spec.get('boundary', 0))
    raise ValueError(f'Unknown run format {fmt}')
