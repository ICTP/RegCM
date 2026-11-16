#!/usr/bin/env python3
#
#    This file is part of ICTP RegCM.
#
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Build the initial profiles and forcing of the coupled single column
# test (scm_coupled) for land boxes: the monthly mean ERA5 state at one
# 6-hourly time, averaged over the land points of a box.
#
# Usage: make_amazon_case.py [--era5 DIR] [--month YYYY-MM] [--hour H]
#                            [--qr K/day] [-o FILE]
#
# Columns: NW Amazon (5S-5N, 75-60W) and S Amazon (12-5S, 65-50W).
# Vertical grid: the 32 MOLOCH levels of the RegCM South America test, on
# the box mean orography; T, q, u, v interpolated in height from the
# ERA5 pressure levels, below the lowest one above the ground a well
# mixed layer (as in make_profiles.py); hydrostatic pressure from the
# box mean surface pressure. Forcing written in the profile file:
#   omega   ERA5 monthly mean (vertical advection in scm_coupled)
#   tten    a constant radiative cooling --qr (default -1.5 K/day) up to
#           200 hPa, going to zero at 100 hPa
#   qten    zero
# Surface fluxes: ERA5 monthly means over the box (scm_coupled applies a
# half sine diurnal cycle). The ERA5 hourly total and convective rain,
# and CAPE, of the box by local hour are written to FILE_era5.txt for
# scm_coupled.py. The box mean profiles at 0, 6, 12 and 18 UTC, the
# relaxation targets of scm_coupled (tgtfile), are written to
# FILE_target.txt.

import argparse
import datetime as dt
import os
import numpy as np
import netCDF4 as nc
from make_profiles import A, B, RD, CP, G, EPS

BOXES = [('NW_Amazon', -5.0, 5.0, -75.0, -60.0),
         ('S_Amazon', -12.0, -5.0, -65.0, -50.0)]


def boxsel(ds, la0, la1, lo0, lo1):
    la = np.asarray(ds['latitude'][:], float)
    lo = np.asarray(ds['longitude'][:], float)
    lo = np.where(lo > 180.0, lo - 360.0, lo)
    j = np.where((la >= la0) & (la <= la1))[0]
    i = np.where((lo >= lo0) & (lo <= lo1))[0]
    return slice(j[0], j[-1] + 1), slice(i[0], i[-1] + 1)


def f8(x):
    return np.asarray(np.ma.filled(x, np.nan), float)


def column(era, month, hour, qr, name, la0, la1, lo0, lo1):
    ym = month.strftime('%Y_%m')
    d6 = lambda v: nc.Dataset('%s/%04d/%s_%s.nc' % (era, month.year, v, ym))
    dm = lambda v: nc.Dataset('%s/monthly/%04d/%s_%s.nc' %
                              (era, month.year, v, ym))
    dh = lambda v: nc.Dataset('%s/hourly/%s_%s.nc' % (era, v, ym))
    lsm = nc.Dataset('%s/fixed/landmask.nc' % era)
    js, is_ = boxsel(lsm, la0, la1, lo0, lo1)
    land = f8(lsm['lsm'][..., js, is_]).squeeze() > 0.5
    oro = f8(nc.Dataset('%s/fixed/orography_m.nc' % era)['z'][..., js, is_]
             ).squeeze()

    def bmean(x):
        return np.nanmean(x[..., land], axis=-1)

    ta = d6('tatm')
    t = nc.num2date(ta['valid_time'][:], ta['valid_time'].units)
    it = [k for k, x in enumerate(t) if x.hour == hour]
    plev = f8(ta['pressure_level'][:]) * 100.0
    prof = {}
    for vn, n in (('tatm', 't'), ('qhum', 'q'), ('geop', 'z'),
                 ('uwnd', 'u'), ('vwnd', 'v')):
        d = d6(vn)
        prof[n] = bmean(f8(d[n][it, :, js, is_]).mean(axis=0))
    T, q, zg = prof['t'], prof['q'], prof['z'] / G
    u, v = prof['u'], prof['v']
    om = bmean(f8(dm('wwnd')['w'][..., js, is_]).squeeze())
    hps = dh('ps')
    th = nc.num2date(hps['valid_time'][:], hps['valid_time'].units)
    ih = [k for k, x in enumerate(th) if x.hour == hour]
    ps = float(bmean(f8(hps['sp'][ih, js, is_]).mean(axis=0)))
    topo = max(float(bmean(oro)), 0.0)
    lh = -float(bmean(f8(dm('lflx')['slhf'][..., js, is_]).squeeze())) \
        / 86400.0
    sh = -float(bmean(f8(dm('sflx')['ishf'][..., js, is_]).squeeze()))
    # ERA5 hourly rain and CAPE by local hour (box centre longitude)
    lon = 0.5 * (lo0 + lo1)
    ref = {}
    for vn, n, fac in (('pr', 'tp', 24000.0), ('prc', 'crr', 86400.0),
                      ('cape', 'cape', 1.0)):
        d = dh(vn)
        x = bmean(f8(d[n][:, js, is_]))
        tt = nc.num2date(d['valid_time'][:], d['valid_time'].units)
        lst = np.array([(y.hour + lon / 15.0) % 24.0 for y in tt])
        # accumulated total rain refers to the hour ending at valid time
        if vn == 'pr':
            lst = (lst - 0.5) % 24.0
        ref[vn] = np.array([np.nanmean(x[np.floor(lst) == h] * fac)
                           for h in range(24)])
    # model levels above the ground, top first; well mixed layer below
    # the lowest ERA5 level above the ground
    zsea = np.array(A) + np.array(B) * topo
    zf = zsea - topo
    zh = np.empty(len(zf) + 1)
    zh[1:-1] = 0.5 * (zf[:-1] + zf[1:]); zh[-1] = 0.0
    zh[0] = zf[0] + 0.5 * (zf[0] - zf[1])
    ok = (plev < ps) & np.isfinite(T)
    zg = zg[ok]; T = T[ok]; q = q[ok]; u = u[ok]; v = v[ok]
    kl = int(np.argmin(zg))
    thl = T[kl] * (1e5 / plev[ok][kl]) ** (RD / CP)
    zs = np.linspace(topo, zg[kl], 6)[:-1]
    psx = plev[ok][kl] * np.exp(G * (zg[kl] - zs) /
                                 (RD * T[kl] * (1 + EPS * q[kl])))
    zg = np.append(zg, zs); T = np.append(T, thl * (psx / 1e5) ** (RD / CP))
    q = np.append(q, np.full(5, q[kl]))
    u = np.append(u, np.full(5, u[kl])); v = np.append(v, np.full(5, v[kl]))
    o = np.argsort(zg)
    Tm = np.interp(zsea, zg[o], T[o]); qm = np.interp(zsea, zg[o], q[o])
    um = np.interp(zsea, zg[o], u[o]); vm = np.interp(zsea, zg[o], v[o])
    tvh = np.interp(zh + topo, zg[o], (T * (1 + EPS * q))[o])
    ph = np.empty(len(zh)); ph[-1] = ps
    for k in range(len(zh) - 2, -1, -1):
        tv = 0.5 * (tvh[k] + tvh[k + 1])
        ph[k] = ph[k + 1] * np.exp(-G * (zh[k] - zh[k + 1]) / (RD * tv))
    pf = np.exp(0.5 * (np.log(ph[:-1]) + np.log(ph[1:])))
    omm = np.interp(np.log(pf), np.log(plev[::-1]), om[::-1])
    tten = np.where(pf >= 2.0e4, qr,
                    qr * np.clip((pf - 1.0e4) / 1.0e4, 0.0, 1.0)) / 86400.0
    when = dt.datetime(month.year, month.month, 1, hour)
    return dict(name=name, lat=0.5 * (la0 + la1), lon=lon, land=1,
                when=when.strftime('%Y%m%d%H'), ps=ps, sh=sh, lh=lh,
                topo=topo, zf=zf, pf=pf, T=Tm, q=qm, u=um, v=vm, om=omm,
                tten=tten, zh=zh, ph=ph, ref=ref)


def main():
    ap = argparse.ArgumentParser(
        description='Coupled single column test case from ERA5 box means')
    ap.add_argument('--era5', default='/home/clima-archive5/ERA5')
    ap.add_argument('--month', default='2023-09')
    ap.add_argument('--hour', type=int, default=12,
                    help='UTC hour of the initial state (default 12)')
    ap.add_argument('--qr', type=float, default=-1.5,
                    help='radiative heating, K/day (default -1.5)')
    ap.add_argument('-o', '--output', default='scm_amazon.txt')
    a = ap.parse_args()
    month = dt.datetime.strptime(a.month, '%Y-%m')
    cols = [column(a.era5, month, a.hour, a.qr, *b) for b in BOXES]
    write(a.output, cols, 'ERA5 %s box means at %02d UTC' %
          (a.month, a.hour))
    tcols = [column(a.era5, month, h, a.qr, *b) for b in BOXES
             for h in (0, 6, 12, 18)]
    write(os.path.splitext(a.output)[0] + '_target.txt', tcols,
          'relaxation targets, ERA5 %s box means at 0, 6, 12, 18 UTC' %
          a.month)
    ref = os.path.splitext(a.output)[0] + '_era5.txt'
    with open(ref, 'w') as f:
        f.write('# ERA5 %s box means by local hour: total rain, '
                'convective rain (mm/day), CAPE (J/kg)\n' % a.month)
        for c in cols:
            for h in range(24):
                f.write('%s %2d %8.3f %8.3f %8.1f\n' % (
                    c['name'], h, c['ref']['pr'][h], c['ref']['prc'][h],
                    c['ref']['cape'][h]))
    for c in cols:
        print('%-10s ps %7.1f topo %5.0f SH %6.1f LH %6.1f ERA5 rain '
              '%5.2f conv %5.2f mm/day CAPE max %5.0f' % (
                  c['name'], c['ps'] / 100, c['topo'], c['sh'], c['lh'],
                  c['ref']['pr'].mean(), c['ref']['prc'].mean(),
                  c['ref']['cape'].max()))


def write(fname, cols, title):
    kz = len(A)
    with open(fname, 'w') as f:
        f.write('# RegCM coupled single column test case, %s, built by '
                'make_amazon_case.py\n' % title)
        f.write('# per column: name lat lon land date ps(Pa) shf lhf '
                '(W/m2, daily mean) zpbl cape conv_rain (unused) '
                'orography(m)\n')
        f.write('# then kz lines top to bottom: z(m above ground) p(Pa) '
                'T(K) q(kg/kg) u v(m/s) omega(Pa/s) tten(K/s) '
                'qten(kg/kg/s)\n')
        f.write('# then kz+1 half levels: z(m above ground) p(Pa)\n')
        f.write('%d %d\n' % (len(cols), kz))
        for c in cols:
            f.write('%s %.2f %.2f %d %s %.2f %.3f %.3f %.1f %.1f %.3f %.1f\n'
                    % (c['name'], c['lat'], c['lon'], c['land'], c['when'],
                       c['ps'], c['sh'], c['lh'], 0.0, 0.0, 0.0, c['topo']))
            for k in range(kz):
                f.write('%10.2f %10.2f %8.3f %12.5e %8.3f %8.3f %12.5e '
                        '%12.5e 0.0\n' % (
                            c['zf'][k], c['pf'][k], c['T'][k], c['q'][k],
                            c['u'][k], c['v'][k], c['om'][k], c['tten'][k]))
            for k in range(kz + 1):
                f.write('%10.2f %10.2f\n' % (c['zh'][k], c['ph'][k]))


if __name__ == '__main__':
    main()
