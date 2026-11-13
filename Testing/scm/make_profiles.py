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
# Build the single column test profiles (scm_profiles.txt) from ERA5.
#
# Usage: make_profiles.py [--era5 DIR] [--date YYYY-MM-DD] [-o FILE]
#
# Twelve columns, six per hemisphere, land and ocean, chosen to cover the
# convective regimes: deep convection over land and ocean, trade wind
# shallow cumulus, stratocumulus, dry land. The deep convection columns
# are taken at the 6-hourly time of the month with the largest ERA5
# convective rain rate at the point, the others at a fixed time of --date.
#
# Vertical grid: the 32 MOLOCH levels of the RegCM South America test
# (z = a + b * orography), heights above the ground. ERA5 6-hourly
# temperature, specific humidity, geopotential and wind on pressure levels
# above the ground are interpolated in height; below the lowest of them a
# well mixed surface layer (constant potential temperature, humidity and
# wind) as after the PBL scheme of a model; the pressure is integrated
# hydrostatically
# from the ERA5 surface pressure. The vertical velocity is the ERA5
# monthly mean omega (only used to find the level of maximum ascent).
#
# Surface fluxes are the ERA5 monthly means. Over land they are scaled by
# a half sine diurnal cycle (peak = pi * daily mean) at the local solar
# time of the profile. The PBL height is diagnosed with a bulk Richardson
# number (critical value 0.25) from the lowest level, with a 0.5 K parcel
# excess over land in daytime.
#
# ERA5 CAPE and convective rain rate at the same time are written as a
# reference (not a truth) for the scheme diagnostics.
#
# The input tendencies are zero, as in RegCM where the cumulus scheme is
# the first physics call of a time step.

import argparse
import datetime as dt
import numpy as np
import netCDF4 as nc

RD, CP, G, LV = 287.04, 1004.64, 9.80665, 2.5008e6
EPS = 0.608

# MOLOCH vertical coordinate of the RegCM South America test (top first)
A = [28251.5756479703, 25047.5736278179, 22197.8752403678, 19663.3010271314,
     17409.0040563733, 15403.9908254211, 13620.6951423589, 12034.5991285594,
     10623.8961313589, 9369.19091238337, 8253.23298952612, 7260.67946639282,
     6377.88408844348, 5592.70962564118, 4894.36100212555, 4273.23687867203,
     3720.79764739849, 3229.44802382694, 2792.43262210149, 2403.74307766461,
     2058.03544045537, 1750.55670289906, 1477.07945254834, 1233.84375093925,
     1017.50543957553, 825.090162319381, 653.952472059271, 501.739459427347,
     366.35840351107, 245.947999800539, 138.852769795099, 43.600300436299]
B = [0.00072479248046875, 0.00638580322265625, 0.0173568725585938,
     0.0332717895507812, 0.0537643432617188, 0.0784683227539062,
     0.107017517089844, 0.139045715332031, 0.174186706542969,
     0.212074279785156, 0.252342224121094, 0.294624328613281,
     0.338554382324219, 0.383766174316406, 0.429893493652344,
     0.476570129394531, 0.523429870605469, 0.570106506347656,
     0.616233825683594, 0.661445617675781, 0.705375671386719,
     0.747657775878906, 0.787925720214844, 0.825813293457031,
     0.860954284667969, 0.892982482910156, 0.921531677246094,
     0.946235656738281, 0.966728210449219, 0.982643127441406,
     0.993614196777344, 0.999275207519531]

# name, latitude, longitude, hour UTC; hour None: the 6-hourly time of the
# month with the largest ERA5 convective rain rate at the point
COLUMNS = [
    ('Sahel_land', 13.0, 2.0, None),
    ('India_land', 23.0, 80.0, None),
    ('GreatPlains_land', 38.0, -98.0, None),
    ('WarmPool_ocean', 8.0, 140.0, None),
    ('California_Sc_ocean', 30.0, -125.0, 12),
    ('NAtl_trades_ocean', 15.0, -55.0, 12),
    ('Amazon_land', -5.0, -60.0, None),
    ('Congo_land', -2.0, 22.0, None),
    ('Australia_land', -25.0, 135.0, 6),
    ('SPCZ_ocean', -10.0, -170.0, None),
    ('SEAtl_Sc_ocean', -15.0, 5.0, 0),
    ('SAtl_trades_ocean', -10.0, -30.0, 12),
]


def wettest(era, lat, lon, month):
    """6-hourly time of the month with the largest convective rain"""
    ym = '%04d_%02d' % (month.year, month.month)
    d = nc.Dataset('%s/hourly/prc_%s.nc' % (era, ym))
    t = nc.num2date(d['valid_time'][:], d['valid_time'].units)
    crr = point(d, 'crr', lat, lon)
    k = [i for i, x in enumerate(t) if x.hour % 6 == 0]
    i = k[int(np.argmax(crr[k]))]
    return dt.datetime(t[i].year, t[i].month, t[i].day, t[i].hour)


def point(ds, var, lat, lon, it=None):
    la = np.asarray(ds['latitude'][:], float)
    lo = np.asarray(ds['longitude'][:], float)
    j = int(np.argmin(np.abs(la - lat)))
    i = int(np.argmin(np.abs(((lo - lon) + 180.0) % 360.0 - 180.0)))
    v = ds[var]
    if it is None:
        x = v[..., j, i]
    else:
        x = v[it, ..., j, i]
    return np.asarray(np.ma.filled(x, np.nan), float).squeeze()


def tindex(ds, when):
    t = nc.num2date(ds['valid_time'][:], ds['valid_time'].units)
    for k, x in enumerate(t):
        if (x.year, x.month, x.day, x.hour) == \
           (when.year, when.month, when.day, when.hour):
            return k
    raise ValueError('time %s not found' % when)


def column(era, name, lat, lon, when):
    ym = '%04d_%02d' % (when.year, when.month)
    d6 = lambda v: nc.Dataset('%s/%04d/%s_%s.nc' % (era, when.year, v, ym))
    dh = lambda v: nc.Dataset('%s/hourly/%s_%s.nc' % (era, v, ym))
    dm = lambda v: nc.Dataset('%s/monthly/%04d/%s_%s.nc' %
                              (era, when.year, v, ym))
    ta = d6('tatm'); it = tindex(ta, when)
    plev = np.asarray(ta['pressure_level'][:], float) * 100.0
    T = point(ta, 't', lat, lon, it)
    q = point(d6('qhum'), 'q', lat, lon, it)
    zg = point(d6('geop'), 'z', lat, lon, it) / G
    u = point(d6('uwnd'), 'u', lat, lon, it)
    v = point(d6('vwnd'), 'v', lat, lon, it)
    om = point(dm('wwnd'), 'w', lat, lon)
    hps = dh('ps'); ih = tindex(hps, when)
    ps = float(point(hps, 'sp', lat, lon, ih))
    cape = float(point(dh('cape'), 'cape', lat, lon, ih))
    crr = float(point(dh('prc'), 'crr', lat, lon, ih)) * 86400.0
    oro = nc.Dataset('%s/fixed/orography_m.nc' % era)
    topo = max(float(point(oro, 'z', lat, lon)), 0.0)
    lsm = nc.Dataset('%s/fixed/landmask.nc' % era)
    land = float(point(lsm, 'lsm', lat, lon)) > 0.5
    lh = -float(point(dm('lflx'), 'slhf', lat, lon)) / 86400.0
    sh = -float(point(dm('sflx'), 'ishf', lat, lon))
    lst = (when.hour + lon / 15.0) % 24.0
    if land:
        fac = np.pi * np.cos(np.pi * (lst - 12.5) / 12.0)
        fac = max(fac, 0.1) if abs(lst - 12.5) < 6.0 else 0.1
        lh *= fac; sh *= fac
    # model levels: heights above the ground, top first
    zsea = np.array(A) + np.array(B) * topo
    zf = zsea - topo
    zh = np.empty(len(zf) + 1)
    zh[1:-1] = 0.5 * (zf[:-1] + zf[1:]); zh[-1] = 0.0
    zh[0] = zf[0] + 0.5 * (zf[0] - zf[1])
    # ERA5 levels above the ground; below the lowest one a well mixed
    # surface layer (constant potential temperature and humidity), as a
    # model column presents itself after the PBL scheme
    ok = plev < ps
    zg = zg[ok]; T = T[ok]; q = q[ok]; u = u[ok]; v = v[ok]
    kl = int(np.argmin(zg))
    th = T[kl] * (1e5 / plev[ok][kl]) ** (RD / CP)
    zs = np.linspace(topo, zg[kl], 6)[:-1]
    psx = plev[ok][kl] * np.exp(G * (zg[kl] - zs) /
                                 (RD * T[kl] * (1 + EPS * q[kl])))
    zg = np.append(zg, zs); T = np.append(T, th * (psx / 1e5) ** (RD / CP))
    q = np.append(q, np.full(5, q[kl]))
    u = np.append(u, np.full(5, u[kl])); v = np.append(v, np.full(5, v[kl]))
    o = np.argsort(zg)
    Tm = np.interp(zsea, zg[o], T[o]); qm = np.interp(zsea, zg[o], q[o])
    um = np.interp(zsea, zg[o], u[o]); vm = np.interp(zsea, zg[o], v[o])
    # hydrostatic pressure from the surface
    tvh = np.interp(zh + topo, zg[o], (T * (1 + EPS * q))[o])
    ph = np.empty(len(zh)); ph[-1] = ps
    for k in range(len(zh) - 2, -1, -1):
        tv = 0.5 * (tvh[k] + tvh[k + 1])
        ph[k] = ph[k + 1] * np.exp(-G * (zh[k] - zh[k + 1]) / (RD * tv))
    pf = np.exp(0.5 * (np.log(ph[:-1]) + np.log(ph[1:])))
    omm = np.interp(np.log(pf), np.log(plev[::-1]), om[::-1])
    # PBL height: bulk Richardson number from the lowest level (critical
    # value 0.25, as in the Holtslag scheme), with a 0.5 K excess of the
    # surface parcel over land in daytime
    thv = Tm * (1 + EPS * qm) * (1e5 / pf) ** (RD / CP)
    th0 = thv[-1] + (0.5 if (land and fac > 0.1) else 0.0)
    zpbl = 100.0
    rib = np.zeros(len(zf))
    for k in range(len(zf) - 2, -1, -1):
        du2 = max((um[k] - um[-1]) ** 2 + (vm[k] - vm[-1]) ** 2, 1.0)
        rib[k] = G * (zf[k] - zf[-1]) * (thv[k] - th0) / (th0 * du2)
        if rib[k] > 0.25:
            zpbl = max(zf[k + 1] + (zf[k] - zf[k + 1]) *
                       (0.25 - rib[k + 1]) / (rib[k] - rib[k + 1]), 100.0)
            break
    return dict(name=name, lat=lat, lon=lon, land=int(land),
                when=when.strftime('%Y%m%d%H'), ps=ps, sh=sh, lh=lh,
                zpbl=zpbl, cape=cape, crr=crr, topo=topo, zf=zf, pf=pf,
                T=Tm, q=qm, u=um, v=vm, om=omm, zh=zh, ph=ph)


def main():
    ap = argparse.ArgumentParser(
        description='Build the RegCM single column test profiles from ERA5')
    ap.add_argument('--era5', default='/home/clima-archive5/ERA5')
    ap.add_argument('--date', default='2023-09-15')
    ap.add_argument('-o', '--output', default='scm_profiles.txt')
    a = ap.parse_args()
    day = dt.datetime.strptime(a.date, '%Y-%m-%d')
    cols = [column(a.era5, n, la, lo,
                   day + dt.timedelta(hours=h) if h is not None
                   else wettest(a.era5, la, lo, day))
            for n, la, lo, h in COLUMNS]
    kz = len(A)
    with open(a.output, 'w') as f:
        f.write('# RegCM single column test profiles from ERA5, built by '
                'make_profiles.py\n')
        f.write('# per column: name lat lon land date ps(Pa) shf lhf '
                '(W/m2, upward) zpbl(m) era5_cape(J/kg) '
                'era5_conv_rain(mm/day) orography(m)\n')
        f.write('# then kz lines top to bottom: z(m above ground) p(Pa) '
                'T(K) q(kg/kg) u v(m/s) omega(Pa/s) tten(K/s) '
                'qten(kg/kg/s)\n')
        f.write('# then kz+1 half levels: z(m above ground) p(Pa)\n')
        f.write('%d %d\n' % (len(cols), kz))
        for c in cols:
            f.write('%s %.2f %.2f %d %s %.2f %.3f %.3f %.1f %.1f %.3f %.1f\n'
                    % (c['name'], c['lat'], c['lon'], c['land'], c['when'],
                       c['ps'], c['sh'], c['lh'], c['zpbl'], c['cape'],
                       c['crr'], c['topo']))
            for k in range(kz):
                f.write('%10.2f %10.2f %8.3f %12.5e %8.3f %8.3f %12.5e '
                        '0.0 0.0\n' % (c['zf'][k], c['pf'][k], c['T'][k],
                                       c['q'][k], c['u'][k], c['v'][k],
                                       c['om'][k]))
            for k in range(kz + 1):
                f.write('%10.2f %10.2f\n' % (c['zh'][k], c['ph'][k]))
    for c in cols:
        print('%-20s %s land %d ps %7.1f SH %6.1f LH %6.1f zpbl %6.0f '
              'ERA5 CAPE %6.0f conv rain %6.2f mm/day' % (
                  c['name'], c['when'], c['land'], c['ps'] / 100, c['sh'],
                  c['lh'],
                  c['zpbl'], c['cape'], c['crr']))


if __name__ == '__main__':
    main()
