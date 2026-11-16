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
# Evaluation of the RegCM ocean surface flux schemes with the stand alone
# driver scm_ocean: Zeng (iocnflux 1, with izengrough 1, 2, 3) and COARE
# 3.5 (iocnflux 2), against the NOAA COARE 3.5 algorithm (Edson et al. 2013)
# run on the same inputs. COARE 3.5 is the bulk algorithm fitted to the
# direct covariance flux measurements; it is the reference here.
#
#   scm_ocean.py --coare DIR curves     transfer coefficients against the
#                                       wind (no air-sea temperature
#                                       difference) and fluxes against the
#                                       air-sea temperature difference
#   scm_ocean.py --coare DIR era5       one month of ERA5 6-hourly states
#                                       (1000 hPa level, daily SST), mean
#                                       fluxes over ocean boxes, also
#                                       against the ERA5 monthly fluxes
#
# DIR is the COARE3.5 Python directory of the NOAA PSL COARE algorithm
# repository (https://github.com/NOAA-PSL/COARE-algorithm, directory
# Python/COARE3.5). The program is looked for in --exe, the current
# directory and the directory of this script. The era5 mode needs
# netCDF4 and the ERA5 archive (--era5).

import argparse
import os
import subprocess
import sys
import tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RD, CP, G, LV = 287.04, 1004.64, 9.80665, 2.5008e6
EP2 = 0.62197

# name, iocnflux, izengrough
CONFIGS = [('Zeng ECMWF (rough 1)', 1, 1), ('Zeng Smith (rough 2)', 1, 2),
           ('Zeng COARE3.0 (rough 3)', 1, 3), ('COARE 3.5', 2, 1)]

BOXES = {'Atl ITCZ 2-10N 45-25W': (2, 10, -45, -25),
         'E Pac ITCZ 4-12N 100-85W': (4, 12, -100, -85),
         'S Atl trades 25-5S 30-10W': (-25, -5, -30, -10),
         'SE Pac 30-10S 100-80W': (-30, -10, -100, -80),
         'S Atl 50-35S 60-20W': (-50, -35, -60, -20),
         'all 50S-15N': (-50, 15, -110, -10)}


def esat(t):
    """saturation vapour pressure over water (Pa), Bolton"""
    return 611.2 * np.exp(17.67 * (t - 273.15) / (t - 29.65))


def find_exe(exe, name):
    for c in (exe, os.path.join(os.getcwd(), name),
              os.path.join(HERE, name)):
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.abspath(c)
    return None


def run_scheme(exe, cfg, s, work):
    """s: dict of sample arrays; returns dict of output arrays"""
    name, flx, rough = cfg
    inp = os.path.join(work, 'in.txt')
    out = os.path.join(work, 'out.txt')
    nml = os.path.join(work, 'scm_ocean.in')
    n = len(s['z'])
    with open(inp, 'w') as f:
        f.write('%d\n' % n)
        np.savetxt(f, np.column_stack([s[k] for k in (
            'z', 'u', 'v', 't', 'qv', 'p', 'ps', 'sst', 'zi', 'swnet',
            'lwnet')]), fmt='%.6e')
    with open(nml, 'w') as f:
        f.write(" &scmparam\n  infile = '%s', outfile = '%s',\n"
                "  iocnflux = %d, izengrough = %d,\n /\n" %
                (inp, out, flx, rough))
    r = subprocess.run([exe, nml], cwd=work, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit('scm_ocean failed for %s:\n%s%s' % (name, r.stdout,
                                                    r.stderr))
    x = np.loadtxt(out, comments='#', skiprows=3)
    # stress from the drag (rho Cd U): BATS does not set taux, tauy
    tau = x[:, 7] * np.hypot(s['u'], s['v'])
    return dict(shf=x[:, 0], lhf=x[:, 1], tau=tau, ustar=x[:, 3],
                u10=x[:, 4])


def run_coare(coare_dir, s):
    """NOAA COARE 3.5 on the same samples (bulk SST, cool skin on);
    the latent heat flux is the evaporation times LV, as scm_ocean"""
    sys.path.insert(0, coare_dir)
    from coare35vn import coare35vn
    n = len(s['z'])
    spd = np.hypot(s['u'], s['v'])
    e = s['qv'] / (EP2 + s['qv']) * s['p']
    rh = np.clip(100.0 * e / esat(s['t']), 1.0, 100.0)
    rs = s['swnet'] / 0.945
    rl = (0.97 * 5.67e-8 * s['sst'] ** 4 - s['lwnet']) / 0.97
    out = {k: np.empty(n) for k in ('shf', 'lhf', 'tau', 'ustar')}
    # the sensor heights must be constants: groups of 1 m
    zr = np.round(s['z'])
    for z in np.unique(zr):
        k = np.where(zr == z)[0]
        a = coare35vn(spd[k], s['t'][k] - 273.15, rh[k],
                      s['sst'][k] - 273.15, P=s['ps'][k] / 100.0,
                      Rs=rs[k], Rl=rl[k], zu=float(z), zt=float(z),
                      zq=float(z), lat=s['lat'][k], zi=s['zi'][k],
                      jcool=1)
        # hsb, hlb: upward positive (cooling of the ocean)
        le = (2.501 - 0.00237 * (s['sst'][k] - 273.15)) * 1e6
        out['shf'][k] = a[:, 2]
        out['lhf'][k] = a[:, 3] / le * LV
        out['tau'][k] = a[:, 1]
        out['ustar'][k] = a[:, 0]
    return out


def samples_curves():
    """moist neutral wind sweep at 10 m and a stability sweep at 40 m"""
    u = np.arange(1.0, 25.01, 1.0)
    sets = {}
    for name, z, du, dts in (('wind', 10.0, u, [0.0]),
                             ('stability', 40.0, np.array([3.0, 7.0, 12.0]),
                              np.arange(-3.0, 5.01, 0.5))):
        uu, dd = np.meshgrid(du, dts, indexing='ij')
        uu = uu.ravel(); dd = dd.ravel(); n = len(uu)
        ps = np.full(n, 101000.0)
        sst = np.full(n, 300.0)
        # air at z: potential temperature referred to ps of sst - dts
        p = ps * np.exp(-G * z / (RD * 298.0))
        tpot = sst - dd
        t = tpot * (p / ps) ** (RD / CP)
        e = 0.8 * esat(t)
        qv = EP2 * e / (p - e)
        sets[name] = dict(z=np.full(n, z), u=uu, v=np.zeros(n), t=t, qv=qv,
                          p=p, ps=ps, sst=sst, zi=np.full(n, 600.0),
                          swnet=np.full(n, 0.0), lwnet=np.full(n, 50.0),
                          lat=np.full(n, 10.0), wind=uu, dts=dd)
    return sets


def coeffs(s, o):
    """Cd, Ce at the level, from the fluxes"""
    spd = np.hypot(s['u'], s['v'])
    rho = s['ps'] / (RD * s['t'])
    q = s['qv'] / (1 + s['qv'])
    qs = 0.98 * EP2 * esat(s['sst']) / (s['ps'] - 0.378 * esat(s['sst']))
    cd = o['tau'] / (rho * spd ** 2)
    ce = o['lhf'] / (rho * LV * spd * (qs - q))
    return cd, ce


def mode_curves(a, exe, work):
    sets = samples_curves()
    res = {}
    for cfg in CONFIGS:
        res[cfg[0]] = {k: run_scheme(exe, cfg, s, work)
                       for k, s in sets.items()}
    if a.coare:
        rname = 'COARE 3.5'
        ref = {k: run_coare(a.coare, s) for k, s in sets.items()}
    else:
        rname = 'Zeng ECMWF'
        ref = res['Zeng ECMWF (rough 1)']
    s = sets['wind']
    wsel = [3, 5, 7, 10, 15, 20]
    idx = [int(np.argmin(abs(s['wind'] - w))) for w in wsel]
    cdr, cer = coeffs(s, ref['wind'])
    print('\nTransfer coefficients at 10 m, no air-sea temperature '
          'difference, RH 80%%, as a ratio to %s' % rname)
    print('%-24s %s' % ('U10 (m/s)', ' '.join('%5d' % w for w in wsel)))
    print('%-24s %s | %s' % (rname + ' Cd*1e3',
                             ' '.join('%5.2f' % (cdr[i] * 1e3) for i in idx),
                             'Ce*1e3 ' + ' '.join('%5.2f' % (cer[i] * 1e3)
                                                 for i in idx)))
    for name, r in res.items():
        cd, ce = coeffs(s, r['wind'])
        print('%-24s %s | %s' % (
            name[:24], ' '.join('%5.2f' % (cd[i] / cdr[i]) for i in idx),
            'Ce     ' + ' '.join('%5.2f' % (ce[i] / cer[i]) for i in idx)))
    s = sets['stability']
    print('\nFluxes at 40 m against the air-sea potential temperature '
          'difference (SST - theta), RH 80%: latent / sensible heat flux, '
          'W/m2')
    for w in (3.0, 7.0, 12.0):
        sel = [i for i in range(len(s['wind'])) if s['wind'][i] == w and
               s['dts'][i] in (-2.0, 0.0, 1.0, 2.0, 4.0)]
        print('U = %4.1f m/s, SST - theta: %s' % (
            w, ' '.join('%10.1f' % s['dts'][i] for i in sel)))
        print('  %-22s %s' % (rname, ' '.join(
            '%5.0f/%4.0f' % (ref['stability']['lhf'][i],
                             ref['stability']['shf'][i]) for i in sel)))
        for name, r in res.items():
            o = r['stability']
            print('  %-22s %s' % (name, ' '.join(
                '%5.0f/%4.0f' % (o['lhf'][i], o['shf'][i]) for i in sel)))


def samples_era5(a):
    import netCDF4 as nc
    E = a.era5
    ym = '2023_09'
    f8 = lambda x: np.asarray(np.ma.filled(x, np.nan), float)
    d = nc.Dataset('%s/2023/tatm_%s.nc' % (E, ym))
    la = f8(d['latitude'][:]); lo = f8(d['longitude'][:])
    lo2 = np.where(lo > 180, lo - 360, lo)
    step = int(round(a.res / 0.25))
    j = np.where((la >= -50) & (la <= 15))[0][::step]
    i = np.where((lo2 >= -110) & (lo2 <= -10))[0][::step]
    # contiguous reads, subsampled afterwards (strided reads of the
    # compressed files are slow and unreliable)
    js, is_ = slice(j[0], j[-1] + 1), slice(i[0], i[-1] + 1)
    sub = lambda x: x[..., ::step, ::step]
    lat2, lon2 = np.meshgrid(la[js][::step], lo2[is_][::step],
                             indexing='ij')
    lsm = nc.Dataset('%s/fixed/landmask.nc' % E)
    sea = sub(f8(lsm['lsm'][..., js, is_]).squeeze()) < 0.01
    t6 = nc.num2date(d['valid_time'][:], d['valid_time'].units)
    times = range(0, len(t6), a.tstep)
    lev = 0   # 1000 hPa, first level of the files
    d6 = {v: nc.Dataset('%s/2023/%s_%s.nc' % (E, v, ym))
          for v in ('tatm', 'qhum', 'geop', 'uwnd', 'vwnd')}
    get6 = lambda v, n, k: sub(f8(d6[v][n][k, lev, js, is_]))
    hps = nc.Dataset('%s/hourly/ps_%s.nc' % (E, ym))
    hssr = nc.Dataset('%s/hourly/ssr_%s.nc' % (E, ym))
    hstrd = nc.Dataset('%s/hourly/strd_%s.nc' % (E, ym))
    th = nc.num2date(hps['valid_time'][:], hps['valid_time'].units)
    hidx = {(x.day, x.hour): k for k, x in enumerate(th)}
    dsst = nc.Dataset('%s/SSTD/sst_%s.nc' % (E, ym))
    cols = {k: [] for k in ('z', 'u', 'v', 't', 'qv', 'p', 'ps', 'sst',
                            'zi', 'swnet', 'lwnet', 'lat', 'lon')}
    for k in times:
        when = t6[k]
        kh = hidx[(when.day, when.hour)]
        ps = sub(f8(hps['sp'][kh, js, is_]))
        sst = sub(f8(dsst['sst'][when.day - 1, js, is_]))
        t = get6('tatm', 't', k)
        q = get6('qhum', 'q', k)
        z = get6('geop', 'z', k) / G
        u = get6('uwnd', 'u', k); v = get6('vwnd', 'v', k)
        sw = sub(f8(hssr['ssr'][kh, js, is_])) / 3600.0
        lwd = sub(f8(hstrd['strd'][kh, js, is_])) / 3600.0
        ok = sea & np.isfinite(sst) & (sst > 272.0) & (z > 15.0) & \
            (ps > 1.0e5)
        lw = 0.97 * (5.67e-8 * sst ** 4 - lwd)
        for n, x in (('z', z), ('u', u), ('v', v), ('t', t),
                     ('qv', q / (1 - q)), ('p', np.full(z.shape, 1.0e5)),
                     ('ps', ps), ('sst', sst),
                     ('zi', np.full(z.shape, 600.0)), ('swnet', sw),
                     ('lwnet', lw), ('lat', lat2), ('lon', lon2)):
            cols[n].append(x[ok])
    s = {n: np.concatenate(x) for n, x in cols.items()}
    # ERA5 monthly mean fluxes at the same points
    em = lambda v, n, fac: sub(f8(nc.Dataset('%s/monthly/2023/%s_%s.nc'
                                             % (E, v, ym))[n][..., js, is_]
                                  ).squeeze()) * fac
    era = {'lhf': em('lflx', 'slhf', -1.0 / 86400.0),
           'shf': em('sflx', 'ishf', -1.0)}
    grid = dict(lat=lat2, lon=lon2, sea=sea)
    return s, era, grid


def mode_era5(a, exe, work):
    s, era, grid = samples_era5(a)
    print('%d samples (%.2f deg, every %d 6-hourly times), level height '
          '%.0f m (%.0f to %.0f)' % (len(s['z']), a.res, a.tstep,
                                    s['z'].mean(), s['z'].min(),
                                    s['z'].max()))
    res = {}
    for cfg in CONFIGS:
        res[cfg[0]] = run_scheme(exe, cfg, s, work)
    rname = 'COARE 3.5 (NOAA)' if a.coare else 'Zeng ECMWF'
    ref = run_coare(a.coare, s) if a.coare else res['Zeng ECMWF (rough 1)']
    # ERA5 monthly fluxes at the sampled points (each point weighted by
    # its number of samples)
    key = {(la, lo): (j, i) for j in range(grid['lat'].shape[0])
           for i in range(grid['lat'].shape[1])
           for la, lo in [(grid['lat'][j, i], grid['lon'][j, i])]}
    jj = np.array([key[(la, lo)][0] for la, lo in zip(s['lat'], s['lon'])])
    ii = np.array([key[(la, lo)][1] for la, lo in zip(s['lat'], s['lon'])])
    e_lhf = era['lhf'][jj, ii]; e_shf = era['shf'][jj, ii]
    for b, (la0, la1, lo0, lo1) in BOXES.items():
        w = (s['lat'] >= la0) & (s['lat'] < la1) & (s['lon'] >= lo0) & \
            (s['lon'] < lo1)
        print('\n%s: LH / SH (W/m2), tau (N/m2), ratio of LH to %s'
              % (b, rname))
        print('  %-24s %6.1f %6.1f' % ('ERA5 monthly', e_lhf[w].mean(),
                                       e_shf[w].mean()))
        print('  %-24s %6.1f %6.1f %6.3f' % (
            rname, ref['lhf'][w].mean(), ref['shf'][w].mean(),
            ref['tau'][w].mean()))
        for name, o in res.items():
            print('  %-24s %6.1f %6.1f %6.3f  %5.2f' % (
                name, o['lhf'][w].mean(), o['shf'][w].mean(),
                o['tau'][w].mean(), o['lhf'][w].mean() /
                ref['lhf'][w].mean()))


def main():
    ap = argparse.ArgumentParser(
        description='Evaluation of the RegCM ocean flux schemes')
    ap.add_argument('mode', choices=('curves', 'era5'))
    ap.add_argument('--coare', help='NOAA COARE-algorithm '
                    'Python/COARE3.5 directory; without it the reference '
                    'is the default Zeng configuration')
    ap.add_argument('--exe', help='scm_ocean program')
    ap.add_argument('--era5', default='/home/clima-archive5/ERA5')
    ap.add_argument('--res', type=float, default=1.0,
                    help='era5: grid step in degrees (default 1)')
    ap.add_argument('--tstep', type=int, default=1,
                    help='era5: use every n-th 6-hourly time')
    a = ap.parse_args()
    exe = find_exe(a.exe, 'scm_ocean')
    if exe is None:
        sys.exit('scm_ocean program not found: build it with "make '
                 'check" in Testing/scm or give --exe')
    work = tempfile.mkdtemp(prefix='scm_ocean_')
    if a.mode == 'curves':
        mode_curves(a, exe, work)
    else:
        mode_era5(a, exe, work)


if __name__ == '__main__':
    main()
