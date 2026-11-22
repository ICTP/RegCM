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
# Build a case of the coupled single column test (scm_coupled) from a
# DEPHY common format SCM driver (https://github.com/GdR-DEPHY/DEPHY-SCM)
# and, optionally, the reference LES profiles of the HIGH-TUNE ensemble
# (Couvreux et al., 2021, JAMES 13, e2020MS002217), averaged on the levels
# of the column.
#
# Usage: make_dephy_case.py --driver X_SCM_driver.nc --name bomex
#                           [--les DIR] [--levels FILE | --dz M]
#                           [--ztop M] [-o DIR]
#
# Writes in DIR (default: this directory) scm_<name>.txt (initial profiles
# and forcing), scm_<name>_target.txt (the geostrophic wind, read by
# scm_coupled with lgeos = .true.; T and q are the initial ones) and, with
# --les, scm_<name>_les.txt.
#
# Vertical grid: the levels of the first column of --levels (default
# scm_itcz.txt, the 32 MOLOCH levels of the South America test over the
# ocean) below ztop (default 10000 m), or a uniform grid of spacing dz.
#
# Forcing (constant in time in BOMEX and RICO, the first record is used):
#   tten  tnta_rad + tnta_adv (K/s), as present in the driver
#   qten  tnqv_adv (kg/kg/s, specific humidity)
#   omega -rho g wa (large scale subsidence)
#   geostrophic wind ug, vg (target file)
# Surface: surface_forcing_temp = "surface_flux": hfss, hfls (constant
# fluxes, run with iocnflux = 0) and ustar (ustarfix); "ts": the SST ts_forc
# (run with the ocean scheme, iocnflux 1 or 2).
#
# LES reference: for every hour of the LES files LES0.nc ... LESn.nc
# (hourly means of the ensemble members), on the column levels below the
# LES top, the ensemble mean, minimum and maximum of T (from theta at the
# pressure of the level), qv, ql (g/kg) and cloud fraction, each the mean
# over the LES levels inside the layer of the column level; per hour the
# liquid water path (g/m2) and the largest cloud fraction of the profile
# (a lower bound of the cloud cover) of each member.
#
# Needs numpy and netCDF4.

import argparse
import glob
import os
import sys

import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
RD, CP, G, P00 = 287.04, 1004.64, 9.80665, 1.0e5
EP1 = 0.608


def model_levels(path, ztop):
    """full and half level heights (top to bottom) of the first column"""
    with open(path) as f:
        lines = [x for x in f.read().splitlines() if not x.startswith('#')]
    nk = int(lines[0].split()[1])
    zfull = np.array([float(lines[2 + k].split()[0]) for k in range(nk)])
    zhalf = np.array([float(lines[2 + nk + k].split()[0])
                      for k in range(nk + 1)])
    k0 = int(np.argmax(zfull < ztop))
    return zfull[k0:], zhalf[k0:]


def uniform_levels(dz, ztop):
    zhalf = np.arange(0.0, ztop + 0.5 * dz, dz)[::-1]
    zfull = 0.5 * (zhalf[:-1] + zhalf[1:])
    return zfull, zhalf


def interp(z, zsrc, f):
    return np.interp(z, zsrc, f)


def read_driver(path):
    d = nc.Dataset(path)
    c = {'case': d.case, 'title': d.title,
         'start': d.start_date.replace('-', '').replace(' ', '')[:10],
         'tforc': d.surface_forcing_temp,
         'z': np.array(d['zh'][0], dtype=float),
         'p': np.array(d['pa'][0], dtype=float),
         'ps': float(d['ps'][0]),
         'lat': float(d['lat'][0]), 'lon': float(d['lon'][0])}
    for v in ('ta', 'qv', 'ua', 'va'):
        c[v] = np.array(d[v][0], dtype=float)
    zf = np.array(d['zh_forc'][0], dtype=float)
    for v in ('ug', 'vg', 'wa', 'tnqv_adv'):
        c[v] = interp(c['z'], zf, np.array(d[v][0], dtype=float))
    c['tten'] = np.zeros_like(c['z'])
    for v in ('tnta_rad', 'tnta_adv'):
        if v in d.variables:
            c['tten'] += interp(c['z'], zf, np.array(d[v][0], dtype=float))
    c['shf'] = float(d['hfss'][0]) if 'hfss' in d.variables else 0.0
    c['lhf'] = float(d['hfls'][0]) if 'hfls' in d.variables else 0.0
    c['ustar'] = float(d['ustar'][0]) if 'ustar' in d.variables else 0.0
    if 'ts_forc' in d.variables:
        c['sst'] = float(d['ts_forc'][0])
    elif 'tskin' in d.variables:
        c['sst'] = float(d['tskin'][0])
    else:
        c['sst'] = float(d['ta'][0, 0])
    return c


def column(c, zfull, zhalf):
    """profile lines of the column, the pressure of the levels"""
    lnp = np.log(c['p'])
    pfull = np.exp(interp(zfull, c['z'], lnp))
    phalf = np.exp(interp(zhalf, c['z'], lnp))
    phalf[-1] = c['ps']
    t = interp(zfull, c['z'], c['ta'])
    q = np.maximum(interp(zfull, c['z'], c['qv']), 1.0e-6)
    rho = pfull / (RD * t * (1.0 + EP1 * q))
    omega = -rho * G * interp(zfull, c['z'], c['wa'])
    prof = {'t': t, 'q': q,
            'u': interp(zfull, c['z'], c['ua']),
            'v': interp(zfull, c['z'], c['va']),
            'ug': interp(zfull, c['z'], c['ug']),
            'vg': interp(zfull, c['z'], c['vg']),
            'omega': omega,
            'tten': interp(zfull, c['z'], c['tten']),
            'qten': interp(zfull, c['z'], c['tnqv_adv'])}
    return pfull, phalf, prof


def write_case(path, name, c, zfull, zhalf, pfull, phalf, prof, uwind,
               vwind, comment, nprof):
    nk = len(zfull)
    with open(path, 'w') as f:
        f.write('# RegCM coupled single column test case, %s: %s, '
                'built by make_dephy_case.py\n' % (c['case'], comment))
        f.write('# per column: name lat lon land date ps(Pa) shf lhf '
                '(W/m2) zpbl cape conv_rain (unused) orography(m) '
                'sst(K)\n')
        f.write('# then kz lines top to bottom: z(m above ground) p(Pa) '
                'T(K) q(kg/kg) u v(m/s) omega(Pa/s) tten(K/s) '
                'qten(kg/kg/s)\n')
        f.write('# then kz+1 half levels: z(m above ground) p(Pa)\n')
        f.write('%d %d\n' % (nprof, nk))
        for _ in range(nprof):
            f.write('%s %.2f %.2f 0 %s %.2f %.3f %.3f 500.0 0.0 0.000 0.0 '
                    '%.2f\n' % (name.upper(), c['lat'], c['lon'],
                                c['start'], c['ps'], c['shf'], c['lhf'],
                                c['sst']))
            for k in range(nk):
                f.write('%10.2f %10.2f %8.3f %12.5e %8.3f %8.3f %12.5e '
                        '%12.5e %12.5e\n' % (
                            zfull[k], pfull[k], prof['t'][k], prof['q'][k],
                            uwind[k], vwind[k], prof['omega'][k],
                            prof['tten'][k], prof['qten'][k]))
            for k in range(nk + 1):
                f.write('%10.2f %10.2f\n' % (zhalf[k], phalf[k]))


def les_reference(path, name, ldir, zfull, zhalf, pfull):
    files = sorted(glob.glob(os.path.join(ldir, 'LES[0-9]*.nc')))
    if not files:
        sys.exit('no LES files LES*.nc in ' + ldir)
    members = []
    for fn in files:
        d = nc.Dataset(fn)
        z = np.ravel(np.array(d['zf'][:], dtype=float))
        tb = np.array(d['time_bnds'][:], dtype=float)
        # full hourly means only (the files end with a short record)
        full = (tb[:, 1] - tb[:, 0]) > 1800.0
        hours = np.rint((tb[full, 1] - tb[0, 0]) / 3600.0).astype(int)
        v = {}
        for n in ('theta', 'qv', 'ql', 'rneb', 'rho'):
            v[n] = np.array(d[n][:], dtype=float).reshape(
                len(tb), -1)[full]
        members.append((z, hours, v))
    ztopl = min(m[0][-1] for m in members)
    keep = zfull < ztopl
    zf, zh, pf = zfull[keep], zhalf[:len(zfull) + 1], pfull[keep]
    zh = zh[-(len(zf) + 1):]
    hours = members[0][1]
    nh = min(len(m[1]) for m in members)
    exner = (pf / P00) ** (RD / CP)

    def layer(z, f):
        out = np.empty(len(zf))
        for k in range(len(zf)):
            sel = (z >= zh[k + 1]) & (z < zh[k])
            out[k] = f[sel].mean() if sel.any() else np.interp(zf[k], z, f)
        return out

    with open(path, 'w') as f:
        f.write('# LES reference for scm_%s.txt: HIGH-TUNE Meso-NH '
                'ensemble (%d members, %s), built by make_dephy_case.py\n'
                % (name, len(files), ldir))
        f.write('# per hour: HOUR h (end of the hourly mean) then lwp '
                '(g/m2) and largest cloud fraction of the profile, mean '
                'min max of the members\n')
        f.write('# then one line per column level below the LES top, top '
                'to bottom: z(m) T(K) qv(g/kg) ql(g/kg) cldfrac, each '
                'mean min max\n')
        f.write('%d %d\n' % (nh, len(zf)))
        for i in range(nh):
            prof, lwp, cc = [], [], []
            for z, _, v in members:
                dz = np.gradient(z)
                lwp.append((v['rho'][i] * v['ql'][i] * dz).sum() * 1.0e3)
                cc.append(v['rneb'][i].max())
                prof.append(np.array([
                    layer(z, v['theta'][i]) * exner,
                    layer(z, v['qv'][i]) * 1.0e3,
                    layer(z, v['ql'][i]) * 1.0e3,
                    layer(z, v['rneb'][i])]))
            prof = np.array(prof)
            f.write('HOUR %d %.3f %.3f %.3f %.4f %.4f %.4f\n' % (
                hours[i], np.mean(lwp), np.min(lwp), np.max(lwp),
                np.mean(cc), np.min(cc), np.max(cc)))
            for k in range(len(zf)):
                f.write('%10.2f' % zf[k])
                for n in range(4):
                    x = prof[:, n, k]
                    f.write(' %.5g %.5g %.5g' % (x.mean(), x.min(), x.max()))
                f.write('\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--driver', required=True, help='DEPHY SCM driver')
    ap.add_argument('--name', required=True, help='case name (bomex)')
    ap.add_argument('--les', help='directory with the LES*.nc files')
    ap.add_argument('--levels', default=os.path.join(HERE, 'scm_itcz.txt'),
                    help='profile file whose first column gives the levels')
    ap.add_argument('--dz', type=float, help='uniform level spacing (m)')
    ap.add_argument('--ztop', type=float, default=10000.0,
                    help='top of the column (m)')
    ap.add_argument('-o', '--outdir', default=HERE)
    a = ap.parse_args()
    if a.dz:
        zfull, zhalf = uniform_levels(a.dz, a.ztop)
        comment = 'uniform levels every %g m to %g m' % (a.dz, a.ztop)
    else:
        zfull, zhalf = model_levels(a.levels, a.ztop)
        comment = 'levels of %s below %g m' % (os.path.basename(a.levels),
                                               a.ztop)
    c = read_driver(a.driver)
    pfull, phalf, prof = column(c, zfull, zhalf)
    base = os.path.join(a.outdir, 'scm_%s' % a.name)
    write_case(base + '.txt', a.name, c, zfull, zhalf, pfull, phalf, prof,
               prof['u'], prof['v'], comment, 1)
    write_case(base + '_target.txt', a.name, c, zfull, zhalf, pfull,
               phalf, prof, prof['ug'], prof['vg'],
               'targets (geostrophic wind)', 4)
    print('%s: %d levels, surface %s, written %s.txt, %s_target.txt' % (
        c['case'], len(zfull), c['tforc'], base, base))
    if c['tforc'] == 'surface_flux':
        print('  run with iocnflux = 0, ustarfix = %g (shf %.2f, lhf %.2f '
              'W/m2)' % (c['ustar'], c['shf'], c['lhf']))
    else:
        print('  run with the ocean scheme, SST %.2f K' % c['sst'])
    if a.les:
        les_reference(base + '_les.txt', a.name, a.les, zfull, zhalf, pfull)
        print('  LES reference %s_les.txt' % base)


if __name__ == '__main__':
    main()
