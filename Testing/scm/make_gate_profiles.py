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
# Build single column test profiles from the GATE Phase III soundings
# (file testgate: 161 three-hourly soundings on 41 pressure levels, from
# 1012 to 70 hPa, with the observed apparent heat source Q1, moisture sink
# Q2, radiative heating QR and the horizontal advection of temperature
# and humidity).
#
# Writes a profile file in the format of scm_profiles.txt (read by the
# scm_* drivers) and a file with the observed Q1, Q2, QR on the same
# levels (read by scm_tiedtke.py --gate).
#
# The 40 levels above the surface are the model levels; the half levels
# are at the middle in log pressure, the lowest one at the surface. The
# input tendencies are the large scale forcing:
#   adv     horizontal advection (ADVT, ADVQ of the file) plus the
#           vertical advection by the observed omega,
#           -omega (dT/dp - R T/(cp p)) and -omega dq/dp
#   advrad  as adv, plus the radiative heating QR (default)
#   none    zero, as in RegCM where cumulus is the first physics call
# The file has no surface fluxes and no boundary layer height: the fluxes
# are constant (--shf, --lhf), the PBL height is a bulk Richardson number
# diagnosis (critical value 0.25) from the lowest level. The observed
# precipitation listed per column is the moisture budget value,
# P = E + (cp/L) int Q2 dp/g, with the assumed evaporation E.
#
# Examples:
#   make_gate_profiles.py                      all 161 soundings
#   make_gate_profiles.py --forcing none -o scm_gate_none.txt
#   make_gate_profiles.py --select 7,33,61,89,101,113,135,150 \
#                         -o scm_gate_check.txt

import argparse
import datetime as dt
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
G = 9.80616
RD = 287.0586
CP = 1004.64
LV = 2.501e6
EPS = 0.6077762876
KAP = RD / CP
# GATE A/B scale array, Phase III from 1974-08-30 00 UTC, 3-hourly
LAT, LON = 8.5, -23.5
START = dt.datetime(1974, 8, 30, 0)


def read_gate(path):
    snd, cur = [], None
    with open(path) as f:
        for line in f:
            if line.startswith('% Soundg'):
                cur = []
                snd.append(cur)
            elif not line.startswith('%') and line.strip():
                cur.append([float(x) for x in line.split()])
    a = np.array(snd)            # (nsound, 41, 12), top to bottom
    names = ['z', 'p', 'T', 'q', 'u', 'v', 'w', 'Q1', 'Q2', 'QR',
             'ADVT', 'ADVQ']
    return {n: a[:, :, i] for i, n in enumerate(names)}


def column(g, n, forcing, shf, lhf):
    z, p = g['z'][n], g['p'][n] * 100.0
    T = g['T'][n] + 273.15
    q = g['q'][n] * 1.0e-3
    w = g['w'][n]
    # vertical advection on the 41 levels (surface included)
    dtdp = np.gradient(T, p)
    dqdp = np.gradient(q, p)
    vadt = -w * (dtdp - KAP * T / p)
    vadq = -w * dqdp
    if forcing == 'none':
        tten = np.zeros_like(T)
        qten = np.zeros_like(q)
    else:
        tten = g['ADVT'][n] / 86400.0 + vadt
        qten = g['ADVQ'][n] + vadq
        if forcing == 'advrad':
            tten = tten + g['QR'][n] / 86400.0
    k = slice(0, -1)             # model levels: all but the surface
    zf, pf = z[k], p[k]
    # half levels: log pressure middle, top extrapolated, bottom surface
    lp = np.log(pf)
    ph = np.empty(len(pf) + 1)
    ph[1:-1] = np.exp(0.5 * (lp[:-1] + lp[1:]))
    ph[0] = np.exp(lp[0] - 0.5 * (lp[1] - lp[0]))
    ph[-1] = p[-1]
    zh = np.empty(len(pf) + 1)
    zh[1:-1] = np.interp(np.log(ph[1:-1]), lp, zf)
    zh[0] = zf[0] + 0.5 * (zf[0] - zf[1])
    zh[-1] = 0.0
    # PBL height: bulk Richardson number from the lowest model level
    qs = q[k]
    thv = T[k] * (1 + EPS * qs) * (1e5 / pf) ** KAP
    u, v = g['u'][n][k], g['v'][n][k]
    zpbl = 100.0
    rib = np.zeros(len(zf))
    for kk in range(len(zf) - 2, -1, -1):
        du2 = max((u[kk] - u[-1]) ** 2 + (v[kk] - v[-1]) ** 2, 1.0)
        rib[kk] = G * (zf[kk] - zf[-1]) * (thv[kk] - thv[-1]) / \
            (thv[-1] * du2)
        if rib[kk] > 0.25:
            zpbl = max(zf[kk + 1] + (zf[kk] - zf[kk + 1]) *
                       (0.25 - rib[kk + 1]) / (rib[kk] - rib[kk + 1]),
                       100.0)
            break
    # moisture budget precipitation [mm/day]
    q2 = g['Q2'][n][k]
    dp = ph[1:] - ph[:-1]
    prec = (lhf / LV + CP / LV * np.sum(q2 / 86400.0 * dp) / G) * 86400.0
    when = START + dt.timedelta(hours=3 * n)
    return dict(name='GATE_%03d' % (n + 1), when=when.strftime('%Y%m%d%H'),
                ps=p[-1], zpbl=zpbl, prec=prec, zf=zf, pf=pf, T=T[k],
                q=qs, u=u, v=v, w=w[k], tten=tten[k], qten=qten[k],
                zh=zh, ph=ph, Q1=g['Q1'][n][k], Q2=q2, QR=g['QR'][n][k],
                shf=shf, lhf=lhf)


def main():
    ap = argparse.ArgumentParser(
        description='Single column test profiles from the GATE soundings')
    ap.add_argument('--gate', default=os.path.join(HERE, 'testgate'))
    ap.add_argument('--forcing', choices=('advrad', 'adv', 'none'),
                    default='advrad')
    ap.add_argument('--shf', type=float, default=10.0,
                    help='sensible heat flux, W/m2 (default 10)')
    ap.add_argument('--lhf', type=float, default=125.0,
                    help='latent heat flux, W/m2 (default 125)')
    ap.add_argument('--select', help='1-based sounding numbers, a,b,c')
    ap.add_argument('-o', '--output', default='scm_gate.txt')
    ap.add_argument('--obs', help='observed Q1, Q2, QR file '
                    '(default: output name with _obs)')
    a = ap.parse_args()
    g = read_gate(a.gate)
    ns = g['p'].shape[0]
    sel = ([int(x) - 1 for x in a.select.split(',')] if a.select
           else range(ns))
    cols = [column(g, n, a.forcing, a.shf, a.lhf) for n in sel]
    kz = len(cols[0]['pf'])
    with open(a.output, 'w') as f:
        f.write('# RegCM single column test profiles from the GATE Phase '
                'III soundings, built by make_gate_profiles.py\n')
        f.write('# forcing %s, shf %.1f lhf %.1f W/m2 (assumed); the '
                'conv_rain column is the budget precipitation\n' %
                (a.forcing, a.shf, a.lhf))
        f.write('# per column: name lat lon land date ps(Pa) shf lhf '
                '(W/m2, upward) zpbl(m) cape(J/kg) conv_rain(mm/day) '
                'orography(m)\n')
        f.write('# then kz lines top to bottom: z(m above ground) p(Pa) '
                'T(K) q(kg/kg) u v(m/s) omega(Pa/s) tten(K/s) '
                'qten(kg/kg/s)\n')
        f.write('# then kz+1 half levels: z(m above ground) p(Pa)\n')
        f.write('%d %d\n' % (len(cols), kz))
        for c in cols:
            f.write('%s %.2f %.2f %d %s %.2f %.3f %.3f %.1f %.1f %.3f %.1f\n'
                    % (c['name'], LAT, LON, 0, c['when'], c['ps'], c['shf'],
                       c['lhf'], c['zpbl'], 0.0, c['prec'], 0.0))
            for k in range(kz):
                f.write('%10.2f %10.2f %8.3f %12.5e %8.3f %8.3f %12.5e '
                        '%12.5e %12.5e\n' % (
                            c['zf'][k], c['pf'][k], c['T'][k], c['q'][k],
                            c['u'][k], c['v'][k], c['w'][k], c['tten'][k],
                            c['qten'][k]))
            for k in range(kz + 1):
                f.write('%10.2f %10.2f\n' % (c['zh'][k], c['ph'][k]))
    obs = a.obs or os.path.splitext(a.output)[0] + '_obs.txt'
    with open(obs, 'w') as f:
        f.write('# GATE observed apparent heat source Q1, moisture sink Q2 '
                'and radiative heating QR (K/day), on the levels of %s\n'
                % os.path.basename(a.output))
        f.write('# per column: name date budget_precipitation(mm/day), '
                'then kz lines top to bottom: p(hPa) Q1 Q2 QR\n')
        f.write('%d %d\n' % (len(cols), kz))
        for c in cols:
            f.write('%s %s %.3f\n' % (c['name'], c['when'], c['prec']))
            for k in range(kz):
                f.write('%8.1f %8.3f %8.3f %8.3f\n' % (
                    c['pf'][k] / 100.0, c['Q1'][k], c['Q2'][k], c['QR'][k]))
    print('%d columns, %d levels, forcing %s: %s, %s' % (
        len(cols), kz, a.forcing, a.output, obs))


if __name__ == '__main__':
    main()
