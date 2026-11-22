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
# Shallow cumulus cases of the coupled single column test (scm_coupled)
# against LES: BOMEX (Siebesma et al., 2003) and RICO (vanZanten et al.,
# 2011), DEPHY forcing (make_dephy_case.py), HIGH-TUNE Meso-NH LES
# ensemble as reference (scm_<case>_les.txt).
#
# Runs the case with its settings (below) plus the ones given on the
# command line and prints, hour by hour, the PBL height, the cumulus cloud
# type, the liquid water path and the largest cloud fraction against the
# LES ensemble (mean and range), then the profiles of T, q, cloud
# condensate and cloud fraction averaged over the evaluation window
# (BOMEX hours 3-6, RICO 20-24, as in the intercomparisons).
#
# Examples:
#   scm_les.py bomex
#   scm_les.py rico --set lmfscv=.false.
#   scm_les.py bomex --set icumulus=6           Kain-Fritsch
#   scm_les.py rico --sweep entshalp=1.0,2.0,3.0
#   scm_les.py bomex --save bomex.out           keep the output
#
# The SCM condensate is the grid mean cloud water + ice of the
# microphysics (qc) and the one seen by the radiation, convective cloud
# included (qrad = cldfrac * in-cloud water / density). The SCM cloud
# fraction is the one seen by the radiation. The program is looked for as
# in scm_coupled.py. Only the standard library is needed.

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scm_coupled as sc  # noqa: E402

RD = 287.04

# Case settings: surface and forcing as prescribed by the case, no
# relaxation, no interactive radiation (the forcing has the radiative
# cooling), winds turned by the Coriolis force around the geostrophic
# wind of the target file.
CASES = {
    'bomex': {'settings': {'ndays': '0.25', 'irad': '0',
                           'taurelax': '1.0e12', 'lgeos': '.true.',
                           'iocnflux': '0', 'ustarfix': '0.28',
                           'nout': '12'},
              'window': (3, 6)},
    'rico': {'settings': {'ndays': '1.0', 'irad': '0',
                          'taurelax': '1.0e12', 'lgeos': '.true.',
                          'iocnflux': '1', 'ocnsw': '0.0', 'ocnlw': '0.0',
                          'nout': '12'},
             'window': (20, 24)},
}


def read_les(path):
    """hour -> {'lwp': (mean, min, max), 'cf': (...), 'levels': [...]}"""
    with open(path) as f:
        lines = [x for x in f.read().splitlines() if not x.startswith('#')]
    nh, nk = (int(x) for x in lines[0].split()[:2])
    les, i = {}, 1
    for _ in range(nh):
        w = lines[i].split()
        h = int(w[1])
        r = {'lwp': tuple(float(x) for x in w[2:5]),
             'cf': tuple(float(x) for x in w[5:8]), 'levels': []}
        for k in range(nk):
            v = [float(x) for x in lines[i + 1 + k].split()]
            r['levels'].append({'z': v[0], 'T': v[1:4], 'q': v[4:7],
                                'ql': v[7:10], 'cf': v[10:13]})
        les[h] = r
        i += nk + 1
    return les


def scm_hourly(steps):
    """hourly values of one column from the output steps"""
    out = {}
    for s in steps:
        h = int(round(s['time']))
        if abs(s['time'] - h) > 1.0e-3:
            continue
        lv = s['levels']
        nk = len(lv['z'])
        lwp = lwr = 0.0
        qrad = []
        for k in range(nk):
            # layer depth from the level heights (top to bottom)
            zt = lv['z'][k - 1] if k > 0 else lv['z'][0]
            zb = lv['z'][k + 1] if k < nk - 1 else 0.0
            dz = 0.5 * (zt - zb) if k > 0 else lv['z'][0] - lv['z'][1]
            rho = lv['p'][k] * 100.0 / (RD * lv['T'][k])
            lwp += lv['qcond'][k] * rho * dz
            lwr += lv['cldf'][k] * lv['lwcr'][k] * dz
            qrad.append(lv['cldf'][k] * lv['lwcr'][k] / rho)
        out[h] = {'zpbl': s['zpbl'], 'ktype': s['ktype'],
                  'kctop': s['kctop'], 'kcbot': s['kcbot'],
                  'rain': s['rainc'] + s['rainl'],
                  'lwp': lwp, 'lwr': lwr, 'cf': max(lv['cldf']),
                  'levels': lv, 'qrad': qrad}
    return out


def report(name, scm, les, window):
    print('%s: hourly, SCM against the LES ensemble mean [min, max]' %
          name.upper())
    print('  h  zpbl(m) ktype  LWP qc  LWP rad   LWP LES              '
          'cfmax   cfmax LES           rain(mm/d)')
    for h in sorted(scm):
        if h == 0:
            continue
        s = scm[h]
        line = '%3d %7.0f %4d %8.1f %8.1f' % (h, s['zpbl'], s['ktype'],
                                             s['lwp'], s['lwr'])
        if h in les:
            r = les[h]
            line += '  %6.1f [%5.1f,%5.1f]  %6.3f  %6.3f [%5.3f,%5.3f]' % (
                r['lwp'] + (s['cf'],) + r['cf'])
        else:
            line += ' ' * 23 + '  %6.3f' % s['cf'] + ' ' * 22
        line += '  %8.3f' % s['rain']
        print(line)
    h0, h1 = window
    hrs = [h for h in range(h0 + 1, h1 + 1) if h in scm and h in les]
    if not hrs:
        return
    print('\n%s: profiles, mean of hours %d-%d (LES layer means on the '
          'SCM levels)' % (name.upper(), h0, h1))
    print('     z(m)   dT(K)   dq(g/kg)   qc  qrad  ql LES [min,max] '
          '(mg/kg)     cf  cf LES [min,max]')
    nl = len(les[hrs[0]]['levels'])
    nk = len(scm[hrs[0]]['levels']['z'])
    for k in range(nl):
        ks = nk - nl + k
        m = len(hrs)

        def sm(f):
            return sum(scm[h]['levels'][f][ks] for h in hrs) / m

        def lm(f, i):
            return sum(les[h]['levels'][k][f][i] for h in hrs) / m

        qrad = sum(scm[h]['qrad'][ks] for h in hrs) / m
        print('%9.0f %7.2f %8.2f %6.1f %5.1f %6.1f [%5.1f,%5.1f]   '
              '%7.3f %6.3f [%5.3f,%5.3f]' % (
                  les[hrs[0]]['levels'][k]['z'], sm('T') - lm('T', 0),
                  sm('q') - lm('q', 0), sm('qcond') * 1.0e3,
                  qrad * 1.0e3, lm('ql', 0) * 1.0e3, lm('ql', 1) * 1.0e3,
                  lm('ql', 2) * 1.0e3, sm('cldf'), lm('cf', 0),
                  lm('cf', 1), lm('cf', 2)))
    print('(dT, dq: SCM minus LES mean; qc, qrad, ql in mg/kg)')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.
                                 RawDescriptionHelpFormatter)
    ap.add_argument('case', choices=sorted(CASES))
    ap.add_argument('--exe', help='scm_coupled program')
    ap.add_argument('--set', action='append', default=[],
                    metavar='NAME=VALUE', help='namelist setting')
    ap.add_argument('--sweep', metavar='NAME=V1,V2,...',
                    help='run once per value')
    ap.add_argument('--save', metavar='FILE', help='keep the output file')
    a = ap.parse_args()
    exe = sc.find_exe(a.exe)
    case = CASES[a.case]
    profiles = os.path.join(HERE, 'scm_%s.txt' % a.case)
    les = read_les(os.path.join(HERE, 'scm_%s_les.txt' % a.case))
    settings = dict(case['settings'])
    for x in a.set:
        k, v = x.split('=', 1)
        settings[k.strip()] = v
    runs = [(None, settings)]
    if a.sweep:
        k, vals = a.sweep.split('=', 1)
        runs = []
        for v in vals.split(','):
            s = dict(settings)
            s[k.strip()] = v
            runs.append(('%s = %s' % (k, v), s))
    for title, s in runs:
        if title:
            print('\n==== %s' % title)
        (header, cols), out, work = sc.run(exe, s, profiles)
        if a.save:
            sc.shutil.copy(out, a.save)
        sc.shutil.rmtree(work, ignore_errors=True)
        for name, steps in cols.items():
            report(name, scm_hourly(steps), les, case['window'])


if __name__ == '__main__':
    main()
