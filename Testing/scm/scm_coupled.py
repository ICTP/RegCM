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
# Python wrapper of the coupled Holtslag + Tiedtke single column test
# (scm_coupled).
#
# Runs the Fortran program on the Amazon case with namelist settings given
# on the command line and prints per column the mean daily rain (from the
# second day on), its spread from day to day, the frequency of the cloud
# types, the largest mean PBL height and the mean diurnal cycle of the
# rain against ERA5 (scm_amazon_era5.txt); sweeps a parameter over a list
# of values, or compares an output with a reference.
#
# Examples:
#   scm_coupled.py                               30 days, RegCM defaults
#   scm_coupled.py --set ndays=10
#   scm_coupled.py --set lmfcapebl=.false.
#   scm_coupled.py --sweep zhnew_fac=0.0,0.25,0.5
#   scm_coupled.py --list                        show the parameters
#   scm_coupled.py --check test.out reference/scm_coupled_amazon.out
#
# Convection switches on and off along the integration, so differences
# of rounding (another compiler or other optimisation flags) grow into
# different rain events: the comparison with the reference holds for the
# build it was made with. The mean daily rain over 30 days is stable to
# about 0.1 mm/day.
#
# The program is looked for in --exe, $SCM_COUPLED, the current directory
# and the directory of this script. Only the standard library is needed.

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

# namelist group of every parameter the program reads
GROUPS = {
    'scmparam': ['dtstep', 'ndays', 'nout', 'dxkm', 'taurelax', 'zrelax',
                 'tgtfile', 'tauwtg', 'zwtg', 'tauwind', 'windfac',
                 'iocnflux', 'izengrough', 'ocnsw', 'ocnlw', 'lnoto'],
    'tiedtkeparam': [
        'entrdd', 'entrpen_lnd', 'entrpen_ocn', 'detrpen_lnd',
        'detrpen_ocn', 'entshalp', 'rcuc_lnd', 'rcuc_ocn', 'rcpec_lnd',
        'rcpec_ocn', 'rhebc_lnd', 'rhebc_ocn', 'rprc_ocn', 'rprc_lnd',
        'revap_lnd', 'revap_ocn', 'lmfpen', 'lmfmid', 'lmfdd', 'lepcld',
        'lmfdudv', 'lmfscv', 'lmfuvdis', 'lmftrac', 'lmfsmooth', 'lmfwstar',
        'lmfcapebl', 'lmfcapeforc', 'rtaucape'],
    'cldparam': ['kfac_deep', 'kfac_shal', 'k2_const', 'iconvcf',
                 'rh0oce', 'rh0land', 'tc0', 'rhmin', 'rhmax',
                 'larcticcorr'],
    'microparam': ['nssopt', 'iautoconv', 'vfqs', 'auto_rate_khair',
                   'auto_rate_kessl', 'auto_rate_klepi', 'rkconv',
                   'skconv', 'rcovpmin', 'rpecons', 'rcldiff'],
    'holtslagparam': ['ricr_ocn', 'ricr_lnd', 'zhnew_fac', 'ifaholtth10',
                      'ifaholt', 'holtth10iter', 'lmoistn2', 'maxprandtl'],
}
GROUP_OF = {v: g for g, vs in GROUPS.items() for v in vs}

STEP_FIELDS = ['time', 'hour', 'shf', 'lhf', 'zpbl', 'kpbl', 'hqfac',
               'ktype', 'kctop', 'kcbot', 'rainc', 'rainl', 'rainc_acc',
               'rainl_acc', 'wspd']
INT_FIELDS = ('kpbl', 'ktype', 'kctop', 'kcbot')
LEVEL_FIELDS = ['z', 'p', 'T', 'q', 'mfu', 'dTcu', 'dTpbl', 'dqcu',
                'dqpbl', 'omega', 'dTmic', 'dqmic', 'qcond', 'cldf']


def find_exe(exe):
    for c in (exe, os.environ.get('SCM_COUPLED'),
              os.path.join(os.getcwd(), 'scm_coupled'),
              os.path.join(HERE, 'scm_coupled')):
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.abspath(c)
    sys.exit('scm_coupled program not found: build it with "make check" '
             'in Testing/scm or give --exe')


def fortran_value(v):
    s = v.strip()
    if s.lower() in ('.true.', '.false.', 't', 'f', 'true', 'false'):
        return '.true.' if s.lower() in ('.true.', 't', 'true') else '.false.'
    try:
        float(s.replace('d', 'e').replace('D', 'e'))
        return s
    except ValueError:
        return "'%s'" % s


def write_namelist(path, settings, infile, outfile):
    groups = {g: {} for g in GROUPS}
    groups['scmparam']['infile'] = "'%s'" % infile
    groups['scmparam']['outfile'] = "'%s'" % outfile
    groups['scmparam']['tgtfile'] = "'%s'" % (
        os.path.splitext(infile)[0] + '_target.txt')
    groups['scmparam']['ndays'] = '30.0'
    for k, v in settings.items():
        if k not in GROUP_OF:
            sys.exit('unknown parameter %s (see --list)' % k)
        groups[GROUP_OF[k]][k] = fortran_value(v)
    with open(path, 'w') as f:
        for g, vals in groups.items():
            f.write(' &%s\n' % g)
            for k, v in vals.items():
                f.write('  %s = %s,\n' % (k, v))
            f.write(' /\n')


def parse_output(path):
    """header text, dict column name -> list of output steps"""
    header, cols = [], {}
    with open(path) as f:
        lines = f.read().splitlines()
    i = 0
    while not lines[i].startswith('NCOL'):
        header.append(lines[i]); i += 1
    nk = int(lines[i].split()[2]); i += 1
    while i < len(lines):
        w = lines[i].split(); i += 1
        s = {}
        for n, x in zip(STEP_FIELDS, w[2:]):
            s[n] = int(x) if n in INT_FIELDS else float(x)
        s['levels'] = {n: [] for n in LEVEL_FIELDS}
        for _ in range(nk):
            for n, x in zip(LEVEL_FIELDS, lines[i].split()):
                s['levels'][n].append(float(x))
            i += 1
        cols.setdefault(w[1], []).append(s)
    return header, cols


def read_era5(path):
    """column name -> {'pr': rain by local hour, 'u10': 10 m wind}"""
    ref = {}
    if not os.path.isfile(path):
        return ref
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            w = line.split()
            r = ref.setdefault(w[0], {'pr': [0.0] * 24, 'u10': None})
            r['pr'][int(w[1])] = float(w[2])
            if len(w) > 5:
                if r['u10'] is None:
                    r['u10'] = [0.0] * 24
                r['u10'][int(w[1])] = float(w[5])
    return ref


def read_profiles(path):
    """column name -> daily mean fluxes, SST, input omega profile"""
    cols = {}
    with open(path) as f:
        lines = [x for x in f.read().splitlines() if not x.startswith('#')]
    ncol, nk = (int(x) for x in lines[0].split()[:2])
    i = 1
    for _ in range(ncol):
        w = lines[i].split()
        c = {'shf': float(w[6]), 'lhf': float(w[7]),
             'sst': float(w[12]) if len(w) > 12 else 0.0,
             'p': [], 'omega': []}
        for k in range(nk):
            v = lines[i + 1 + k].split()
            c['p'].append(float(v[1]) / 100.0)
            c['omega'].append(float(v[6]))
        cols[w[0]] = c
        i += 2 * nk + 2
    return cols


def run(exe, settings, profiles, keep=None):
    work = os.path.abspath(keep or tempfile.mkdtemp(prefix='scm_coupled_'))
    os.makedirs(work, exist_ok=True)
    out = os.path.join(work, 'scm_coupled.out')
    nml = os.path.join(work, 'scm_coupled.in')
    write_namelist(nml, settings, os.path.abspath(profiles), out)
    r = subprocess.run([exe, nml], cwd=work, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.isfile(out):
        sys.exit('scm_coupled failed:\n' + r.stdout + r.stderr)
    return parse_output(out), out, work


def mean(x):
    return sum(x) / len(x) if x else float('nan')


def summary(cols, era5, prof=None, title=None):
    if title:
        print(title)
    for name, steps in cols.items():
        sel = [s for s in steps if s['time'] > 24.0]
        if not sel:
            sel = steps
        nday = int(max(s['time'] for s in sel) // 24.0 + 0.999)
        daily = []
        for d in range(1, nday + 1):
            x = [s['rainc'] + s['rainl'] for s in sel
                 if 24.0 * d < s['time'] <= 24.0 * (d + 1)]
            if x:
                daily.append(mean(x))
        rain = mean(daily)
        sd = math.sqrt(mean([(x - rain) ** 2 for x in daily]))
        freq = [100.0 * mean([1.0 if s['ktype'] == k else 0.0
                              for s in sel]) for k in (1, 2, 3)]
        diur = [[s['rainc'] + s['rainl'] for s in sel
                 if int(s['hour']) % 24 == h] for h in range(24)]
        zpbl = max(mean([s['zpbl'] for s in sel if int(s['hour']) == h])
                   for h in range(24) if diur[h])
        print('%-10s rain %5.2f mm/day (+- %4.2f, large scale %4.2f), '
              'ERA5 %s; deep %2.0f%% shallow %2.0f%% mid %2.0f%%; max '
              'mean zpbl %4.0f m' % (
                  name, rain, sd, mean([s['rainl'] for s in sel]),
                  '%4.2f' % mean(era5[name]['pr']) if name in era5 else '-',
                  freq[0], freq[1], freq[2], zpbl))
        print('  local hour ' + ' '.join('%5d' % h for h in range(0, 24, 2)))
        print('  SCM        ' + ' '.join(
            '%5.1f' % mean(diur[h] + diur[h + 1]) for h in range(0, 24, 2)))
        if name in era5:
            print('  ERA5       ' + ' '.join(
                '%5.1f' % mean(era5[name]['pr'][h:h + 2])
                for h in range(0, 24, 2)))
        pc = (prof or {}).get(name)
        if 'wspd' in sel[0]:
            u10 = era5.get(name, {}).get('u10')
            print('  surface    SH %5.1f LH %6.1f W/m2, lowest level wind '
                  '%5.2f m/s; ERA5 SH %s LH %s U10 %s' % (
                      mean([s['shf'] for s in sel]),
                      mean([s['lhf'] for s in sel]),
                      mean([s['wspd'] for s in sel]),
                      '%5.1f' % pc['shf'] if pc else '-',
                      '%6.1f' % pc['lhf'] if pc else '-',
                      '%5.2f' % mean(u10) if u10 else '-'))
        if 'omega' in sel[0]['levels'] and sel[0]['levels']['omega']:
            p = sel[0]['levels']['p']
            out = []
            for pl in (850.0, 700.0, 500.0, 300.0):
                k = min(range(len(p)), key=lambda i: abs(p[i] - pl))
                w = mean([s['levels']['omega'][k] for s in sel])
                ref = '-'
                if pc:
                    kk = min(range(len(pc['p'])),
                             key=lambda i: abs(pc['p'][i] - pl))
                    ref = '%6.3f' % pc['omega'][kk]
                out.append('%3.0f: %6.3f (%s)' % (pl, w, ref))
            print('  omega Pa/s (input) ' + '  '.join(out))


def compare(test, ref, rtol, atol_frac):
    """compare two outputs: integers exactly, reals within tolerance
    (rtol relative, plus atol_frac times the largest magnitude of the
    field in the reference step); returns the list of failures"""
    _, a = parse_output(test)
    _, b = parse_output(ref)
    fails = []
    if list(a) != list(b) or any(len(a[n]) != len(b[n]) for n in a):
        return ['different columns or steps']
    for name in a:
        for sa, sb in zip(a[name], b[name]):
            where = '%s %.2f h' % (name, sb['time'])
            for n in STEP_FIELDS:
                if n not in sa or n not in sb:
                    continue
                u, v = sa[n], sb[n]
                if n in INT_FIELDS:
                    if u != v:
                        fails.append('%s %s %d != %d' % (where, n, u, v))
                elif abs(u - v) > rtol * abs(v) + atol_frac * abs(v) + \
                        1e-12:
                    fails.append('%s %s %g != %g' % (where, n, u, v))
            for n in LEVEL_FIELDS:
                xa, xb = sa['levels'][n], sb['levels'][n]
                if not xa or not xb:
                    continue
                scale = max(abs(x) for x in xb)
                for k, (u, v) in enumerate(zip(xa, xb)):
                    if abs(u - v) > rtol * abs(v) + atol_frac * scale + \
                            1e-30:
                        fails.append('%s level %d %s %g != %g' %
                                     (where, k + 1, n, u, v))
                        break
    return fails


def main():
    ap = argparse.ArgumentParser(
        description='Coupled Holtslag + Tiedtke single column test (RegCM)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Parameters: ' + ', '.join(sorted(GROUP_OF)))
    ap.add_argument('--exe', help='scm_coupled program')
    ap.add_argument('--profiles', default=os.path.join(HERE,
                                                       'scm_amazon.txt'))
    ap.add_argument('--era5', default=None,
                    help='ERA5 rain by local hour (default: profiles '
                         'name with _era5)')
    ap.add_argument('--set', action='append', default=[],
                    metavar='NAME=VALUE', help='namelist setting')
    ap.add_argument('--sweep', metavar='NAME=V1,V2,...',
                    help='one run per value, summary per value')
    ap.add_argument('--list', action='store_true',
                    help='print the parameters and their values')
    ap.add_argument('--save', metavar='FILE', help='keep the output file')
    ap.add_argument('--keep', metavar='DIR', help='working directory')
    ap.add_argument('--compare', metavar='REF',
                    help='compare the run with a reference output')
    ap.add_argument('--check', nargs=2, metavar=('TEST', 'REF'),
                    help='only compare two output files')
    ap.add_argument('--rtol', type=float,
                    default=float(os.environ.get('SCM_RTOL', '1e-4')))
    ap.add_argument('--atol', type=float,
                    default=float(os.environ.get('SCM_ATOL', '1e-5')),
                    help='absolute tolerance as a fraction of the '
                         'largest value of the field in the step')
    a = ap.parse_args()

    if a.check:
        fails = compare(a.check[0], a.check[1], a.rtol, a.atol)
        for x in fails[:20]:
            print('DIFF: %s' % x)
        print('FAIL' if fails else 'PASS', ': %d differences' % len(fails))
        sys.exit(1 if fails else 0)

    exe = find_exe(a.exe)
    era5 = read_era5(a.era5 or os.path.splitext(a.profiles)[0] +
                     '_era5.txt')
    prof = read_profiles(a.profiles)
    settings = dict(s.split('=', 1) for s in a.set)
    if a.sweep:
        name, values = a.sweep.split('=', 1)
        for v in values.split(','):
            s = dict(settings); s[name] = v
            (_, cols), _, work = run(exe, s, a.profiles)
            summary(cols, era5, prof, '\n== %s = %s' % (name, v))
            shutil.rmtree(work, ignore_errors=True)
        return
    (header, cols), out, work = run(exe, settings, a.profiles, a.keep)
    if a.list:
        print('\n'.join(header))
    else:
        summary(cols, era5, prof)
    status = 0
    if a.compare:
        fails = compare(out, a.compare, a.rtol, a.atol)
        for x in fails[:20]:
            print('DIFF: %s' % x)
        print('compare with %s: %s (%d differences)' % (
            a.compare, 'FAIL' if fails else 'PASS', len(fails)))
        status = 1 if fails else 0
    if a.save:
        shutil.copy(out, a.save)
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    sys.exit(status)


if __name__ == '__main__':
    main()
