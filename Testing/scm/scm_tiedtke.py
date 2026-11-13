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
# Python wrapper of the Tiedtke single column test (scm_tiedtke).
#
# Runs the Fortran program on the test profiles with namelist settings
# given on the command line, prints a summary per column, optionally
# plots the profiles, sweeps a parameter over a list of values, or
# compares an output with a reference.
#
# Examples:
#   scm_tiedtke.py                              run with the RegCM defaults
#   scm_tiedtke.py --list                       show the parameters
#   scm_tiedtke.py --set entrpen_ocn=1.0e-3 --set lmfwstar=.true.
#   scm_tiedtke.py --sweep entshalp=1,2,3       one run per value
#   scm_tiedtke.py --plot tiedtke.png           profiles (needs matplotlib)
#   scm_tiedtke.py --save run.out --compare reference/scm_tiedtke.out
#   scm_tiedtke.py --check test.out reference/scm_tiedtke.out
#
# The program is looked for in --exe, $SCM_TIEDTKE, the current directory
# and the directory of this script. Only the standard library is needed,
# except for --plot.

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

# namelist group of every parameter the program reads
GROUPS = {
    'scmparam': ['dtcum', 'dxkm'],
    'tiedtkeparam': ['entrdd', 'entrpen_lnd', 'entrpen_ocn', 'detrpen_lnd',
                     'detrpen_ocn', 'entshalp', 'rcuc_lnd', 'rcuc_ocn',
                     'rcpec_lnd', 'rcpec_ocn', 'rhebc_lnd', 'rhebc_ocn',
                     'rprc_ocn', 'rprc_lnd', 'revap_lnd', 'revap_ocn',
                     'lmfpen', 'lmfmid', 'lmfdd', 'lepcld', 'lmfdudv',
                     'lmfscv', 'lmfuvdis', 'lmftrac', 'lmfsmooth',
                     'lmfwstar', 'lmfcapebl'],
    'cldparam': ['kfac_deep', 'kfac_shal', 'k2_const'],
}
GROUP_OF = {v: g for g, vs in GROUPS.items() for v in vs}

LEVEL_FIELDS = ['z', 'p', 'mfu', 'dT', 'dq', 'dqc', 'dqi', 'du', 'dv',
                'cldfrc', 'qdetr']
COLUMN_FIELDS = ['land', 'ktype', 'kctop', 'kcbot', 'rain', 'snow', 'heat',
                 'moist', 'era5_cape', 'era5_crr']
KTYPE = {0: '-', 1: 'deep', 2: 'shallow', 3: 'midlevel'}


def find_exe(exe):
    for c in (exe, os.environ.get('SCM_TIEDTKE'),
              os.path.join(os.getcwd(), 'scm_tiedtke'),
              os.path.join(HERE, 'scm_tiedtke')):
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.abspath(c)
    sys.exit('scm_tiedtke program not found: build it with "make check" '
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
    """header namelist text, list of columns (dict with levels)"""
    header, cols = [], []
    with open(path) as f:
        lines = f.read().splitlines()
    i = 0
    while not lines[i].startswith('NCOL'):
        header.append(lines[i]); i += 1
    ncol, nk = map(int, lines[i].split()[1:3]); i += 1
    for _ in range(ncol):
        w = lines[i].split(); i += 1
        c = {'name': w[1]}
        for n, x in zip(COLUMN_FIELDS, w[2:]):
            c[n] = int(x) if n in ('land', 'ktype', 'kctop', 'kcbot') \
                else float(x)
        c['levels'] = {n: [] for n in LEVEL_FIELDS}
        for _ in range(nk):
            for n, x in zip(LEVEL_FIELDS, lines[i].split()):
                c['levels'][n].append(float(x))
            i += 1
        cols.append(c)
    return header, cols


def run(exe, settings, profiles, keep=None):
    work = keep or tempfile.mkdtemp(prefix='scm_tiedtke_')
    os.makedirs(work, exist_ok=True)
    out = os.path.join(work, 'scm_tiedtke.out')
    nml = os.path.join(work, 'scm_tiedtke.in')
    write_namelist(nml, settings, os.path.abspath(profiles), out)
    r = subprocess.run([exe, nml], cwd=work, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.isfile(out):
        sys.exit('scm_tiedtke failed:\n' + r.stdout + r.stderr)
    res = parse_output(out)
    return res, out, work


def cloud_height(c, k):
    if k <= 0:
        return '    -'
    return '%5.0f' % c['levels']['z'][k - 1]


def summary(cols, title=None):
    if title:
        print(title)
    print('%-20s %4s %-8s %5s %5s %8s %9s %9s %8s %7s %6s' % (
        'column', 'land', 'type', 'base', 'top', 'mfu_max', 'rain',
        'heating', 'moisten', 'ERA5', 'ERA5'))
    print('%-20s %4s %-8s %5s %5s %8s %9s %9s %8s %7s %6s' % (
        '', '', '', '(m)', '(m)', 'kg/m2/s', 'mm/day', 'W/m2', 'W/m2',
        'CAPE', 'crr'))
    for c in cols:
        print('%-20s %4d %-8s %5s %5s %8.4f %9.3f %9.1f %8.1f %7.0f %6.1f'
              % (c['name'], c['land'], KTYPE.get(c['ktype'], '?'),
                 cloud_height(c, c['kcbot']), cloud_height(c, c['kctop']),
                 max(c['levels']['mfu']), c['rain'] + c['snow'], c['heat'],
                 c['moist'], c['era5_cape'], c['era5_crr']))


def compare(test, ref, rtol, atol_frac):
    """compare two outputs: integers exactly, reals within tolerance
    (rtol relative, plus atol_frac times the largest magnitude of the
    field in the reference column); returns the list of failures"""
    _, a = parse_output(test)
    _, b = parse_output(ref)
    fails = []
    if [c['name'] for c in a] != [c['name'] for c in b]:
        return ['different columns']
    for ca, cb in zip(a, b):
        for n in ('ktype', 'kctop', 'kcbot'):
            if ca[n] != cb[n]:
                fails.append('%s %s %d != %d' % (ca['name'], n, ca[n], cb[n]))
        for n in ('rain', 'snow', 'heat', 'moist'):
            if abs(ca[n] - cb[n]) > rtol * abs(cb[n]) + 1e-12 + \
               atol_frac * abs(cb[n]):
                fails.append('%s %s %g != %g' % (ca['name'], n, ca[n], cb[n]))
        for n in LEVEL_FIELDS:
            xa, xb = ca['levels'][n], cb['levels'][n]
            scale = max(abs(x) for x in xb)
            for k, (u, v) in enumerate(zip(xa, xb)):
                if abs(u - v) > rtol * abs(v) + atol_frac * scale + 1e-30:
                    fails.append('%s level %d %s %g != %g' %
                                 (ca['name'], k + 1, n, u, v))
                    break
    return fails


def plot(cols, png):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n = len(cols)
    nx = 4
    ny = (n + nx - 1) // nx
    fig, axs = plt.subplots(ny, nx, figsize=(4 * nx, 3.6 * ny),
                            squeeze=False)
    for ax, c in zip(axs.flat, cols):
        z = [x / 1000.0 for x in c['levels']['z']]
        ax.plot(c['levels']['dT'], z, label='dT (K/day)')
        ax.plot([2.5 * x for x in c['levels']['dq']], z,
                label='L/cp dq (K/day)')
        ax2 = ax.twiny()
        ax2.plot(c['levels']['mfu'], z, 'k--', lw=0.8, label='mfu')
        ax2.set_xlabel('mfu (kg/m2/s)', fontsize=7)
        ax2.tick_params(labelsize=7)
        ax.axvline(0, color='grey', lw=0.5)
        ax.set_ylim(0, 16)
        ax.set_title('%s: %s, %.2f mm/day' % (
            c['name'], KTYPE.get(c['ktype'], '?'), c['rain'] + c['snow']),
            fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_ylabel('z (km)', fontsize=7)
    for ax in list(axs.flat)[n:]:
        ax.axis('off')
    axs.flat[0].legend(fontsize=7, loc='upper right')
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    print('plot written to %s' % png)


def main():
    ap = argparse.ArgumentParser(
        description='Tiedtke single column test (RegCM)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Parameters: ' + ', '.join(sorted(GROUP_OF)))
    ap.add_argument('--exe', help='scm_tiedtke program')
    ap.add_argument('--profiles', default=os.path.join(HERE,
                                                       'scm_profiles.txt'))
    ap.add_argument('--set', action='append', default=[],
                    metavar='NAME=VALUE', help='namelist setting')
    ap.add_argument('--sweep', metavar='NAME=V1,V2,...',
                    help='one run per value, summary table per value')
    ap.add_argument('--list', action='store_true',
                    help='print the parameters and their values')
    ap.add_argument('--plot', metavar='PNG', help='plot the profiles')
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
                         'largest value of the field in the column')
    a = ap.parse_args()

    if a.check:
        fails = compare(a.check[0], a.check[1], a.rtol, a.atol)
        for x in fails[:20]:
            print('DIFF: %s' % x)
        print('FAIL' if fails else 'PASS', ': %d differences' % len(fails))
        sys.exit(1 if fails else 0)

    exe = find_exe(a.exe)
    settings = dict(s.split('=', 1) for s in a.set)
    if a.sweep:
        name, values = a.sweep.split('=', 1)
        for v in values.split(','):
            s = dict(settings); s[name] = v
            (_, cols), _, work = run(exe, s, a.profiles)
            summary(cols, '\n== %s = %s' % (name, v))
            shutil.rmtree(work, ignore_errors=True)
        return
    (header, cols), out, work = run(exe, settings, a.profiles, a.keep)
    if a.list:
        print('\n'.join(header))
    else:
        summary(cols)
    if a.plot:
        plot(cols, a.plot)
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
