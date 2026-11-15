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
# Python wrapper of the Holtslag single column test (scm_holtslag).
#
# Runs the Fortran program on the test profiles with namelist settings
# given on the command line, prints a summary per column, optionally
# plots the profiles, sweeps a parameter over a list of values, or
# compares an output with a reference.
#
# Examples:
#   scm_holtslag.py                             one call, RegCM defaults
#   scm_holtslag.py --set nsteps=90             6 hours with dtpbl = 240 s
#   scm_holtslag.py --list                      show the parameters
#   scm_holtslag.py --set lmoistn2=.true. --set maxprandtl=4.0
#   scm_holtslag.py --sweep ricr_ocn=0.25,0.275,0.3
#   scm_holtslag.py --plot holtslag.png         profiles (needs matplotlib)
#   scm_holtslag.py --save run.out --compare reference/scm_holtslag.out
#   scm_holtslag.py --check test.out reference/scm_holtslag.out
#
# The program is looked for in --exe, $SCM_HOLTSLAG, the current directory
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
    'scmparam': ['dtpbl', 'nsteps'],
    'holtslagparam': ['ricr_ocn', 'ricr_lnd', 'zhnew_fac', 'ifaholtth10',
                      'ifaholt', 'holtth10iter', 'lmoistn2', 'maxprandtl'],
}
GROUP_OF = {v: g for g, vs in GROUPS.items() for v in vs}

LEVEL_FIELDS = ['z', 'p', 'dT', 'dq', 'dqc', 'du', 'dv', 'Tdiff', 'qdiff',
                'udiff', 'vdiff']
COLUMN_FIELDS = ['land', 'kpbl', 'zpbl', 'ustar', 'heat', 'moist', 'shf',
                 'lhf', 'era5_zpbl']
INT_FIELDS = ('land', 'kpbl')


def find_exe(exe):
    for c in (exe, os.environ.get('SCM_HOLTSLAG'),
              os.path.join(os.getcwd(), 'scm_holtslag'),
              os.path.join(HERE, 'scm_holtslag')):
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.abspath(c)
    sys.exit('scm_holtslag program not found: build it with "make check" '
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
            c[n] = int(x) if n in INT_FIELDS else float(x)
        c['levels'] = {n: [] for n in LEVEL_FIELDS}
        for _ in range(nk):
            for n, x in zip(LEVEL_FIELDS, lines[i].split()):
                c['levels'][n].append(float(x))
            i += 1
        cols.append(c)
    return header, cols


def run(exe, settings, profiles, keep=None):
    work = os.path.abspath(keep or tempfile.mkdtemp(prefix='scm_holtslag_'))
    os.makedirs(work, exist_ok=True)
    out = os.path.join(work, 'scm_holtslag.out')
    nml = os.path.join(work, 'scm_holtslag.in')
    write_namelist(nml, settings, os.path.abspath(profiles), out)
    r = subprocess.run([exe, nml], cwd=work, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.isfile(out):
        sys.exit('scm_holtslag failed:\n' + r.stdout + r.stderr)
    res = parse_output(out)
    return res, out, work


def summary(cols, title=None):
    if title:
        print(title)
    print('%-20s %4s %5s %7s %7s %6s %8s %8s %8s %8s' % (
        'column', 'land', 'kpbl', 'zpbl', 'ERA5', 'ustar', 'heating',
        'shf', 'moisten', 'lhf'))
    print('%-20s %4s %5s %7s %7s %6s %8s %8s %8s %8s' % (
        '', '', '', '(m)', 'zpbl', '(m/s)', 'W/m2', 'W/m2', 'W/m2', 'W/m2'))
    for c in cols:
        print('%-20s %4d %5d %7.0f %7.0f %6.3f %8.1f %8.1f %8.1f %8.1f'
              % (c['name'], c['land'], c['kpbl'], c['zpbl'],
                 c['era5_zpbl'], c['ustar'], c['heat'], c['shf'],
                 c['moist'], c['lhf']))


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
        for n in ('kpbl',):
            if ca[n] != cb[n]:
                fails.append('%s %s %d != %d' % (ca['name'], n, ca[n], cb[n]))
        for n in ('zpbl', 'ustar', 'heat', 'moist'):
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


def plot(cols, png, diff=False):
    """first call tendencies, or the changes after nsteps if diff"""
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
        if diff:
            ax.plot(c['levels']['Tdiff'], z, label='T-T0 (K)')
            ax.plot([2.5 * x for x in c['levels']['qdiff']], z,
                    label='L/cp (q-q0) (K)')
            ax.plot(c['levels']['udiff'], z, ':', label='u-u0 (m/s)')
        else:
            ax.plot(c['levels']['dT'], z, label='dT (K/day)')
            ax.plot([2.5 * x for x in c['levels']['dq']], z,
                    label='L/cp dq (K/day)')
            ax.plot(c['levels']['du'], z, ':', label='du (m/s/day)')
        ax.axhline(c['zpbl'] / 1000.0, color='k', lw=0.8, ls='--')
        ax.axhline(c['era5_zpbl'] / 1000.0, color='grey', lw=0.8, ls=':')
        ax.axvline(0, color='grey', lw=0.5)
        ax.set_ylim(0, max(3.0, 1.5 * c['zpbl'] / 1000.0))
        ax.set_title('%s: zpbl %.0f m (ERA5 %.0f)' % (
            c['name'], c['zpbl'], c['era5_zpbl']), fontsize=8)
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
        description='Holtslag single column test (RegCM)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Parameters: ' + ', '.join(sorted(GROUP_OF)))
    ap.add_argument('--exe', help='scm_holtslag program')
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
        plot(cols, a.plot, diff=int(settings.get('nsteps', '1')) > 1)
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
