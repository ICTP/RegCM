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
# Compact digest of a single column test output, kept as the reference of
# "make check" in place of the full output (reference/<case>.dig).
#
# One line per quantity and column: the column quantities, and for each
# level field its pressure weighted column mean and its largest absolute
# value. For the coupled test (time series) the values are sampled at up to
# eight output times evenly spaced to the end of the run, and the final T
# and q profiles are added. Values are written with five significant
# digits.
#
# Usage:
#   scm_digest.py make KIND OUTPUT DIGEST     write the digest of an output
#   scm_digest.py check KIND OUTPUT DIGEST    compare an output with a digest
# KIND is tiedtke, holtslag or coupled. The check compares integers exactly
# and reals within rtol relative plus atol times the largest magnitude on
# the line (SCM_RTOL, SCM_ATOL, as the full output comparison).

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

NSAMPLE = 8
PROFILES = ('T', 'q')


def fmt(x):
    if isinstance(x, int):
        return '%d' % x
    return '%.5g' % x


def colstats(p, f):
    """pressure weighted mean and largest absolute value of a profile"""
    n = len(p)
    if n == 0:
        return 0.0, 0.0
    if n == 1:
        return f[0], abs(f[0])
    w = []
    for k in range(n):
        a = p[k - 1] if k > 0 else p[k]
        b = p[k + 1] if k < n - 1 else p[k]
        w.append(abs(b - a) * 0.5 if 0 < k < n - 1 else abs(b - a))
    sw = sum(w)
    mean = sum(x * y for x, y in zip(w, f)) / sw if sw > 0 else 0.0
    return mean, max(abs(x) for x in f)


def digest_single(mod, path):
    """tiedtke, holtslag: one call per column"""
    _, cols = mod.parse_output(path)
    ints = getattr(mod, 'INT_FIELDS', ('land', 'ktype', 'kctop', 'kcbot'))
    lines = ['COLUMNS ' + ' '.join(c['name'] for c in cols)]
    for n in mod.COLUMN_FIELDS:
        vals = []
        for c in cols:
            v = c.get(n, 0)
            vals.append(int(v) if n in ints else float(v))
        lines.append('%s %s' % (n, ' '.join(fmt(v) for v in vals)))
    for n in mod.LEVEL_FIELDS:
        if n in ('z', 'p'):
            continue
        st = [colstats(c['levels']['p'], c['levels'][n]) for c in cols]
        lines.append('%s:mean %s' % (n, ' '.join(fmt(s[0]) for s in st)))
        lines.append('%s:max %s' % (n, ' '.join(fmt(s[1]) for s in st)))
    return lines


def digest_coupled(mod, path):
    _, cols = mod.parse_output(path)
    lines = []
    for name, steps in cols.items():
        ns = len(steps)
        m = min(NSAMPLE, ns)
        idx = sorted(set(round((i + 1) * (ns - 1) / m) for i in range(m)))
        sel = [steps[i] for i in idx]
        lines.append('COLUMN %s' % name)
        for n in mod.STEP_FIELDS:
            if n not in sel[0]:
                continue
            lines.append('%s %s' % (n, ' '.join(fmt(s[n]) for s in sel)))
        for n in mod.LEVEL_FIELDS:
            if n in ('z', 'p'):
                continue
            st = [colstats(s['levels']['p'], s['levels'][n]) for s in sel]
            lines.append('%s:mean %s' % (n, ' '.join(fmt(x[0]) for x in st)))
            lines.append('%s:max %s' % (n, ' '.join(fmt(x[1]) for x in st)))
        for n in PROFILES:
            lines.append('%s:final %s' % (
                n, ' '.join(fmt(x) for x in steps[-1]['levels'][n])))
    return lines


def module(kind):
    if kind == 'tiedtke':
        import scm_tiedtke as mod
    elif kind == 'holtslag':
        import scm_holtslag as mod
    elif kind == 'coupled':
        import scm_coupled as mod
    else:
        sys.exit('unknown kind ' + kind)
    return mod


def digest(kind, path):
    mod = module(kind)
    if kind == 'coupled':
        return digest_coupled(mod, path)
    return digest_single(mod, path)


def parse_value(x):
    try:
        return int(x)
    except ValueError:
        return float(x)


def compare(test, ref, rtol, atol):
    """compare two digests (lists of lines); returns the failures"""
    fails = []
    if len(test) != len(ref):
        return ['different number of lines %d != %d' % (len(test), len(ref))]
    col = ''
    for a, b in zip(test, ref):
        wa, wb = a.split(), b.split()
        if wa[0] in ('COLUMN', 'COLUMNS'):
            if wa != wb:
                fails.append('different columns: %s / %s' % (a, b))
            col = ' '.join(wb[1:2])
            continue
        if wa[0] != wb[0] or len(wa) != len(wb):
            fails.append('different line: %s / %s' % (wa[0], wb[0]))
            continue
        va = [parse_value(x) for x in wa[1:]]
        vb = [parse_value(x) for x in wb[1:]]
        scale = max(abs(x) for x in vb) if vb else 0.0
        for i, (u, v) in enumerate(zip(va, vb)):
            if isinstance(v, int) and isinstance(u, int):
                bad = u != v
            else:
                bad = abs(u - v) > rtol * abs(v) + atol * scale + 1e-30
            if bad:
                fails.append('%s %s [%d] %s != %s' % (col, wb[0], i + 1,
                                                       wa[i + 1], wb[i + 1]))
                break
    return fails


def main():
    if len(sys.argv) != 5 or sys.argv[1] not in ('make', 'check'):
        sys.exit('usage: scm_digest.py make|check KIND OUTPUT DIGEST')
    what, kind, out, dig = sys.argv[1:]
    lines = digest(kind, out)
    if what == 'make':
        with open(dig, 'w') as f:
            f.write('# digest of %s (scm_digest.py, %s)\n' %
                    (os.path.basename(out), kind))
            f.write('\n'.join(lines) + '\n')
        return
    with open(dig) as f:
        ref = [x.rstrip('\n') for x in f if not x.startswith('#')]
    rtol = float(os.environ.get('SCM_RTOL', '1e-4'))
    atol = float(os.environ.get('SCM_ATOL', '1e-5'))
    fails = compare(lines, ref, rtol, atol)
    for x in fails[:20]:
        print('DIFF: ' + x)
    print('%s : %d differences' % ('FAIL' if fails else 'PASS', len(fails)))
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
