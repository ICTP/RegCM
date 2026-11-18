#
#    This file is part of ICTP RegCM.
#
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
"""Catalog builder: finds the reference datasets known to rcmeval
(known.yaml) under a directory, or in a list of files of a remote machine,
writes the site catalog, and reports what is available, what is missing
and where to download it.

  python -m rcmeval.builder --scan /data/OBS -o mysite.yaml
  python -m rcmeval.builder --list files.txt --root /data/OBS -o mysite.yaml

A list of files is the output of 'cd /data/OBS; find . -name "*.nc*"'.
"""

import argparse
import os
import re
import sys
import yaml

from .units import TARGET

HERE = os.path.dirname(os.path.abspath(__file__))


def template_regex(t):
    """Regular expression of a file template: {year}, {month:02d} and *."""
    out = ''
    i = 0
    seen = set()
    while i < len(t):
        if t.startswith('{year}', i):
            out += r'(?P=year)' if 'year' in seen else r'(?P<year>\d{4})'
            seen.add('year')
            i += 6
        elif t.startswith('{month:02d}', i):
            out += r'(?P=month)' if 'month' in seen else r'(?P<month>\d{2})'
            seen.add('month')
            i += 11
        elif t[i] == '*':
            out += r'[^/]*'
            i += 1
        else:
            out += re.escape(t[i])
            i += 1
    return re.compile(r'^(?P<root>(?:.*/)?)' + out + '$')


def literal(t):
    """Longest literal piece of the file name of a template, a cheap
    filter before the regular expression."""
    base = t.split('/')[-1]
    pieces = re.split(r'\{[^}]*\}|\*', base)
    return max(pieces, key=len)


def list_files(scan=None, listfile=None):
    """Files under the scan directories (absolute paths), or of the list
    (relative to the directory of the list)."""
    files = []
    if scan:
        for top in scan:
            for dp, dn, fn in os.walk(top, followlinks=True):
                dn[:] = [d for d in dn if not d.startswith('.')]
                for f in fn:
                    if f.endswith(('.nc', '.nc4')):
                        files.append(os.path.join(os.path.abspath(dp), f))
    else:
        with open(listfile) as f:
            for line in f:
                p = line.strip()
                if p.endswith(('.nc', '.nc4')):
                    files.append(os.path.normpath(p))
    return files


def match(files, known):
    """For every known dataset the layout and directory with the most
    variables found: {name: (root, layout, {var: (years, [files])})}."""
    found = {}
    for name, ds in known.items():
        best = None
        for li, layout in enumerate(ds['layouts']):
            hits = {}
            for var, spec in layout.items():
                rx = template_regex(spec['file'])
                lit = literal(spec['file'])
                for f in files:
                    if lit not in f:
                        continue
                    m = rx.match(f)
                    if m:
                        root = m.group('root').rstrip('/')
                        y = m.groupdict().get('year')
                        h = hits.setdefault(root, {}).setdefault(
                            var, (set(), []))
                        if y:
                            h[0].add(int(y))
                        h[1].append(f)
            words = [w.lower() for w in ds.get('dirs', [name.split('_')[0]])]
            verified = all(spec.get('verified', True)
                           for spec in layout.values())
            for root, vh in hits.items():
                named = any(w in root.lower() for w in words)
                # more variables, a named directory, a verified layout,
                # the latest version (file names), more files
                latest = max(f for v in vh.values() for f in v[1])
                latest = os.path.basename(latest)
                score = (len(vh), named, verified, latest,
                         sum(len(v[1]) for v in vh.values()))
                if best is None or score > best[0]:
                    best = (score, root, li, vh)
        if best:
            found[name] = (best[1], best[2], best[3])
    return found


def years_text(ys):
    if not ys:
        return 'single file'
    ys = sorted(ys)
    gaps = len(ys) != ys[-1] - ys[0] + 1
    return f'{ys[0]}-{ys[-1]}' + (' (with gaps)' if gaps else '')


def top_dir(f, scan_roots):
    """First directory of f below its scan root."""
    for r in scan_roots or []:
        r = os.path.abspath(r)
        if f.startswith(r + '/'):
            f = f[len(r) + 1:]
            break
    return f.split('/')[0] if '/' in f else '.'


def build(files, prefix, known, want=None, scan_roots=None):
    found = match(files, known)
    catalog = {}
    lines = []
    used_tops = set()
    for name, (root, li, vh) in sorted(found.items()):
        ds = known[name]
        layout = ds['layouts'][li]
        entry = {'root': os.path.join(prefix, root) if prefix else
                 (root or '/'),
                 'frequency': ds.get('frequency', 'mon')}
        if ds.get('surface'):
            entry['surface'] = ds['surface']
        variables = {}
        for var, (ys, fl) in vh.items():
            spec = {k: v for k, v in layout[var].items() if k != 'verified'}
            if '*' in spec['file'] and '{' not in spec['file']:
                # versions in the name: the latest one
                spec['file'] = os.path.relpath(sorted(fl)[-1], root) \
                    if root else sorted(fl)[-1]
            variables[var] = spec
            used_tops.add(top_dir(fl[0], scan_roots))
        entry['variables'] = variables
        catalog[name] = entry
        yrs = set().union(*(v[0] for v in vh.values()))
        miss = [v for v in layout if v not in vh]
        unver = [v for v in vh if layout[v].get('verified') is False]
        lines.append(f'FOUND   {name:11s} {years_text(yrs):22s} '
                     f'{", ".join(sorted(vh))}  [{entry["root"]}]')
        if miss:
            lines.append(f'        {"":11s} missing variables: '
                         f'{", ".join(miss)}')
        if unver:
            lines.append(f'        {"":11s} check the units of: '
                         f'{", ".join(unver)} (layout not verified)')
    for name, ds in known.items():
        if name not in found:
            lines.append(f'ABSENT  {name:11s} {ds["description"]}; '
                         f'{ds["coverage"]}')
            lines.append(f'        {"":11s} get it from: {ds["access"]}')
    tops = {}
    for f in files:
        t = top_dir(f, scan_roots)
        tops[t] = tops.get(t, 0) + 1
    other = sorted((t for t in tops if t not in used_tops and t != '.'),
                   key=lambda t: -tops[t])
    if other:
        lines.append('')
        lines.append('Directories not recognized (add them to the site '
                     'catalog by hand, or to known.yaml):')
        for t in other:
            lines.append(f'        {t} ({tops[t]} files)')
    if want:
        lines.append('')
        lines.append('References by variable (found / downloadable):')
        for v in want:
            have = [n for n in catalog if v in catalog[n]['variables']]
            can = [n for n, ds in known.items() if n not in catalog and
                   any(v in lay for lay in ds['layouts'])]
            lines.append(f'  {v:8s} found: {", ".join(have) or "-":30s} '
                         f'downloadable: {", ".join(can) or "-"}')
    return catalog, lines


def main(argv=None):
    ap = argparse.ArgumentParser(prog='rcmeval.builder',
                                 description=__doc__.split('\n\n')[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--scan', action='append',
                   help='directory with the observations (repeatable)')
    g.add_argument('--list', help='file list (find . -name "*.nc*")')
    ap.add_argument('--root', default='',
                    help='directory of the list on its machine')
    ap.add_argument('--known', default=os.path.join(HERE, 'known.yaml'))
    ap.add_argument('--want', default=None,
                    help='variables of interest (comma separated), '
                    'default all those of the known datasets')
    ap.add_argument('-o', '--output', help='catalog to write')
    a = ap.parse_args(argv)
    with open(a.known) as f:
        known = yaml.safe_load(f)['datasets']
    files = list_files(a.scan, a.list)
    prefix = '' if a.scan else a.root
    if a.want:
        want = a.want.split(',')
    else:
        want = sorted({v for ds in known.values() for lay in ds['layouts']
                       for v in lay}, key=lambda v: list(TARGET).index(v)
                      if v in TARGET else 99)
    catalog, lines = build(files, prefix, known, want, a.scan)
    text = '\n'.join(lines)
    print(text)
    if a.output:
        with open(a.output, 'w') as f:
            f.write('# Catalog written by rcmeval.builder from '
                    f'{", ".join(a.scan) if a.scan else a.list}\n#\n')
            f.write(''.join('# ' + ln + '\n' for ln in lines))
            yaml.safe_dump({'datasets': catalog}, f, sort_keys=False,
                           width=100)
        print(f'\nWritten {a.output}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
