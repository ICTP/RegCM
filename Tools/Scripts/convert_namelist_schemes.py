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
# Convert the scheme selection of a RegCM namelist (&physicsparam) from the
# integer codes to the scheme names (Share/mod_scheme_names.F90):
#
#   ibltyp   -> pbl_scheme         icldfrac -> cloud_fraction
#   icup_lnd -> cumulus_land       irrtm    -> radiation_scheme
#   icup_ocn -> cumulus_ocean      icumcloud   removed
#   ipptls   -> microphysics
#
# and removes the parameters no more read: iconvlwp (&physicsparam), tc0
# and cllwcv (&cldparam).
#
# The SUBEX microphysics (ipptls = 1) is removed: it becomes noto, with a
# warning. The Kuo (1) and Betts-Miller (3) cumulus are removed: the
# conversion stops. The &subexparam group, no more read, is kept.
#
# Usage:
#   convert_namelist_schemes.py regcm.in            print the converted file
#   convert_namelist_schemes.py -i *.in             convert in place

import argparse
import re
import sys

NAMES = {
    'ibltyp': ('pbl_scheme', {0: 'none', 1: 'holtslag', 2: 'uw', 3: 'gfs',
                              4: 'myj', 5: 'shin-hong'}),
    'icup_lnd': ('cumulus_land', None),
    'icup_ocn': ('cumulus_ocean', None),
    'ipptls': ('microphysics', {0: 'none', 1: 'noto', 2: 'noto', 3: 'wsm5',
                                4: 'wsm7', 5: 'wdm7'}),
    'icldfrac': ('cloud_fraction', {0: 'sundqvist', 1: 'xu-randall',
                                    2: 'thompson', 3: 'gultepe-isaac',
                                    4: 'teixeira', 5: 'tompkins',
                                    6: 'echam5', 7: 'digital'}),
    'irrtm': ('radiation_scheme', {0: 'ccm', 1: 'rrtmg'}),
}
CUMULUS = {0: 'none', -1: 'shallow', 2: 'grell', 4: 'emanuel', 5: 'tiedtke',
           6: 'kain-fritsch'}
REMOVED_CUMULUS = {1: 'Kuo', 3: 'Betts-Miller'}

LINE = re.compile(r'^(\s*)(\w+)(\s*=\s*)([-+]?\d+)(\s*,?)(.*)$')


def convert(text, path):
    out, ingroup, incld = [], False, False
    for line in text.splitlines(keepends=True):
        s = line.strip().lower()
        if s.startswith('&'):
            ingroup = s.startswith('&physicsparam')
            incld = s.startswith('&cldparam')
        elif s.startswith('/'):
            ingroup = False
            incld = False
        k = s.split('=')[0].strip()
        if ((ingroup and k == 'iconvlwp') or
                (incld and k in ('tc0', 'cllwcv'))) and '=' in s:
            continue
        m = LINE.match(line.rstrip('\n')) if ingroup else None
        if not m or m.group(2).lower() not in NAMES and \
                m.group(2).lower() != 'icumcloud':
            out.append(line)
            continue
        indent, key, _, val, comma, rest = m.groups()
        key, code = key.lower(), int(val)
        if key == 'icumcloud':
            continue
        name, table = NAMES[key]
        if table is None:
            if code in REMOVED_CUMULUS:
                sys.exit('%s: %s = %d, the %s cumulus scheme is removed '
                         '(Main/historical): choose another one' %
                         (path, key, code, REMOVED_CUMULUS[code]))
            table = CUMULUS
        if code not in table:
            sys.exit('%s: %s = %d is not a valid code' % (path, key, code))
        if key == 'ipptls' and code == 1:
            print('%s: the SUBEX microphysics (ipptls = 1) is removed, '
                  'converted to noto' % path, file=sys.stderr)
        newline = "%s%s = '%s'%s" % (indent, name, table[code],
                                    comma or ',')
        # keep a trailing comment, dropping the description of the codes
        rest = rest.strip()
        if rest.startswith('!'):
            newline += ' ' + rest
        out.append(newline + '\n')
    return ''.join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.
                                 RawDescriptionHelpFormatter)
    ap.add_argument('-i', '--in-place', action='store_true')
    ap.add_argument('files', nargs='+')
    a = ap.parse_args()
    for path in a.files:
        with open(path) as f:
            new = convert(f.read(), path)
        if a.in_place:
            with open(path, 'w') as f:
                f.write(new)
        else:
            sys.stdout.write(new)


if __name__ == '__main__':
    main()
