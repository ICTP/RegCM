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
# Compare two sets of mpp_regression dumps record by record.
#
# Usage: compare_mpp_regression.py <reference_dir> <test_dir>
#
# Each directory holds one sub-directory per run (np<N>_c<case>) with one
# res_<rank>.bin file per processor. Every record is a 48 character label,
# the element count and the element size, followed by the raw data.
# Exit status is 0 only if every record is bit for bit identical.

import glob
import os
import struct
import sys

LABEL = 48


def records(path):
    data = open(path, 'rb').read()
    # The byte order depends on the compiler conversion flags: the first
    # record holds 12 integers of 4 bytes.
    endian = '>' if struct.unpack('>ii', data[LABEL:LABEL+8]) == (12, 4) else '<'
    pos = 0
    out = []
    while pos < len(data):
        label = data[pos:pos+LABEL].decode().strip()
        count, size = struct.unpack(endian + 'ii', data[pos+LABEL:pos+LABEL+8])
        pos += LABEL + 8
        out.append((label, data[pos:pos+count*size]))
        pos += count * size
    return out


def main():
    if len(sys.argv) != 3:
        sys.exit('Usage: compare_mpp_regression.py <reference_dir> <test_dir>')
    ref, new = sys.argv[1:]
    files = sorted(glob.glob(os.path.join(ref, '*', 'res_*.bin')))
    if not files:
        sys.exit('No dumps found in ' + ref)
    nrec = 0
    bad = 0
    for fref in files:
        rel = os.path.relpath(fref, ref)
        fnew = os.path.join(new, rel)
        if not os.path.exists(fnew):
            print('MISSING', rel)
            bad += 1
            continue
        a, b = records(fref), records(fnew)
        if len(a) != len(b):
            print('RECORD COUNT', rel, len(a), len(b))
            bad += 1
            continue
        for (la, da), (lb, db) in zip(a, b):
            nrec += 1
            if la != lb or da != db:
                print('DIFF', rel, la)
                bad += 1
    print('%d records compared, %d differences' % (nrec, bad))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
