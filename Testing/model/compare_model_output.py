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
# Compare the output files of two RegCM runs.
#
# Usage: compare_model_output.py <reference_dir> <test_dir>
#
# Every floating point variable of the ATM, SRF, RAD, STS, SHF, LAK, OPT and
# CHE files found in the reference directory is compared bit for bit with
# the same file in the test directory. For each differing variable the
# maximum absolute difference and the first differing output time step are
# reported. Exit status is 0 only if every variable is identical.

import glob
import os
import sys

import numpy as np
from netCDF4 import Dataset

STREAMS = ('ATM', 'SRF', 'RAD', 'STS', 'SHF', 'LAK', 'OPT', 'CHE')


def output_files(directory):
    files = []
    for stream in STREAMS:
        files += glob.glob(os.path.join(directory, '*_%s.*.nc' % stream))
    return sorted(os.path.basename(f) for f in files)


def first_difference(a, b):
    """Index of the first time step where a and b differ"""
    for it in range(a.shape[0]):
        if not np.array_equal(a[it], b[it], equal_nan=True):
            return it
    return None


def compare_file(reffile, testfile):
    ndiff = 0
    nvar = 0
    with Dataset(reffile) as ref, Dataset(testfile) as test:
        ref.set_auto_mask(False)
        test.set_auto_mask(False)
        for name, var in ref.variables.items():
            if var.dtype.kind != 'f':
                continue
            nvar += 1
            if name not in test.variables:
                print('  %-12s missing in the test file' % name)
                ndiff += 1
                continue
            a = var[:]
            b = test.variables[name][:]
            if a.shape != b.shape:
                print('  %-12s shape %s instead of %s' %
                      (name, b.shape, a.shape))
                ndiff += 1
                continue
            if np.array_equal(a, b, equal_nan=True):
                continue
            ndiff += 1
            diff = np.abs(a.astype(np.float64) - b.astype(np.float64))
            message = '  %-12s max |diff| %.6g' % (name, np.nanmax(diff))
            if 'time' in var.dimensions and var.dimensions[0] == 'time':
                it = first_difference(a, b)
                message += ', first differing time step %d of %d' % \
                    (it + 1, a.shape[0])
            print(message)
    return nvar, ndiff


def main():
    if len(sys.argv) != 3:
        print('Usage: %s <reference_dir> <test_dir>' % sys.argv[0])
        return 2
    refdir, testdir = sys.argv[1:]
    files = output_files(refdir)
    if not files:
        print('No RegCM output files in %s' % refdir)
        return 2
    total_var = 0
    total_diff = 0
    for f in files:
        testfile = os.path.join(testdir, f)
        if not os.path.exists(testfile):
            print('%s: missing in %s' % (f, testdir))
            total_diff += 1
            continue
        nvar, ndiff = compare_file(os.path.join(refdir, f), testfile)
        total_var += nvar
        total_diff += ndiff
        print('%s: %d variables, %s' %
              (f, nvar, 'identical' if ndiff == 0 else '%d differ' % ndiff))
    if total_diff == 0:
        print('IDENTICAL: %d variables in %d files' % (total_var, len(files)))
        return 0
    print('DIFFERENT: %d differences' % total_diff)
    return 1


if __name__ == '__main__':
    sys.exit(main())
