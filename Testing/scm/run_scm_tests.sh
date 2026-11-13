#!/bin/bash
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
# Run the single column tests on the test profiles with the default
# settings. Fails if a run does not complete. If python3 is available the
# output is compared with the reference in reference/ (tolerances from
# SCM_RTOL and SCM_ATOL, see scm_tiedtke.py).
#
# Environment:
#   SCM_OUTDIR     where the outputs are written    (default: ./scm_output)

here=$(cd $(dirname $0) && pwd)
outdir=${SCM_OUTDIR:-$PWD/scm_output}
status=0

for scheme in tiedtke; do
  prog=$PWD/scm_$scheme
  [ -x $prog ] || prog=$here/scm_$scheme
  d=$outdir/$scheme
  rm -rf $d
  mkdir -p $d
  cat > $d/scm_$scheme.in << EOF
 &scmparam
  infile = '$here/scm_profiles.txt',
  outfile = '$d/scm_$scheme.out',
 /
EOF
  if (cd $d && $prog scm_$scheme.in > log 2>&1) && \
     grep -q "columns written" $d/log; then
    echo "PASS: scm_$scheme run"
  else
    echo "FAIL: scm_$scheme run, see $d/log"
    status=1
    continue
  fi
  ref=$here/reference/scm_$scheme.out
  if command -v python3 > /dev/null 2>&1; then
    if python3 $here/scm_$scheme.py --check $d/scm_$scheme.out $ref; then
      echo "PASS: scm_$scheme reference"
    else
      echo "FAIL: scm_$scheme reference"
      status=1
    fi
  else
    echo "SKIP: scm_$scheme reference (no python3)"
  fi
done
exit $status
