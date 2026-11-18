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

# Test cases, scheme:profile file:name. The output reference/<name>.out
# is the reference of the case. The coupled cases run one day, with the
# relaxation targets of their profiles; the _kf_ ones with the
# Kain-Fritsch cumulus (icumulus = 6) instead of Tiedtke.
cases="tiedtke:scm_profiles.txt:scm_tiedtke
       holtslag:scm_profiles.txt:scm_holtslag
       tiedtke:scm_gate_check.txt:scm_tiedtke_gate
       coupled:scm_amazon.txt:scm_coupled_amazon
       coupled:scm_itcz.txt:scm_coupled_itcz
       coupled:scm_amazon.txt:scm_coupled_kf_amazon
       coupled:scm_itcz.txt:scm_coupled_kf_itcz"

for case in $cases; do
  scheme=${case%%:*}
  rest=${case#*:}
  profiles=${rest%%:*}
  name=${rest#*:}
  prog=$PWD/scm_$scheme
  [ -x $prog ] || prog=$here/scm_$scheme
  d=$outdir/$name
  rm -rf $d
  mkdir -p $d
  extra=
  if [ $scheme = coupled ]; then
    extra="tgtfile = '$here/${profiles%.txt}_target.txt', ndays = 1.0,"
    # the ITCZ case also runs the WTG omega and the wind relaxation
    case $name in
      *itcz) extra="$extra tauwtg = 10800.0, tauwind = 10800.0," ;;
    esac
    case $name in
      *_kf_*) extra="$extra icumulus = 6," ;;
    esac
  fi
  cat > $d/$name.in << EOF
 &scmparam
  infile = '$here/$profiles',
  outfile = '$d/$name.out',
  $extra
 /
EOF
  if (cd $d && $prog $name.in > log 2>&1) && \
     grep -q "columns written" $d/log; then
    echo "PASS: $name run"
  else
    echo "FAIL: $name run, see $d/log"
    status=1
    continue
  fi
  ref=$here/reference/$name.out
  if command -v python3 > /dev/null 2>&1; then
    if python3 $here/scm_$scheme.py --check $d/$name.out $ref; then
      echo "PASS: $name reference"
    else
      echo "FAIL: $name reference"
      status=1
    fi
  else
    echo "SKIP: $name reference (no python3)"
  fi
done
exit $status
