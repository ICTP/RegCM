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
# Run mpp_regression for every case on several processor counts. Fails if a
# run does not complete. If MPP_REFERENCE points to the output directory of
# a reference build, the dumps are also compared to it record by record.
#
# Environment:
#   MPIEXEC        MPI launcher                     (default: mpirun)
#   MPIEXEC_FLAGS  extra launcher flags             (default: none)
#   MPP_NPROCS     processor counts to run          (default: "1 2 4 6")
#   MPP_OUTDIR     where the dumps are written      (default: ./mpp_output)
#   MPP_REFERENCE  dumps of a reference build       (default: no comparison)

here=$(cd $(dirname $0) && pwd)
prog=$PWD/mpp_regression
[ -x $prog ] || prog=$here/mpp_regression
mpiexec=${MPIEXEC:-mpirun}
nprocs=${MPP_NPROCS:-"1 2 4 6"}
outdir=${MPP_OUTDIR:-$PWD/mpp_output}

status=0
for np in $nprocs; do
  for c in 0 1 2 3; do
    d=$outdir/np${np}_c${c}
    rm -rf $d
    mkdir -p $d
    if $mpiexec $MPIEXEC_FLAGS -np $np $prog $d $c > $d/log 2>&1 && \
       grep -q "mpp_regression case *$c *done" $d/log; then
      echo "PASS: np=$np case=$c"
    else
      echo "FAIL: np=$np case=$c (see $d/log)"
      status=1
    fi
  done
done

if [ -n "$MPP_REFERENCE" ]; then
  python3 $here/compare_mpp_regression.py $MPP_REFERENCE $outdir || status=1
fi
exit $status
