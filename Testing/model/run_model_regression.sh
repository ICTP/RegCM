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
# Run RegCM namelists and compare their output with reference runs.
#
# Usage: run_model_regression.sh [options] namelist [namelist ...]
#
#   -d dir   run directory, where the namelists refer to their input
#            (default: current directory)
#   -x exe   model executable (default: bin/regcmMPICLM45 or bin/regcmMPI
#            in the run directory)
#   -n np    number of MPI processes (default: 6)
#   -r dir   reference directory (default: <run directory>/reference)
#   -o dir   test output directory (default: <run directory>/regression)
#   -s       save the runs as the new reference instead of comparing
#   -c       compare the existing test output, do not run the model
#
# Each namelist is run with its dirout set to <test dir>/<namelist name>
# and compared with <reference dir>/<namelist name>. The restart files are
# not kept in the reference. MPIEXEC and MPIEXEC_FLAGS select the MPI
# launcher. The exit status is 0 only if every run completes and is bit for
# bit identical to its reference.

here=$(cd $(dirname $0) && pwd)
rundir=$PWD
exe=""
nproc=6
refdir=""
outdir=""
save=no
runmodel=yes

while getopts "d:x:n:r:o:sc" opt; do
  case $opt in
    d) rundir=$(cd $OPTARG && pwd) ;;
    x) exe=$OPTARG ;;
    n) nproc=$OPTARG ;;
    r) refdir=$OPTARG ;;
    o) outdir=$OPTARG ;;
    s) save=yes ;;
    c) runmodel=no ;;
    *) sed -n '16,32p' $0 | cut -c3- ; exit 2 ;;
  esac
done
shift $((OPTIND-1))
if [ $# -lt 1 ]; then
  sed -n '16,32p' $0 | cut -c3-
  exit 2
fi

refdir=${refdir:-$rundir/reference}
outdir=${outdir:-$rundir/regression}
if [ -z "$exe" ]; then
  for e in $rundir/bin/regcmMPICLM45 $rundir/bin/regcmMPI; do
    [ -x $e ] && exe=$e && break
  done
fi
mpiexec=${MPIEXEC:-mpirun}
mkdir -p $refdir $outdir
refdir=$(cd $refdir && pwd)
outdir=$(cd $outdir && pwd)

status=0
for namelist in "$@"; do
  [ -f $namelist ] || namelist=$rundir/$namelist
  name=$(basename $namelist .in)
  test=$outdir/$name
  if [ $runmodel = yes ]; then
    if [ ! -x "$exe" ]; then
      echo "No model executable found: use -x"
      exit 2
    fi
    rm -rf $test
    mkdir -p $test
    sed "s|^ *dirout *=.*| dirout = '$test',|" $namelist > $test/$name.in
    start=$(date +%s)
    ( cd $rundir && $mpiexec $MPIEXEC_FLAGS -np $nproc $exe \
        $test/$name.in > $test/$name.log 2>&1 )
    seconds=$(( $(date +%s) - start ))
    if ! grep -a -q "simulation successfully reached end" $test/$name.log
    then
      echo "FAIL: $name did not complete (see $test/$name.log)"
      status=1
      continue
    fi
    echo "$name: completed in $seconds s"
  fi
  if [ $save = yes ]; then
    rm -rf $refdir/$name
    mkdir -p $refdir/$name
    for stream in ATM SRF RAD STS SHF LAK OPT CHE; do
      cp $test/*_$stream.*.nc $refdir/$name/ 2> /dev/null
    done
    echo "$name: saved as reference in $refdir/$name"
  else
    echo "$name: comparing with $refdir/$name"
    python3 $here/compare_model_output.py $refdir/$name $test || status=1
  fi
done
exit $status
