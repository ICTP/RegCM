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
# digest of the output (scm_digest.py) is compared with the reference digest
# reference/<case>.dig (tolerances from SCM_RTOL and SCM_ATOL).
#
# Environment:
#   SCM_OUTDIR     where the outputs are written    (default: ./scm_output)
#   SCM_FULLREF    directory of full reference outputs <case>.out: if set,
#                  the full outputs are also compared with them
#   SCM_RRTMDIR    directory of the RRTMG tables rrtmg_sw.nc and
#                  rrtmg_lw.nc (RegCM_Data/RRTM): without it the case
#                  with RRTMG is skipped
#   SCM_UPDATE     if set to 1, write the reference digests (and the full
#                  outputs to SCM_FULLREF if set) from this run instead of
#                  comparing

here=$(cd $(dirname $0) && pwd)
outdir=${SCM_OUTDIR:-$PWD/scm_output}
status=0

# Test cases, scheme:profile file:name. The output reference/<name>.out
# is the reference of the case. The coupled cases run one day, with the
# relaxation targets of their profiles; the _kf_, _em_ and _grell_ ones
# with the Kain-Fritsch (icumulus = 6), Emanuel (4) and Grell (2) cumulus
# instead of Tiedtke. The BOMEX (6
# hours) and RICO (one day) shallow cumulus cases and the FIRE (one day)
# stratocumulus case run with the settings of scm_les.py: DEPHY forcing,
# geostrophic wind of the target file, no relaxation; interactive
# radiation only in FIRE. FIRE also runs with the sequential split
# (isplit = 1, IFS order). BOMEX also runs with the UW, MYJ and Shin-Hong
# boundary layer schemes (ipbl = 2, 4, 5), FIRE with the WSM5 and WDM7
# microphysics (imicro = 3, 5) and the Amazon with WSM7 (imicro = 4);
# FIRE with RRTMG (irrtm = 1) if SCM_RRTMDIR gives its tables.
cases="tiedtke:scm_profiles.txt:scm_tiedtke
       holtslag:scm_profiles.txt:scm_holtslag
       tiedtke:scm_gate_check.txt:scm_tiedtke_gate
       coupled:scm_amazon.txt:scm_coupled_amazon
       coupled:scm_itcz.txt:scm_coupled_itcz
       coupled:scm_amazon.txt:scm_coupled_kf_amazon
       coupled:scm_itcz.txt:scm_coupled_kf_itcz
       coupled:scm_amazon.txt:scm_coupled_em_amazon
       coupled:scm_itcz.txt:scm_coupled_em_itcz
       coupled:scm_amazon.txt:scm_coupled_grell_amazon
       coupled:scm_itcz.txt:scm_coupled_grell_itcz
       coupled:scm_bomex.txt:scm_coupled_bomex
       coupled:scm_rico.txt:scm_coupled_rico
       coupled:scm_fire.txt:scm_coupled_fire
       coupled:scm_fire.txt:scm_coupled_fire_seq
       coupled:scm_bomex.txt:scm_coupled_bomex_uw
       coupled:scm_bomex.txt:scm_coupled_bomex_myj
       coupled:scm_bomex.txt:scm_coupled_bomex_shinhong
       coupled:scm_fire.txt:scm_coupled_fire_wsm5
       coupled:scm_fire.txt:scm_coupled_fire_wdm7
       coupled:scm_amazon.txt:scm_coupled_wsm7_amazon
       coupled:scm_fire.txt:scm_coupled_fire_rrtmg"

for case in $cases; do
  scheme=${case%%:*}
  rest=${case#*:}
  profiles=${rest%%:*}
  name=${rest#*:}
  prog=$PWD/scm_$scheme
  [ -x $prog ] || prog=$here/scm_$scheme
  d=$outdir/$name
  case $name in
    *_rrtmg)
      if [ -z "$SCM_RRTMDIR" ] || [ ! -f "$SCM_RRTMDIR/rrtmg_sw.nc" ]; then
        echo "SKIP: $name (no RRTMG tables, set SCM_RRTMDIR)"
        continue
      fi ;;
  esac
  rm -rf $d
  mkdir -p $d
  extra=
  if [ $scheme = coupled ]; then
    extra="tgtfile = '$here/${profiles%.txt}_target.txt',"
    # the ITCZ case also runs the WTG omega and the wind relaxation
    case $name in
      *itcz) extra="$extra ndays = 1.0, tauwtg = 10800.0, tauwind = 10800.0," ;;
      *bomex*) extra="$extra ndays = 0.25, irad = 0, taurelax = 1.0e12,
        lgeos = .true., iocnflux = 0, ustarfix = 0.28," ;;
      *rico) extra="$extra ndays = 1.0, irad = 0, taurelax = 1.0e12,
        lgeos = .true., iocnflux = 1, ocnsw = 0.0, ocnlw = 0.0," ;;
      *fire|*fire_wsm5|*fire_wdm7|*fire_rrtmg) extra="$extra ndays = 1.0, irad = 1,
        ladvt = .true.,
        taurelax = 1.0e12, lgeos = .true., iocnflux = 1," ;;
      *fire_seq) extra="$extra ndays = 1.0, irad = 1, ladvt = .true.,
        taurelax = 1.0e12, lgeos = .true., iocnflux = 1, isplit = 1,
        csplit = 'rad pbl cu mic'," ;;
      *) extra="$extra ndays = 1.0," ;;
    esac
    case $name in
      *_kf_*) extra="$extra icumulus = 6," ;;
      *_em_*) extra="$extra icumulus = 4," ;;
      *_grell_*) extra="$extra icumulus = 2," ;;
      *_uw) extra="$extra ipbl = 2," ;;
      *_myj) extra="$extra ipbl = 4," ;;
      *_shinhong) extra="$extra ipbl = 5," ;;
      *_wsm5) extra="$extra imicro = 3," ;;
      *_wsm7_*) extra="$extra imicro = 4," ;;
      *_wdm7) extra="$extra imicro = 5," ;;
      *_rrtmg) extra="$extra irrtm = 1, rrtmdir = '$SCM_RRTMDIR'," ;;
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
  ref=$here/reference/$name.dig
  if ! command -v python3 > /dev/null 2>&1; then
    echo "SKIP: $name reference (no python3)"
    continue
  fi
  if [ "${SCM_UPDATE:-0}" = 1 ]; then
    python3 $here/scm_digest.py make $scheme $d/$name.out $ref && \
      echo "UPDATED: $name reference digest"
    if [ -n "$SCM_FULLREF" ]; then
      mkdir -p $SCM_FULLREF && cp $d/$name.out $SCM_FULLREF/$name.out
    fi
    continue
  fi
  if python3 $here/scm_digest.py check $scheme $d/$name.out $ref; then
    echo "PASS: $name reference"
  else
    echo "FAIL: $name reference"
    status=1
  fi
  if [ -n "$SCM_FULLREF" ] && [ -f $SCM_FULLREF/$name.out ]; then
    if python3 $here/scm_$scheme.py --check $d/$name.out \
        $SCM_FULLREF/$name.out; then
      echo "PASS: $name full reference"
    else
      echo "FAIL: $name full reference"
      status=1
    fi
  fi
done
exit $status
