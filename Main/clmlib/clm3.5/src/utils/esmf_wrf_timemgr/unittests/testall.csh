#!/bin/csh
#
# Build and run Test1.exe and compare results with known-good output
#

foreach type ( GREGORIAN NOLEAP )
  echo "Run for $type calendar"
  if ( $type == "GREGORIAN" )then
     set correctfile = "Test1.out.correct"
  else
     set correctfile = "Test1.out.correct.noleap"
  endif
  gmake superclean
  set arg=""
  if ( $type == "GREGORIAN" ) set arg="WRFTEST=TRUE"
  gmake tests $arg
  set ok = $status
  if ( $ok != 0 )then
    echo "ERROR problem building tests"
    exit 9999
  endif
  echo "Run test suite..."
  ./Test1.exe >&! Test1.out
  set ok = $status
  if ( $ok != 0 )then
    echo "ERROR problem with testing!"
    exit 9999
  endif
  grep FAIL Test1.out
  set ok = $status
  if ( $ok != 1 )then
    echo "ERROR problem with testing, some tests fail!"
    exit 9999
  endif
  which xxdiff >& /dev/null
  set ok = $status
  if ( $ok == 0 ) then
    xxdiff $correctfile Test1.out
    set diffok = $status
  else
    diff $correctfile Test1.out
    set diffok = $status
  endif
  if ( $diffok != 0 )then
    echo "ERROR problem with testing!"
    exit 9999
  endif
  echo "$type Testing Passed"
  echo "Results from tests"
  cat Test1.out
  gmake superclean
end
echo "Testing Passed\!\!\!"
