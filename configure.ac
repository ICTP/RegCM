#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#
#    ICTP RegCM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

AC_INIT([RegCM], [5.0], [regcm-developer@eforge.escience-lab.org])

AC_SUBST([PACKAGE_URL], ["https://eforge.escience-lab.org/gf/project/regcm"])

LT_INIT

AM_INIT_AUTOMAKE([foreign])

AC_PROG_CC
AC_PROG_FC

AC_LANG([Fortran])

AC_ARG_WITH([netcdf],
  AC_HELP_STRING([--with-netcdf],
            [Path to NetCDF installation (default: NETCDF environment)]),
  [NC_PREFIX="$withval"], [NC_PREFIX="${NETCDF}"]
)

AC_ARG_WITH([hdf5],
  AC_HELP_STRING([--with-hdf5],
            [Path to HDF5 installation (default: HDF5 environment)]),
  [HDF5_PREFIX="$withval"], [HDF5_PREFIX="${HDF5}"]
)

AC_ARG_WITH([szip],
  AC_HELP_STRING([--with-szip],
            [Path to SZIP installation (default: SZIP environment)]),
  [SZIP_PREFIX="$withval"], [SZIP_PREFIX="${SZIP}"]
)

PARALLEL="1"
AC_ARG_ENABLE(mpiserial,
  AC_HELP_STRING([--enable-mpiserial],
                 [Use the included MPI replacement library for single processor]),
  [DONOTCHECKMPI="DONOTCHECKMPI"; PARALLEL="0"]
)
if eval "test x$PARALLEL = x1"; then
  AM_CPPFLAGS="-DMPP1 $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

AC_ARG_ENABLE(selectedrealkind,
  AC_HELP_STRING([--enable-selectedrealkind],
                 [define single precision and double precision numbers using the selected_real_kind function. Default uses the kind inquiry function.]),
                 [SRKDEF="SELECTEDREALKIND"]
)

AC_ARG_ENABLE(sequence,
  AC_HELP_STRING([--enable-sequence],
                 [Modify MCT types to make them contiguous in memory.]),
  [SRKDEF="SEQUENCE"],
)

USEBAND="0"
AC_ARG_ENABLE(band,
  AC_HELP_STRING([--enable-band],
                 [Supply this option if you plan on using tropical band option.]),
  [USEBAND="1"],
)
AM_CONDITIONAL(DO_COMPILE_BAND, [test x"$USEBAND" = x1])
if eval "test x$USEBAND = x1"; then
  AM_CPPFLAGS="-DBAND $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

USECLM="0"
AC_ARG_ENABLE(clm,
  AC_HELP_STRING([--enable-clm],
                 [Supply this option if you plan on using CLM option.]),
  [USECLM="1"],
)
AM_CONDITIONAL(DO_COMPILE_LIBCLM, [test x"$USECLM" = x1])
if eval "test x"$USECLM" = x1"; then
  if eval "test x"$PARALLEL" = x0"; then
    AC_MSG_ERROR([CLM enabled code works only with MPI version])
  fi
  AC_MSG_WARN([Compiling CLM enabled code])
  AM_CPPFLAGS="-DCLM $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

AC_ARG_ENABLE(nc4-gzip,
  AC_HELP_STRING([--enable-nc4-gzip],
                 [Supply this option if you plan on using NetCDF4/HDF5 output compression option.]),
  [NETCDF4_HDF5="1"],[NETCDF4_HDF5="0"]
)
if eval "test x$NETCDF4_HDF5 = x1"; then
  AM_CPPFLAGS="-DNETCDF4_HDF5 $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

DEBUG="0"
AC_ARG_ENABLE(debug,
  AC_HELP_STRING([--enable-debug],
                 [Enable debugging flags and per processor log file]),
  [DEBUG="1"]
)
if eval "test x$DEBUG = x1"; then
  AM_CPPFLAGS="-DDEBUG $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

if eval "test -d $NC_PREFIX/include"; then
  NC_INCLUDES="-I$NC_PREFIX/include"
else
  NC_PREFIX="/usr"
  NC_INCLUDES="-I/usr/include"
fi

if eval "test -d $NC_PREFIX/lib"; then
  NC_LDFLAGS="-L$NC_PREFIX/lib"
  NC_LIBS="-lnetcdf"
else
  NC_PREFIX="/usr"
  NC_LDFLAGS="-L/usr/lib"
  NC_LIBS="-lnetcdf"
fi
AC_SUBST(NC_PREFIX)

AC_LANG_PUSH([C])
AC_C_BIGENDIAN
AC_LANG_POP([C])
__AC_FC_NAME_MANGLING

case $ac_cv_fc_mangling in
  "lower case, underscore, no extra underscore")
    name_mangling="FORTRAN_UNDERSCORE_";;
  "lower case, no underscore")
    name_mangling="FORTRAN_SAME";;
  "upper case, no underscore")
    name_mangling="FORTRAN_CAPS_";;
  "lower case, double underscore")
    name_mangling="FORTRAN_GNUF2C";;
  *)
    name_mangling="FORTRAN_MANGLING_UNKNOWN"
    AC_MSG_WARN([UNKNOWN FORTRAN NAME-MANGLING SCHEME]);;
esac

RR_PATH_NETCDF
RR_PATH_NETCDF_F90

test "$osname"=NONE && osname=`uname -s`
OSNAME=$osname
AC_SUBST(OSNAME)

test "$machinename"=NONE && machinename=`uname -m`
fullhostname=`hostname -f`

# DEFINE VARIABLES ACCORDING TO OS AND COMPILER

if echo $osname | grep -i aix >/dev/null 2>&1; then
   SYSDEF="AIX"
fi
if echo $osname | grep -i darwin >/dev/null 2>&1; then
   SYSDEF="DARWIN"
fi
if echo $osname | grep -i unix_system_v >/dev/null 2>&1; then
   SYSDEF="UNIXSYSTEMV"
fi
if echo $osname | grep -i irix >/dev/null 2>&1; then
   SYSDEF="IRIX"
fi
if echo $osname | grep -i irix64 >/dev/null 2>&1; then
   SYSDEF="IRIX64"
fi
if echo $osname | grep -i linux >/dev/null 2>&1; then
   SYSDEF="LINUX"
fi
if echo $osname | grep -i osf1 >/dev/null 2>&1; then
   SYSDEF="OSF1"
fi
if echo $osname | grep -i super >/dev/null 2>&1; then
   SYSDEF="SUPERUX"
fi
if echo $osname | grep -i sun >/dev/null 2>&1; then
   SYSDEF="SUNOS"
fi
if echo $osname | grep -i t3e >/dev/null 2>&1; then
   SYSDEF="T3E"
fi
if echo $osname | grep -i unicos >/dev/null 2>&1; then
   SYSDEF="UNICOS"
fi
if test -z "$SYSDEF"; then
   AC_MSG_WARN([OPERATING SYSTEM UNKNOWN])
   SYSDEF="UNKNOWNOS"
fi

AM_CONDITIONAL(OSNAME_AIX, [test x"$osname" = xAIX])
AM_CONDITIONAL(OSNAME_LINUX, [test x"$osname" = xLinux])

AM_CONDITIONAL(COMPILER_IBM, [test x"$FC" = xxlf90])
AM_CONDITIONAL(COMPILER_GNU, [test x"$FC" = xgfortran])
AM_CONDITIONAL(COMPILER_G95, [test x"$FC" = xg95])
AM_CONDITIONAL(COMPILER_SUN, [test x"$FC" = xsunf95])
AM_CONDITIONAL(COMPILER_PGI, [test x"$FC" = xpgf90])
AM_CONDITIONAL(COMPILER_INTEL, [test x"$FC" = xifort])

AM_CONDITIONAL(REGCM_PARALLEL, [test x$PARALLEL = x1])

# Force Big endian I/O for known compilers
if eval "test x"$FC" = xgfortran"; then
  FCFLAGS="-fconvert=big-endian $FCFLAGS"
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
fi
if eval "test x"$FC" = xifort"; then
  FCFLAGS="-convert big_endian -assume byterecl -heap-arrays $FCFLAGS"
  CPRDEF="INTEL"
  AC_SUBST(FCFLAGS)
fi
if eval "test x"$FC" = xg95"; then
  FCFLAGS="-fendian=big -fmultiple-save $FCFLAGS"
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
fi
if eval "test x"$FC" = xpgf90"; then
  FCFLAGS="-byteswapio $FCFLAGS"
  CPRDEF="PGI"
  AC_SUBST(FCFLAGS)
fi
if eval "test x"$FC" = xxlf90"; then
  CPPFLAGS="-DIBM $CPPFLAGS"
  CPRDEF="XLF"
  AC_SUBST(CPPFLAGS)
fi
if eval "test x"$FC" = xsun95"; then
  CPRDEF="WORKSHOP"
fi

# Set debug/production flags
if eval "test x$DEBUG = x1"; then
  if eval "test x"$FC" = xg95"; then
    FCFLAGS="-O0 -g -Wall -Wextra -fbounds-check -ftrace=full $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xgfortran"; then
    FCFLAGS="-O0 -g -ffpe-trap=zero,invalid,overflow,underflow -fbacktrace \
            -fdump-core -Wall -pedantic $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xifort"; then
    FCFLAGS="-O0 -g \
             -check arg_temp_created,bounds,format,output_conversion,uninit \
             -traceback -warn all -fp-model strict $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xxlf90"; then
    FCFLAGS="-O0 -qsuppress=1501-245 -bnoquiet -qcheck -qffltrap=nanq,overflow \
             -qmaxmem=-1 -qzerosize $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xsunf95"; then
    FCFLAGS="-g -O0 -C -w3 -ansi -xfilebyteorder=big16:%all \
             -D__DATE__="'`date +%D`'" -D__TIME__="'`date +%T`'" $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xpgf90"; then
    FCFLAGS="-Mbounds -Mcache_align -Mchkfpstk -Mchkstk -Mdalign -Mdclchk \
             -Mdepchk -Miomutex -Mrecursive -Msave -Ktrap=fp -O0 -g $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
else
  if eval "test x"$FC" = xg95"; then
    FCFLAGS="-O3 $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xgfortran"; then
    FCFLAGS="-O3 $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xifort"; then
    FCFLAGS="-O3 -fp-model precise $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xxlf90"; then
    FCFLAGS="-qsuppress=1501-245 -bnoquiet -O3 -qstrict \
             -qmaxmem=-1 -qzerosize $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xsunf95"; then
    FCFLAGS="-O3 -xfilebyteorder=big16:%all \
             -D__DATE__="'`date +%D`'" -D__TIME__="'`date +%T`'" $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$FC" = xpgf90"; then
    FCFLAGS="-fastsse $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
fi

AM_CONDITIONAL(MPISERIAL, [test x"$PARALLEL" = x1])
if eval "test x$PARALLEL = x1"; then
  ACX_MPI
  FC=$MPIFC
else
  MPIF90=$FC
  MPIFC=$FC
  FORT_SIZE="real4double8"
  MPIPATH="${PWD}/external"
  CPPFLAGS="-I$MPIPATH $CPPFLAGS"
  LDFLAGS="-L$MPIPATH $LDFLAGS"
  LIBS="$LIBS -lmpiser"
  AC_SUBST(CPPFLAGS)
  AC_SUBST(LDFLAGS)
fi

MPEUCFLAGS="-DSYS$SYSDEF -DCPR$CPRDEF -D$name_mangling"
MPEUFCFLAGS="-DSYS$SYSDEF -DCPR$CPRDEF"
if test -n "$SRKDEF"; then
  MPEUCFLAGS="$MPEUCFLAGS -D$SRKDEF"
  MPEUFCFLAGS="$MPEUFCFLAGS -D$SRKDEF"
fi
AC_SUBST(MPEUCFLAGS)
AC_SUBST(MPEUFCFLAGS)

F90FLAGS="$FCFLAGS"
AC_SUBST(F90FLAGS)

REGCM_ROOT=${PWD}
AC_SUBST(REGCM_ROOT)
bindir=${REGCM_ROOT}/Bin
AC_SUBST(bindir)

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile PreProc/Terrain/Makefile PreProc/ICBC/Makefile \
                 PreProc/CLM/Makefile PostProc/Makefile Main/netlib/Makefile \
                 Main/clmlib/Makefile Main/Makefile external/Makefile])

AC_OUTPUT
