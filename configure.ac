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

AC_INIT([RegCM], [4.5], [ggiulian@ictp.it])

AC_SUBST([PACKAGE_URL], ["https://gforge.ictp.it/gf/project/regcm"])

#AC_PREREQ([2.64])

LT_INIT

AM_INIT_AUTOMAKE([foreign])

AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_FC

AC_LANG([Fortran])

AX_PROG_NF_CONFIG([NC_PREFIX=`nf-config --prefix`; NFAUTO="yes"],
                  [NC_PREFIX="${NETCDF}"; NFAUTO="no"])

AX_PROG_NC_CONFIG([NC_PREFIX=`nc-config --prefix`; NCAUTO="yes"],
                  [NC_PREFIX="${NETCDF}"; NCAUTO="no"])

AC_ARG_WITH([netcdf],
  AC_HELP_STRING([--with-netcdf],
            [Path to NetCDF installation (default: NETCDF environment)]),
  [NC_PREFIX="$withval"; NFAUTO="no"; NCAUTO="no"], [FAKE="no"]
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

AC_ARG_WITH([oasis-include],
  AC_HELP_STRING([--with-oasis-include],
                 [Path to OASIS installation (default: do not use)]),
  [USEOASIS="1"; OASIS_INCLUDE="$withval"], [USEOASIS="0"]
)

PARALLEL="1"
AC_ARG_ENABLE(mpiserial,
  AC_HELP_STRING([--enable-mpiserial],
                 [Use the included MPI replacement library for single processor]),
  [DONOTCHECKMPI="DONOTCHECKMPI"; PARALLEL="0"]
)

SP_REAL="0"
AC_ARG_ENABLE(singleprecision,
  AC_HELP_STRING([--enable-singleprecision],
                 [Default real type is real4 single precision.]),
                 [SP_REAL="1"]
)
if eval "test x"$SP_REAL" = x1"; then
  AM_CPPFLAGS="\$(DEFINE)SINGLE_PRECISION_REAL $AM_CPPFLAGS"
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

USECLM="0"
USEMEGAN="0"
AC_ARG_ENABLE(clm,
  AC_HELP_STRING([--enable-clm],
                 [Supply this option if you plan on using CLM option.]),
  [USECLM="1"],
)
AM_CONDITIONAL(DO_COMPILE_LIBCLM, [test x"$USECLM" = x1])

AC_ARG_ENABLE(megan,
  AC_HELP_STRING([--enable-megan],
                 [Supply this option if you plan on using CLM option with MEGAN.]),
  [USEMEGAN="1"],
)

if eval "test x"$USECLM" = x1"; then
  AC_MSG_WARN([Compiling CLM enabled code])
  AM_CPPFLAGS="\$(DEFINE)CLM $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
  if eval "test x"$USEMEGAN" = x1"; then
    AC_MSG_WARN([Compiling CLM with MEGAN enabled code])
    AM_CPPFLAGS="\$(DEFINE)VOC $AM_CPPFLAGS"
    AC_SUBST(AM_CPPFLAGS)
  fi
fi

USECLM45="0"
AC_ARG_ENABLE(clm45,
  AC_HELP_STRING([--enable-clm45],
                 [Supply this option if you plan on using CLM45 option.]),
  [USECLM45="1"],
)
AM_CONDITIONAL(DO_COMPILE_LIBCLM45, [test x"$USECLM45" = x1])

if eval "test x"$USECLM45" = x1"; then
  AC_MSG_WARN([Compiling CLM45 enabled code])
  AM_CPPFLAGS="\$(DEFINE)CLM45 $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

USECPL="0"
AC_ARG_ENABLE(cpl,
  AC_HELP_STRING([--enable-cpl],
                 [Supply this option if you plan on using REGESM.]),
  [USECPL="1"],
)
AM_CONDITIONAL(DO_COMPILE_CPL, [test x"$USECPL" = x1])

if eval "test x"$USECPL" = x1"; then
  AM_CPPFLAGS="\$(DEFINE)CPL $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

if eval "test x"$USEOASIS" = x1"; then
  AM_CPPFLAGS="\$(DEFINE)OASIS -I$OASIS_INCLUDE $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

AC_ARG_ENABLE(nc4-gzip,
  AC_HELP_STRING([--enable-nc4-gzip],
                 [Supply this option if you plan on using NetCDF4/HDF5 output compression option.]),
  [NETCDF4_COMPRESS="1"],[NETCDF4_COMPRESS="0"]
)
if eval "test x$NETCDF4_COMPRESS = x1"; then
  AM_CPPFLAGS="\$(DEFINE)NETCDF4_COMPRESS $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

AC_ARG_ENABLE(nc4,
  AC_HELP_STRING([--enable-nc4],
                 [Supply this option if you plan on using NetCDF4/HDF5 option.]),
  [NETCDF4_HDF5="1"],[NETCDF4_HDF5="0"]
)
if eval "test x$NETCDF4_HDF5 = x1"; then
  AM_CPPFLAGS="\$(DEFINE)NETCDF4_HDF5 $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

DEBUG="0"
AC_ARG_ENABLE(debug,
  AC_HELP_STRING([--enable-debug],
                 [Enable debugging flags and per processor log file]),
  [DEBUG="1"]
)
if eval "test x$DEBUG = x1"; then
  AM_CPPFLAGS="\$(DEFINE)DEBUG $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

NC_INCLUDES=""
NC_LDFLAGS=""
NC_LIBS=""

if eval "test x$NFAUTO = xyes"; then
  if eval "test x$NCAUTO = xyes"; then
    NC_INCLUDES="`nc-config --cflags` `nf-config --fflags`"
    NC_LDFLAGS=""
    NC_LIBS="`nc-config --libs` `nf-config --flibs`"
  else
    NC_INCLUDES="`nf-config --fflags`"
    NC_LDFLAGS=""
    NC_LIBS="`nf-config --flibs`"
  fi
else
  if eval "test x$NCAUTO = xno"; then
    NC_INCLUDES="-I$NC_PREFIX/include"
    NC_LDFLAGS="-L$NC_PREFIX/lib"
    NC_LIBS="-lnetcdf"
  else
    NCVER=`nc-config --version | cut -d " " -f 2`
    if eval "test x$NCVER = x4.0.1"; then
      NC_INCLUDES="`nc-config --fflags`"
    else
      # new netcdf library has this
      NC_INCLUDES="-I`nc-config --includedir`"
    fi
    NC_LDFLAGS=""
    NC_LIBS="`nc-config --libs`"
  fi
fi

AC_SUBST(NC_PREFIX)
AC_SUBST(NC_INCLUDES)
AC_SUBST(NC_LIBS)
AC_SUBST(NC_LDFLAGS)

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
  "lower case, no underscore, no extra underscore")
    name_mangling="FORTRAN_SAME";;
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

AM_CONDITIONAL(OSNAME_AIX, [ test -z "${osname##AIX}" ] )
AM_CONDITIONAL(OSNAME_LINUX, [ test -z "${osname##Linux}" ] )

# Cut away flags for this test
FCPROG=`echo $FC | cut -d " " -f 1`
AC_SUBST(FCPROG)

AM_CONDITIONAL(COMPILER_IBM, [ test -z "${FCPROG##xxlf200*}" -o \
                                    -z "${FCPROG##xbgxlf200*}" ] )
AM_CONDITIONAL(COMPILER_CRAY, [ test -z "${FCPROG##xftn*}" ] )
AM_CONDITIONAL(COMPILER_GNU, [ test -z "${FCPROG##gfortran*}" -o \
                                    -z "${FCPROG##\/gfortran*}" -o \
                                    -z "${FCPROG##scorep}" ] )
AM_CONDITIONAL(COMPILER_G95, [ test -z "${FCPROG##g95*}" ] )
AM_CONDITIONAL(COMPILER_SUN, [ test -z "${FCPROG##sunf9*}" ] )
AM_CONDITIONAL(COMPILER_PGI, [ test -z "${FCPROG##pgf9*}" -o \
                                    -z "${FCPROG##pgfortran*}" ] )
AM_CONDITIONAL(COMPILER_INTEL, [ test -z "${FCPROG##ifort*}" -o \
                                      -z "${FCPROG##impif*}" ] )
AM_CONDITIONAL(REGCM_PARALLEL, [ test x"$PARALLEL" = "x1"])

# Force Big endian I/O for known compilers
if test -z "$COMPILER_GNU_TRUE"; then
  # Following check seems not to work on SuSE, need investigation
#  GFORTRAN_VERSION="`gfortran --version | head -1 | cut -d " " -f 4`"
#  GFORTRAN_MAJOR="`echo $GFORTRAN_VERSION | cut -d "." -f 1`"
#  GFORTRAN_MINOR="`echo $GFORTRAN_VERSION | cut -d "." -f 2`"
#  echo "gfortran is version $GFORTRAN_VERSION"
#  if eval "test "$GFORTRAN_MAJOR" -lt 4"; then
#    AC_MSG_WARN([We need at least gfortran 4.6.x.....])
#    exit 1
#  else
#    if eval "test "$GFORTRAN_MINOR" -lt 6"; then
#      AC_MSG_WARN([We need at least gfortran 4.6.x.....])
#      exit 1
#    fi
#  fi
  FCFLAGS="-fconvert=big-endian -fno-range-check -fPIC"
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_INTEL_TRUE"; then
  FCFLAGS="-convert big_endian -assume byterecl \
           -heap-arrays -fPIC"
  CPRDEF="INTEL"
  AC_SUBST(FCFLAGS)
  LIBS=`echo $LIBS | sed -e 's/-Wl,--as-needed//g'`
  AC_SUBST(LIBS)
fi
if test -z "$COMPILER_G95_TRUE"; then
  FCFLAGS="-fendian=big -fmultiple-save"
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_PGI_TRUE"; then
  FCFLAGS="-byteswapio"
  CPRDEF="PGI"
  AC_SUBST(FCFLAGS)
fi
if test -z "$COMPILER_IBM_TRUE"; then
  FCFLAGS="-WF,-DIBM"
  CPRDEF="XLF"
  AC_SUBST(FCFLAGS)
fi
if test -x "$COMPILER_SUN_TRUE"; then
  CPRDEF="WORKSHOP"
fi
if test -z "$COMPILER_CRAY_TRUE"; then
  FCFLAGS=""
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
fi

AC_SUBST(CPRDEF)

# Set debug/production flags
if eval "test x$DEBUG = x1"; then
  if test -z "$COMPILER_G95_TRUE"; then
    FCFLAGS="-O0 -g -Wall -Wextra -fbounds-check -ftrace=full $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if test -z "$COMPILER_GNU_TRUE"; then
    FCFLAGS="-O0 -g -pg -ffpe-trap=zero,invalid,overflow,underflow -fbacktrace \
            -fdump-core -Wall -pedantic -Wconversion -Wsurprising -Waliasing \
            -Wampersand -Warray-bounds -Wcharacter-truncation -Wuninitialized \
            -Wline-truncation -Wintrinsics-std -Wno-tabs -Wunderflow \
            -Wunused-parameter -Wno-align-commons $FCFLAGS \
            -fno-range-check -Wimplicit-interface -Wimplicit-procedure \
            -Wintrinsics-std -Wintrinsic-shadow \
            -fcheck=array-temps,bounds,do,mem,recursion"
    AC_SUBST(FCFLAGS)
  fi
  if test -z "$COMPILER_INTEL_TRUE"; then
    FCFLAGS="-O0 -g \
             -check arg_temp_created,bounds,format,output_conversion,uninit \
             -traceback -warn all -fp-model source -ftz \
             -fimf-arch-consistency=true $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if test -z "$COMPILER_IBM_TRUE"; then
    FCFLAGS="-O0 -qsuppress=1501-245 -qcheck -qflttrap=nanq,overflow \
             -qxlf2003=polymorphic -qmaxmem=-1 -qzerosize $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if test -z "$COMPILER_PGI_TRUE"; then
    FCFLAGS="-Mbounds -Mcache_align -Mchkfpstk -Mchkstk -Mdalign -Mdclchk \
             -Mdepchk -Miomutex -Mrecursive -Msave -Ktrap=denorm -O0 -g \
             -traceback $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if test -z "$COMPILER_CRAY_TRUE"; then
    FCFLAGS="-g $(FCFLAGS)"
    AC_SUBST(FCFLAGS)
    CPRDEF="GNU"
  fi
else
  if test -z "$COMPILER_G95_TRUE"; then
    FCFLAGS="-O3 $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if test -z "$COMPILER_GNU_TRUE"; then
    if test -z "${osname##Linux}"; then
      FCFLAGS="-O3 -finline-functions -mtune=native -march=native $FCFLAGS"
    else
      # MacOSX may choke on using AVX extensions...
      FCFLAGS="-O3 $FCFLAGS"
    fi
    AC_SUBST(FCFLAGS)
  fi
  if test -z "$COMPILER_INTEL_TRUE"; then
    AR=xiar
    FCFLAGS="-O3 -finline-functions -fltconsistency -mieee-fp \
             -traceback $FCFLAGS"
    AC_SUBST(FCFLAGS)
    AC_SUBST(AR)
  fi
  if test -z "$COMPILER_IBM_TRUE"; then
    FCFLAGS="-qsuppress=1501-245 -O3 -qstrict -qxlf2003=polymorphic \
             -qmaxmem=-1 -qzerosize $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if test -z "$COMPILER_SUN_TRUE"; then
    FCFLAGS="-O3 -xfilebyteorder=big16:%all \
             -D__DATE__="'`date +%D`'" -D__TIME__="'`date +%T`'" $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if test -z "$COMPILER_PGI_TRUE"; then
    FCFLAGS="-fastsse $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  if test -z "$COMPILER_CRAY_TRUE"; then
    FCFLAGS="$(FCFLAGS)"
    AC_SUBST(FCFLAGS)
    CPRDEF="GNU"
  fi
fi

MCTPATH="${PWD}/external"
if test -z "$COMPILER_SUN_TRUE"; then
  MCTFLAGS="-M$MCTPATH"
  AM_CPPFLAGS="$MCTFLAGS $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
else
  MCTFLAGS="-I$MCTPATH"
  AM_CPPFLAGS="$MCTFLAGS $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi
AM_LDFLAGS="-L$MCTPATH $AM_LDFLAGS"
AC_SUBST(AM_LDFLAGS)
if eval "test x"$PARALLEL" = x1"; then
  ACX_MPI
else
  FORT_SIZE="real4double8"
  AC_SUBST(FORT_SIZE)
  MPIFC=$FC
  AC_SUBST(MPIFC)
  LIBS="$LIBS -lmpiser"
  AC_SUBST(LIBS)
  AM_CPPFLAGS="\$(DEFINE)MPI_SERIAL $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

MPEUCFLAGS="\$(DEFINE)SYS$SYSDEF \$(DEFINE)CPR$CPRDEF \$(DEFINE)$name_mangling"
MPEUFCFLAGS="\$(DEFINE)SYS$SYSDEF \$(DEFINE)CPR$CPRDEF"
if test -n "$SRKDEF"; then
  MPEUCFLAGS="$MPEUCFLAGS \$(DEFINE)$SRKDEF"
  MPEUFCFLAGS="$MPEUFCFLAGS \$(DEFINE)$SRKDEF"
fi
AC_SUBST(MPEUCFLAGS)
AC_SUBST(MPEUFCFLAGS)

F90FLAGS="$FCFLAGS"
AC_SUBST(F90FLAGS)

REGCM_ROOT=${PWD}
AC_SUBST(REGCM_ROOT)
AC_PREFIX_DEFAULT(${PWD})

if test -z "$COMPILER_SUN_TRUE"; then
  CPPFLAGS="$AM_CPPFLAGS -M$REGCM_ROOT/Share"
  AC_SUBST(CPPFLAGS)
else
  CPPFLAGS="$AM_CPPFLAGS -I$REGCM_ROOT/Share"
  AC_SUBST(CPPFLAGS)
fi
LIBS="-L$REGCM_ROOT/Share -lrcmlib $LIBS"
SYSLIBS="$LIBS"
AC_SUBST(SYSLIBS)

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile external/Makefile Share/Makefile PreProc/Makefile \
                 PreProc/Terrain/Makefile PreProc/ICBC/Makefile \
                 PreProc/CH_EMISSIONS/Makefile \
                 PreProc/CH_EMISSIONS/interpolation/Makefile \
                 PreProc/CH_EMISSIONS/interpolation/tools/Makefile \
                 PreProc/CH_EMISSIONS/interpolation/tools/REGCM_grid/Makefile \
                 PreProc/CLM/Makefile PreProc/CLM45/Makefile PostProc/Makefile \
                 Main/mpplib/Makefile Main/netlib/Makefile \
                 Main/Makefile Main/batslib/Makefile Main/clmlib/Makefile \
                 Main/ocnlib/Makefile Main/clmlib/clm4.5/Makefile \
                 Main/chemlib/Makefile Main/chemlib/GAS_CBMZ_NEW/Makefile \
                 Main/cumlib/Makefile Main/radlib/Makefile \
                 Main/radlib/RRTMG_SW/Makefile Main/pbllib/Makefile \
                 Tools/Programs/CheckSun/Makefile \
                 Tools/Programs/RegCM_read/Makefile \
                 Tools/Programs/timeseries/Makefile])

AC_OUTPUT
