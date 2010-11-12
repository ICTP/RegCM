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
            [Path to NetCDF installation (default: NETCDF environment)],
  [NC_PREFIX="$withval"], [NC_PREFIX="${NETCDF}"]
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

RR_SVN_DEF
RR_PATH_NETCDF
RR_PATH_NETCDF_F90

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile PreProc/Terrain/Makefile PreProc/ICBC/Makefile \
                 PreProc/CLM/Makefile PostProc/Makefile Main/netlib/Makefile \
                 Main/clmlib/Makefile Main/Makefile])

test "$osname"=NONE && osname=`uname -s`
test "$machinename"=NONE && machinename=`uname -m`
fullhostname=`hostname -f`

AM_CONDITIONAL(OSNAME_AIX, [test x"$osname" = xAIX])
AM_CONDITIONAL(COMPILER_GFORTRAN, [test x"$FC" = xgfortran])
AM_CONDITIONAL(COMPILER_SUN, [test x"$FC" = xsunf95])
AM_CONDITIONAL(REGCM_PARALLEL, [test x$PARALLEL = x1])

if eval "test x$PARALLEL = x1"; then
  ACX_MPI
  FC=$MPIFC
fi

REGCM_ROOT=${PWD}
AC_SUBST(REGCM_ROOT)
bindir=${REGCM_ROOT}/Bin
AC_SUBST(bindir)

AC_OUTPUT
