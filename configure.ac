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

AC_INIT([RegCM], [5.0.0], [ggiulian@ictp.it])

AC_SUBST([PACKAGE_URL], ["https://github.com/ICTP/RegCM"])

#AC_PREREQ([2.64])

LT_INIT

AM_INIT_AUTOMAKE([foreign])

AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_FC

AC_LANG([Fortran])

AC_ARG_WITH([netcdf],
  AC_HELP_STRING([--with-netcdf],
            [Path to NetCDF installation (default: NETCDF environment)]),
  [NC_PREFIX="$withval"; NF_PREFIX="$withval"; NFAUTO="no"; NCAUTO="no"],
  [FAKE="no"]
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

AC_ARG_WITH([oasis-path],
  AC_HELP_STRING([--with-oasis-path],
		 [Path to OASIS3-MCT installation (require --enable-oasis) [default: ~/oasis3-mct/compile_oa3-mct]]),
  [OASIS_PREFIX="$withval"], [OASIS_PREFIX="$HOME/oasis3-mct/compile_oa3-mct"]
)

AC_ARG_WITH([oasis-comm],
  AC_HELP_STRING([--with-oasis-comm],
		 [Communication technique used in OASIS3-MCT (require --enable-oasis) [default: MPI1]]),
  [OASIS_CHAN="$withval"], [OASIS_CHAN="MPI1"]
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

USE_KNL="0"
AC_ARG_ENABLE(knl,
  AC_HELP_STRING([--enable-knl],
                 [Compile with Xeon Phi flags.]),
                 [USE_KNL="1"]
)
if eval "test x"$USE_KNL" = x1"; then
  AM_CPPFLAGS="\$(DEFINE)KNL $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

USE_ZEN3="0"
AC_ARG_ENABLE(zen3,
  AC_HELP_STRING([--enable-zen3],
                 [Compile with Zen3 processor optimization]),
                 [USE_ZEN3="1"]
)
if eval "test x"$USE_ZEN3" = x1"; then
  AM_CPPFLAGS="\$(DEFINE)ZEN3 $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

USE_BDW="0"
AC_ARG_ENABLE(bdw,
  AC_HELP_STRING([--enable-bdw],
                 [Compile with AVX2 flags.]),
                 [USE_BDW="1"]
)
if eval "test x"$USE_BDW" = x1"; then
  AM_CPPFLAGS="\$(DEFINE)BDW $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

USE_NHL="0"
AC_ARG_ENABLE(nhl,
  AC_HELP_STRING([--enable-nhl],
                 [Compile with SSE 4.2 flags.]),
                 [USE_NHL="1"]
)
if eval "test x"$USE_NHL" = x1"; then
  AM_CPPFLAGS="\$(DEFINE)NHL $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

USE_SKL="0"
AC_ARG_ENABLE(skl,
  AC_HELP_STRING([--enable-skl],
                 [Compile with Sky Lake flags.]),
                 [USE_SKL="1"]
)
if eval "test x"$USE_SKL" = x1"; then
  AM_CPPFLAGS="\$(DEFINE)SKL $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

USE_OPENACC="0"
# TODO: make this configurable
PGI_GPU_ARCH="cc80"

AC_ARG_ENABLE(openacc,
  AC_HELP_STRING([--enable-openacc],
                 [Compile OpenACC directives (only for PGI fortran compilation)]),
                 [USE_OPENACC="1"]
)
if eval "test x"$USE_OPENACC" = x1"; then
  AM_CPPFLAGS="\$(DEFINE)OPENACC $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

AC_ARG_ENABLE(openacc-debug,
  AC_HELP_STRING([--enable-openacc-debug],
                 [Compile OpenACC directives with autocompare (only for PGI fortran compilation)]),
                 [USE_OPENACC="2"]
)
if eval "test x"$USE_OPENACC" = x2"; then
  AM_CPPFLAGS="\$(DEFINE)OPENACC $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

AC_ARG_ENABLE(openacc-managed,
  AC_HELP_STRING([--enable-openacc-managed],
                 [Compile OpenACC directives with managed memory (only for PGI fortran compilation)]),
                 [USE_OPENACC="3"]
)
if eval "test x"$USE_OPENACC" = x3"; then
  AM_CPPFLAGS="\$(DEFINE)OPENACC $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

AC_ARG_ENABLE(openacc-stdpar,
  AC_HELP_STRING([--enable-openacc-stdpar],
                 [Compile OpenACC directives with stdpar (only for PGI fortran compilation)]),
                 [USE_OPENACC="4"]
)
if eval "test x"$USE_OPENACC" = x4"; then
  AM_CPPFLAGS="\$(DEFINE)OPENACC $AM_CPPFLAGS -DSTDPAR"
  AC_SUBST(AM_CPPFLAGS)
fi

USE_PROF="0"
AC_ARG_ENABLE(profile,
  AC_HELP_STRING([--enable-profile],
                 [Compile with profile flags.]),
                 [USE_PROF="1"]
)
if eval "test x"$USE_PROF" = x1"; then
  AM_CPPFLAGS="\$(DEFINE)PROF $AM_CPPFLAGS"
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

USEOASIS="0"
AC_ARG_ENABLE(oasis,
  AC_HELP_STRING([--enable-oasis],
                 [Supply this option if you plan on using OASIS option.]),
  [USEOASIS="1"]
)
AM_CONDITIONAL(DO_COMPILE_OASIS, [test x"$USEOASIS" = x1])

if eval "test x"$USEOASIS" = x1"; then
  AC_MSG_WARN([Compiling OASIS enabled code])
  AC_MSG_WARN([OASIS path: $OASIS_PREFIX])
  AC_MSG_WARN([OASIS communication technique: $OASIS_CHAN])
  AM_CPPFLAGS="\$(DEFINE)OASIS $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

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

PNETCDF_IN_NETCDF=0
AC_ARG_ENABLE(parallel-nc,
  AC_HELP_STRING([--enable-parallel-nc],
                 [Supply this option if you plan on using NetCDF4 parallel option based upon Parallel Netcdf.]),
  [PNETCDF_IN_NETCDF="1"],[PNETCDF_IN_NETCDF="0"]
)
if eval "test x$PNETCDF_IN_NETCDF = x1"; then
  AM_CPPFLAGS="\$(DEFINE)PNETCDF_IN_NETCDF $AM_CPPFLAGS"
  AC_SUBST(AM_CPPFLAGS)
fi

PNETCDF="0"
AC_ARG_ENABLE(pnetcdf,
  AC_HELP_STRING([--enable-pnetcdf],
                 [Use Parallel NetCDF from NU-ANL instead of netcdf.]),
  [PNETCDF="1"],[PNETCDF="0"]
)
if eval "test x$PNETCDF = x1"; then
  NULLOPTION=$( echo $AM_CPPFLAGS | sed -e 's/\$(DEFINE)PNETCDF_IN_NETCDF//')
  AM_CPPFLAGS="\$(DEFINE)PNETCDF $NULLOPTION"
  AC_SUBST(AM_CPPFLAGS)
fi

AC_ARG_ENABLE(nc4-filters,
  AC_HELP_STRING([--enable-nc4-filters],
		 [Supply this option if you plan on using NetCDF4/HDF5 output filtering option.]),
  [NETCDF4_COMPRESS="1"],[NETCDF4_COMPRESS="0"]
  )
if eval "test x$NETCDF4_COMPRESS = x1"; then
  NULLOPTION=$( echo $AM_CPPFLAGS | sed -e 's/\$(DEFINE)PNETCDF_IN_NETCDF//' \
                -e 's/\$(DEFINE)PNETCDF//')
  AM_CPPFLAGS="\$(DEFINE)NETCDF4_COMPRESS $NULLOPTION"
  AC_SUBST(AM_CPPFLAGS)
fi

NETCDF_CDF5=0
AC_ARG_ENABLE(cdf5,
  AC_HELP_STRING([--enable-cdf5],
                 [Supply this option if you plan on using NetCDF/CDF5 option.]),
  [NETCDF_CDF5="1"],[NETCDF_CDF5="0"]
)

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

AX_PROG_NF_CONFIG([NF_PREFIX=`nf-config --prefix`; NFAUTO="yes"],
                  [NF_PREFIX="${NETCDF}"; NFAUTO="no"])

AX_PROG_NC_CONFIG([NC_PREFIX=`nc-config --prefix`; NCAUTO="yes"],
                  [NC_PREFIX="${NETCDF}"; NCAUTO="no"])

if eval "test x$NFAUTO = xyes"; then
  if eval "test x$NCAUTO = xyes"; then
    NC_INCLUDES="`nc-config --cflags` `nf-config --fflags`"
    NC_LDFLAGS="-L`nf-config --prefix`/lib -L`nc-config --prefix`/lib"
    NC_LIBS="`nc-config --libs` `nf-config --flibs`"
  else
    NC_INCLUDES="-I$NC_PREFIX/include `nf-config --fflags`"
    NC_LDFLAGS="-L$NC_PREFIX/lib -L`nc-config --prefix`/lib"
    NC_LIBS="-L$NC_PREFIX/lib `nf-config --flibs`"
  fi
else
  if eval "test x$NCAUTO = xno"; then
    NC_INCLUDES="-I$NC_PREFIX/include -I$NF_PREFIX/include"
    NC_LDFLAGS="-L$NC_PREFIX/lib -L$NF_PREFIX/lib"
    NC_LIBS="-lnetcdff -lnetcdf"
  else
    NCVER=`nc-config --version | cut -d " " -f 2`
    if eval "test x$NCVER = x4.0.1"; then
      NC_INCLUDES="`nc-config --fflags`"
    else
      # new netcdf library has this
      NC_INCLUDES="-I`nc-config --includedir`"
    fi
    NC_LDFLAGS="-L`nc-config --prefix`/lib"
    NC_LIBS="`nc-config --libs`"
  fi
fi

AC_SUBST(NC_PREFIX)
AC_SUBST(NF_PREFIX)
AC_SUBST(NC_INCLUDES)
AC_SUBST(NC_LDFLAGS)
AC_SUBST(NC_LIBS)

if eval "test x$PNETCDF = x1"; then
  AX_PROG_PNETCDF_CONFIG([NCP_PREFIX=`pnetcdf-config --prefix`; NCPAUTO="yes"],
                         [NCP_PREFIX="${NETCDF}"; NCPAUTO="no"])
  if eval "test x$NCPAUTO = xyes"; then
    CC=`pnetcdf-config --cc`
    AC_SUBST(CC)
    NCVER=`pnetcdf-config --version | cut -d " " -f 2`
    NC_INCLUDES="$NC_INCLUDES -I`pnetcdf-config --includedir`"
    NC_LDFLAGS="$NC_LDFLAGS -L`pnetcdf-config --libdir`"
    NC_LIBS="$NC_LIBS -L`pnetcdf-config --libdir` -lpnetcdf"
  else
    NC_INCLUDES="$NC_INCLUDES -I$NCP_PREFIX/include"
    NC_LDFLAGS="$NC_LDFLAGS -L$NCP_PREFIX/lib"
    NC_LIBS="$NC_LIBS -L$NCP_PREFIX/lib -lpnetcdf"
  fi
fi

AC_SUBST(NC_INCLUDES)
AC_SUBST(NC_LDFLAGS)
AC_SUBST(NC_LIBS)

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
FCPROG=`basename $(echo $FC | cut -d " " -f 1)`
AC_SUBST(FCPROG)

AM_CONDITIONAL(COMPILER_IBM, [ test -z "${FCPROG##xlf200*}" -o \
                                    -z "${FCPROG##bgxlf200*}" ] )
AM_CONDITIONAL(COMPILER_CRAY, [ test -z "${FCPROG##ftn*}" ] )
AM_CONDITIONAL(COMPILER_GNU, [ test -z "${FCPROG##gfortran*}" -o \
                                    -z "${FCPROG##\/gfortran*}" -o \
                                    -z "${FCPROG##scorep}" ] )
AM_CONDITIONAL(COMPILER_G95, [ test -z "${FCPROG##g95*}" ] )
AM_CONDITIONAL(COMPILER_SUN, [ test -z "${FCPROG##sunf9*}" ] )
AM_CONDITIONAL(COMPILER_PGI, [ test -z "${FCPROG##pgf9*}" -o \
                                    -z "${FCPROG##nvf9*}" -o \
                                    -z "${FCPROG##pgfortran*}" -o \
                                    -z "${FCPROG##nvfortran*}" ] )
AM_CONDITIONAL(COMPILER_INTEL, [ test -z "${FCPROG##ifort*}" -o \
                                      -z "${FCPROG##impif*}" -o \
                                      -z "${FCPROG##mpiifort*}" ] )
AM_CONDITIONAL(COMPILER_INTEL_LLVM, [ test -z "${FCPROG##ifx*}" -o \
                                           -z "${FCPROG##mpiifx*}" ] )
AM_CONDITIONAL(REGCM_PARALLEL, [ test x"$PARALLEL" = "x1"])

DEFINE="-D"
AC_SUBST(DEFINE)

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
  FCFLAGS="$FCFLAGS -fconvert=big-endian -fno-range-check -fPIC"
  OMPFLAGS="-fopenmp"
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
  AC_SUBST(OMPFLAGS)
  if eval "test x"$USE_KNL" = x1"; then
    FCFLAGS="$FCFLAGS -march=knl"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$USE_ZEN3" = x1"; then
    FCFLAGS="$FCFLAGS -march=znver3"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$USE_SKL" = x1"; then
    FCFLAGS="$FCFLAGS -march=skylake-avx512"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$USE_BDW" = x1"; then
    FCFLAGS="$FCFLAGS -march=broadwell"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$USE_NHL" = x1"; then
    FCFLAGS="$FCFLAGS -march=nehalem"
    AC_SUBST(FCFLAGS)
  fi
fi
if test -z "$COMPILER_INTEL_LLVM_TRUE"; then
  OMPFLAGS="-qopenmp -qopt-report"
  FCFLAGS="$FCFLAGS -convert big_endian -assume byterecl \
           -heap-arrays -fPIC -DSTDPAR $OMPFLAGS"
  CPRDEF="INTEL"
  AC_SUBST(FCFLAGS)
  AC_SUBST(OMPFLAGS)
  LIBS=`echo $LIBS | sed -e 's/-Wl,--as-needed//g'`
  AC_SUBST(LIBS)
fi
if test -z "$COMPILER_INTEL_TRUE"; then
  FCFLAGS="$FCFLAGS -convert big_endian -assume byterecl \
           -heap-arrays -fPIC"
  IFORTVER=`ifort --version | grep IFORT | cut -d " " -f 3 | cut -d "." -f 1`
  if test $IFORTVER -ge 17; then
    OMPFLAGS="-qopenmp"
  else
    OMPFLAGS="-openmp"
  fi
  CPRDEF="INTEL"
  AC_SUBST(FCFLAGS)
  AC_SUBST(OMPFLAGS)
  LIBS=`echo $LIBS | sed -e 's/-Wl,--as-needed//g'`
  AC_SUBST(LIBS)
  if eval "test x"$USE_KNL" = x1"; then
    FCFLAGS="$FCFLAGS -xMIC-AVX512"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$USE_ZEN3" = x1"; then
    FCFLAGS="$FCFLAGS -axCORE-AVX2"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$USE_SKL" = x1"; then
    FCFLAGS="$FCFLAGS -xCORE-AVX512"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$USE_BDW" = x1"; then
    FCFLAGS="$FCFLAGS -xCORE-AVX2"
    AC_SUBST(FCFLAGS)
  fi
  if eval "test x"$USE_NHL" = x1"; then
    FCFLAGS="$FCFLAGS -xSSE4.2"
    AC_SUBST(FCFLAGS)
  fi
fi
if test -z "$COMPILER_G95_TRUE"; then
  FCFLAGS="$FCFLAGS -fendian=big -fmultiple-save"
  OMPFLAGS=""
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
  AC_SUBST(OMPFLAGS)
fi
if test -z "$COMPILER_PGI_TRUE"; then
  FCFLAGS="$FCFLAGS -byteswapio"
  OMPFLAGS="-mp"
  CPRDEF="PGI"
  AC_SUBST(FCFLAGS)
  AC_SUBST(OMPFLAGS)
fi
if test -z "$COMPILER_IBM_TRUE"; then
  DEFINE="-WF,-D"
  FCFLAGS="$FCFLAGS -WF,-DIBM"
  OMPFLAGS="-qsmp=omp"
  CPRDEF="XLF"
  AC_SUBST(FCFLAGS)
  AC_SUBST(OMPFLAGS)
  AC_SUBST(DEFINE)
fi
if test -x "$COMPILER_SUN_TRUE"; then
  CPRDEF="WORKSHOP"
  OMPFLAGS="-xopenmp=parallel"
  AC_SUBST(OMPFLAGS)
fi
if test -z "$COMPILER_CRAY_TRUE"; then
  FCFLAGS="$FCFLAGS"
  CPRDEF="GNU"
  AC_SUBST(FCFLAGS)
fi

AC_SUBST(CPRDEF)

# Set debug/production flags
if eval "test x$DEBUG = x1"; then
  if test -z "$COMPILER_G95_TRUE"; then
    FCFLAGS="$FCFLAGS -O0 -g -Wall -Wextra -fbounds-check -ftrace=full"
    NOOPTFC="$FCFLAGS"
    AC_SUBST(FCFLAGS)
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_GNU_TRUE"; then
    if test -z "${osname##Linux}"; then
      FCFLAGS="-O0 -g -pg -fbacktrace \
            -fdump-core -Wall -pedantic -Wconversion -Wsurprising -Waliasing \
            -Wampersand -Warray-bounds -Wcharacter-truncation -Wuninitialized \
            -Wline-truncation -Wintrinsics-std -Wno-tabs -Wunderflow \
            -Wunused-parameter -Wno-align-commons $FCFLAGS \
            -fno-range-check -fcheck=bounds -Wimplicit-interface \
            -Wimplicit-procedure -Wintrinsic-shadow \
            -fcheck=array-temps,bounds,do,mem,recursion,pointer"
            # -ffpe-trap=zero,invalid,overflow,underflow
    else
      # Newver Macos versions have a missing library problem for -pg
      FCFLAGS="-O0 -g -fbacktrace \
            -fdump-core -Wall -pedantic -Wconversion -Wsurprising -Waliasing \
            -Wampersand -Warray-bounds -Wcharacter-truncation -Wuninitialized \
            -Wline-truncation -Wintrinsics-std -Wno-tabs -Wunderflow \
            -Wunused-parameter -Wno-align-commons $FCFLAGS \
            -fno-range-check -Wimplicit-interface -Wimplicit-procedure \
            -Wintrinsic-shadow -fcheck=array-temps,bounds,do,mem,recursion"
            # -ffpe-trap=zero,invalid,overflow,underflow
    fi
    NOOPTFC="$FCFLAGS"
    AC_SUBST(FCFLAGS)
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_INTEL_LLVM_TRUE"; then
    FCFLAGS="-O0 -fpe0 -g -check all -traceback -warn all -ftz \
             -fimf-arch-consistency=true -fltconsistency -mieee-fp $FCFLAGS"
    NOOPTFC="$FCFLAGS"
    AC_SUBST(FCFLAGS)
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_INTEL_TRUE"; then
    FCFLAGS="-O0 -fpe0 -g -check all -traceback -warn all -ftz \
             -fimf-arch-consistency=true -fltconsistency -mieee-fp $FCFLAGS"
    NOOPTFC="$FCFLAGS"
    AC_SUBST(FCFLAGS)
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_IBM_TRUE"; then
    FCFLAGS="-O0 -qsuppress=1501-245 -qcheck -qflttrap=nanq,overflow \
             -qxlf2003=polymorphic -qmaxmem=-1 -qzerosize $FCFLAGS"
    NOOPTFC="$FCFLAGS"
    AC_SUBST(FCFLAGS)
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_PGI_TRUE"; then
    FCFLAGS="-Mbounds -Mcache_align -Mchkstk -Mdalign \
             -Miomutex -Mrecursive -Ktrap=divz -O0 -g \
             -gopt -traceback $FCFLAGS"
    if test "test x$USE_OPENACC = x1"; then
      FCFLAGS="$FCFLAGS -acc=gpu -gpu=$PGI_GPU_ARCH,lineinfo -Minfo=accel"
    fi
    if test "test x$USE_OPENACC = x4"; then
      FCFLAGS="$FCFLAGS -acc=gpu \
	      -gpu=$PGI_GPU_ARCH,managed,ptxinfo,lineinfo -Mnodepchk -Minfo=all"
    fi
    NOOPTFC="$FCFLAGS"
    AC_SUBST(FCFLAGS)
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_CRAY_TRUE"; then
    FCFLAGS="-g $(FCFLAGS)"
    NOOPTFC="$FCFLAGS"
    AC_SUBST(FCFLAGS)
    AC_SUBST(NOOPTFC)
    CPRDEF="GNU"
  fi
else
  if test -z "$COMPILER_G95_TRUE"; then
    NOOPTFC="$FCFLAGS"
    FCFLAGS="-O3 $FCFLAGS"
    AC_SUBST(FCFLAGS)
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_GNU_TRUE"; then
    # The GNU compiler works with -O2
    NOOPTFC="-O2 -finline-functions $FCFLAGS" # -fno-var-tracking-assignments
    if test -z "${osname##Linux}"; then
      FCFLAGS="-O3 -finline-functions -mtune=native -ftree-vectorize \
               -march=native $FCFLAGS" # -ftree-parallelize-loops=2"
      if eval "test x"$USE_PROF" = x1"; then
        FCFLAGS="-g -pg $FCFLAGS"
      fi
    else
      # MacOSX may choke on using AVX extensions...
      FCFLAGS="-O3 -finline-functions $FCFLAGS"
      # Newer versions of MacOS have problems with pg options
      if eval "test x"$USE_PROF" = x1"; then
        FCFLAGS="-g $FCFLAGS"
      fi
    fi
    AC_SUBST(FCFLAGS)
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_INTEL_LLVM_TRUE"; then
    AR=xiar
    AC_SUBST(AR)
    NOOPTFC="-ftz -O0 -fp-model source $FCFLAGS"
    if eval "test x"$USE_PROF" = x1"; then
      FCFLAGS="-ftz -O2 -g -p -fp-model source $FCFLAGS"
      AC_SUBST(FCFLAGS)
    else
	    FCFLAGS="-ftz -O3 -finline-functions -fp-model source $FCFLAGS"
      AC_SUBST(FCFLAGS)
    fi
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_INTEL_TRUE"; then
    AR=xiar
    AC_SUBST(AR)
    NOOPTFC="-ftz -O0 -fp-model source $FCFLAGS"
    if eval "test x"$USE_PROF" = x1"; then
      FCFLAGS="-ftz -O2 -g -p -fp-model source $FCFLAGS"
      AC_SUBST(FCFLAGS)
    else
      FCFLAGS="-ftz -ip -O3 -finline-functions -fp-model source $FCFLAGS"
#      FCFLAGS="-DF2008 -ftz -fpe0 -O2 -g -traceback -fp-stack-check \
#	      -check all -fimf-arch-consistency=true -fltconsistency \
#	      -mieee-fp $FCFLAGS"
      AC_SUBST(FCFLAGS)
    fi
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_IBM_TRUE"; then
    NOOPTFC="$FCFLAGS"
    FCFLAGS="-qsuppress=1501-245 -O3 -qstrict -qxlf2003=polymorphic \
             -qmaxmem=-1 -qzerosize $FCFLAGS"
    AC_SUBST(FCFLAGS)
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_SUN_TRUE"; then
    NOOPTFC="-xfilebyteorder=big16:%all \
	    -D__DATE__="'`date +%D`'" -D__TIME__="'`date +%T`'" $FCFLAGS"
    FCFLAGS="-O3 -xfilebyteorder=big16:%all \
             -D__DATE__="'`date +%D`'" -D__TIME__="'`date +%T`'" $FCFLAGS"
    AC_SUBST(FCFLAGS)
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_PGI_TRUE"; then
    NOOPTFC="-O1 $FCFLAGS"
    FCFLAGS="-O3 -fast -fastsse $FCFLAGS"
    if eval "test x"$USE_OPENACC" = x1"; then
      FCFLAGS="$FCFLAGS -acc=gpu -gpu=$PGI_GPU_ARCH,lineinfo -Minfo=accel"
      NOOPTFC="$NOOPTFC -acc=gpu -gpu=$PGI_GPU_ARCH,lineinfo -Minfo=accel"
    fi
    if eval "test x"$USE_OPENACC" = x2"; then
      FCFLAGS="$NOOPTFC -acc=gpu \
          -gpu="${PGI_GPU_ARCH}",autocompare -Minfo=accel"
      NOOPTFC="$NOOPTFC -acc=gpu \
          -gpu="${PGI_GPU_ARCH}",autocompare -Minfo=accel"
    fi
    if eval "test x"$USE_OPENACC" = x3"; then
      FCFLAGS="$FCFLAGS -acc=gpu \
          -gpu=$PGI_GPU_ARCH,lineinfo,managed -Minfo=accel"
      NOOPTFC="$NOOPTFC -acc=gpu \
          -gpu=$PGI_GPU_ARCH,lineinfo,managed -Minfo=accel"
    fi
    if eval "test x"$USE_OPENACC" = x4"; then
      FCFLAGS="$FCFLAGS -stdpar=gpu,multicore \
          -gpu=$PGI_GPU_ARCH,lineinfo,managed -Minfo=accel"
      NOOPTFC="$NOOPTFC -stdpar=gpu,multicore \
          -gpu=$PGI_GPU_ARCH,lineinfo,managed -Minfo=accel"
    fi
    AC_SUBST(FCFLAGS)
    AC_SUBST(NOOPTFC)
  fi
  if test -z "$COMPILER_CRAY_TRUE"; then
    NOOPTFC="$FCFLAGS"
    FCFLAGS="$(FCFLAGS)"
    AC_SUBST(FCFLAGS)
    CPRDEF="GNU"
    AC_SUBST(NOOPTFC)
  fi
fi

RR_PATH_NETCDF
RR_PATH_NETCDF_F90

if eval "test x$PNETCDF = x1"; then
  if eval "test x$PNETCDF_IN_NETCDF = x1"; then
    RR_NETCDF4
  else
    RR_PATH_PNETCDF
    if eval "test x$NETCDF_CDF5 = x1"; then
      RR_CDF5
    fi
  fi
else
  if eval "test x$NETCDF_CDF5 = x1"; then
    RR_CDF5
  else
    RR_NETCDF4
  fi
fi

AM_CONDITIONAL(RRTM_NC, [test x"$netcdf" = xyes])

RCM_FC_CHECK_IEEE_ARITHMETIC
RCM_FC_CHECK_QUAD_PRECISION

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
  RCM_MPI_CHECK_MPI3
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

if eval "test x"$USEOASIS" = x1"; then
  if test -z "$COMPILER_SUN_TRUE"; then
    AM_CPPFLAGS="$AM_CPPFLAGS -M$OASIS_PREFIX/build/lib/psmile.$OASIS_CHAN"
    AM_CPPFLAGS="$AM_CPPFLAGS -M$OASIS_PREFIX/build/lib/scrip"
    AM_CPPFLAGS="$AM_CPPFLAGS -M$OASIS_PREFIX/build/lib/mct"
    AC_SUBST(AM_CPPFLAGS)
  else
    AM_CPPFLAGS="$AM_CPPFLAGS -I$OASIS_PREFIX/build/lib/psmile.$OASIS_CHAN"
    AM_CPPFLAGS="$AM_CPPFLAGS -I$OASIS_PREFIX/build/lib/scrip"
    AM_CPPFLAGS="$AM_CPPFLAGS -I$OASIS_PREFIX/build/lib/mct"
    AC_SUBST(AM_CPPFLAGS)
  fi
  OASIS3_LIBS="-L$OASIS_PREFIX/lib -lpsmile.$OASIS_CHAN -lmct -lmpeu -lscrip"
  LIBS="$LIBS $OASIS3_LIBS"
  AC_SUBST(LIBS)
fi

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
                 PreProc/Terrain/Makefile PreProc/ICBC/cmip6/Makefile \
                 PreProc/ICBC/Makefile PreProc/CH_EMISSIONS/Makefile \
                 PreProc/CH_EMISSIONS/interpolation/Makefile \
                 PreProc/CH_EMISSIONS/interpolation/tools/Makefile \
                 PreProc/CH_EMISSIONS/interpolation/tools/REGCM_grid/Makefile \
                 PreProc/CLM/Makefile PreProc/CLM45/Makefile PostProc/Makefile \
                 Main/mpplib/Makefile Main/netlib/Makefile \
                 Main/Makefile Main/batslib/Makefile Main/clmlib/Makefile \
                 Main/ocnlib/Makefile Main/clmlib/clm4.5/Makefile \
                 Main/chemlib/Makefile Main/chemlib/GAS_CBMZ_NEW/Makefile \
                 Main/cumlib/Makefile Main/microlib/Makefile \
                 Main/radlib/Makefile Main/radlib/RRTMG/Makefile \
                 Main/pbllib/Makefile Main/cloudlib/Makefile \
                 Main/oasislib/Makefile \
                 Tools/Programs/CheckSun/Makefile \
                 Tools/Programs/RegCM_read/Makefile \
                 Tools/Programs/timeseries/Makefile])

AC_OUTPUT
