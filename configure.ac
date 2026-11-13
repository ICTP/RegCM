#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

AC_INIT([RegCM], [5.0.0], [ggiulian@ictp.it])

# Whether the user gave FCFLAGS: AC_PROG_FC otherwise sets a default (-g,
# plus -O2 for a GNU compiler) that some vendor branches below drop
rcm_user_fcflags=${FCFLAGS+yes}

AC_SUBST([PACKAGE_URL], ["https://github.com/ICTP/RegCM"])

#AC_PREREQ([2.64])

AM_INIT_AUTOMAKE([foreign])
AM_PROG_AR

LT_INIT

AC_PROG_CC
AC_PROG_FC

AC_PROG_RANLIB

AC_LANG([Fortran])

AC_ARG_WITH([netcdf],
  AS_HELP_STRING([--with-netcdf],
            [Path to NetCDF installation (default: NETCDF environment)]),
  [NC_PREFIX="$withval"; NF_PREFIX="$withval"; NFAUTO="no"; NCAUTO="no"],
  [FAKE="no"]
)

AC_ARG_WITH([hdf5],
  AS_HELP_STRING([--with-hdf5],
            [Path to HDF5 installation (default: HDF5 environment)]),
  [HDF5_PREFIX="$withval"], [HDF5_PREFIX="${HDF5}"]
)

AC_ARG_WITH([szip],
  AS_HELP_STRING([--with-szip],
            [Path to SZIP installation (default: SZIP environment)]),
  [SZIP_PREFIX="$withval"], [SZIP_PREFIX="${SZIP}"]
)

AC_ARG_WITH([oasis-path],
  AS_HELP_STRING([--with-oasis-path],
		 [Path to OASIS3-MCT installation (require --enable-oasis) [default: OASIS_PREFIX]]),
  [OASIS_PREFIX="$withval"], [OASIS_PREFIX="${OASIS_PREFIX}"]
)

AC_ARG_WITH([oasis-comm],
  AS_HELP_STRING([--with-oasis-comm],
		 [Communication technique used in OASIS3-MCT (require --enable-oasis) [default: MPI1]]),
  [OASIS_CHAN="$withval"], [OASIS_CHAN="MPI1"]
)

SP_REAL="0"
AC_ARG_ENABLE(singleprecision,
  AS_HELP_STRING([--enable-singleprecision],
                 [Default real type is real4 single precision.]),
                 [AS_IF([test "x$enableval" != xno], [SP_REAL="1"])]
)
if test "x$SP_REAL" = x1; then
  AM_CPPFLAGS="-DSINGLE_PRECISION_REAL $AM_CPPFLAGS"
fi

USE_KNL="0"
AC_ARG_ENABLE(knl,
  AS_HELP_STRING([--enable-knl],
                 [Compile with Xeon Phi flags.]),
                 [AS_IF([test "x$enableval" != xno], [USE_KNL="1"])]
)

USE_ZEN3="0"
AC_ARG_ENABLE(zen3,
  AS_HELP_STRING([--enable-zen3],
                 [Compile with Zen3 processor optimization]),
                 [AS_IF([test "x$enableval" != xno], [USE_ZEN3="1"])]
)

USE_BDW="0"
AC_ARG_ENABLE(bdw,
  AS_HELP_STRING([--enable-bdw],
                 [Compile with AVX2 flags.]),
                 [AS_IF([test "x$enableval" != xno], [USE_BDW="1"])]
)

USE_NHL="0"
AC_ARG_ENABLE(nhl,
  AS_HELP_STRING([--enable-nhl],
                 [Compile with SSE 4.2 flags.]),
                 [AS_IF([test "x$enableval" != xno], [USE_NHL="1"])]
)

USE_SKL="0"
AC_ARG_ENABLE(skl,
  AS_HELP_STRING([--enable-skl],
                 [Compile with Sky Lake flags.]),
                 [AS_IF([test "x$enableval" != xno], [USE_SKL="1"])]
)

USE_OPENACC="0"
# TODO: make this configurable
PGI_GPU_ARCH="ccnative"

AC_ARG_ENABLE(openacc,
  AS_HELP_STRING([--enable-openacc],
                 [Compile OpenACC directives (NVIDIA nvfortran only)]),
                 [AS_IF([test "x$enableval" != xno], [USE_OPENACC="1"])]
)

AC_ARG_ENABLE(openacc-debug,
  AS_HELP_STRING([--enable-openacc-debug],
                 [Compile OpenACC directives with autocompare (NVIDIA nvfortran only)]),
                 [AS_IF([test "x$enableval" != xno], [USE_OPENACC="2"])]
)

AC_ARG_ENABLE(openacc-managed,
  AS_HELP_STRING([--enable-openacc-managed],
                 [Compile OpenACC directives with managed memory (NVIDIA nvfortran only)]),
                 [AS_IF([test "x$enableval" != xno], [USE_OPENACC="3"])]
)

AC_ARG_ENABLE(openacc-stdpar,
  AS_HELP_STRING([--enable-openacc-stdpar],
                 [Compile OpenACC directives with stdpar (NVIDIA nvfortran only)]),
                 [AS_IF([test "x$enableval" != xno], [USE_OPENACC="4"])]
)
if test "x$USE_OPENACC" != x0; then
  AM_CPPFLAGS="-DOPENACC $AM_CPPFLAGS"
fi
if test "x$USE_OPENACC" = x4; then
  AM_CPPFLAGS="$AM_CPPFLAGS -DSTDPAR"
fi

USE_PROF="0"
AC_ARG_ENABLE(profile,
  AS_HELP_STRING([--enable-profile],
                 [Compile with profile flags.]),
                 [AS_IF([test "x$enableval" != xno], [USE_PROF="1"])]
)

USEOASIS="0"
AC_ARG_ENABLE(oasis,
  AS_HELP_STRING([--enable-oasis],
                 [Supply this option if you plan on using OASIS option.]),
  [AS_IF([test "x$enableval" != xno], [USEOASIS="1"])]
)
AM_CONDITIONAL(DO_COMPILE_OASIS, [test x"$USEOASIS" = x1])

if test "x$USEOASIS" = x1; then
  AC_MSG_WARN([Compiling OASIS enabled code])
  AC_MSG_WARN([OASIS path: $OASIS_PREFIX])
  AC_MSG_WARN([OASIS communication technique: $OASIS_CHAN])
  AM_CPPFLAGS="-DOASIS $AM_CPPFLAGS"
fi

USEECLM="0"
AC_ARG_ENABLE(eclm,
  AS_HELP_STRING([--enable-eclm],
                 [Supply this option if you plan on using ECLM option.]),
  [AS_IF([test "x$enableval" != xno], [USEECLM="1"])]
)
AM_CONDITIONAL(DO_COMPILE_ECLM, [test x"$USEECLM" = x1])

if test "x$USEECLM" = x1; then
  USEOASIS=1
  AM_CONDITIONAL(DO_COMPILE_OASIS, [test x"$USEOASIS" = x1])
  AC_MSG_WARN([Compiling ECLM and OASIS enabled code])
  AC_MSG_WARN([OASIS path: $OASIS_PREFIX])
  AC_MSG_WARN([OASIS communication technique: $OASIS_CHAN])
  AM_CPPFLAGS="-DOASIS -DECLM $AM_CPPFLAGS"
fi

USECLM45="0"
AC_ARG_ENABLE(clm45,
  AS_HELP_STRING([--enable-clm45],
                 [Supply this option if you plan on using CLM45 option.]),
  [AS_IF([test "x$enableval" != xno], [USECLM45="1"])],
)
AM_CONDITIONAL(DO_COMPILE_LIBCLM45, [test x"$USECLM45" = x1])

if test "x$USECLM45" = x1; then
  AC_MSG_WARN([Compiling CLM45 enabled code])
  AM_CPPFLAGS="-DCLM45 $AM_CPPFLAGS"
fi

USECPL="0"
AC_ARG_ENABLE(cpl,
  AS_HELP_STRING([--enable-cpl],
                 [Supply this option if you plan on using REGESM.]),
  [AS_IF([test "x$enableval" != xno], [USECPL="1"])],
)
AM_CONDITIONAL(DO_COMPILE_CPL, [test x"$USECPL" = x1])

if test "x$USECPL" = x1; then
  AM_CPPFLAGS="-DCPL $AM_CPPFLAGS"
fi

PNETCDF_IN_NETCDF=0
AC_ARG_ENABLE(parallel-nc,
  AS_HELP_STRING([--enable-parallel-nc],
                 [Supply this option if you plan on using NetCDF4 parallel option based upon Parallel Netcdf.]),
  [AS_IF([test "x$enableval" != xno], [PNETCDF_IN_NETCDF="1"])]
)
if test "x$PNETCDF_IN_NETCDF" = x1; then
  AM_CPPFLAGS="-DPNETCDF_IN_NETCDF $AM_CPPFLAGS"
fi

PNETCDF="0"
AC_ARG_ENABLE(pnetcdf,
  AS_HELP_STRING([--enable-pnetcdf],
                 [Use Parallel NetCDF from NU-ANL instead of netcdf.]),
  [AS_IF([test "x$enableval" != xno], [PNETCDF="1"])]
)
if test "x$PNETCDF" = x1; then
  NULLOPTION=$( echo $AM_CPPFLAGS | sed -e 's/-DPNETCDF_IN_NETCDF//')
  AM_CPPFLAGS="-DPNETCDF $NULLOPTION"
fi

NETCDF4_COMPRESS="0"
AC_ARG_ENABLE(nc4-filters,
  AS_HELP_STRING([--enable-nc4-filters],
		 [Supply this option if you plan on using NetCDF4/HDF5 output filtering option.]),
  [AS_IF([test "x$enableval" != xno], [NETCDF4_COMPRESS="1"])]
  )
if test "x$NETCDF4_COMPRESS" = x1; then
  NULLOPTION=$( echo $AM_CPPFLAGS | sed -e 's/-DPNETCDF_IN_NETCDF//' \
                -e 's/-DPNETCDF//')
  AM_CPPFLAGS="-DNETCDF4_COMPRESS $NULLOPTION"
fi

NETCDF_CDF5=0
AC_ARG_ENABLE(cdf5,
  AS_HELP_STRING([--enable-cdf5],
                 [Supply this option if you plan on using NetCDF/CDF5 option.]),
  [AS_IF([test "x$enableval" != xno], [NETCDF_CDF5="1"])]
)

ASYNC_NETCDF="0"
AC_ARG_ENABLE(async-netcdf,
  AS_HELP_STRING([--enable-async-netcdf],
                 [Enable pthread-backed async NetCDF output and ICBC prefetching.]),
  [AS_IF([test "x$enableval" != xno], [ASYNC_NETCDF="1"])]
)
AM_CONDITIONAL(DO_ASYNC_NETCDF, [test x"$ASYNC_NETCDF" = x1])
if test "x$ASYNC_NETCDF" = x1; then
  AM_CPPFLAGS="-DASYNC_NETCDF $AM_CPPFLAGS"
fi

REPRODUCIBLE="0"
AC_ARG_ENABLE(reproducible,
  AS_HELP_STRING([--enable-reproducible],
                 [GNU Fortran: no glibc vector math library (libmvec) and
                  no fused multiply-add contraction, so that results do
                  not depend on the MPI domain decomposition. Slower.
                  [default: no]]),
  [AS_IF([test "x$enableval" != xno], [REPRODUCIBLE="1"])]
)

DEBUG="0"
AC_ARG_ENABLE(debug,
  AS_HELP_STRING([--enable-debug],
                 [Enable debugging flags and per processor log file]),
  [AS_IF([test "x$enableval" != xno], [DEBUG="1"])]
)
if test "x$DEBUG" = x1; then
  AM_CPPFLAGS="-DDEBUG $AM_CPPFLAGS"
fi

NC_INCLUDES=""
NC_LDFLAGS=""
NC_LIBS=""

AX_PROG_NF_CONFIG([NF_PREFIX=`nf-config --prefix`; NFAUTO="yes"],
                  [NF_PREFIX="${NETCDF}"; NFAUTO="no"])

AX_PROG_NC_CONFIG([NC_PREFIX=`nc-config --prefix`; NCAUTO="yes"],
                  [NC_PREFIX="${NETCDF}"; NCAUTO="no"])

NC_CC="$CC"
if test "x$NFAUTO" = xyes; then
  if test "x$NCAUTO" = xyes; then
    NC_CC="`nc-config --cc`"
    NC_INCLUDES="`nc-config --cflags` `nf-config --fflags`"
    NC_LDFLAGS="-L`nf-config --prefix`/lib -L`nc-config --prefix`/lib"
    NC_LIBS="`nc-config --libs` `nf-config --flibs`"
  else
    NC_INCLUDES="-I$NC_PREFIX/include `nf-config --fflags`"
    NC_LDFLAGS="-L$NC_PREFIX/lib -L`nc-config --prefix`/lib"
    NC_LIBS="-L$NC_PREFIX/lib `nf-config --flibs`"
  fi
else
  if test "x$NCAUTO" = xno; then
    NC_INCLUDES="-I$NC_PREFIX/include -I$NF_PREFIX/include"
    NC_LDFLAGS="-L$NC_PREFIX/lib -L$NF_PREFIX/lib"
    NC_LIBS="-lnetcdff -lnetcdf"
  else
    NCVER=`nc-config --version | cut -d " " -f 2`
    if test "x$NCVER" = x4.0.1; then
      NC_INCLUDES="`nc-config --fflags`"
    else
      # new netcdf library has this
      NC_INCLUDES="-I`nc-config --includedir`"
    fi
    NC_LDFLAGS="-L`nc-config --prefix`/lib"
    NC_LIBS="`nc-config --libs`"
  fi
fi

if test "x$PNETCDF" = x1; then
  AX_PROG_PNETCDF_CONFIG([NCP_PREFIX=`pnetcdf-config --prefix`; NCPAUTO="yes"],
                         [NCP_PREFIX="${NETCDF}"; NCPAUTO="no"])
  if test "x$NCPAUTO" = xyes; then
    CC=`pnetcdf-config --cc`
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

# Fortran compiler vendor, from the compiler's own version banner
RCM_FC_VENDOR([$FC], [fc_vendor])
AC_MSG_CHECKING([for the Fortran compiler vendor])
AC_MSG_RESULT([$fc_vendor])
if test "$fc_vendor" = unknown; then
  AC_MSG_WARN([Fortran compiler $FC not recognised: no vendor flags
               (big-endian unformatted I/O, OpenMP, optimisation) are set])
fi
AM_CONDITIONAL(COMPILER_GNU, [test "$fc_vendor" = gnu])
AM_CONDITIONAL(COMPILER_INTEL, [test "$fc_vendor" = intel])
AM_CONDITIONAL(COMPILER_INTEL_LLVM, [test "$fc_vendor" = intel_llvm])
AM_CONDITIONAL(COMPILER_PGI, [test "$fc_vendor" = nvhpc])
AM_CONDITIONAL(COMPILER_NAG, [test "$fc_vendor" = nag])
AM_CONDITIONAL(COMPILER_IBM, [test "$fc_vendor" = ibm])
AM_CONDITIONAL(COMPILER_CRAY, [test "$fc_vendor" = cray])

# Base flags: big-endian unformatted I/O, OpenMP, target architecture
case $fc_vendor in
  gnu)
    # The default -g -O2 would follow, and override, the -O3 or -O0 below
    test "x$rcm_user_fcflags" = xyes || FCFLAGS=""
    FCFLAGS="$FCFLAGS -fconvert=big-endian -fno-range-check -fPIC"
    OMPFLAGS="-fopenmp"
    test "x$USE_KNL" = x1 && FCFLAGS="$FCFLAGS -march=knl"
    test "x$USE_ZEN3" = x1 && FCFLAGS="$FCFLAGS -march=znver3"
    test "x$USE_SKL" = x1 && FCFLAGS="$FCFLAGS -march=skylake-avx512"
    test "x$USE_BDW" = x1 && FCFLAGS="$FCFLAGS -march=broadwell"
    test "x$USE_NHL" = x1 && FCFLAGS="$FCFLAGS -march=nehalem"
    ;;
  nag)
    FCFLAGS="$FCFLAGS -convert=big-endian -pic -f2018"
    OMPFLAGS="-openmp"
    ;;
  intel_llvm)
    if test "x$DEBUG" = x1; then
      OMPFLAGS=""
    else
      OMPFLAGS="-fiopenmp -fno-openmp-target-do-concurrent"
    fi
    test "x$rcm_user_fcflags" = xyes || FCFLAGS=""
    FCFLAGS="$FCFLAGS -convert big_endian -assume byterecl \
             -heap-arrays -qopt-report -fPIC"
    LIBS=`echo $LIBS | sed -e 's/-Wl,--as-needed//g'`
    test "x$USE_ZEN3" = x1 && FCFLAGS="$FCFLAGS -axCORE-AVX2"
    test "x$USE_SKL" = x1 && FCFLAGS="$FCFLAGS -xCORE-AVX512"
    test "x$USE_BDW" = x1 && FCFLAGS="$FCFLAGS -xCORE-AVX2"
    test "x$USE_NHL" = x1 && FCFLAGS="$FCFLAGS -xSSE4.2"
    ;;
  intel)
    test "x$rcm_user_fcflags" = xyes || FCFLAGS=""
    FCFLAGS="$FCFLAGS -convert big_endian -assume byterecl -heap-arrays -fPIC"
    if test -n "$fc_vendor_major" && test "$fc_vendor_major" -lt 17; then
      OMPFLAGS="-openmp"
    else
      OMPFLAGS="-qopenmp"
    fi
    LIBS=`echo $LIBS | sed -e 's/-Wl,--as-needed//g'`
    test "x$USE_KNL" = x1 && FCFLAGS="$FCFLAGS -xMIC-AVX512"
    test "x$USE_ZEN3" = x1 && FCFLAGS="$FCFLAGS -axCORE-AVX2"
    test "x$USE_SKL" = x1 && FCFLAGS="$FCFLAGS -xCORE-AVX512"
    test "x$USE_BDW" = x1 && FCFLAGS="$FCFLAGS -xCORE-AVX2"
    test "x$USE_NHL" = x1 && FCFLAGS="$FCFLAGS -xSSE4.2"
    ;;
  nvhpc)
    test "x$rcm_user_fcflags" = xyes || FCFLAGS=""
    FCFLAGS="$FCFLAGS -byteswapio"
    OMPFLAGS="-mp"
    ;;
  ibm)
    FCFLAGS="$FCFLAGS -DIBM"
    OMPFLAGS="-qsmp=omp"
    ;;
  cray)
    # Drop the autoconf defaults: -g limits the CCE optimisation, and a
    # later -O2 would override the optimisation level set below
    FCFLAGS=`echo " $FCFLAGS " | sed -e 's/ -g / /g' -e 's/ -g / /g' \
             -e 's/ -O2 / /g' -e 's/  */ /g' -e 's/^ //' -e 's/ $//'`
    FCFLAGS="$FCFLAGS -h byteswapio"
    OMPFLAGS="-fopenmp"
    ;;
esac

# Only the GNU optimisation and debug flags depend on the OS: Linux or
# not (macOS).
osname=`uname -s`

# Debug or production flags
if test "x$DEBUG" = x1; then
  case $fc_vendor in
    nag)
      FCFLAGS="-O0 -g -pg -mtrace=all -gline -ieee=nonstd -nan \
               -Warn=allocation -C=all $FCFLAGS"
      NOOPTFC="$FCFLAGS"
      ;;
    gnu)
      if test -z "${osname##Linux}"; then
        FCFLAGS="-O0 -g -pg -fbacktrace \
              -fdump-core -Wall -pedantic -Wconversion -Wsurprising -Waliasing \
              -Wampersand -Warray-bounds -Wcharacter-truncation -Wuninitialized \
              -Wline-truncation -Wintrinsics-std -Wno-tabs -Wunderflow \
              -Wunused-parameter -Wno-align-commons $FCFLAGS \
              -fno-range-check -fcheck=bounds -Wimplicit-interface \
              -Wimplicit-procedure -Wintrinsic-shadow \
              -fcheck=array-temps,bounds,do,mem,recursion,pointer \
              -ffpe-trap=zero,invalid,overflow,underflow"
      else
        # Newer macOS versions have a missing library problem for -pg
        FCFLAGS="-O0 -g -fbacktrace \
              -fdump-core -Wall -pedantic -Wconversion -Wsurprising -Waliasing \
              -Wampersand -Warray-bounds -Wcharacter-truncation -Wuninitialized \
              -Wline-truncation -Wintrinsics-std -Wno-tabs -Wunderflow \
              -Wunused-parameter -Wno-align-commons $FCFLAGS \
              -fno-range-check -Wimplicit-interface -Wimplicit-procedure \
              -Wintrinsic-shadow -fcheck=array-temps,bounds,do,mem,recursion"
      fi
      NOOPTFC="$FCFLAGS"
      ;;
    intel_llvm)
      FCFLAGS="-O0 -g -traceback -warn all -check bounds,contiguous \
               -mieee-fp -fpe0 $FCFLAGS"
      NOOPTFC="$FCFLAGS"
      ;;
    intel)
      FCFLAGS="-O0 -fpe0 -g -check all -traceback -warn all -ftz \
               -fimf-arch-consistency=true -fltconsistency -mieee-fp $FCFLAGS"
      NOOPTFC="$FCFLAGS"
      ;;
    ibm)
      FCFLAGS="-O0 -qsuppress=1501-245 -qcheck -qflttrap=nanq,overflow \
               -qxlf2003=polymorphic -qmaxmem=-1 -qzerosize $FCFLAGS"
      NOOPTFC="$FCFLAGS"
      ;;
    nvhpc)
      FCFLAGS="-Mbounds -Mcache_align -Mchkstk -Mdalign \
               -Miomutex -Mrecursive -Ktrap=divz -O0 -g \
               -gopt -traceback $FCFLAGS"
      NOOPTFC="$FCFLAGS"
      ;;
    cray)
      FCFLAGS="-O0 -g $FCFLAGS"
      NOOPTFC="$FCFLAGS"
      ;;
  esac
else
  case $fc_vendor in
    nag)
      NOOPTFC="-O2 $FCFLAGS"
      FCFLAGS="-O3 $FCFLAGS"
      ;;
    gnu)
      NOOPTFC="-O2 -finline-functions $FCFLAGS"
      if test -z "${osname##Linux}"; then
        FCFLAGS="-O3 -finline-functions -mtune=native -ftree-vectorize \
                 -march=native $FCFLAGS"
        test "x$USE_PROF" = x1 && FCFLAGS="-g -pg $FCFLAGS"
      else
        # macOS may choke on AVX extensions and on -pg
        FCFLAGS="-O3 -finline-functions $FCFLAGS"
        test "x$USE_PROF" = x1 && FCFLAGS="-g $FCFLAGS"
      fi
      ;;
    intel_llvm)
      NOOPTFC="-ftz -O0 -fp-model source $FCFLAGS"
      if test "x$USE_PROF" = x1; then
        FCFLAGS="-ftz -O2 -g -p -fp-model source $FCFLAGS"
      else
        FCFLAGS="-ftz -O3 -flto -finline-functions -fp-model source $FCFLAGS"
      fi
      ;;
    intel)
      NOOPTFC="-ftz -O0 -fp-model source $FCFLAGS"
      if test "x$USE_PROF" = x1; then
        FCFLAGS="-ftz -O2 -g -p -fp-model source $FCFLAGS"
      else
        FCFLAGS="-ftz -ip -O3 -finline-functions -fp-model source $FCFLAGS"
      fi
      ;;
    ibm)
      NOOPTFC="$FCFLAGS"
      FCFLAGS="-qsuppress=1501-245 -O3 -qstrict -qxlf2003=polymorphic \
               -qmaxmem=-1 -qzerosize $FCFLAGS"
      ;;
    nvhpc)
      NOOPTFC="-O1 $FCFLAGS"
      FCFLAGS="-O3 -fast -fastsse $FCFLAGS"
      ;;
    cray)
      NOOPTFC="$FCFLAGS"
      FCFLAGS="-O3 $FCFLAGS"
      ;;
  esac
fi

# Which points of a loop fall in its vector body and which in the scalar
# remainder depends on the subdomain size, so the two must give the same
# result to the last bit:
# - Ubuntu and Debian gfortran pre-include math-vector-fortran.h, which
#   makes vectorised loops call libmvec, whose functions differ from the
#   scalar ones in the last bit. Leaving out the standard include
#   directories drops the header; the intrinsic modules are then found
#   with an explicit path.
# - The compiler may contract a*b+c into a fused multiply-add in one part
#   of the loop and not in the other.
if test "x$REPRODUCIBLE" = x1; then
  if test "$fc_vendor" = gnu; then
    GFC_FINCLUDE=`$FC -print-file-name=finclude`
    if test -d "$GFC_FINCLUDE"; then
      REPROFLAGS="-nostdinc -fintrinsic-modules-path $GFC_FINCLUDE"
    else
      REPROFLAGS="-nostdinc"
    fi
    REPROFLAGS="$REPROFLAGS -ffp-contract=off"
    AC_MSG_NOTICE([reproducible build: adding $REPROFLAGS])
    FCFLAGS="$FCFLAGS $REPROFLAGS"
    NOOPTFC="$NOOPTFC $REPROFLAGS"
  else
    AC_MSG_WARN([--enable-reproducible has an effect only with GNU Fortran])
  fi
fi

if test "x$USE_OPENACC" != x0; then
  if test "$fc_vendor" = nvhpc; then
    case $USE_OPENACC in
      1) ACCFLAGS="-acc=gpu -cuda -gpu=$PGI_GPU_ARCH,lineinfo -Minfo=accel" ;;
      2) ACCFLAGS="-acc=gpu -cuda -gpu=$PGI_GPU_ARCH,autocompare -Minfo=accel"
         # autocompare is run without optimisation
         FCFLAGS="$NOOPTFC" ;;
      3) ACCFLAGS="-acc=gpu -cuda \
                   -gpu=$PGI_GPU_ARCH,lineinfo,mem:managed -Minfo=accel" ;;
      4) ACCFLAGS="-stdpar=gpu -cuda \
                   -gpu=$PGI_GPU_ARCH,lineinfo,mem:unified -Minfo=accel" ;;
    esac
    FCFLAGS="$FCFLAGS $ACCFLAGS"
    NOOPTFC="$NOOPTFC $ACCFLAGS"
    LDFLAGS="$LDFLAGS -traceback -cudalib=nvtx,cusparse,curand"
  else
    AC_MSG_WARN([OpenACC options have an effect only with NVIDIA nvfortran])
  fi
fi

RR_PATH_NETCDF
RR_PATH_NETCDF_F90

if test "x$PNETCDF" = x1; then
  if test "x$PNETCDF_IN_NETCDF" = x1; then
    RR_NETCDF4
  else
    RR_PATH_PNETCDF
    if test "x$NETCDF_CDF5" = x1; then
      RR_CDF5
    fi
  fi
else
  if test "x$NETCDF_CDF5" = x1; then
    RR_CDF5
  else
    RR_NETCDF4
  fi
fi

AM_CONDITIONAL(RRTM_NC, [test x"$netcdf" = xyes])

RCM_FC_CHECK_IEEE_ARITHMETIC
RCM_FC_CHECK_QUAD_PRECISION

EXTPATH="${PWD}/external"
AM_CPPFLAGS="-I$EXTPATH $AM_CPPFLAGS"
AM_LDFLAGS="-L$EXTPATH $AM_LDFLAGS"
ACX_MPI
RCM_MPI_CHECK_MPI3

# The model is compiled with MPIFC, the flags were chosen for FC
RCM_FC_VENDOR([$MPIFC], [mpifc_vendor])
if test "$mpifc_vendor" != "$fc_vendor"; then
  AC_MSG_WARN([MPIFC=$MPIFC is $mpifc_vendor but FC=$FC is $fc_vendor:
               the compiler flags are those for $fc_vendor])
fi

if test "x$USEOASIS" = x1; then
  AM_CPPFLAGS="$AM_CPPFLAGS -I$OASIS_PREFIX/include"
  OASIS3_LIBS="-L$OASIS_PREFIX/lib -lpsmile.$OASIS_CHAN -lmct -lmpeu -lscrip"
  LIBS="$LIBS $OASIS3_LIBS"
fi

F90FLAGS="$FCFLAGS"

REGCM_ROOT=${PWD}
AC_PREFIX_DEFAULT(${PWD})

CPPFLAGS="$AM_CPPFLAGS -I$REGCM_ROOT/Share"
if test "x$ASYNC_NETCDF" = x1; then
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([pthread.h], [], [
    AC_MSG_ERROR([pthread.h is required for --enable-async-netcdf])
  ])
  CFLAGS="$CFLAGS -pthread"
  LIBS="$LIBS -pthread"
  AC_LANG_POP([C])
fi
LIBS="-L$REGCM_ROOT/Share -lrcmlib $LIBS"
SYSLIBS="$LIBS"

AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_LDFLAGS)
AC_SUBST(CC)
AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(FCFLAGS)
AC_SUBST(F90FLAGS)
AC_SUBST(NOOPTFC)
AC_SUBST(OMPFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(LIBS)
AC_SUBST(SYSLIBS)
AC_SUBST(REGCM_ROOT)
AC_SUBST(NC_CC)
AC_SUBST(NC_PREFIX)
AC_SUBST(NF_PREFIX)
AC_SUBST(NC_INCLUDES)
AC_SUBST(NC_LDFLAGS)
AC_SUBST(NC_LIBS)

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile external/Makefile Share/Makefile PreProc/Makefile \
                 PreProc/Terrain/Makefile PreProc/ICBC/cmip6/Makefile \
                 PreProc/ICBC/cmip7/Makefile PreProc/ICBC/Makefile \
                 PreProc/CH_EMISSIONS/Makefile \
                 PreProc/CH_EMISSIONS/interpolation/Makefile \
                 PreProc/CH_EMISSIONS/interpolation/tools/Makefile \
                 PreProc/CH_EMISSIONS/interpolation/tools/REGCM_grid/Makefile \
                 PreProc/CLM45/Makefile PostProc/Makefile \
                 Main/mpplib/Makefile Main/netlib/Makefile \
                 Main/Makefile Main/batslib/Makefile Main/clmlib/Makefile \
                 Main/ocnlib/Makefile Main/clmlib/clm4.5/Makefile \
                 Main/chemlib/Makefile Main/chemlib/GAS_CBMZ_NEW/Makefile \
                 Main/cumlib/Makefile Main/microlib/Makefile \
                 Main/radlib/Makefile Main/radlib/RRTMG/Makefile \
                 Main/pbllib/Makefile Main/cloudlib/Makefile \
                 Main/oasislib/Makefile \
                 Testing/Makefile Testing/mpplib/Makefile \
                 Testing/share/Makefile \
                 Tools/Programs/CheckSun/Makefile \
                 Tools/Programs/RegCM_read/Makefile \
                 Tools/Programs/timeseries/Makefile])

AC_OUTPUT
