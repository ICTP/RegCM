dnl
dnl autoconf macro for detecting NetCDF
dnl

AC_DEFUN([AX_PROG_NC_CONFIG], [
  AC_REQUIRE([AC_PROG_EGREP])

  AC_CACHE_CHECK([if nc-config program is present],[ax_cv_prog_nc_config],[
  AS_IF([nc-config --version 2>/dev/null | egrep -q '^netCDF '],
        [ax_cv_prog_nc_config=yes], [ax_cv_prog_nc_config=no])
      ])
  AS_IF([test "$ax_cv_prog_nc_config" = "yes"], [[$1]], [[$2]])
])

AC_DEFUN([AX_PROG_NF_CONFIG], [
  AC_REQUIRE([AC_PROG_EGREP])

  AC_CACHE_CHECK([if nf-config program is present],[ax_cv_prog_nf_config],[
  AS_IF([nf-config --version 2>/dev/null | egrep -q '.*'],
        [ax_cv_prog_nf_config=yes], [ax_cv_prog_nf_config=no])
      ])
  AS_IF([test "$ax_cv_prog_nf_config" = "yes"], [[$1]], [[$2]])
])

AC_DEFUN([RR_PATH_NETCDF],[

  AC_CHECKING([for NetCDF])

  save_CPPFLAGS="$CPPFLAGS"
  save_LDFLAGS="$LDFLAGS"

  CPPFLAGS="$CPPFLAGS $NC_INCLUDES"
  AMDEPFLAGS="$AMDEPFLAGS $NC_INCLUDES"
  LIBS="$LIBS $NC_LIBS"
  LDFLAGS="$LDFLAGS $NC_LDFLAGS"

  AC_SUBST([AMDEPFLAGS])

  netcdf=no

  AC_LANG_PUSH([C])
  AC_CHECKING([for netcdf.h])
  AC_CHECK_HEADER([netcdf.h],
                  [netcdf=yes], [netcdf=no])

  if test "x$netcdf" = xno; then
      AC_MSG_ERROR([NetCDF include not found])
  fi

  AC_CHECKING([for libnetcdf.a])
  AC_CHECK_LIB([netcdf], [nc_close],
               [netcdf=yes], [netcdf=no])
  if test "x$netcdf" = xno; then
    AC_CHECKING([if we need to link jpeg library for HDF4])
    LDFLAGS="$LDFLAGS $NC_LDFLAGS"
    LIBS="$LIBS -ljpeg"
    AC_CHECK_LIB([netcdf], [nc_enddef],
    [netcdf=yes], [netcdf=no])
    if test "x$netcdf" = xno; then
      AC_CHECKING([if we need to link hdf5 library])
      NC_LDFLAGS="$NC_LDFLAGS -L$HDF5_PREFIX/lib"
      LDFLAGS="$LDFLAGS $NC_LDFLAGS"
      LIBS="$LIBS -lhdf5 -lhdf5_hl"
      AC_CHECK_LIB([netcdf], [nc_sync],
                   [netcdf=yes], [netcdf=no])
      if test "x$netcdf" = xno; then
        AC_CHECKING([if we need to link szlib library])
        NC_LDFLAGS="$NC_LDFLAGS -L$SZIP_PREFIX/lib"
        LDFLAGS="$LDFLAGS $NC_LDFLAGS"
        LIBS="$LIBS -lsz"
        AC_CHECK_LIB([netcdf], [nc_inq_libvers],
                     [netcdf=yes], [netcdf=no])
      fi
    fi
    if test "x$netcdf" = xno; then
      AC_MSG_ERROR([NetCDF library not found])
    fi
  fi

# Put them back to how they used to be and set the AM versions
# The AM versions must be substituted explicitly

  CPPFLAGS="$save_CPPFLAGS"
  LDFLAGS="$save_LDFLAGS"
  AM_CPPFLAGS="$NC_INCLUDES $AM_CPPFLAGS"
  AM_LDFLAGS="$NC_LDFLAGS $AM_LDFLAGS"
  AC_SUBST([AM_CPPFLAGS])
  AC_SUBST([AM_LDFLAGS])
  AC_LANG_POP([C])

# Netcdf Fortran interface can be placed in a separate libnetcdff
#
  LDFLAGS="$LDFLAGS $NC_LDFLAGS"
  AC_CHECKING([for libnetcdff.a])
  AC_CHECK_LIB([netcdf], [nf_close],
               [netcdf=yes], [netcdf=no])
  if test "x$netcdf" = xno; then
    LIBS="-lnetcdff $LIBS"
    AC_CHECKING([for libnetcdff.a])
    AC_CHECK_LIB([netcdff], [nf_close],
                 [netcdf=yes], [netcdf=no])
    if test "x$netcdf" = xno; then
      AC_MSG_ERROR([NetCDF library not found])
    fi
  fi
  LDFLAGS="$save_LDFLAGS"

])

dnl
dnl autoconf macro for detecting NetCDF module file
dnl 

AC_DEFUN([RR_PATH_NETCDF_F90],[

  AC_CHECKING([for NetCDF module file])
  save_FCFLAGS="$FCFLAGS"

  if test -z "$NC_INCLUDES"; then
    for flag in "-I" "-M" "-p"; do
      FCFLAGS="$flag$NC_PREFIX/include $save_FCFLAGS"
      AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([[ ]],
                         [[      use netcdf]])],
                         [netcdf=yes; NC_FCFLAGS=$flag],
                         [netcdf=no])
      if test "x$netcdf" = xyes; then
        break
      fi
    done
    if test "x$netcdf" = xno; then
      AC_MSG_ERROR([NetCDF module not found])
    fi

    FCFLAGS="$save_FCFLAGS"
    AM_CPPFLAGS="$NC_FCFLAGS$NC_PREFIX/include $AM_CPPFLAGS"
  else
    FCFLAGS="$NC_INCLUDES $save_FCFLAGS"
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([[ ]],
                       [[      use netcdf]])],
                       [netcdf=yes],
                       [netcdf=no])

    if test "x$netcdf" = xno; then
      AC_MSG_ERROR([NetCDF module not found])
    fi

    FCFLAGS="$save_FCFLAGS"
    AM_CPPFLAGS="$NC_INCLUDES $AM_CPPFLAGS"
  fi

  AC_SUBST([AM_CPPFLAGS])
])

AC_DEFUN([RR_CDF5],[
  AC_CHECKING([for NetCDF CDF5])
  save_FCFLAGS="$FCFLAGS"

  if test -z "$NC_INCLUDES"; then
    for flag in "-I" "-M" "-p"; do
      FCFLAGS="$flag$NC_PREFIX/include $save_FCFLAGS"
      AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([[ ]],
                         [[
         use netcdf
         implicit none
         integer :: imode
         imode = nf90_cdf5]])],
                         [cdf5=yes],
                         [cdf5=no])
      if test "x$cdf5" = xyes; then
        AM_CPPFLAGS="-DNETCDF_CDF5 $AM_CPPFLAGS"
        break
      fi
    done
    FCFLAGS="$save_FCFLAGS"
  else
    FCFLAGS="$NC_INCLUDES $save_FCFLAGS"
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([[ ]],
                       [[
         use netcdf
         implicit none
         integer :: imode
         imode = nf90_cdf5]])],
                       [cdf5=yes],
                       [cdf5=no])

    if test "x$cdf5" = xyes; then
      AM_CPPFLAGS="-DNETCDF_CDF5 $AM_CPPFLAGS"
    fi
    FCFLAGS="$save_FCFLAGS"
  fi

  AC_SUBST([AM_CPPFLAGS])
])

dnl @synopsis ACX_MPI([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl @summary figure out how to compile/link code with MPI
dnl
dnl This macro tries to find out how to compile programs that use MPI
dnl (Message Passing Interface), a standard API for parallel process
dnl communication (see http://www-unix.mcs.anl.gov/mpi/)
dnl
dnl On success, it sets the MPICC, MPICXX, or MPIF77 output variable to
dnl the name of the MPI compiler, depending upon the current language.
dnl (This may just be $CC/$CXX/$F77, but is more often something like
dnl mpicc/mpiCC/mpif77.) It also sets MPILIBS to any libraries that are
dnl needed for linking MPI (e.g. -lmpi, if a special
dnl MPICC/MPICXX/MPIF77 was not found).
dnl
dnl If you want to compile everything with MPI, you should set:
dnl
dnl     CC="$MPICC" #OR# CXX="$MPICXX" #OR# F77="$MPIF77"
dnl     LIBS="$MPILIBS $LIBS"
dnl
dnl NOTE: The above assumes that you will use $CC (or whatever) for
dnl linking as well as for compiling. (This is the default for automake
dnl and most Makefiles.)
dnl
dnl The user can force a particular library/compiler by setting the
dnl MPICC/MPICXX/MPIF77 and/or MPILIBS environment variables.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if an MPI
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands to
dnl run it if it is not found. If ACTION-IF-FOUND is not specified, the
dnl default action will define HAVE_MPI.
dnl
dnl @category InstalledPackages
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl @author Julian Cummings <cummings@cacr.caltech.edu>
dnl @version 2006-10-13
dnl @license GPLWithACException

AC_DEFUN([ACX_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

AC_LANG_CASE([C], [
	AC_REQUIRE([AC_PROG_CC])
	AC_ARG_VAR(MPICC,[MPI C compiler command])
	AC_CHECK_PROGS(MPICC, mpicc hcc mpxlc_r mpxlc mpcc cmpicc, $CC)
	acx_mpi_save_CC="$CC"
	CC="$MPICC"
	AC_SUBST(MPICC)
],
[C++], [
	AC_REQUIRE([AC_PROG_CXX])
	AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
	AC_CHECK_PROGS(MPICXX, mpic++ mpicxx mpiCC hcp mpxlC_r mpxlC mpCC cmpic++, $CXX)
	acx_mpi_save_CXX="$CXX"
	CXX="$MPICXX"
	AC_SUBST(MPICXX)
],
[Fortran 77], [
	AC_REQUIRE([AC_PROG_F77])
	AC_ARG_VAR(MPIF77,[MPI Fortran 77 compiler command])
	AC_CHECK_PROGS(MPIF77, mpif77 hf77 mpxlf mpf77 mpif90 mpf90 mpixlf90 mpixlf95 mpixlf_r cmpifc cmpif90c, $F77)
	acx_mpi_save_F77="$F77"
	F77="$MPIF77"
	AC_SUBST(MPIF77)
],
[Fortran], [
	AC_REQUIRE([AC_PROG_FC])
	AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
	AC_CHECK_PROGS(MPIFC, mpifort mpiifort mpixfl2003 mpixlf2008 mpf90 cmpifc cmpif90c hf90 mpif90, $FC)
	acx_mpi_save_FC="$FC"
	FC="$MPIFC"
	AC_SUBST(MPIFC)
])

if test x = x"$MPILIBS"; then
	AC_LANG_CASE([C], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[C++], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[Fortran 77], [AC_MSG_CHECKING([for MPI_Init])
			AC_LINK_IFELSE([AC_LANG_PROGRAM([],[      call MPI_Init])],[MPILIBS=" "
				AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])],
		[Fortran], [AC_MSG_CHECKING([for MPI_Init])
			AC_LINK_IFELSE([AC_LANG_PROGRAM([],[      call MPI_Init])],[MPILIBS=" "
				AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])])
fi
AC_LANG_CASE([Fortran 77], [
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
	fi
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpich, MPI_Init, [MPILIBS="-lfmpich"])
	fi
],
[Fortran], [
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
	fi
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(mpichf90, MPI_Init, [MPILIBS="-lmpichf90"])
	fi
])
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
fi
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
fi

dnl We have to use AC_TRY_COMPILE and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
AC_LANG_CASE([C], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[C++], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[Fortran 77], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpif.h])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[      include 'mpif.h'])],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[Fortran], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpif.h])
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[      include 'mpif.h'])],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi])

AC_LANG_CASE([C], [CC="$acx_mpi_save_CC"],
	[C++], [CXX="$acx_mpi_save_CXX"],
	[Fortran 77], [F77="$acx_mpi_save_F77"],
	[Fortran], [FC="$acx_mpi_save_FC"])

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])dnl ACX_MPI

dnl
dnl autoconf macro for detecting pNetCDF
dnl

AC_DEFUN([AX_PROG_PNETCDF_CONFIG], [
  AC_REQUIRE([AC_PROG_EGREP])

  AC_CACHE_CHECK([if pnetcdf-config program is present],[ax_cv_prog_pnetcdf_config],[
  AS_IF([pnetcdf-config --version 2>/dev/null | egrep -q '^PnetCDF '],
        [ax_cv_prog_pnetcdf_config=yes], [ax_cv_prog_pnetcdf_config=no])
      ])
  AS_IF([test "$ax_cv_prog_pnetcdf_config" = "yes"], [[$1]], [[$2]])
])

AC_DEFUN([RR_PATH_PNETCDF],[

  AC_CHECKING([for Parallel NetCDF])

  save_CPPFLAGS="$CPPFLAGS"
  save_FCFLAGS="$FCFLAGS"
  save_LDFLAGS="$LDFLAGS"

  CPPFLAGS="$CPPFLAGS $NC_INCLUDES"
  FCFLAGS="$CPPFLAGS $NC_INCLUDES"
  AMDEPFLAGS="$AMDEPFLAGS $NC_INCLUDES"
  LIBS="$LIBS $NC_LIBS"
  LDFLAGS="$LDFLAGS $NC_LDFLAGS"

  AC_SUBST([AMDEPFLAGS])

  pnetcdf=no

  AC_LANG_PUSH([C])
  AC_CHECKING([for pnetcdf.h])
  AC_CHECK_HEADER([pnetcdf.h],
                  [pnetcdf=yes], [pnetcdf=no])

  if test "x$pnetcdf" = xno; then
      AC_MSG_ERROR([Parallel NetCDF include not found])
  fi

  AC_CHECKING([for libpnetcdf.a])
  AC_CHECK_LIB([pnetcdf], [ncmpi_close],
               [pnetcdf=yes], [pnetcdf=no])
  if test "x$pnetcdf" = xno; then
    AC_MSG_ERROR([Parallel NetCDF library not found])
  fi

  AC_LANG_POP([C])

  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([[ ]],
                      [[
         use pnetcdf
         implicit none
         integer :: imode
         imode = nf90_cdf5]])],
                     [cdf5=yes],
                     [cdf5=no])
  if test "x$cdf5" = xyes; then
    AM_CPPFLAGS="-DNETCDF_CDF5 $AM_CPPFLAGS"
  fi

# Put them back to how they used to be and set the AM versions
# The AM versions must be substituted explicitly

  CPPFLAGS="$save_CPPFLAGS"
  FCFLAGS="$save_FCFLAGS"
  LDFLAGS="$save_LDFLAGS"
  AM_CPPFLAGS="$NC_INCLUDES $AM_CPPFLAGS"
  AM_LDFLAGS="$NC_LDFLAGS $AM_LDFLAGS"
  AC_SUBST([AM_CPPFLAGS])
  AC_SUBST([AM_LDFLAGS])
])

AC_DEFUN([RCM_FC_CHECK_IEEE_ARITHMETIC],[
  # Init
  fc_has_ieee_arithmetic="no"
  AC_MSG_CHECKING([whether the Fortran compiler supports IEEE_ARITHMETIC])
  # Try to compile a piece of code that uses the module.
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
      use, intrinsic :: ieee_arithmetic
      real :: val
      if (ieee_is_nan(val)) then  ! NaN
        write(*,*)"Hello NAN"
      end if
    ]])], [fc_has_ieee_arithmetic="yes"])
  AC_LANG_POP([Fortran])
  if test "${fc_has_ieee_arithmetic}" = "yes"; then
    AC_DEFINE([HAVE_FC_IEEE_ARITHMETIC],1, 
      [Define to 1 if your Fortran compiler supports IEEE_ARITHMETIC module.])
    FCFLAGS="-DF2008 $FCFLAGS"
    AC_SUBST(FCFLAGS)
  fi
  AC_MSG_RESULT([${fc_has_ieee_arithmetic}])
]) # RCM_FC_CHECK_IEEE_ARITHMETIC
