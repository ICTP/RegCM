!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_realkinds

!program test
!  use, intrinsic :: iso_fortran_env
!  use, intrinsic :: ieee_arithmetic
!  integer(int32) i
!  real(real32) x
!  integer(int64) ii
!  real(real64) xx
!
!  x = ieee_value(x, ieee_quiet_nan)
!  i = transfer(x,i)
!  write(*,'(a,i20,a)') '#define __SYSTEM_NAN_32__', i,'_int32'
!  x = ieee_value(x, ieee_positive_inf)
!  i = transfer(x,i)
!  write(*,'(a,i20,a)') '#define __SYSTEM_INF_32__', i,'_int32'
!  xx = ieee_value(xx, ieee_quiet_nan)
!  ii = transfer(xx,ii)
!  write(*,'(a,i20,a)') '#define __SYSTEM_NAN_64__', ii,'_int64'
!  xx = ieee_value(xx, ieee_positive_inf)
!  ii = transfer(xx,ii)
!  write(*,'(a,i20,a)') '#define __SYSTEM_INF_64__', ii,'_int64'
!end program test

#ifdef F2008
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: ieee_arithmetic
#endif

#define __SYSTEM_NAN_32__            -4194304_int32
#define __SYSTEM_INF_32__          2139095040_int32
#define __SYSTEM_NAN_64__   -2251799813685248_int64
#define __SYSTEM_INF_64__ 9218868437227405312_int64

  implicit none

  public

#ifdef F2008
  integer, parameter :: rk4  = real32
  integer, parameter :: rk8  = real64
  integer, parameter :: rk16 = real128
#ifdef NAGFOR
#ifdef SINGLE_PRECISION_REAL
  integer, parameter :: rkx = rk4
  real(rk4) :: nan
  real(rk4) :: inf
#else
  integer, parameter :: rkx = rk8
  real(rk8) :: nan
  real(rk8) :: inf
#endif
#else
#ifdef SINGLE_PRECISION_REAL
  integer, parameter :: rkx = rk4
  real(rk4), parameter :: nan = transfer(__SYSTEM_NAN_32__, 1._real32)
  real(rk4), parameter :: inf = transfer(__SYSTEM_INF_32__, 1._real32)
#else
  integer, parameter :: rkx = rk8
  real(rk8), parameter :: nan = transfer(__SYSTEM_NAN_64__, 1._real64)
  real(rk8), parameter :: inf = transfer(__SYSTEM_INF_64__, 1._real64)
#endif
#endif
#else
  ! Kind helpers as suggested by
  !   Metcalf, M., J. Reid, and M. Cohen
  !   Fortran 95/2003 Explained 2004, p. 71
  integer, parameter :: rk4  = kind(1.0)
  integer, parameter :: rk8  = selected_real_kind(2*precision(1.0_rk4))
  integer, parameter :: rk16 = selected_real_kind(2*precision(1.0_rk8))
  real(rk4), parameter :: inf_r4 = O'07760000000'
  real(rk4), parameter :: nan_r4 = O'07770000000'
  real(rk8), parameter :: inf_r8 = O'0777600000000000000000'
  real(rk8), parameter :: nan_r8 = O'0777700000000000000000'
#ifdef SINGLE_PRECISION_REAL
  integer, parameter :: rkx = rk4
  real(rk4), parameter :: inf = inf_r4
  real(rk4), parameter :: nan = nan_r4
#else
  integer, parameter :: rkx = rk8
  real(rk8), parameter :: inf = inf_r8
  real(rk8), parameter :: nan = nan_r8
#endif
#endif

  interface is_nan
    module procedure is_nan_single
    module procedure is_nan_double
  end interface

  interface is_inf
    module procedure is_inf_single
    module procedure is_inf_double
  end interface

  contains

#ifdef F2008

#ifdef NAGFOR
  subroutine init_realkinds( )
    implicit none
#ifdef SINGLE_PRECISION_REAL
    nan = 1.0_rk4/0.0_rk4
    inf = 2.0_rk4 * huge(1.0_rk4)
#else
    nan = tiny(1.0_rk8)
    inf = huge(1.0_rk8)
#endif
  end subroutine init_realkinds
#endif

  logical elemental function is_nan_double(x)
    implicit none
    real(rk8), intent(in) :: x
    is_nan_double = (ieee_class(x) == ieee_quiet_nan .or. &
                     ieee_class(x) == ieee_signaling_nan)
  end function is_nan_double

  logical elemental function is_inf_double(x)
    implicit none
    real(rk8), intent(in) :: x
    is_inf_double = (ieee_class(x) == ieee_negative_inf .or. &
                     ieee_class(x) == ieee_positive_inf)
  end function is_inf_double

  logical elemental function is_nan_single(x)
    implicit none
    real(rk4), intent(in) :: x
    is_nan_single = (ieee_class(x) == ieee_quiet_nan .or. &
                     ieee_class(x) == ieee_signaling_nan)
  end function is_nan_single

  logical elemental function is_inf_single(x)
    implicit none
    real(rk4), intent(in) :: x
    is_inf_single = (ieee_class(x) == ieee_negative_inf .or. &
                     ieee_class(x) == ieee_positive_inf)
  end function is_inf_single

#else

  logical elemental function is_nan_double(x)
    implicit none
    real(rk8), intent(in) :: x
    is_nan_double = ( (x /= x) .or. ((x > 0.0D0) .eqv. (x <= 0.0D0)) )
  end function is_nan_double

  logical elemental function is_inf_double(x)
    implicit none
    real(rk8), intent(in) :: x
    is_inf_double = ( x > huge(x) )
  end function is_inf_double

  logical elemental function is_nan_single(x)
    implicit none
    real(rk4), intent(in) :: x
    is_nan_single = ( (x /= x) .or. ((x > 0.0) .eqv. (x <= 0.0)) )
  end function is_nan_single

  logical elemental function is_inf_single(x)
    implicit none
    real(rk4), intent(in) :: x
    is_inf_single = ( x > huge(x) )
  end function is_inf_single

#endif

end module mod_realkinds

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
