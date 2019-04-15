!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_realkinds

#ifdef F2008
  use , intrinsic :: iso_fortran_env
  use , intrinsic :: ieee_arithmetic
#endif

  implicit none

  public

#ifdef F2008
  integer , parameter :: rk4  = REAL32
  integer , parameter :: rk8  = REAL64
  integer , parameter :: rk16 = REAL128
#else
  ! Kind helpers as suggested by
  !   Metcalf, M., J. Reid, and M. Cohen
  !   Fortran 95/2003 Explained 2004, p. 71
  integer , parameter :: rk4  = kind(1.0)
  integer , parameter :: rk8  = selected_real_kind(2*precision(1.0_rk4))
  integer , parameter :: rk16 = selected_real_kind(2*precision(1.0_rk8))
#endif

#ifdef SINGLE_PRECISION_REAL
  integer , parameter :: rkx = rk4
#ifdef __PGI
  ! quiet nan for portland group compilers
  real(rk4), parameter :: inf = O'07760000000'
  real(rk4), parameter :: nan = O'07770000000'
#else
  ! signaling nan otherwise
  real(rk4), parameter :: inf = O'07760000000'
  real(rk4), parameter :: nan = O'07761000000'
#endif
#else
  integer , parameter :: rkx = rk8
#ifdef __PGI
  ! quiet nan for portland group compilers
  real(rk8), parameter :: inf = O'0777600000000000000000'
  real(rk8), parameter :: nan = O'0777700000000000000000'
#else
  ! signaling nan otherwise
  real(rk8), parameter :: inf = O'0777600000000000000000'
  real(rk8), parameter :: nan = O'0777610000000000000000'
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

  logical elemental function is_nan_double(x)
    implicit none
    real(rk8) , intent(in) :: x
    is_nan_double = (ieee_class(x) == ieee_quiet_nan .or. &
                     ieee_class(x) == ieee_signaling_nan)
  end function is_nan_double

  logical elemental function is_inf_double(x)
    implicit none
    real(rk8) , intent(in) :: x
    is_inf_double = (ieee_class(x) == ieee_negative_inf .or. &
                     ieee_class(x) == ieee_positive_inf)
  end function is_inf_double

  logical elemental function is_nan_single(x)
    implicit none
    real(rk4) , intent(in) :: x
    is_nan_single = (ieee_class(x) == ieee_quiet_nan .or. &
                     ieee_class(x) == ieee_signaling_nan)
  end function is_nan_single

  logical elemental function is_inf_single(x)
    implicit none
    real(rk4) , intent(in) :: x
    is_inf_single = (ieee_class(x) == ieee_negative_inf .or. &
                     ieee_class(x) == ieee_positive_inf)
  end function is_inf_single

#else

  logical elemental function is_nan_double(x)
    implicit none
    real(rk8) , intent(in) :: x
    is_nan_double = ( (x /= x) .or. ((x > 0.0D0) .eqv. (x <= 0.0D0)) )
  end function is_nan_double

  logical elemental function is_inf_double(x)
    implicit none
    real(rk8) , intent(in) :: x
    is_inf_double = ( x > huge(x) )
  end function is_inf_double

  logical elemental function is_nan_single(x)
    implicit none
    real(rk4) , intent(in) :: x
    is_nan_single = ( (x /= x) .or. ((x > 0.0) .eqv. (x <= 0.0)) )
  end function is_nan_single

  logical elemental function is_inf_single(x)
    implicit none
    real(rk4) , intent(in) :: x
    is_inf_single = ( x > huge(x) )
  end function is_inf_single

#endif

end module mod_realkinds

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
