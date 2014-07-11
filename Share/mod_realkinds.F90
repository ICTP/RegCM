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
!
  public

  ! Kind helpers
  integer , parameter :: rk16 = selected_real_kind(P=27,R=2400)
  integer , parameter :: rk8 = selected_real_kind(P=13,R=300)
  integer , parameter :: rk4 = selected_real_kind(P= 6,R=37)

#ifdef __PGI
  ! quiet nan for portland group compilers
  real(rk8), parameter :: inf = O'0777600000000000000000'
  real(rk8), parameter :: nan = O'0777700000000000000000'
#else
  ! signaling nan otherwise
  real(rk8), parameter :: inf = O'0777600000000000000000'
  real(rk8), parameter :: nan = O'0777610000000000000000'
#endif

  contains

  logical elemental function is_nan(x)
    implicit none
    real(rk8) , intent(in) :: x
    is_nan = ( (x /= x) .or. ((x > 0.0D0) .eqv. (x <= 0.0D0)) )
  end function is_nan

  logical elemental function is_inf(x)
    implicit none
    real(rk8) , intent(in) :: x
    is_inf = ( x > huge(x) )
  end function is_inf

end module mod_realkinds
