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

module mod_intkinds
#ifdef F2008
  use , intrinsic :: iso_fortran_env
#endif

  implicit none

  public

#define __SYSTEM_INTMAX_32__          2147483647_int32
#define __SYSTEM_INTMAX_64__ 9223372036854775807_int64

  ! Kind helpers
  integer , parameter :: ik8 = selected_int_kind(R=18)
  integer , parameter :: ik4 = selected_int_kind(R=9)
  integer , parameter :: ik2 = selected_int_kind(R=4)
  integer , parameter :: ik1 = selected_int_kind(R=2)

#ifdef F2008
  integer(ik4),  parameter :: bigint = __SYSTEM_INTMAX_32__
#else
  integer(ik4),  parameter :: bigint = 2147483647
#endif

end module mod_intkinds
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
