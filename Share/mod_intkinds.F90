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

module mod_intkinds
#ifdef F2008
  use, intrinsic :: iso_fortran_env
#endif

  implicit none

  public

#define __SYSTEM_INTMAX_32__          2147483647_int32
#define __SYSTEM_INTMAX_64__ 9223372036854775807_int64

  ! Kind helpers

#ifdef F2008
  integer, parameter :: ik8 = int64
  integer, parameter :: ik4 = int32
  integer, parameter :: ik2 = int16
  integer, parameter :: ik1 = int8
  integer(ik4),  parameter :: bigint = __SYSTEM_INTMAX_32__
#else
  integer, parameter :: ik8 = selected_int_kind(R=18)
  integer, parameter :: ik4 = selected_int_kind(R=9)
  integer, parameter :: ik2 = selected_int_kind(R=4)
  integer, parameter :: ik1 = selected_int_kind(R=2)
  integer(ik4),  parameter :: bigint = 2147483647
#endif

end module mod_intkinds
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
