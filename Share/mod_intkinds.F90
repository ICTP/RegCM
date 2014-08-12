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

  implicit none

  public

  ! Kind helpers
  integer , parameter :: ik8 = selected_int_kind(R=18)
  integer , parameter :: ik4 = selected_int_kind(R=9)
  integer , parameter :: ik2 = selected_int_kind(R=4)
  integer , parameter :: ik1 = selected_int_kind(R=2)

#ifdef __PGI
  ! quiet nan for portland group compilers
  integer(ik4),  parameter :: bigint = O'17777777777'
#else
  ! signaling nan otherwise
  integer(ik4),  parameter :: bigint = O'17777777777'
#endif

end module mod_intkinds
