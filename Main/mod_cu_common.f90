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

module mod_cu_common
!
! Storage and constants related to cumulus convection schemes.
!
  use mod_dynparam
  use mod_memutil

  implicit none
!
  real(8) :: clfrcv ! Cloud fractional cover for convective precip
  real(8) :: cllwcv ! Cloud liquid water content for convective precip.

  integer , pointer , dimension(:) :: icon ! Precip points counter
!
  contains
!
  subroutine allocate_mod_cu_common
  implicit none
  call getmem1d(icon,1,jxp,'mod_cu_common:icon')
  end subroutine allocate_mod_cu_common
!
end module mod_cu_common
