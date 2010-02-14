!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_date
      implicit none
!
! COMMON /CTLTIME/
!
      integer :: ldatez , lday , lhour , lmonth , lyear , nnbase ,      &
               & nnnnnn
!
! COMMON /ICONTR/
!
      integer :: nnnchk , nnnend , nstart , nstrt0
!
! COMMON /NOWADAYS/
!
      integer :: nmonth , nyear
!
! COMMON /TIME10/
!
      integer :: mmrec , ndate0 , ndate1
!
! COMMON /UNIFY/
!
      integer :: idate0 , idate1 , idate2 , idatee , idates , idatex

      end module mod_date
