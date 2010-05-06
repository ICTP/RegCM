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

      module mod_outrad
      use mod_dynparam
      implicit none

      real(4) ,allocatable, dimension(:,:,:) :: frad2d
      real(4) ,allocatable, dimension(:,:,:,:) :: frad3d

      contains

      subroutine allocate_mod_outrad
        implicit none
#ifdef MPP1
        allocate(frad2d(jxp,iym2,nrad2d))
        allocate(frad3d(jxp,iym2,kz,nrad3d))
#else
        allocate(frad2d(jxm2,iym2,nrad2d))
        allocate(frad3d(jxm2,iym2,kz,nrad3d))
#endif

        end subroutine allocate_mod_outrad

      end module mod_outrad
