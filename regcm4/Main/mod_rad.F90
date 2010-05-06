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

      module mod_rad

      use mod_dynparam

      implicit none

      real(8) , allocatable, dimension(:,:) :: cldfra , cldlwc
      real(8) , allocatable, dimension(:,:,:) :: heatrt
      real(8) , allocatable, dimension(:,:,:) :: o3prof

      contains 

        subroutine allocate_mod_rad
        implicit none

        allocate(cldfra(iym1,kz))
        allocate(cldlwc(iym1,kz))
#ifdef MPP1
        allocate(heatrt(iym1,kz,jxp))
        allocate(o3prof(iym1,kzp1,jxp))
#else
        allocate(heatrt(iym1,kz,jxm1))
        allocate(o3prof(iym1,kzp1,jxm1))
#endif 
        end subroutine  allocate_mod_rad

      end module mod_rad
