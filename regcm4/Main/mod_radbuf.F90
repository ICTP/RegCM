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

      module mod_radbuf

      use mod_regcm_param

      implicit none

      ! absnxt  - Nearest layer absorptivities
      ! abstot  - Non-adjacent layer absorptivites
      ! emstot  - Total emissivity

#ifdef MPP1
      real(8) , allocatable, dimension(:,:,:,:)  :: absnxt
      real(8),  allocatable, dimension(:,:,:,:)  :: abstot
      real(8) , allocatable, dimension(:,:,:) :: emstot
#else
      real(8) , dimension(iym1,kz,4,jxm1) :: absnxt
      real(8) , dimension(iym1,kzp1,kz + 1,jxm1) :: abstot
      real(8) , dimension(iym1,kzp1,jxm1) :: emstot
#endif

contains 
	subroutine allocate_mod_radbuf 
	
#ifdef MPP1
	allocate(absnxt(iym1,kz,4,jxp))
	allocate(abstot(iym1,kzp1,kz + 1,jxp))
	allocate(emstot(iym1,kzp1,jxp))
#else
#endif 
	end subroutine allocate_mod_radbuf 

      end module mod_radbuf

