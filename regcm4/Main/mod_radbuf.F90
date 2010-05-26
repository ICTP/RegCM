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

      use mod_dynparam

      implicit none

      ! absnxt  - Nearest layer absorptivities
      ! abstot  - Non-adjacent layer absorptivites
      ! emstot  - Total emissivity

      real(8) , allocatable, dimension(:,:,:,:)  :: absnxt,absnxt0
      real(8),  allocatable, dimension(:,:,:,:)  :: abstot,abstot0
      real(8) , allocatable, dimension(:,:,:) :: emstot,emstot0
      real(8), allocatable, dimension(:,:,:,:):: xuinpl
      contains 

        subroutine allocate_mod_radbuf 
        implicit none        
#ifdef MPP1
        allocate(absnxt(iym1,kz,4,jxp))
        allocate(abstot(iym1,kzp1,kzp1,jxp))
        allocate(emstot(iym1,kzp1,jxp))
        allocate(absnxt0(iym1,kz,4,jxp))
        allocate(abstot0(iym1,kzp1,kzp1,jxp))
        allocate(emstot0(iym1,kzp1,jxp))
        allocate(xuinpl(iym1,kzp1,4,jxp))
#else
        allocate(absnxt(iym1,kz,4,jxm1))
        allocate(abstot(iym1,kzp1,kzp1,jxm1))
        allocate(emstot(iym1,kzp1,jxm1))
        allocate(absnxt0(iym1,kz,4,jxm1))
        allocate(abstot0(iym1,kzp1,kzp1,jxm1))
        allocate(emstot0(iym1,kzp1,jxm1))
        allocate(xuinpl(iym1,kzp1,4,jxm1))       

#endif 
        end subroutine allocate_mod_radbuf 

      end module mod_radbuf

