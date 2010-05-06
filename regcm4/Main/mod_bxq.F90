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

      module mod_bxq

      use mod_dynparam

      implicit none
!
      real(8) ,allocatable, dimension(:,:,:) :: ddsum
      real(8) ,allocatable, dimension(:,:,:,:) :: deld
      real(8) ,allocatable, dimension(:,:,:,:) :: delh
      real(8) ,allocatable, dimension(:,:,:) :: dhsum
      real(8) ,allocatable, dimension(:,:) :: psdot
      real(8) ,allocatable, dimension(:,:,:) :: work
      real(8) ,allocatable, dimension(:,:) :: uu , vv
      real(8) ,allocatable, dimension(:,:,:) :: uuu , vvv

      contains 

        subroutine allocate_mod_bxq
        implicit none
#ifdef MPP1
        allocate(ddsum(iy,jxp,nsplit))
        allocate(deld(iy,jxp,nsplit,3))
        allocate(delh(iy,0:jxp,nsplit,3))
        allocate(dhsum(iy,0:jxp,nsplit))
        allocate(psdot(iy,jxp))
        allocate(work(iy,jxp,3))
        allocate(uu(iy,jxp+1))
        allocate(vv(iy,jxp+1))
        allocate(uuu(iy,kz,jxp+1))
        allocate(vvv(iy,kz,jxp+1))
#else
        allocate(ddsum(iy,jx,nsplit))
        allocate(deld(iy,jx,nsplit,3))
        allocate(delh(iy,jx,nsplit,3))
        allocate(dhsum(iy,jx,nsplit))
        allocate(psdot(iy,jx))
        allocate(work(iy,jx,3))
        allocate(uu(iy,jx))
        allocate(vv(iy,jx))
        allocate(uuu(iy,kz,jx))
        allocate(vvv(iy,kz,jx))
#endif 
        end subroutine allocate_mod_bxq

      end module mod_bxq
