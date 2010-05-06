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

      module mod_pbldim

      use mod_dynparam

      implicit none
!
      real(8) , allocatable , dimension(:,:) :: zq
      real(8) , allocatable ,  dimension(:,:) :: rhox2d
      real(8) , allocatable , dimension(:,:,:) :: dzq , thvx , thx3d
      real(8) , allocatable , dimension(:,:,:) :: za

      contains

      subroutine allocate_mod_pbldim
        
#ifdef MPP1
        allocate(dzq(iy,kz,jxp))
        allocate(thvx(iy,kz,jxp))
        allocate(thx3d(iy,kz,jxp))
        allocate(za(iy,kz,jxp))
        allocate(rhox2d(iy,jxp))
#else 
        allocate(dzq(iy,kz,jx))
        allocate(thvx(iy,kz,jx))
        allocate(thx3d(iy,kz,jx))
        allocate(za(iy,kz,jx))
        allocate(rhox2d(iy,jx))
#endif 
        allocate(zq(iy,kzp1))

      end subroutine allocate_mod_pbldim

      end module mod_pbldim
