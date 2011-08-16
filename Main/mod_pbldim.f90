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

!
! Storage parameters and constants related to
!     the boundary layer
!
      use mod_constants
      use mod_dynparam
      use mod_memutil
!
      private

      real(8) , pointer , dimension(:,:) :: zq
      real(8) , pointer ,  dimension(:,:) :: rhox2d
      real(8) , pointer , dimension(:,:,:) :: dzq , thvx
      real(8) , pointer , dimension(:,:,:) :: za
      real(8) , pointer , dimension(:,:) :: kpbl
!
      public :: allocate_mod_pbldim
      public :: zq
      public :: za , dzq , thvx
      public :: kpbl
      public :: rhox2d
!
      contains
!
      subroutine allocate_mod_pbldim
        implicit none
!        
        call getmem2d(zq,1,iy,1,kzp1,'pbldim:zq')
        call getmem3d(dzq,1,iy,1,kz,1,jxp,'pbldim:dzq')
        call getmem3d(thvx,1,iy,1,kz,1,jxp,'pbldim:thvx')
        call getmem3d(za,1,iy,1,kz,1,jxp,'pbldim:za')
        call getmem2d(rhox2d,1,iy,1,jxp,'pbldim:rhox2d')
        call getmem2d(kpbl,1,iy,1,jxp,'pbldim:kpbl')

      end subroutine allocate_mod_pbldim
!
      end module mod_pbldim
