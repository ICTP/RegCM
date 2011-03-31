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
      real(8) , pointer , dimension(:,:,:) :: dzq , thvx , thx3d
      real(8) , pointer , dimension(:,:,:) :: za
!
      public :: allocate_mod_pbldim
      public :: zq
      public :: za , dzq , thvx , thx3d
      public :: rhox2d
!
      contains
!
      subroutine allocate_mod_pbldim(lmpi)
        implicit none
        logical , intent(in) :: lmpi
!        
        call getmem2d(zq,iy,kzp1,'pbldim:zq')
        if (lmpi) then
          call getmem3d(dzq,iy,kz,jxp,'pbldim:dzq')
          call getmem3d(thvx,iy,kz,jxp,'pbldim:thvx')
          call getmem3d(thx3d,iy,kz,jxp,'pbldim:thx3d')
          call getmem3d(za,iy,kz,jxp,'pbldim:za')
          call getmem2d(rhox2d,iy,jxp,'pbldim:rhox2d')
        else
          call getmem3d(dzq,iy,kz,jx,'pbldim:dzq')
          call getmem3d(thvx,iy,kz,jx,'pbldim:thvx')
          call getmem3d(thx3d,iy,kz,jx,'pbldim:thx3d')
          call getmem3d(za,iy,kz,jx,'pbldim:za')
          call getmem2d(rhox2d,iy,jx,'pbldim:rhox2d')
        end if
      end subroutine allocate_mod_pbldim
!
      end module mod_pbldim
