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

      use mod_regcm_param

      implicit none
!
      real(8) , dimension(iy,kzp1) :: zq

#ifdef MPP1
      real(8) ,allocatable, dimension(:,:,:) :: dzq , thvx , thx3d , za
      real(8) ,allocatable,  dimension(:,:) :: rhox2d
#else
      real(8) , dimension(iy,kz,jx) :: dzq , thvx , thx3d , za
      real(8) , dimension(iy,jx) :: rhox2d
#endif
contains
     subroutine allocate_mod_pbldim
	
	allocate(dzq(iy,kz,jxp))
	allocate(thvx(iy,kz,jxp))
	allocate(thx3d(iy,kz,jxp))
	allocate(za(iy,kz,jxp))

	allocate(rhox2d(iy,jxp))

     end subroutine allocate_mod_pbldim
      end module mod_pbldim
