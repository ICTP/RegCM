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

      module mod_blh_tmp
      use mod_regcm_param

      implicit none

#ifdef MPP1
      real(8) ,allocatable, dimension(:,:,:) :: cgh , kvc , kvh , kvm , kvq
      real(8) ,allocatable, dimension(:,:) :: hfxv , obklen , th10 , ustr ,      &
                                   & xhfx , xqfx
#else
      real(8) , dimension(iy,kz,jxm1) :: cgh , kvc , kvh , kvm , kvq
      real(8) , dimension(iy,jx) :: hfxv , obklen , th10 , ustr ,       &
                                   & xhfx , xqfx
#endif

contains 

      subroutine allocate_mod_blh_tmp

      allocate(cgh(iy,kz,jxp))	
      allocate(kvc(iy,kz,jxp))	
      allocate(kvh(iy,kz,jxp))	
      allocate(kvm(iy,kz,jxp))	
      allocate(kvq(iy,kz,jxp))

	
      allocate(hfxv(iy,jxp))	
      allocate(obklen(iy,jxp))	
      allocate(th10(iy,jxp))	
      allocate(ustr(iy,jxp))	
      allocate(xhfx(iy,jxp))	
      allocate(xqfx(iy,jxp))	


      end  subroutine allocate_mod_blh_tmp

      end module mod_blh_tmp
