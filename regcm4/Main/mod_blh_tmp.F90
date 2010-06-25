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
      use mod_dynparam

      implicit none

      real(8) ,allocatable, dimension(:,:,:) :: cgh , kvc , kvh , kvm , &
                                             &  kvq
      real(8) ,allocatable, dimension(:,:) :: hfxv , obklen , th10 ,    &
                                             & ustr , xhfx , xqfx

      contains 

      subroutine allocate_mod_blh_tmp
      implicit none
#ifdef MPP1
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

#else 

#ifdef BAND
      allocate(cgh(iy,kz,jx))        
      allocate(kvc(iy,kz,jx))        
      allocate(kvh(iy,kz,jx))        
      allocate(kvm(iy,kz,jx))        
      allocate(kvq(iy,kz,jx))
#else
      allocate(cgh(iy,kz,jxm1))        
      allocate(kvc(iy,kz,jxm1))        
      allocate(kvh(iy,kz,jxm1))        
      allocate(kvm(iy,kz,jxm1))        
      allocate(kvq(iy,kz,jxm1))
#endif

      allocate(hfxv(iy,jx))        
      allocate(obklen(iy,jx))        
      allocate(th10(iy,jx))        
      allocate(ustr(iy,jx))        
      allocate(xhfx(iy,jx))        
      allocate(xqfx(iy,jx))        

#endif 

      end  subroutine allocate_mod_blh_tmp

      end module mod_blh_tmp
