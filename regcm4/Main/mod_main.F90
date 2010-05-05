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

      module mod_main

      use mod_regcm_param

      implicit none
!
#ifdef MPP1
      real(8) , allocatable, dimension(:,:) :: cldefi , f , hfx , hgfact , htsd , &
                                   & qfx , rainc , rainnc , tgb , tgbb ,&
                                   & xlat , xlong , zpbl
      real(8) , allocatable, dimension(:,:) :: ht
      real(8) , allocatable, dimension(:,:) :: msfd , msfx , pdotb , psa ,   &
           & rainp
      real(8) , allocatable, dimension(:,:) :: psb
      real(8) , allocatable, dimension(:,:,:) :: qca , qcb , qva , qvb ,    &
           & ta , tb , ua , ub , va , vb
      real(8) , allocatable, dimension(:,:) :: satbrt , tga
      real(8) , allocatable, dimension(:,:,:) :: snowc
      real(8) , allocatable, dimension(:,:,:) :: so4 , tbase
      real(8) , allocatable, dimension(:,:) :: uvdrag
#else
      real(8) , dimension(iy,jx) :: cldefi , f , hfx , hgfact , ht ,    &
                                  & htsd , msfd , msfx , pdotb , psa ,  &
                                  & psb , qfx , rainc , rainnc ,        &
                                  & satbrt , tga , tgb , tgbb , xlat ,  &
                                  & xlong , zpbl
      real(8) , dimension(iy,kz,jx) :: qca , qcb , qva , qvb , so4 ,    &
                                     & ta , tb , tbase , ua , ub , va , &
                                     & vb
      real(8) , dimension(nnsg,iy,jx) :: snowc
      real(8) , dimension(iy,jx) :: uvdrag
#endif

contains 
	subroutine allocate_mod_main

	allocate(cldefi(iy,jxp) )
	allocate(f(iy,jxp) )
	allocate(hfx(iy,jxp) )
	allocate(hgfact(iy,jxp) )
	allocate(htsd(iy,jxp) )
        allocate(qfx(iy,jxp) )
	allocate(rainc(iy,jxp) )
	allocate(rainnc(iy,jxp) )
	allocate(tgb(iy,jxp) )
	allocate(tgbb(iy,jxp) )
	allocate(xlat(iy,jxp) )
	allocate(xlong(iy,jxp) )
	allocate(zpbl(iy,jxp) )
        
        allocate(ht(iy,0:jxp+1) )

	allocate(msfd(iy,-1:jxp+2))
	allocate(msfx(iy,-1:jxp+2))
	allocate(pdotb(iy,-1:jxp+2))
	allocate(psa(iy,-1:jxp+2))
	allocate(rainp(iy,-1:jxp+2)) 

	allocate(psb(iy,0:jxp+2))

	allocate(qca(iy,kz,-1:jxp+2) )
	allocate(qcb(iy,kz,-1:jxp+2) )
	allocate(qva(iy,kz,-1:jxp+2) )
	allocate(qvb(iy,kz,-1:jxp+2) )
        allocate(ta(iy,kz,-1:jxp+2) )
	allocate(tb(iy,kz,-1:jxp+2) )
	allocate(ua(iy,kz,-1:jxp+2) )
	allocate(ub(iy,kz,-1:jxp+2) )
	allocate(va(iy,kz,-1:jxp+2) )
	allocate(vb(iy,kz,-1:jxp+2) )


        allocate(satbrt(iy,jxp+1)) 
	allocate(tga(iy,jxp+1)) 
        allocate(snowc(nnsg,iy,jxp)) 
        allocate(so4(iy,kz,jxp))
	allocate(tbase(iy,kz,jxp)) 
       
        allocate(uvdrag(iy,0:jxp))
	end subroutine allocate_mod_main 

 end module mod_main
