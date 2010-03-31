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
      real(8) , dimension(iy,jxp) :: cldefi , f , hfx , hgfact , htsd , &
                                   & qfx , rainc , rainnc , tgb , tgbb ,&
                                   & xlat , xlong , zpbl
      real(8) , dimension(iy,0:jxp+1) :: ht
      real(8) , dimension(iy,-1:jxp+2) :: msfd , msfx , pdotb , psa ,   &
           & rainp
      real(8) , dimension(iy,0:jxp+2) :: psb
      real(8) , dimension(iy,kz,-1:jxp+2) :: qca , qcb , qva , qvb ,    &
           & ta , tb , ua , ub , va , vb
      real(8) , dimension(iy,jxp+1) :: satbrt , tga
      real(8) , dimension(nnsg,iy,jxp) :: snowc
      real(8) , dimension(iy,kz,jxp) :: so4 , tbase
      real(8) , dimension(iy,0:jxp) :: uvdrag
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

      end module mod_main
