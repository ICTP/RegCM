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

      module mod_pmoist

      use mod_dynparam

      implicit none
!
      real(8) :: caccr , cevap , clfrcv , clfrcvmax , cllwcv , conf ,   &
               & dtauc , edtmax , edtmaxo , edtmaxx , edtmin , edtmino ,&
               & edtminx , fcmax , gulland , guloce , htmax , htmin ,   &
               & mincld , pbcmax , qck10 , qck1land , qck1oce , qcth ,  &
               & qdcrit , rh0land , rh0oce , rhmax , shrmax , tc0 ,     &
               & shrmin , skbmax
      real(8) , allocatable , dimension(:,:) :: cbmf2d , cgul ,         &
               & dtauc2d , edtmax2d , edtmaxo2d , edtmaxx2d ,           &
                                   & edtmin2d , edtmino2d , edtminx2d , &
                                   & htmax2d , htmin2d , mincld2d ,     &
                                   & pbcmax2d , qck1 , rh0 , shrmax2d , &
                                   & shrmin2d
      real(8) , allocatable , dimension(:,:,:) :: fcc , rsheat , rswat
      real(8) , allocatable , dimension(:) :: qwght
      real(8) , allocatable , dimension(:,:,:) :: twght , vqflx
      integer ,allocatable, dimension(:) :: icon
      integer ,allocatable, dimension(:,:) :: kbmax2d
      integer :: kbmax

      contains

      subroutine allocate_mod_pmoist
      implicit none
#ifdef MPP1
      allocate(cbmf2d(iy,jxp))
      allocate(cgul(iy,jxp))
      allocate(dtauc2d(iy,jxp))
      allocate(edtmax2d(iy,jxp))
      allocate(edtmaxo2d(iy,jxp))
      allocate(edtmaxx2d(iy,jxp))
      allocate(edtmin2d(iy,jxp))
      allocate(edtmino2d(iy,jxp))
      allocate(edtminx2d(iy,jxp))
      allocate(htmax2d(iy,jxp))
      allocate(htmin2d(iy,jxp))
      allocate(mincld2d(iy,jxp))
      allocate(pbcmax2d(iy,jxp))
      allocate(qck1(iy,jxp))
      allocate(rh0(iy,jxp))
      allocate(shrmax2d(iy,jxp))
      allocate(shrmin2d(iy,jxp))
      allocate(fcc(iy,kz,jxp))
      allocate(rsheat(iy,kz,jxp))
      allocate(rswat(iy,kz,jxp))
      allocate(icon(jxp))
      allocate(kbmax2d(iy,jxp))
#else
      allocate(cbmf2d(iy,jx))
      allocate(cgul(iy,jx))
      allocate(dtauc2d(iy,jx))
      allocate(edtmax2d(iy,jx))
      allocate(edtmaxo2d(iy,jx))
      allocate(edtmaxx2d(iy,jx))
      allocate(edtmin2d(iy,jx))
      allocate(edtmino2d(iy,jx))
      allocate(edtminx2d(iy,jx))
      allocate(htmax2d(iy,jx))
      allocate(htmin2d(iy,jx))
      allocate(mincld2d(iy,jx))
      allocate(pbcmax2d(iy,jx))
      allocate(qck1(iy,jx))
      allocate(rh0(iy,jx))
      allocate(shrmax2d(iy,jx))
      allocate(shrmin2d(iy,jx))
      allocate(fcc(iy,kz,jx))
      allocate(rsheat(iy,kz,jx))
      allocate(rswat(iy,kz,jx))
      allocate(icon(jx))
      allocate(kbmax2d(iy,jx))
#endif
      allocate(qwght(kz))
      allocate(twght(kz,5:kz,1:kz-3))
      allocate(vqflx(kz,5:kz,1:kz-3))
      end subroutine allocate_mod_pmoist

      end module mod_pmoist
