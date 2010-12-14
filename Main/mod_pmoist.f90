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
!
! Storage, parameters and constants related to
!     moisture calculations.
!
      use mod_dynparam

      implicit none
!
      real(8) :: caccr , cevap , clfrcv , clfrcvmax , cllwcv , conf ,   &
               & dtauc , edtmax , edtmaxo , edtmaxx , edtmin , edtmino ,&
               & edtminx , fcmax , gulland , guloce , htmax , htmin ,   &
               & mincld , pbcmax , qck10 , qck1land , qck1oce , qcth ,  &
               & rh0land , rh0oce , rhmax , shrmax , tc0 , shrmin ,     &
               & skbmax, shrmax_ocn, shrmin_ocn, edtmax_ocn, edtmin_ocn,&
               & edtmaxo_ocn, edtmino_ocn, edtmaxx_ocn, edtminx_ocn

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
!
      contains
!
      subroutine allocate_mod_pmoist(lmpi)
      implicit none
      logical , intent(in) :: lmpi
      integer :: nj

      if (lmpi) then
        nj = jxp
      else
        nj = jx
      end if

      allocate(cbmf2d(iy,nj))
      allocate(cgul(iy,nj))
      allocate(dtauc2d(iy,nj))
      allocate(edtmax2d(iy,nj))
      allocate(edtmaxo2d(iy,nj))
      allocate(edtmaxx2d(iy,nj))
      allocate(edtmin2d(iy,nj))
      allocate(edtmino2d(iy,nj))
      allocate(edtminx2d(iy,nj))
      allocate(htmax2d(iy,nj))
      allocate(htmin2d(iy,nj))
      allocate(mincld2d(iy,nj))
      allocate(pbcmax2d(iy,nj))
      allocate(qck1(iy,nj))
      allocate(rh0(iy,nj))
      allocate(shrmax2d(iy,nj))
      allocate(shrmin2d(iy,nj))
      allocate(fcc(iy,kz,nj))
      allocate(rsheat(iy,kz,nj))
      allocate(rswat(iy,kz,nj))
      cbmf2d = 0.0D0
      cgul = 0.0D0
      dtauc2d = 0.0D0
      edtmax2d = 0.0D0
      edtmaxo2d = 0.0D0
      edtmaxx2d = 0.0D0
      edtmin2d = 0.0D0
      edtmino2d = 0.0D0
      edtminx2d = 0.0D0
      htmax2d = 0.0D0
      htmin2d = 0.0D0
      mincld2d = 0.0D0
      pbcmax2d = 0.0D0
      qck1 = 0.0D0
      rh0 = 0.0D0
      shrmax2d = 0.0D0
      shrmin2d = 0.0D0
      fcc = 0.0D0
      rsheat = 0.0D0
      rswat = 0.0D0
      allocate(icon(nj))
      allocate(kbmax2d(iy,nj))
      icon = 0
      kbmax2d = 0
      allocate(qwght(kz))
      allocate(twght(kz,5:kz,1:kz-3))
      allocate(vqflx(kz,5:kz,1:kz-3))
      qwght = 0.0D0
      twght = 0.0D0
      vqflx = 0.0D0
!
      end subroutine allocate_mod_pmoist
!
      end module mod_pmoist
