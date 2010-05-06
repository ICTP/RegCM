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

      use mod_regcm_param

      implicit none
!
      real(8) :: caccr , cevap , clfrcv , clfrcvmax , cllwcv , conf ,   &
               & dtauc , edtmax , edtmaxo , edtmaxx , edtmin , edtmino ,&
               & edtminx , fcmax , gulland , guloce , htmax , htmin ,   &
               & mincld , pbcmax , qck10 , qck1land , qck1oce , qcth ,  &
               & qdcrit , rh0land , rh0oce , rhmax , shrmax , tc0 ,     &
               & shrmin , skbmax
#ifdef MPP1
      real(8) ,allocatable, dimension(:,:) :: cbmf2d , cgul , dtauc2d ,          &
                                   & edtmax2d , edtmaxo2d , edtmaxx2d , &
                                   & edtmin2d , edtmino2d , edtminx2d , &
                                   & htmax2d , htmin2d , mincld2d ,     &
                                   & pbcmax2d , qck1 , rh0 , shrmax2d , &
                                   & shrmin2d
      real(8) ,allocatable, dimension(:,:,:) :: fcc , rsheat , rswat
#else
      real(8) , dimension(iy,jx) :: cbmf2d , cgul , dtauc2d , edtmax2d ,&
                                  & edtmaxo2d , edtmaxx2d , edtmin2d ,  &
                                  & edtmino2d , edtminx2d , htmax2d ,   &
                                  & htmin2d , mincld2d , pbcmax2d ,     &
                                  & qck1 , rh0 , shrmax2d , shrmin2d
      real(8) , dimension(iy,kz,jx) :: fcc , rsheat , rswat
#endif
      real(8) , dimension(kz) :: qwght
      real(8) , dimension(kz,5:kz,1:kz-3) :: twght , vqflx
!
#ifdef MPP1
      integer ,allocatable, dimension(:) :: icon
      integer ,allocatable, dimension(:,:) :: kbmax2d
#else
      integer , dimension(jx) :: icon
      integer , dimension(iy,jx) :: kbmax2d
#endif
      integer :: kbmax

contains

     subroutine allocate_mod_pmoist

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
#endif

     end subroutine allocate_mod_pmoist
      end module mod_pmoist
