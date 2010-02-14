!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_pmoist

      use mod_regcm_param

      implicit none
!
! COMMON /PMOIST/
!
      real(8) :: caccr , cevap , clfrcv , clfrcvmax , cllwcv , conf ,   &
               & dtauc , edtmax , edtmaxo , edtmaxx , edtmin , edtmino ,&
               & edtminx , fcmax , gulland , guloce , htmax , htmin ,   &
               & mincld , pbcmax , qck10 , qck1land , qck1oce , qcth ,  &
               & qdcrit , rh0land , rh0oce , rhmax , rv , shrmax ,      &
               & shrmin , skbmax , tc0 
#ifdef MPP1
      real(8) , dimension(ix,jxp) :: cbmf2d , cgul , dtauc2d ,          &
                                   & edtmax2d , edtmaxo2d , edtmaxx2d , &
                                   & edtmin2d , edtmino2d , edtminx2d , &
                                   & htmax2d , htmin2d , mincld2d ,     &
                                   & pbcmax2d , qck1 , rh0 , shrmax2d , &
                                   & shrmin2d
      real(8) , dimension(ix,kx,jxp) :: fcc , rsheat , rswat
#else
      real(8) , dimension(ix,jx) :: cbmf2d , cgul , dtauc2d , edtmax2d ,&
                                  & edtmaxo2d , edtmaxx2d , edtmin2d ,  &
                                  & edtmino2d , edtminx2d , htmax2d ,   &
                                  & htmin2d , mincld2d , pbcmax2d ,     &
                                  & qck1 , rh0 , shrmax2d , shrmin2d
      real(8) , dimension(ix,kx,jx) :: fcc , rsheat , rswat
#endif
      real(8) , dimension(kx) :: qwght
      real(8) , dimension(kx,5:kx,1:kx-3) :: twght , vqflx
!
! COMMON /PMOISTINT/
!
#ifdef MPP1
      integer , dimension(jxp) :: icon
      integer , dimension(ix,jxp) :: kbmax2d
#else
      integer , dimension(jx) :: icon
      integer , dimension(ix,jx) :: kbmax2d
#endif
      integer :: kbmax

      end module mod_pmoist
