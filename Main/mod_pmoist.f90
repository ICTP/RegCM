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
  use mod_constants
  use mod_dynparam
  use mod_memutil

  implicit none
!
  real(8) :: clfrcv , clfrcvmax , cllwcv , conf ,   &
             dtauc , fcmax , gulland , guloce , htmax , htmin ,     &
             mincld , pbcmax , qck10 , qck1land , qck1oce ,         &
             rh0land , rh0oce , tc0 , skbmax

  real(8) , pointer , dimension(:,:) :: cbmf2d , dtauc2d ,        &
        edtmax2d , edtmaxo2d , edtmaxx2d , edtmin2d , edtmino2d , &
        edtminx2d , htmax2d , htmin2d , mincld2d , pbcmax2d ,     &
        rh0 , shrmax2d , shrmin2d
  real(8) , pointer , dimension(:,:,:) :: rsheat , rswat
  real(8) , pointer , dimension(:) :: qwght
  real(8) , pointer , dimension(:,:,:) :: twght , vqflx
  integer , pointer , dimension(:) :: icon
  integer , pointer , dimension(:,:) :: kbmax2d
  integer :: kbmax
!
  contains
!
  subroutine allocate_mod_pmoist
  implicit none

  call getmem2d(cbmf2d,1,iy,1,jxp,'mod_pmoist:cbmf2d')
  call getmem2d(dtauc2d,1,iy,1,jxp,'mod_pmoist:dtauc2d')
  call getmem2d(edtmax2d,1,iy,1,jxp,'mod_pmoist:edtmax2d')
  call getmem2d(edtmaxo2d,1,iy,1,jxp,'mod_pmoist:edtmaxo2d')
  call getmem2d(edtmaxx2d,1,iy,1,jxp,'mod_pmoist:edtmaxx2d')
  call getmem2d(edtmin2d,1,iy,1,jxp,'mod_pmoist:edtmin2d')
  call getmem2d(edtmino2d,1,iy,1,jxp,'mod_pmoist:edtmino2d')
  call getmem2d(edtminx2d,1,iy,1,jxp,'mod_pmoist:edtminx2d')
  call getmem2d(htmax2d,1,iy,1,jxp,'mod_pmoist:htmax2d')
  call getmem2d(htmin2d,1,iy,1,jxp,'mod_pmoist:htmin2d')
  call getmem2d(mincld2d,1,iy,1,jxp,'mod_pmoist:mincld2d')
  call getmem2d(pbcmax2d,1,iy,1,jxp,'mod_pmoist:pbcmax2d')
  call getmem2d(rh0,1,iy,1,jxp,'mod_pmoist:rh0')
  call getmem2d(shrmax2d,1,iy,1,jxp,'mod_pmoist:shrmax2d')
  call getmem2d(shrmin2d,1,iy,1,jxp,'mod_pmoist:shrmin2d')
  call getmem3d(rsheat,1,iy,1,kz,1,jxp,'mod_pmoist:rsheat')
  call getmem3d(rswat,1,iy,1,kz,1,jxp,'mod_pmoist:rswat')
  call getmem1d(icon,1,jxp,'mod_pmoist:icon')
  call getmem2d(kbmax2d,1,iy,1,jxp,'mod_pmoist:kbmax2d')
  call getmem1d(qwght,1,kz,'mod_pmoist:qwght')
  call getmem3d(twght,1,kz,5,kz,1,kz-3,'mod_pmoist:twght')
  call getmem3d(vqflx,1,kz,5,kz,1,kz-3,'mod_pmoist:vqflx')
!
  end subroutine allocate_mod_pmoist
!
end module mod_pmoist
