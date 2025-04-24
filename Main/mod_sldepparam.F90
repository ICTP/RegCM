!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_sldepparam
  !
  ! for the departure point calculation
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil

  implicit none

  public

  ! GTD the departe point location (for the interpolation) for dot points
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: xndp_d => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: xnnm1dp_d => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: xnnm2dp_d => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: xnnp1dp_d => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: yndp_d => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: ynnm1dp_d => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: ynnm2dp_d => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: ynnp1dp_d => null( )

  ! GTD the departe point location (for the interpolation) for cross points
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: xndp_x => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: xnnm1dp_x => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: xnnm2dp_x => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: xnnp1dp_x => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: yndp_x => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: ynnm1dp_x => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: ynnm2dp_x => null( )
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: ynnp1dp_x => null( )

  ! GTD the weighting coefficient of interpolation) for dot points
  real(rkx), pointer, contiguous, dimension(:,:,:) :: alfdp_d, alfp1dp_d => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: alfm1dp_d => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: alfm2dp_d => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: betdp_d => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: betp1dp_d => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: betm1dp_d => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: betm2dp_d => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: alffbl_d => null( )

  ! GTD the weighting coefficient of interpolation) for dot points
  real(rkx), pointer, contiguous, dimension(:,:,:) :: alfdp_x => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: alfp1dp_x => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: alfm1dp_x => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: alfm2dp_x => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: betdp_x => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: betp1dp_x => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: betm1dp_x => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: betm2dp_x => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: alffbl_x => null( )

  ! GTD the advective velocity for the dot and cross points
  real(rkx), pointer, contiguous, dimension(:,:,:) :: vadvy_d => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uadvx_d => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: vadvy_x => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uadvx_x => null( )

  ! GTD the advective velocity near the arrival point for mcgregor calculation
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uadxp1_d => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uadxm1_d => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: vadyp1_d => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: vadym1_d => null( )

  ! GTD the advective velocity near the arrival point for mcgregor calculation
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uadxp1_x => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uadxm1_x => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: vadyp1_x => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: vadym1_x => null( )

  contains

  subroutine allocate_mod_sldepparam
    implicit none
    call getmem3d(xndp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:xndp_x')
    call getmem3d(xnnm1dp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:xnnm1dp_x')
    call getmem3d(xnnm2dp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:xnnm2dp_x')
    call getmem3d(xnnp1dp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:xnnp1dp_x')
    call getmem3d(yndp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:yndp_x')
    call getmem3d(ynnm1dp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:ynnm1dp_x')
    call getmem3d(ynnm2dp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:ynnm2dp_x')
    call getmem3d(ynnp1dp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:ynnp1dp_x')
    call getmem3d(alfdp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:alfdp_x')
    call getmem3d(alfm1dp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:alfm1dp_x')
    call getmem3d(alfm2dp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:alfm2dp_x')
    call getmem3d(alfp1dp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:alfp1dp_x')
    call getmem3d(betdp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:betdp_x')
    call getmem3d(betm1dp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:betm1dp_x')
    call getmem3d(betm2dp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:betm2dp_x')
    call getmem3d(betp1dp_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:betp1dp_x')
    call getmem3d(alffbl_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:alffbl_x')
    call getmem3d(uadvx_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:uadvx_x')
    call getmem3d(vadvy_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:vadvy_x')
    call getmem3d(uadxp1_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:uadxp1_x')
    call getmem3d(uadxm1_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:uadxm1_x')
    call getmem3d(vadyp1_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:vadyp1_x')
    call getmem3d(vadym1_x,jci1,jci2,ici1,ici2,1,kz,'sldepparam:vadym1_x')
    call getmem3d(xndp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:xndp_d')
    call getmem3d(xnnm1dp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:xnnm1dp_d')
    call getmem3d(xnnm2dp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:xnnm2dp_d')
    call getmem3d(xnnp1dp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:xnnp1dp_d')
    call getmem3d(yndp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:yndp_d')
    call getmem3d(ynnm1dp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:ynnm1dp_d')
    call getmem3d(ynnm2dp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:ynnm2dp_d')
    call getmem3d(ynnp1dp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:ynnp1dp_d')
    call getmem3d(alfdp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:alfdp_d')
    call getmem3d(alfm1dp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:alfm1dp_d')
    call getmem3d(alfm2dp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:alfm2dp_d')
    call getmem3d(alfp1dp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:alfp1dp_d')
    call getmem3d(betdp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:betdp_d')
    call getmem3d(betm1dp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:betm1dp_d')
    call getmem3d(betm2dp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:betm2dp_d')
    call getmem3d(betp1dp_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:betp1dp_d')
    call getmem3d(alffbl_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:alffbl_d')
    call getmem3d(uadvx_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:uadvx_d')
    call getmem3d(vadvy_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:vadvy_d')
    call getmem3d(uadxp1_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:uadxp1_d')
    call getmem3d(uadxm1_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:uadxm1_d')
    call getmem3d(vadyp1_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:vadyp1_d')
    call getmem3d(vadym1_d,jdi1,jdi2,idi1,idi2,1,kz,'sldepparam:vadym1_d')
  end subroutine allocate_mod_sldepparam

end module mod_sldepparam
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
