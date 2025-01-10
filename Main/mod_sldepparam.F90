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
  integer(ik4), pointer , dimension(:,:,:) :: xndp_d => null( )
  integer(ik4), pointer , dimension(:,:,:) :: xnnm1dp_d => null( )
  integer(ik4), pointer , dimension(:,:,:) :: xnnm2dp_d => null( )
  integer(ik4), pointer , dimension(:,:,:) :: xnnp1dp_d => null( )
  integer(ik4), pointer , dimension(:,:,:) :: yndp_d => null( )
  integer(ik4), pointer , dimension(:,:,:) :: ynnm1dp_d => null( )
  integer(ik4), pointer , dimension(:,:,:) :: ynnm2dp_d => null( )
  integer(ik4), pointer , dimension(:,:,:) :: ynnp1dp_d => null( )

  ! GTD the departe point location (for the interpolation) for cross points
  integer(ik4) , pointer , dimension(:,:,:) :: xndp_x => null( )
  integer(ik4) , pointer , dimension(:,:,:) :: xnnm1dp_x => null( )
  integer(ik4) , pointer , dimension(:,:,:) :: xnnm2dp_x => null( )
  integer(ik4) , pointer , dimension(:,:,:) :: xnnp1dp_x => null( )
  integer(ik4) , pointer , dimension(:,:,:) :: yndp_x => null( )
  integer(ik4) , pointer , dimension(:,:,:) :: ynnm1dp_x => null( )
  integer(ik4) , pointer , dimension(:,:,:) :: ynnm2dp_x => null( )
  integer(ik4) , pointer , dimension(:,:,:) :: ynnp1dp_x => null( )

  ! GTD the weighting coefficient of interpolation) for dot points
  real(rkx) , pointer , dimension(:,:,:) :: alfdp_d , alfp1dp_d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: alfm1dp_d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: alfm2dp_d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: betdp_d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: betp1dp_d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: betm1dp_d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: betm2dp_d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: alffbl_d => null( )

  ! GTD the weighting coefficient of interpolation) for dot points
  real(rkx) , pointer , dimension(:,:,:) :: alfdp_x => null( )
  real(rkx) , pointer , dimension(:,:,:) :: alfp1dp_x => null( )
  real(rkx) , pointer , dimension(:,:,:) :: alfm1dp_x => null( )
  real(rkx) , pointer , dimension(:,:,:) :: alfm2dp_x => null( )
  real(rkx) , pointer , dimension(:,:,:) :: betdp_x => null( )
  real(rkx) , pointer , dimension(:,:,:) :: betp1dp_x => null( )
  real(rkx) , pointer , dimension(:,:,:) :: betm1dp_x => null( )
  real(rkx) , pointer , dimension(:,:,:) :: betm2dp_x => null( )
  real(rkx) , pointer , dimension(:,:,:) :: alffbl_x => null( )

  ! GTD the advective velocity for the dot and cross points
  real(rkx) , pointer , dimension(:,:,:) :: vadvy_d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: uadvx_d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: vadvy_x => null( )
  real(rkx) , pointer , dimension(:,:,:) :: uadvx_x => null( )

  ! GTD the advective velocity near the arrival point for mcgregor calculation
  real(rkx) , pointer , dimension(:,:,:) :: uadxp1_d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: uadxm1_d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: vadyp1_d => null( )
  real(rkx) , pointer , dimension(:,:,:) :: vadym1_d => null( )

  ! GTD the advective velocity near the arrival point for mcgregor calculation
  real(rkx) , pointer , dimension(:,:,:) :: uadxp1_x => null( )
  real(rkx) , pointer , dimension(:,:,:) :: uadxm1_x => null( )
  real(rkx) , pointer , dimension(:,:,:) :: vadyp1_x => null( )
  real(rkx) , pointer , dimension(:,:,:) :: vadym1_x => null( )

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
