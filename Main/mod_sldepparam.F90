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
  integer(ik4), pointer , dimension(:,:,:) :: xndp_d , xnnm1dp_d , &
       xnnm2dp_d , xnnp1dp_d , yndp_d , ynnm1dp_d , ynnm2dp_d , ynnp1dp_d

  ! GTD the departe point location (for the interpolation) for cross points
  integer(ik4) , pointer , dimension(:,:,:) :: xndp_x , xnnm1dp_x , &
       xnnm2dp_x , xnnp1dp_x , yndp_x , ynnm1dp_x , ynnm2dp_x , ynnp1dp_x

  ! GTD the weighting coefficient of interpolation) for dot points
  real(rk8) , pointer , dimension(:,:,:) :: alfdp_d , alfp1dp_d , &
       alfm1dp_d , alfm2dp_d , betdp_d , betp1dp_d , betm1dp_d ,  &
       betm2dp_d , alffbl_d

  ! GTD the weighting coefficient of interpolation) for dot points
  real(rk8) , pointer , dimension(:,:,:) :: alfdp_x , alfp1dp_x , &
       alfm1dp_x , alfm2dp_x , betdp_x , betp1dp_x , betm1dp_x ,  &
       betm2dp_x , alffbl_x

  ! GTD the advective velocity for the dot and cross points
  real(rk8) , pointer , dimension(:,:,:) :: vadvy_d , uadvx_d
  real(rk8) , pointer , dimension(:,:,:) :: vadvy_x , uadvx_x

  ! GTD the advective velocity near the arrival point for mcgregor calculation
  real(rk8) , pointer , dimension(:,:,:) :: uadxp1_d , uadxm1_d
  real(rk8) , pointer , dimension(:,:,:) :: vadyp1_d , vadym1_d

  ! GTD the advective velocity near the arrival point for mcgregor calculation
  real(rk8) , pointer , dimension(:,:,:) :: uadxp1_x , uadxm1_x
  real(rk8) , pointer , dimension(:,:,:) :: vadyp1_x , vadym1_x

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
