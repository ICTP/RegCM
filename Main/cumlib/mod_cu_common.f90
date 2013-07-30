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

module mod_cu_common
!
! Storage and constants related to cumulus convection schemes.
!
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_runparams

  public
!
  real(rk8) :: cevapu ! Raindrop evap rate coef [[(kg m-2 s-1)-1/2]/s]

  integer(ik4) , pointer , dimension(:,:) :: cucontrol ! which scheme to use

  ! Grell shared parameters for tracers removal
  integer(ik4) , pointer , dimension(:,:) :: icumtop , icumbot
!
  real(rk8) , pointer , dimension(:,:) :: sfhgt    ! mddom%ht
  real(rk8) , pointer , dimension(:,:,:) :: hgt    ! za
  real(rk8) , pointer , dimension(:,:,:) :: ptatm  ! atm1%t
  real(rk8) , pointer , dimension(:,:,:) :: puatm  ! atm1%u
  real(rk8) , pointer , dimension(:,:,:) :: pvatm  ! atm1%v
  real(rk8) , pointer , dimension(:,:,:,:) :: pvqxtm ! atm1%qx
  real(rk8) , pointer , dimension(:,:,:) :: tas    ! atms%tb3d
  real(rk8) , pointer , dimension(:,:,:) :: uas    ! atms%ubx3d
  real(rk8) , pointer , dimension(:,:,:) :: vas    ! atms%vbx3d
  real(rk8) , pointer , dimension(:,:,:) :: pas    ! atms%pb3d
  real(rk8) , pointer , dimension(:,:,:) :: qsas   ! atms%qsb3d
  real(rk8) , pointer , dimension(:,:,:,:) :: qxas ! atms%qxb3d
  real(rk8) , pointer , dimension(:,:,:,:) :: chias   ! atms%chib3d
  real(rk8) , pointer , dimension(:,:,:) :: tten   ! aten%t
  real(rk8) , pointer , dimension(:,:,:) :: uten   ! aten%u
  real(rk8) , pointer , dimension(:,:,:) :: vten   ! aten%v
  real(rk8) , pointer , dimension(:,:,:,:) :: qxten    ! aten%qx
  real(rk8) , pointer , dimension(:,:,:,:) :: tchiten  ! chiten 
  real(rk8) , pointer , dimension(:,:) :: psfcps   ! sfs%psa
  real(rk8) , pointer , dimension(:,:) :: sfcps    ! sfs%psb
  real(rk8) , pointer , dimension(:,:) :: rainc    ! sfs%rainc
  real(rk8) , pointer , dimension(:,:,:) :: convpr ! prec rate ( used in chem) 
  real(rk8) , pointer , dimension(:,:) :: qfx      ! sfs%qfx
  real(rk8) , pointer , dimension(:,:,:) :: svv    ! qdot
  real(rk8) , pointer , dimension(:,:) :: lmpcpc   ! pptc
  integer(ik4) , pointer , dimension(:,:) :: lmask    ! ldmsk
  real(rk8) , pointer , dimension(:,:,:) :: rcldlwc  ! rcldlwc 
  real(rk8) , pointer , dimension(:,:,:) :: rcldfra  ! rcldfra
  integer(ik4) , pointer , dimension(:,:) :: rktrop  ! ktrop

  integer(ik4) :: total_precip_points

  real(rk8) , dimension(10) :: cld_profile
  real(rk8) , parameter :: maxcloud_dp =  85.0D0 ! In cb

  contains
!
  subroutine allocate_mod_cu_common(ichem)
    implicit none
    integer(ik4) , intent(in) :: ichem
    integer :: k
    if ( icup == 99 .or. icup == 98) then
      call getmem2d(cucontrol,jci1,jci2,ici1,ici2,'mod_cu_common:cucontrol')
    end if
    call getmem2d(icumbot,jci1,jci2,ici1,ici2,'mod_cu_common:icumbot')
    call getmem2d(icumtop,jci1,jci2,ici1,ici2,'mod_cu_common:icumtop')
    if ( ichem == 1 ) then
      call getmem3d(convpr,jci1,jci2,ici1,ici2,1,kz,'mod_cu_common:convpr')
    end if
    cld_profile(1) = 3.0D0/dble(kz)
    cld_profile(2) = 2.0D0/dble(kz)
    do k = 3 , 9
      cld_profile(k) = 1.00D0/dble(kz)
    end do
    cld_profile(10) = 2.0D0/dble(kz)
  end subroutine allocate_mod_cu_common
!
  subroutine model_cumulus_cloud
    implicit none
    real(rk8) :: akclth , tcel , scalep , scalef
    integer(ik4):: i , j , k , ktop , kbot , kclth , ikh

    rcldfra(:,:,:) = d_zero
    rcldlwc(:,:,:) = d_zero

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          tcel = tas(j,i,k)-tzero
          ! Temperature dependance for convective cloud water content
          ! in g/m3 (Lemus et al., 1997)
          if ( tcel < -50D0 ) then
            rcldlwc(j,i,k) = 0.001D0
          else
            rcldlwc(j,i,k) = 0.127D+00 + 6.78D-03 * tcel + &
                             1.29D-04 * tcel**2 + 8.36D-07 * tcel**3
            if ( rcldlwc(j,i,k) > 0.300D0 ) rcldlwc(j,i,k) = 0.300D0
            if ( rcldlwc(j,i,k) < 0.001D0 ) rcldlwc(j,i,k) = 0.001D0
          end if
        end do
      end do
    end do

    scalef = (d_one-clfrcv)

    if ( icumcloud == 1 ) then
      iloop1: &
      do i = ici1 , ici2
        jloop1: &
        do j = jci1 , jci2
          ! The regcm model is top to bottom
          ktop = icumtop(j,i)
          kbot = icumbot(j,i)
          kclth = kbot - ktop + 1
          if ( kclth < 2 ) cycle jloop1
          akclth = d_one/dble(kclth)
          do k = ktop , kbot
            rcldfra(j,i,k) = d_one - scalef**akclth
          end do
        end do jloop1
      end do iloop1
    else if ( icumcloud == 2 ) then
      iloop3: &
      do i = ici1 , ici2
        jloop3: &
        do j = jci1 , jci2
          ktop = icumtop(j,i)
          kbot = icumbot(j,i)
          kclth = kbot - ktop + 1
          if ( kclth < 2 ) cycle jloop3
          scalep = min((pas(j,i,kbot)-pas(j,i,ktop))/maxcloud_dp,d_one)
          do k = ktop , kbot
            ikh = max(1,min(10,int((dble(k-ktop+1)/dble(kclth))*d_10)))
            rcldfra(j,i,k) = cld_profile(ikh)*clfrcv*scalep
          end do
        end do jloop3
      end do iloop3
    end if

    where ( rcldfra < 1.0D-4 )
      rcldlwc = 0.001D0
    end where

  end subroutine model_cumulus_cloud
!
end module mod_cu_common
