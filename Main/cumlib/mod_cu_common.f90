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
  use mod_dynparam
  use mod_memutil
  use mod_realkinds

  public
!

  real(dp) :: clfrcv ! Cloud fractional cover for convective precip
  real(dp) :: cllwcv ! Cloud liquid water content for convective precip.
  real(dp) :: cevapu ! Raindrop evap rate coef [[(kg m-2 s-1)-1/2]/s]

  integer , pointer , dimension(:,:) :: cucontrol ! which scheme to use

  ! Grell shared parameters for tracers removal
  integer , pointer , dimension(:,:) :: icumtop , icumbot , icumdwd
!
  real(dp) , pointer , dimension(:,:) :: sfhgt    ! mddom%ht
  real(dp) , pointer , dimension(:,:,:) :: hgt    ! za
  real(dp) , pointer , dimension(:,:,:) :: ptatm  ! atm1%t
  real(dp) , pointer , dimension(:,:,:) :: puatm  ! atm1%u
  real(dp) , pointer , dimension(:,:,:) :: pvatm  ! atm1%v
  real(dp) , pointer , dimension(:,:,:) :: pvqvtm ! atm1%qv
  real(dp) , pointer , dimension(:,:,:) :: tas    ! atms%tb3d
  real(dp) , pointer , dimension(:,:,:) :: uas    ! atms%ubx3d
  real(dp) , pointer , dimension(:,:,:) :: vas    ! atms%vbx3d
  real(dp) , pointer , dimension(:,:,:) :: pas    ! atms%pb3d
  real(dp) , pointer , dimension(:,:,:) :: qsas   ! atms%qsb3d
  real(dp) , pointer , dimension(:,:,:) :: qcas   ! atms%qcb3d
  real(dp) , pointer , dimension(:,:,:) :: qvas   ! atms%qvb3d
  real(dp) , pointer , dimension(:,:,:,:) :: chias   ! atms%chib3d
  real(dp) , pointer , dimension(:,:,:) :: tten   ! aten%t
  real(dp) , pointer , dimension(:,:,:) :: uten   ! aten%u
  real(dp) , pointer , dimension(:,:,:) :: vten   ! aten%v
  real(dp) , pointer , dimension(:,:,:) :: qvten  ! aten%qv
  real(dp) , pointer , dimension(:,:,:) :: qcten  ! aten%qc
  real(dp) , pointer , dimension(:,:,:,:) :: tchiten  !chiten 
  real(dp) , pointer , dimension(:,:) :: psfcps   ! sfs%psa
  real(dp) , pointer , dimension(:,:) :: sfcps    ! sfs%psb
  real(dp) , pointer , dimension(:,:) :: rainc    ! sfs%rainc
  real(dp) , pointer , dimension(:,:) :: qfx      ! sfs%qfx
  real(dp) , pointer , dimension(:,:,:) :: svv    ! qdot
  real(dp) , pointer , dimension(:,:) :: lmpcpc   ! pptc
  integer , pointer , dimension(:,:) :: lmask    ! ldmsk
  real(dp) , pointer , dimension(:,:,:) :: rcldlwc  ! rcldlwc 
  real(dp) , pointer , dimension(:,:,:) :: rcldfra  ! rcldfra

  real(dp) , pointer , dimension(:) :: flev , hlev , dflev , wlev
                                    ! sigma, a,     dsigma, qcon
  real(dp) :: dtmdl
  real(dp) :: dtcum , aprdiv ! dtsec , d_one/dble(ntsrf)

  logical :: lchem
  integer :: icup
  integer :: igcc
  integer :: total_precip_points

  data lchem /.false./

  contains
!
  subroutine allocate_mod_cu_common
    implicit none
    if ( icup == 99 .or. icup == 98) then
      call getmem2d(cucontrol,jci1,jci2,ici1,ici2,'mod_cu_common:cucontrol')
    end if
    call getmem2d(icumbot,jci1,jci2,ici1,ici2,'mod_cu_common:icumbot')
    call getmem2d(icumtop,jci1,jci2,ici1,ici2,'mod_cu_common:icumtop')
    call getmem2d(icumdwd,jci1,jci2,ici1,ici2,'mod_cu_common:icumtop')
  end subroutine allocate_mod_cu_common
!
end module mod_cu_common
