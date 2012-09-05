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
  real(rk8) :: clfrcv ! Cloud fractional cover for convective precip
  real(rk8) :: cllwcv ! Cloud liquid water content for convective precip.
  real(rk8) :: cevapu ! Raindrop evap rate coef [[(kg m-2 s-1)-1/2]/s]

  integer(ik4) , pointer , dimension(:,:) :: cucontrol ! which scheme to use

  ! Grell shared parameters for tracers removal
  integer(ik4) , pointer , dimension(:,:) :: icumtop , icumbot , icumdwd

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

  real(rk8) , pointer , dimension(:) :: flev , hlev , dflev , wlev
                                    ! sigma, a,     dsigma, qcon
  real(rk8) :: dtmdl
  real(rk8) :: dtcum , aprdiv ! dtsec , d_one/dble(ntsrf)

  logical :: lchem
  integer(ik4) :: total_precip_points

  data lchem /.false./

  contains
!
  subroutine allocate_mod_cu_common(ichem)
    implicit none
    integer(ik4) , intent(in) :: ichem
    if ( icup == 99 .or. icup == 98) then
      call getmem2d(cucontrol,jci1,jci2,ici1,ici2,'mod_cu_common:cucontrol')
    end if
    call getmem2d(icumbot,jci1,jci2,ici1,ici2,'mod_cu_common:icumbot')
    call getmem2d(icumtop,jci1,jci2,ici1,ici2,'mod_cu_common:icumtop')
    call getmem2d(icumdwd,jci1,jci2,ici1,ici2,'mod_cu_common:icumtop')
    if ( ichem == 1 ) then
      call getmem3d(convpr,jci1,jci2,ici1,ici2,1,kz,'mod_cu_common:convpr')
    end if
  end subroutine allocate_mod_cu_common
!
end module mod_cu_common
