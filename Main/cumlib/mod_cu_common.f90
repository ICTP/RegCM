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

  public
!

  real(8) :: clfrcv ! Cloud fractional cover for convective precip
  real(8) :: cllwcv ! Cloud liquid water content for convective precip.
  real(8) :: cevapu ! Raindrop evap rate coef [[(kg m-2 s-1)-1/2]/s]

  integer , pointer , dimension(:) :: icon ! Precip points counter
  integer , pointer , dimension(:,:) :: cucontrol ! which scheme to use
  integer , pointer , dimension(:,:) :: icumtop , icumbot
!
  real(8) , pointer , dimension(:,:) :: sfhgt    ! mddom%ht
  real(8) , pointer , dimension(:,:) :: psfcps   ! sps1%ps
  real(8) , pointer , dimension(:,:) :: sfcps    ! sps2%ps
  real(8) , pointer , dimension(:,:,:) :: hgt    ! za
  real(8) , pointer , dimension(:,:,:) :: ptatm  ! atm1%t
  real(8) , pointer , dimension(:,:,:) :: puatm  ! atm1%u
  real(8) , pointer , dimension(:,:,:) :: pvatm  ! atm1%v
  real(8) , pointer , dimension(:,:,:) :: pvqvtm ! atm1%qv
  real(8) , pointer , dimension(:,:,:) :: tas    ! atms%tb3d
  real(8) , pointer , dimension(:,:,:) :: uas    ! atms%ubx3d
  real(8) , pointer , dimension(:,:,:) :: vas    ! atms%vbx3d
  real(8) , pointer , dimension(:,:,:) :: pas    ! atms%pb3d
  real(8) , pointer , dimension(:,:,:) :: qsas   ! atms%qsb3d
  real(8) , pointer , dimension(:,:,:) :: qcas   ! atms%qcb3d
  real(8) , pointer , dimension(:,:,:) :: qvas   ! atms%qvb3d
  real(8) , pointer , dimension(:,:,:) :: tten   ! aten%t
  real(8) , pointer , dimension(:,:,:) :: uten   ! aten%u
  real(8) , pointer , dimension(:,:,:) :: vten   ! aten%v
  real(8) , pointer , dimension(:,:,:) :: qvten  ! aten%qv
  real(8) , pointer , dimension(:,:,:) :: qcten  ! aten%qc
  real(8) , pointer , dimension(:,:) :: rainc    ! sfsta%rainc
  real(8) , pointer , dimension(:,:) :: qfx      ! sfsta%qfx
  real(8) , pointer , dimension(:,:,:) :: svv    ! qdot
  real(8) , pointer , dimension(:,:) :: lmpcpc   ! pptc
  integer , pointer , dimension(:,:) :: lmask    ! ldmsk
  real(8) , pointer , dimension(:,:) :: rcldlwc  ! rcldlwc 
  real(8) , pointer , dimension(:,:) :: rcldfra  ! rcldfra

  real(8) , pointer , dimension(:) :: flev , hlev , dflev , wlev
                                    ! sigma, a,     dsigma, qcon

  real(8) :: dtmdl
  real(8) :: dtcum , aprdiv ! dtsec , d_one/dble(ntsrf)
  real(8) :: dtcum2 ! 2*dtcum

  logical :: lchem
  integer :: icup
  integer :: igcc

  data lchem /.false./

  contains
!
  subroutine allocate_mod_cu_common
    implicit none
    call getmem1d(icon,1,jxp,'mod_cu_common:icon')
    if ( icup == 99 .or. icup == 98) then
      call getmem2d(cucontrol,1,iy,1,jxp,'mod_cu_common:cucontrol')
    end if
    call getmem2d(icumbot,1,iy,1,jxp,'mod_cu_common:icumbot')
    call getmem2d(icumtop,1,iy,1,jxp,'mod_cu_common:icumtop')
  end subroutine allocate_mod_cu_common
!
end module mod_cu_common
