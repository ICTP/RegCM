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

module mod_ocn_albedo

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_runparams
  use mod_service
  use mod_ocn_internal

  implicit none

  private
  !
  ! Solar flux partitioned at wavelength of 0.7micr
  !
  real(rkx) , parameter :: fsol1 = 0.5_rkx
  real(rkx) , parameter :: fsol2 = 0.5_rkx
  !
  ! Short and long wave albedo for new snow
  !
  real(rkx) , parameter :: snal0 = 0.95_rkx
  real(rkx) , parameter :: snal1 = 0.65_rkx
  !
  ! Short and long wave albedo for sea ice
  !
  real(rkx) , parameter :: sical0 = 0.6_rkx
  real(rkx) , parameter :: sical1 = 0.4_rkx

  public :: ocn_albedo

  contains
  !
  ! Albedo calculates fragmented albedos (direct and diffuse) in
  ! wavelength regions split at 0.7um.
  !
  subroutine ocn_albedo
    implicit none
!
    real(rkx) :: age , albg , albgl , albgld , albgs , albgsd , &
                 cf1 , cff , conn , cons ,       &
                 czeta , czf , dfalbl , dfalbs , dralbl , dralbs , sl , &
                 sl2 , sli , tdiff , tdiffs
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'ocn_albedo'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    do i = iocnbeg , iocnend
      czeta = czenith(i)
      !
      !================================================================
      !       2.   get albedo
      !================================================================
      !
      if ( mask(i) == 1 .or. mask(i) == 3 ) then
        ! ocean albedo depends on zenith angle
        if ( czeta >= d_zero ) then
          ! albedo independent of wavelength
          albg = 0.05_rkx/(czeta+0.15_rkx)
          albgs = albg
          albgl = albg
          albgsd = 0.08_rkx
          albgld = 0.08_rkx
        else
          albg = 0.05_rkx
          albgs = albg
          albgl = albg
          albgsd = 0.08_rkx
          albgld = 0.08_rkx
        end if
      else if ( mask(i) == 2 .or. mask(i) == 4 ) then
        ! Ice over ocean or lake
        tdiffs = sts(i) - tzero
        tdiff = max(tdiffs,d_zero)
        tdiffs = min(tdiff,20.0_rkx)
        albgl = sical1 - 1.1e-2_rkx*tdiffs
        albgs = sical0 - 2.45e-2_rkx*tdiffs
        albg = fsol1*albgs + fsol2*albgl
        albgsd = albgs
        albgld = albgl
        if ( sncv(i) > d_zero ) then
          ! Snow albedo depends on snow-age, zenith angle, and thickness
          ! of snow. snow albedoes for visible and ir solar rad visible
          ! albedo depends on snow age
          ! age gives reduction of visible rad snow albedo due to age
          cons = 0.2_rkx
          conn = 0.5_rkx
          age = (d_one-d_one/(d_one+snag(i)))
          ! sl helps control albedo zenith dependence
          sl = d_two
          sli = d_one/sl
          sl2 = d_two*sl
          ! snal0= new snow albedo for vis rad, sol zen le 6
          ! snal1= new snow albedo for long-wave rad
          dfalbs = snal0*(d_one-cons*age)
          ! czf corrects albedo of new snow for solar zenith
          cf1 = ((d_one+sli)/(d_one+sl2*czeta)-sli)
          cff = max(cf1,d_zero)
          czf = 0.4_rkx*cff*(d_one-dfalbs)
          dralbs = dfalbs + czf
          dfalbl = snal1*(d_one-conn*age)
          czf = 0.4_rkx*cff*(d_one-dfalbl)
          dralbl = dfalbl + czf
          !------------------------------------------
          ! 4.1  compute albedo for snow on sea ice
          !------------------------------------------
          albgs = (d_one-scvk(i))*albgs + dralbs*scvk(i)
          albgl = (d_one-scvk(i))*albgl + dralbl*scvk(i)
          albgsd = (d_one-scvk(i))*albgsd + dfalbs*scvk(i)
          albgld = (d_one-scvk(i))*albgld + dfalbl*scvk(i)
        end if
      else
        albgs = d_zero
        albgl = d_zero
        albgsd = d_zero
        albgld = d_zero
      end if
      !
      ! not part of albedo in the ccm
      !
      swdiral(i) = albgs
      lwdiral(i) = albgl
      swdifal(i) = albgsd
      lwdifal(i) = albgld
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine ocn_albedo

end module mod_ocn_albedo
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
