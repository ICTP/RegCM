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
    real(rkx) :: age , albg , albgl , albgld , albgs , albgsd , &
                 cf1 , cff , conn , cons , czeta , czf , sl2 ,  &
                 dfalbl , dfalbs , dralbl , dralbs , sl , sli , &
                 tdiff , tdiffs , scrat , scvk , wfac , onemc , &
                 w0 , wspd
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
        if ( iwhitecap == 1 ) then
          wspd = um10(i)
          ! Monahan and O'Muircheartaigh [1980]
          ! Fraction of whitecapping function of windspeed
          wfac = 2.95e-6_rkx * wspd**3.52
          ! Ocean albedo depends on zenith angle.
          ! Solar zenith dependence from Briegleb et al., [1986]
          albg = 0.026_rkx / (czeta**1.7_rkx + 0.065_rkx) + &
                 0.15_rkx * (czeta-0.1_rkx)*(czeta-0.5_rkx)*(czeta-1.0_rkx)
          if ( czeta > 0.01_rkx .and. czeta < 0.12_rkx ) then
            ! Katsaros et al [1985] , reduction by wind waves at low angles.
            albg = max(0.05_rkx,albg*(d_one - 0.0036_rkx*wspd))
          end if
          ! Koepke [1984] - Increase by whitecapping
          albg = albg + 0.22_rkx * wfac
          albgs = albg
          albgl = albg
          albgsd = 0.05_rkx + 0.11 * wfac
          albgld = 0.05_rkx
        else if ( iwhitecap == 2 ) then
          wspd = um10(i)
          onemc = 1.0_rkx-czeta
          w0 = 180.0_rkx*(czeta**3)*onemc**2
          wfac = 3.84_rkx * 1.0e-6_rkx * wspd**3.41_rkx
          if ( wspd > w0 ) then
            albg = 0.021_rkx + 0.0421_rkx*onemc**2 + &
                   0.128_rkx*onemc**3 - 0.04_rkx*onemc**6 + &
                   (4.0_rkx / (5.68_rkx+wspd-w0) + (0.074_rkx * &
                   onemc)/(1.0_rkx+3.0_rkx*(wspd-w0))) * onemc**6
          else
            albg = (1.0_rkx + (5.4_rkx*czeta**2*onemc**2*wspd * &
                               (wspd-1.1_rkx*w0)**2)/w0**3) * &
                   (0.021_rkx + 0.0421_rkx*onemc**2 + &
                    0.128_rkx*onemc**3 - 0.04*onemc**6 + &
                    (4.0_rkx/5.68_rkx + 0.074_rkx*onemc)*onemc**6)
          end if
          albg = (1.0_rkx-wfac)*albg + 0.3_rkx*wfac
          albgs = albg
          albgl = 0.05_rkx/(czeta+0.15_rkx)
          albgsd = 0.022_rkx*(1.0_rkx+0.55_rkx*exp(-(wspd/7.0_rkx)**2) + &
                   1.45_rkx*exp(-(wspd/40.0_rkx)**2))
          albgsd = (1.0_rkx-wfac)*albgsd + 0.3_rkx*wfac
          albgld = 0.08_rkx
        else
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
        end if
      else
        ! Ice over ocean or lake
        tdiffs = t2m(i) - icetriggert
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
          scrat = sncv(i)*0.01_rkx/(d_one+d_three*age)
          scvk = scrat/(0.1_rkx+scrat)
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
          albgs = (d_one-scvk)*albgs + dralbs*scvk
          albgl = (d_one-scvk)*albgl + dralbl*scvk
          albgsd = (d_one-scvk)*albgsd + dfalbs*scvk
          albgld = (d_one-scvk)*albgld + dfalbl*scvk
        end if
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
