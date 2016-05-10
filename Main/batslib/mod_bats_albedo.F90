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

module mod_bats_albedo

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_runparams
  use mod_service
  use mod_bats_param
  use mod_bats_internal
  use mod_bats_drag

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

  logical :: ldesseas = .false.

  public :: albedo
  public :: ldesseas

  contains
!
! Albedo calculates fragmented albedos (direct and diffuse) in
! wavelength regions split at 0.7um.
!
! CM hands albedos to radiation package which computes
! rswf(i) = net solar absorbed over full grid square
! vegswab(j,i) = vegetation absorbed (full solar spectrum)
! solar(j,i) = shortwave  solar incident
!
! Here these are calculated at the end of albedo - they use only
! direct albedos for now
!
! in both versions :  lftemp uses vegswab
! tgrund uses vegswab & rswf(i) to get
! ground absorbed solar
! photosynthesis uses solar - see subrouts
! stomat and co2 (carbon)
!
! For sea, sea-ice veg albedos are not set these albedos are not
! treated as arrays here
!
! (depuv/10.0)= the ratio of upper soil layer to total root depth
! Used to compute "wet" for soil albedo
!
  subroutine albedo
    implicit none
    real(rkx) :: age , albg , albgl , albgld , albgs , albgsd , albl ,  &
                 albld , albs , albsd , alwet , cf1 , cff ,     &
                 conn , cons , czeta , czf , dfalbl , dfalbs , dralbl , &
                 dralbs , sl , sl2 , sli , wet
    ! real(rkx) :: sfac
    integer(ik4) :: kolour , i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'albedo'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! =================================================================
    ! 1. set initial parameters
    ! =================================================================
    !
    !
    ! Desert seasonal albedo
    ! Works for Sahara desert and generally northern emisphere
    ! In souther emisphere only some points have this class
    !
    if ( ldesseas ) then
      if ( xmonth == 1 .or. xmonth == 2 .or. xmonth == 12 ) then
        solour(1) = 0.12_rkx
      endif
      if ( xmonth == 3 .or. xmonth == 4 .or. xmonth == 5 ) then
        solour(1) = 0.15_rkx
      endif
      if ( xmonth == 6 .or. xmonth == 7 .or. xmonth == 8) then
        solour(1) = 0.18_rkx
      endif
      if ( xmonth == 9 .or. xmonth == 10 .or. xmonth == 11) then
        solour(1) = 0.15_rkx
      endif
    end if
    !
    ! In depth, wt is frac of grid square covered by snow;
    ! depends on average snow depth, vegetation, etc.
    !
    call depth
    !
    ! 1.2  set default vegetation and albedo
    !
    do i = ilndbeg , ilndend
      czeta = czenith(i)
      swal(i) = d_zero
      lwal(i) = d_zero
      albld = d_zero
      albsd = d_zero
      !
      !================================================================
      !       2.   get albedo
      !================================================================
      !
      ! sfac = d_one - aseas(i)
      !
      ! ccm tests here on land mask for veg and soils data
      ! reduces albedo at low temps !!!!!
      ! should respond to moisture too (commented out) (pat, 27 oct 86)
      ! lncl(i) = lncl(iveg1(i)) - seasf(iveg1(i)) * sfac
      !
      albs = albvgs(lveg(i))
      albl = albvgl(lveg(i))
      !---------------------------------------------------------------
      if ( (lveg(i) < 12) .or. (lveg(i) > 15) ) then
        ! 2.1  bare soil albedos
        !      (soil albedo depends on moisture)
        kolour = kolsol(lveg(i))
        wet = ssw(i)/depuv(lveg(i))
        alwet = max((11.0_rkx-40.0_rkx*wet),d_zero)*0.01_rkx
        alwet = min(alwet,solour(kolour))
        albg = solour(kolour) + alwet
        albgs = albg
        albgl = d_two*albg
        ! higher nir albedos set diffuse albedo
        albgld = albgl
        albgsd = albgs
        albsd = albs
        albld = albl
        swal(i) = albs
        lwal(i) = albl
      else if ( lveg(i) == 12 ) then
        ! 2.2   permanent ice sheet
        albgs = 0.8_rkx
        albgsd = 0.8_rkx
        albgl = 0.55_rkx
        albgld = 0.55_rkx
      else
        ! 2.3  inland water, swamps, rice paddies etc.
        albg = 0.05_rkx/(czeta+0.15_rkx)
        albgs = albg
        albgsd = albg
        albgl = albg
        albgld = albg
      end if
      ! ==============================================================
      ! 4.  correct for snow cover
      ! ==============================================================
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
        if ( lncl(i) > 0.001_rkx ) then
          ! effective albedo over vegetation with snow
          albl = (d_one-wt(i))*albl + dralbl*wt(i)
          albld = (d_one-wt(i))*albld + dfalbl*wt(i)
          albs = (d_one-wt(i))*albs + dralbs*wt(i)
          albsd = (d_one-wt(i))*albsd + dfalbs*wt(i)
        end if
        !----------------------------------------------------------------
        !         4.1  compute albedo for snow on bare ground
        !----------------------------------------------------------------
        albgs = (d_one-scvk(i))*albgs + dralbs*scvk(i)
        albgl = (d_one-scvk(i))*albgl + dralbl*scvk(i)
        albgsd = (d_one-scvk(i))*albgsd + dfalbs*scvk(i)
        albgld = (d_one-scvk(i))*albgld + dfalbl*scvk(i)
      end if
      !
      ! not part of albedo in the ccm
      !
      swdiral(i) = (d_one-lncl(i))*albgs + lncl(i)*albs
      lwdiral(i) = (d_one-lncl(i))*albgl + lncl(i)*albl
      swdifal(i) = (d_one-lncl(i))*albgsd + lncl(i)*albsd
      lwdifal(i) = (d_one-lncl(i))*albgld + lncl(i)*albld
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine albedo

end module mod_bats_albedo
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
