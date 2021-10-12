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

module mod_bats_bndry

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_service
  use mod_bats_leaftemp
  use mod_bats_param
  use mod_bats_drag
  use mod_bats_internal
  use mod_constants

  implicit none

  private

  public :: soilbc , bndry

  real(rkx) , parameter :: rainsnowtemp = 2.2_rkx
  real(rkx) , parameter :: xnu = twopi/secpd

  contains

  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  !
  ! this subrout overwrites many of the soil constants
  ! as a function of location(jlon,jlat)
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  !
  subroutine soilbc
    implicit none
    real(rkx) :: ck , dmax , dmin , dmnor , phi0 , tweak1
    integer(ik4) :: i
    !
    ! ================================================================
    ! new soils data as a fn of texture make porosity, soil suction,
    ! hydraul conduc, wilting frac variables rather than consts
    ! relfc is the ratio of field capacity to saturated water content,
    ! defined so the rate of gravitational drainage at field
    ! capacity is assumed to be 2 mm/day (baver et al., 1972)
    ! ===============================================================
    !
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'soilbc'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    do i = ilndbeg , ilndend
      freza(lveg(i)) = 0.15_rkx*deprv(lveg(i))
      frezu(lveg(i)) = 0.15_rkx*depuv(lveg(i))
      texrat(i) = skrat(ltex(i))
      porsl(i) = xmopor(ltex(i))
      xkmx(i) = xmohyd(ltex(i))
      bsw(i) = bee(ltex(i))
      bfc(i) = 5.8_rkx - bsw(i)*(0.8_rkx+0.12_rkx*(bsw(i)-d_four) * &
               log10(1.0e2_rkx*xkmx(i)))
      phi0 = xmosuc(ltex(i))
      dmax = bsw(i)*phi0*xkmx(i)/porsl(i)
      dmin = 1.0e-3_rkx
      dmnor = 1550.0_rkx*dmin/dmax
      tweak1 = (bsw(i)*(bsw(i)-6.0_rkx)+10.3_rkx) / &
               (bsw(i)*bsw(i)+40.0_rkx*bsw(i))
      ck = (d_one+dmnor)*tweak1*0.23_rkx/0.02356_rkx
      evmx0(i) = 1.02_rkx*dmax*ck/sqrt(depuv(lveg(i))*deprv(lveg(i)))
      gwmx0(i) = depuv(lveg(i))*porsl(i)
      gwmx1(i) = deprv(lveg(i))*porsl(i)
      gwmx2(i) = deptv(lveg(i))*porsl(i)
      wiltr(i) = xmowil(ltex(i))
      ! force irrigated crop to be at field capacity
      relfc(i) = xmofc(ltex(i))
      ! Imported Lara Kuepper's Irrigated Crop modification from RegCM3
      ! see Kueppers et al. (2008)
      ! relaw is between field capacity and wilting point
      ! relaw(i) = 0.75_rkx*(xmofc(ltex(i))-xmowil(ltex(i)))+xmowil(ltex(i))
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine soilbc
  !
  !=======================================================================
  !l  based on: bats version 1e          copyright 18 august 1989
  !=======================================================================
  !
  !     This is the main routine when interfacing with a
  !     meteorological model
  !
  !                f l o w   d i a g r a m   f o r   b n d r y
  !
  !                bndry ===>    drag ===> dragdn ===> depth
  !                             satur
  !                            vcover
  !                              drip
  !                            lftemp ======================> stomat
  !                               co2 ===> carbon             frawat
  !                            tgrund                           root
  !                              snow                          satur
  !                             water                         lfdrag
  !                                                           condch
  !                                                           condcq
  !                                                            deriv
  !  **  type1  = crop
  !  **  type2  = short grass
  !  **  type3  = evergreen needle leaf tree
  !  **  type4  = deciduous needle leaf tree
  !  **  type5  = deciduous broadleaf tree
  !  **  type6  = evergreen brodaleaf tree
  !  **  type7  = tall grass
  !  **  type8  = desert
  !  **  type9  = tundra
  !  **  type10 = irrig crop
  !  **  type11 = semi-desert
  !  **  type12 = ice
  !  **  type13 = bog or marsh
  !  **  type14 = inland water
  !  **  type15 = sea
  !  **  type16 = evgr shrub
  !  **  type17 = decid shrub
  !  **  type18 = mixed tree
  !
  !  ** note: water and soil parameters are in mm
  !
  subroutine bndry
    implicit none
    real(rkx) :: qsatd , rai
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'bndry'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    !=======================================================================
    !   1.   initialize
    !=======================================================================

    do i = ilndbeg , ilndend
      htvp(i) = wlh(tgrd(i))
      flnet(i) = d_zero
      fevpg(i) = d_zero
      fseng(i) = d_zero
      vegt(i) = d_zero
      efpr(i) = d_zero
      etr(i) = d_zero
      pw(i) = d_zero
      evapw(i) = d_zero
      ps(i) = d_zero
      evaps(i) = d_zero
      ! switch between rain and snow /tm is ref. temp set= anemom temp - 2.2
      tm(i) = sts(i) - rainsnowtemp
      ! soil moisture ratio (to max) as used in subrouts tgrund,
      ! water, and root (called by lftemp): watu=upper, watr=root,
      ! watt=total
      watu(i) = ssw(i)/gwmx0(i)
      watr(i) = rsw(i)/gwmx1(i)
      watt(i) = tsw(i)/gwmx2(i)
      watu(i) = min(watu(i),d_one)
      watr(i) = min(watr(i),d_one)
      watt(i) = min(watt(i),d_one)
      watr(i) = max(watr(i),minwrat)
      watu(i) = max(watu(i),minwrat)
      watt(i) = max(watt(i),minwrat)
    end do
    !
    !=======================================================================
    !   2.  calculate transfer coefficients at layer 1
    !=======================================================================
    ! 2.1  sigf corrected in drag: remove snow-covered veg
    call dragc
    ! 2.2  get saturation vapor pressure of soil surface
    call satur(qgrd,tgrd,sfcp)

    !=======================================================================
    ! 3.   bare land
    !=======================================================================
    ! 3.1  get derivative of fluxes with repect to tg
    !
    do i = ilndbeg , ilndend
      if ( sigf(i) <= minsigf ) then
        qsatd = pfqsdt(tgrd(i),sfcp(i))
        qsatd = qsatd * gwet(i)
        ! call bats_qsdt(tgrd(i),qgrd(i),qsatd)
        ! qsatd = qsatd * gwet(i)
        rai = cdrx(i)*vspda(i)*rhs(i)
        cgrnds(i) = rai*cpd
        cgrndl(i) = rai*qsatd
        cgrnd(i) = cgrnds(i) + cgrndl(i)*htvp(i)
        ! 3.2  sensible and latent fluxes using soil temperatures
        ! from previous time step
        delq(i) = (qs(i)-qgrd(i)) * gwet(i)
        delt(i) = sts(i) - tgrd(i)
        evpr(i) = -rai*delq(i)
        sent(i) = -cgrnds(i)*delt(i)
        if ( abs(sent(i)) < dlowval ) sent(i) = d_zero
        if ( abs(evpr(i)) < dlowval ) evpr(i) = d_zero
        ! 3.3  fluxes to subrout tgrund (evap is in kg/m**2/s)
        fseng(i) = sent(i)
        fevpg(i) = evpr(i)
        ! 3.4  equate canopy to air, for temp, wind over bare grnd;
        ! needed as factors of sigf(=0) in subr water (uaf) and
        ! subr drag (tlef(i) carried over to next tstep).
        tlef(i) = sts(i)
        uaf(i) = vspda(i)
      end if
    end do
    do i = ilndbeg , ilndend
      if ( sigf(i) > minsigf ) then   !  check each point
        !================================================================
        ! 4.   vegetation
        !================================================================
        ! 4.1  add precipitation to leaf water
        ldew(i) = ldew(i)+dtbat*sigf(i)*prcp(i)
        ldew(i) = max(ldew(i),d_zero)
      else
        ldew(i) = d_zero
      end if
    end do
    !
    ! 4.2  distribute excess leaf water to soil
    call vcover
    call drip
    !
    ! 4.3  calculate canopy temperature, soil and total fluxes,
    !      and leaf water change by evapotranspiration
    call lftemp
    !
    ! 4.4  calculate carbon sources and sinks
    !call co2
    !
    ! 5.2  calculate soil temp and surface hydrology
    call tgrund
    call snow
    call water

    !=======================================================================
    !   6.   linkage to meteorological model
    !=======================================================================

    do i = ilndbeg , ilndend
      ! 6.1  rate of momentum transfer per velocity
      drag(i) = cdrx(i)*vspda(i)*rhs(i)
      if ( sigf(i) <= minsigf ) taf(i) = tgrd(i)
      ! 6.2  parameters for temperature difference at anemometer level
      ! zdelt(i) = zdelt(i)*delt(i)
      ! 6.4  evaporative flux, accounting for sublimation
      ! evprr(i) = htvp(i)*(evpr(i)-fevpg) + htvp(i)*fevpg
      ! 6.5  nondimensional equivalent bucket capacity for comparisons
      !      with bucket models; usually 1 or less, except where
      !      saturated (then around 2)
      ! rh2ox(i) = (watr-wiltr)/(relfc-wiltr)
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
    contains

#include <wlh.inc>
#include <pfesat.inc>
#include <pfdesatdt.inc>
#include <pqderiv.inc>

  end subroutine bndry
  !
  !=======================================================================
  ! VCOVER
  !     Provides leaf and stem area parameters;
  !     depends on climate through subsoil temperatures.
  !=======================================================================
  !
  subroutine vcover
    implicit none
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'vcover'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    do i = ilndbeg , ilndend
      if ( sigf(i) > minsigf ) then
        xlai(i) = xla(lveg(i))
        xlai(i) = xlai(i) + (xlai0(lveg(i))-xlai(i))*(d_one-aseas(i))
        rlai(i) = xlai(i) + sai(lveg(i))
        xlsai(i) = xlai(i) + sai(lveg(i))
        vegt(i) = sigf(i)*xlsai(i)
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine vcover
  !
  !=======================================================================
  !  DRIP
  !     Excess leaf water is put into rain or snow;
  !     leaf water is reset to its maximum value.
  !=======================================================================
  !
  subroutine drip
    implicit none
    integer(ik4) :: i
    real(rkx) :: xrun
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'drip'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    do i = ilndbeg , ilndend
      etrrun(i) = d_zero
      sdrop(i) = d_zero
      if ( sigf(i) > minsigf ) then
        ! xrun = leaf drip ; sdrop = snow drop off foliage
        xrun = ldew(i) - dewmax*vegt(i)
        ! test on maximum value of dew
        if ( xrun > d_zero ) then
          etrrun(i) = xrun
          ldew(i) = dewmax*vegt(i)
        end if
        ! below freezing excess leaf water falls as snow
        if ( (xrun > d_zero) .and. (tm(i) < tzero) ) then
          etrrun(i) = d_zero
          sdrop(i) = xrun
        end if
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine drip
  !
  !=======================================================================
  ! WATER
  !     update soil moisture and runoff
  !
  !     new algorithms for three soil layers (dickinson & kennedy 8-88)
  !     calculate fluxes through air, surface layer, and root layer faces
  !
  !          b = b of clapp and hornberger
  !       est0 = soil water flux, out top
  !      gwatr = net input of water to the soil surface
  !       ircp = leaf interception
  !     wflux1 = soil water flux, 10 cm
  !     wflux2 = soil water flux, 1 m
  !     rsubss = soil water flux by grav. drainage thru 10 cm interface
  !     rsubsr = soil water flux by grav. drainage thru 1 m interface
  !     rsubst = soil water flux by grav. drainage thru 10 m interface
  !     trnof(i) = total runoff mm
  !     srnof(i) = surface runoff mm
  !
  !     xkmxr and wflux1 determine flow thru upper/root soil interface
  !     evmxt, xkmx1, and xkmx2 determine flow thru lower interfaces
  !
  !     veg type 10 "irrigated crop" is irrigated through reducing
  !          the runoff (srnof), i.e., by adding a negative number
  !          if the land isn't at least 60% saturated.
  !     veg type 13 and 14 are water covered (lake, swamp, rice paddy);
  !          negative runoff keeps this land saturated.
  !=======================================================================
  !
  subroutine water
    implicit none
    real(rkx) :: b , bfac , bfac2 , delwat , est0 , evmax , evmxr , &
               evmxt , rap , vakb , wtg2c , xxkb , gwatr , rsubsr , &
               rsubss , xkmx1 , xkmx2 , xkmxr , wata , b1 , b2 , b3
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'water'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !=======================================================================
    !   1.   define soil water fluxes
    !=======================================================================
    !
    do i = ilndbeg , ilndend
      !
      ! 1.1  reduce infiltration for frozen ground
      !
      if ( tgbrd(i) > tzero ) then
        xkmxr = xkmx(i)
      else
        xkmxr = d_zero
      end if
      !
      ! 1.11 permafrost or ice sheet
      !
      if ( lveg(i) == 9 .or. lveg(i) == 12 ) then
        xkmx1 = d_zero
        xkmx2 = d_zero
      else
        xkmx1 = xkmx(i)
        xkmx2 = drain
      end if
      !
      ! 1.2  diffusive fluxes
      !
      evmxr = evmx0(i)*xkmxr/xkmx(i)
      evmxt = evmx0(i)*xkmx1/xkmx(i)
      b = bsw(i)
      bfac = watr(i)**(d_three+bfc(i)) * watu(i)**(b-bfc(i)-d_one)
      bfac2 = watt(i)**(d_two+bfc(i)) * watr(i)**(b-bfc(i))
      wfluxc(i) = evmxr*(depuv(lveg(i))/deprv(lveg(i)))**0.4_rkx*bfac
      wflux1(i) = wfluxc(i)*(watr(i)-watu(i))
      wflux2(i) = evmxt*sqrt(depuv(lveg(i))/deprv(lveg(i)))* &
                  bfac2*(watt(i)-watr(i))
      !
      ! 1.3  gravitational drainage
      !
      b1 = b+d_half
      b2 = b+2.5_rkx
      b3 = d_two*b+d_three
      rsubss = xkmxr*(watr(i)**b1)*(watu(i)**b2)
      rsubsr = xkmx1*(watt(i)**b1)*(watr(i)**b2)
      rsubst(i) = max(d_zero,xkmx2*(watt(i)**b3))
      !
      ! 1.32 bog
      !
      if ( lveg(i) == 13 ) then
        rsubst(i) = d_zero
        rsubss = d_zero
        rsubsr = d_zero
      end if
      !
      ! 1.4  fluxes through internal surfaces
      !
      wflux1(i) = wflux1(i) - rsubss
      wflux2(i) = wflux2(i) - rsubsr
    end do
    !
    ! 1.5  net flux at air interface
    !
    do i = ilndbeg , ilndend
      gwatr = pw(i) + sm(i) + etrrun(i)/dtbat - evapw(i)
      !
      !=================================================================
      ! 2.   define runoff terms
      !=================================================================
      !
      ! 2.1  surface runoff
      !
      wata = d_half*(watu(i)+watr(i))
      !
      ! 2.11 increase surface runoff over frozen ground
      !
      srnof(i) = d_zero
      if ( tgrd(i) < tzero ) then
        srnof(i) = min(d_one,wata**1) * max(d_zero,gwatr)
      else
        srnof(i) = min(wata**4,1.0_rkx) * max(0.0_rkx,gwatr)
      end if
      !
      ! 2.12 irrigate cropland
      !
      if ( lveg(i) == 10 .and. watr(i) < relfc(i) ) then
        srnof(i) = srnof(i) + min(d_zero,(rsw(i)-relfc(i)*gwmx1(i))/dtbat)
      end if
      ! Imported Lara Kuepper's Irrigated Crop modification from RegCM3
      ! see Kueppers et al. (2008)
      ! if ( lveg(i) == 10 .and. watr(i) < relaw(i) ) then
      !   srnof(i) = srnof(i) + (rsw(i)-relaw(i)*gwmx1(i))/dtbat
      ! end if
      !
      ! 2.13 saturate swamp or rice paddy
      !
      if ( lveg(i) == 13 ) then
        srnof(i) = srnof(i) + min(d_zero,(rsw(i)-relfc(i)*gwmx1(i))/dtbat)
      end if
      !
      ! 2.2  total runoff
      !
      srnof(i) = max(d_zero,srnof(i))
      trnof(i) = max(srnof(i)+rsubst(i),d_zero)
      !
      !=================================================================
      !         3.   increment soil moisture
      !=================================================================
      !
      ! 3.1  update top layer with implicit treatment of flux from below
      !
      gwatr = gwatr - efpr(i)*etr(i)
      ssw(i) = ssw(i) + dtbat * (max(gwatr-srnof(i),d_zero) + wflux1(i))
      ssw(i) = ssw(i) / (d_one + dtbat * wfluxc(i)/gwmx0(i))
      ssw(i) = max(ssw(i),gwmx0(i)*minwrat)
      !
      ! 3.2  update root zone
      !
      rsw(i) = rsw(i) + dtbat * (max(gwatr-srnof(i),d_zero) + wflux2(i))
      rsw(i) = max(rsw(i),gwmx1(i)*minwrat)
      !
      ! 3.3  update total water
      !
      tsw(i) = tsw(i) + dtbat * (max(gwatr-srnof(i),d_zero) - rsubst(i))
      tsw(i) = max(tsw(i),gwmx2(i)*minwrat)
    end do

    do i = ilndbeg , ilndend
      !
      !=================================================================
      !  4.   check whether soil water exceeds maximum capacity or
      !       becomes negative (should rarely or never happen)
      !=================================================================
      !
      ! 4.1  surface water assumed to move downward into soil
      !
      if ( ssw(i) > gwmx0(i) ) then
        delwat = ssw(i) - gwmx0(i)
        ssw(i) = gwmx0(i)
        rsw(i) = rsw(i) + delwat
      end if
      !
      ! 4.2  excess root layer water assumed to move downward
      !
      if ( rsw(i) > gwmx1(i) ) then
        delwat = rsw(i) - gwmx1(i)
        rsw(i) = gwmx1(i)
        tsw(i) = tsw(i) + delwat
      end if
      !
      ! 4.3  excess total water assumed to go to subsurface runoff
      !
      if ( tsw(i) > gwmx2(i) ) then
        delwat = tsw(i) - gwmx2(i)
        tsw(i) = gwmx2(i)
        trnof(i) = trnof(i) + delwat/dtbat
      end if
    end do
    !
    !=======================================================================
    !   7.   calculate potential evaporation and use it to determine
    !        wetness factor, allowing for snow being saturated
    !=======================================================================
    !
    do i = ilndbeg , ilndend
      xxkb = min(rough(lveg(i)),d_one)
      vakb = (d_one-sigf(i))*vspda(i) + sigf(i) * &
             (xxkb*uaf(i)+(d_one-xxkb)*vspda(i))
      wtg2c = (d_one-sigf(i))*cdrx(i)*vakb
      rap = rhs(i) * (csoilc*uaf(i)*sigf(i)*(qgrd(i)+delq(i)-qs(i)) + &
                      wtg2c*(qgrd(i)-qs(i)))
      bfac = watr(i)**(d_three+bfc(i)) * watu(i)**(bsw(i)-bfc(i)-d_one)
      est0 = evmx0(i)*bfac*watu(i)
      evmax = max(est0,d_zero)
      gwet(i) = min(d_one,evmax/max(1.0e-14_rkx,rap))
      gwet(i) = scvk(i) + gwet(i)*(d_one-scvk(i))
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine water
  !
  !=======================================================================
  ! SNOW
  !     update snow cover and snow age
  !
  !     three-part if block:
  !       if snow cover < 0, then snow cover and snow age = 0
  !       if antarctica, snow age = 0 (katabatic winds keep snow fresh)
  !       if elsewhere, snow age follows given formulae
  !
  !        ps = snow precipitation rate
  !     evaps = moisture flux from ground to atmosphere
  !        sm = snow melt rate
  !     sdrop = snow fallen from vegetation
  !
  !     aging of snow consists of three factors:
  !           age1: snow crystal growth
  !           age2: surface melting
  !           age3: accumulation  of other particles, soot, etc., which
  !                      is small in southern hemisphere
  !
  !=======================================================================
  !
  subroutine snow
    implicit none
    real(rkx) :: age1 , age2 , arg , arg2 , dela , dela0 , &
                 dels , tage , sold! , fsnts , tpw
    integer(ik4) :: i
    real(rkx) , parameter :: age3 = 0.3_rkx
    !real(rkx) , parameter :: ax = -48.23_rkx
    !real(rkx) , parameter :: bx = 0.75_rkx
    !real(rkx) , parameter :: cx = 1.16_rkx
    !real(rkx) , parameter :: dx = 1.02_rkx
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'snow'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    !=======================================================================
    !   1.   partition soil evaporation and precipitation
    !        between water and snow
    !=======================================================================
    !
    do i = ilndbeg , ilndend
      evaps(i) = scvk(i) * fevpg(i)
      evapw(i) = fevpg(i) - evaps(i)
      ! tm is temperature of precipitation
      ! Aiguo Dai
      ! Temperature and pressure dependence of the rain-snow phase
      ! transition over land and ocean
      ! GRL, VOL. 35, L12802, doi:10.1029/2008GL033295, 2008
      !tpw = prcp(i)*(d_one-sigf(i))
      !if ( tm(i)-tzero < 6.0_rkx ) then
      !  if ( tm(i)-tzero < -6.0_rkx ) then
      !    fsnts = 1.0_rkx
      !  else
      !    fsnts = ax * (tanh(bx*(tm(i)-tzero-cx))-dx)
      !  end if
      !else
      !  fsnts = d_zero
      !end if
      !ps(i) = tpw*fsnts
      !pw(i) = tpw-ps(i)
      if ( tm(i) >= tzero ) then
        pw(i) = prcp(i)*(d_one-sigf(i))
        ps(i) = d_zero
      else ! snowing
        pw(i) = d_zero
        ps(i) = prcp(i)*(d_one-sigf(i))
      end if
    end do
    !
    !=======================================================================
    !   2.   update snow cover
    !=======================================================================
    !
    do i = ilndbeg , ilndend
      sold = sncv(i)
      sncv(i) = sncv(i) + dtbat*(ps(i)-evaps(i)-sm(i)) + sdrop(i)
      if ( sncv(i) < dlowval ) then
        sncv(i) = d_zero
        snag(i) = d_zero
      end if
      ! snow cover except for antarctica
      !==================================================================
      !  3.   increment non-dimensional "age" of snow;
      !       10 mm snow restores surface to that of new snow.
      !==================================================================
      if ( sncv(i) > d_zero ) then
        arg = 5.0e3_rkx*(d_one/tzero-d_one/tgrd(i))
        age1 = exp(arg)
        arg2 = max(min(d_zero,d_10*arg),-25.0_rkx)
        age2 = exp(arg2)
        tage = age1 + age2 + age3
        dela0 = 1.0e-6_rkx*dtbat
        dela = dela0*tage
        dels = d_r10*max(d_zero,sncv(i)-sold)
        snag(i) = (snag(i)+dela)*(d_one-dels)
        if ( snag(i) < dlowval ) snag(i) = d_zero
      end if
      ! antarctica
      if ( sncv(i) > 800.0_rkx ) snag(i) = d_zero
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine snow
  !
  !=======================================================================
  ! TGRUND
  !
  !     Present version of diurnal and seasonal force restore (red 7-88)
  !     based on Dickinson (1988) force restore paper in j. climate.
  !     in particular, shows that a 0.1 m or thicker surface layer will
  !     dominate thermal conductivity for surface temperature over the
  !     diurnal cycle, and the first 1 m gives appropriate conductivity
  !     over annual cycle.
  !
  !     for snow cover, use weighted average of soil and snow properties.
  !     asymptotes to snow values for snow depth large compared to
  !     daily and annual temperature waves, respectively.
  !
  !     the term fct2 provides latent heat for the freezing or
  !     thawing of ground over temperatures between -4 and 0 deg c;
  !     this range may be changed if the 0.25 in the arithmetic fcn
  !     fct1=1/range and the range for subsoil temperature freezing
  !     are correspondingly changed.
  !
  !      secpd = seconds in a day
  !       wlhv = latent heat of vaporization of water
  !       wlhs = latent heat of sublimation of water
  !       wlhf = latent heat of freezing of water
  !       htvp =  wlhv or = wlhs if snow or tg<zero
  !     depsnw = depth of snow
  !     xdtime = nondimensional diurnal time step
  !     dtimea = nondimensional annual time step
  !     fsk(x) = soil conductivity fcn (x = v(h2o)/v(soil+pores))
  !     fsc(x) = soil heat capacity fcn (x = v(h2o)/v(soil+pores))
  !         hs = net energy flux into the surface
  !         sg = solar flux absorbed by bare ground
  !     skd,ska= diurnal and annual thermal diffusivities
  !=======================================================================
  !
  subroutine tgrund
    implicit none
    real(rkx) :: bcoefd , bcoefs , c31 , c3t , c41 , c4t , cder , depr , &
             depu , xdt2 , xdtime , dtimea , froze2 , frozen , rscss ,   &
             tbef , tg , tinc , wtas , wtax , wtd , wtds , hs , depann , &
             depdiu , rscsa , rscsd , ska , skd , sks
    real(rkx) :: dtbat2 , rdtbat2 , xlexp , xnua , swtrta , swtrtd
    integer(ik4) :: i
    real(rkx) , parameter :: xkperi = 1.4e-6_rkx
    real(rkx) , parameter :: t3 = 271.0_rkx ! permafrost temperature
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'tgrund'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    dtbat2  = dtbat*d_two
    rdtbat2 = d_one/dtbat2

    !=======================================================================
    ! 1.   define thermal conductivity, heat capacity,
    !      and other force restore parameters
    !=======================================================================

    xnua = xnu/dayspy
    xdtime = dtbat*xnu
    dtimea = dtbat*xnua
    xdt2 = d_half*xdtime

    do i = ilndbeg , ilndend
      ! 1.1  frozen ground values using 44 m2/yr for frozen ground
      !      thermal diffusion coefficient, based on the values of
      !      50 and 38 quoted by osterkamp; ice contribution to
      !      specific heat only 0.49 that of water
      swtrtd = watu(i)*porsl(i)
      if ( tgrd(i) < tzero ) then
        frozen = 0.85_rkx*min(d_one,d_rfour*(tzero-tgrd(i)))
        skd = xkperi
        rscsd = fsc(swtrtd*(d_one-0.51_rkx*frozen))
      else
        skd = fsk(swtrtd)*texrat(i)
        rscsd = fsc(swtrtd)
      end if
      swtrta = watr(i)*porsl(i)
      if ( tgbrd(i) < tzero ) then
        froze2 = 0.85_rkx*min(d_one,d_rfour*(tzero-tgbrd(i)))
        ska = xkperi
        rscsa = fsc(swtrta*(d_one-0.51_rkx*froze2))
      else
        ska = fsk(swtrta)*texrat(i)
        rscsa = fsc(swtrta)
      end if
      ! 1.2  correct for snow cover, if significant
      depdiu = sqrt(d_two*skd/xnu)
      bcoef(i) = xdtime*depdiu/(rscsd*skd)
      sks = d_zero
      if ( scrat(i) > 0.001_rkx ) then
        xlexp = d_two*scrat(i)/depdiu
        ! Graziano : Limit exponential argument
        if ( xlexp < 25.0_rkx ) then
          wtd = exp(-xlexp)
        else
          wtd = d_zero
        end if
        rscss = csnw*rhosw(i)
        sks = 7.0e-7_rkx*cws*rhosw(i)
        bcoefs = sqrt(d_two*sks/xnu)/(rscss*sks)
        wtds = (d_one-wtd)*scvk(i)
        bcoefd = sqrt(d_two*skd/xnu)/(rscsd*skd)
        bcoef(i) = xdtime*(wtds*bcoefs+(d_one-wtds)*bcoefd)
        depdiu = wtds*sqrt(d_two*sks/xnu) + (d_one-wtds)*depdiu
      end if
      depann = sqrt(d_two*ska/xnua)
      if ( scrat(i) > 0.02_rkx ) then
        xlexp = d_two*scrat(i)/depann
        if ( xlexp < 25.0_rkx ) then
          wtax = exp(-xlexp)
        else
          wtax = d_zero
        end if
        wtas = (d_one-wtax)*scvk(i)
        depann = wtas*sqrt(d_two*sks/xnua) + (d_one-wtas)*depann
      end if
      deprat(i) = depann/depdiu
      !=================================================================
      !         2.   collect force restore terms
      !=================================================================
      cc(i) = d_one
      fct2(i) = d_zero
      !
      ! 2.1  add freezing thermal inertia
      if ( (tgrd(i) < tzero) .and. (tgrd(i) > (tzero-d_four)) ) then
        depu = depuv(lveg(i))*d_r1000
        cc(i) = d_one + max(ssw(i)-frezu(lveg(i)),d_zero)*fct1(depu*rscsd)
      end if
      if ( (tgbrd(i) < tzero) .and. (tgbrd(i) > (tzero-d_four)) ) then
        depr = deprv(lveg(i))*d_r1000
        fct2(i) = max(rsw(i)-freza(lveg(i)),d_zero) * fct1(depr*rscsa)
      end if
      ! 2.2  large thermal inertial for permanent ice cap
      if ( lveg(i) == 12 ) fct2(i) = d_1000*fct2(i)
      ! 2.3  collect energy flux terms
      rnet(i) = swflx(i) - sigf(i)*(abswveg(i)-flnet(i)) - &
                (d_one-sigf(i))*(lwflx(i) - sigf(i)*flneto(i))
      hs = rnet(i) - fseng(i) - fevpg(i)*htvp(i)
      bb(i) = bcoef(i)*hs + xdtime*tgbrd(i)
      ! 2.4  add in snowmelt (melt enough snow to reach freezing temp)
      sm(i) = d_zero
      if ( sncv(i) > d_zero ) then
        cder = bcoef(i)*cgrnd(i)
        sm(i) = (bb(i) + (cc(i)-xdt2+cder)*tgrd(i) - &
                         (cc(i)+xdt2+cder)*tzero) / (bcoef(i)*wlhf)
        ! snow melt always between 0 and total snow
        sm(i) = max(d_zero,min(sm(i),sncv(i)*d_two*rdtbat2))
        bb(i) = bb(i) - bcoef(i)*wlhf*sm(i)
      end if
    end do

    !=======================================================================
    !   3.   update soil temperatures
    !=======================================================================
    !   3.1  update surface soil temperature
    do i = ilndbeg , ilndend
      tbef = tgrd(i)
      cder = bcoef(i)*cgrnd(i)
      tg = (bb(i)+(cc(i)-xdt2+cder)*tgrd(i)) / (cc(i)+xdt2+cder)
      tgrd(i) = tg
      ! 3.2  put brakes on large temperature excursions
      tgrd(i) = min(tbef+d_10,tgrd(i))
      tgrd(i) = max(tbef-d_10,tgrd(i))
      ! 3.3  correct fluxes to present soil temperature
      tinc = tgrd(i) - tbef
      sent(i) = sent(i) + tinc*cgrnds(i)
      evpr(i) = evpr(i) + tinc*cgrndl(i)
      if ( abs(sent(i)) < dlowval ) sent(i) = d_zero
      if ( abs(evpr(i)) < dlowval ) evpr(i) = d_zero
      fseng(i) = fseng(i) + tinc*cgrnds(i)
      fevpg(i) = fevpg(i) + tinc*cgrndl(i)
      !
      ! 3.5  couple to deep temperature in permafrost
      ! 3.6  update subsoil temperature
      if ( lveg(i) == 9 .or. lveg(i) == 12 ) then
        c31 = d_half*dtimea*(d_one+deprat(i))
        c41 = dtimea*deprat(i)
        tgbrd(i) = ((d_one-c31+fct2(i)) * tgbrd(i)+c41*tgrd(i) + &
                   dtimea*t3)/(d_one+c31+fct2(i))
      else
        c3t = d_half*dtimea*deprat(i)
        c4t = dtimea*deprat(i)
        tgbrd(i) = ((d_one-c3t+fct2(i))*tgbrd(i)+c4t*tgrd(i)) / &
                (d_one+c3t+fct2(i))
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  contains

    pure real(rkx) function fsk(x)
      implicit none
      real(rkx) , intent(in) :: x
      fsk = (2.9e-7_rkx*x+4.0e-9_rkx) / &
           (((d_one-0.6_rkx*x)*x+0.09_rkx)*(0.23_rkx+x))
    end function fsk
    pure real(rkx) function fsc(x)
      implicit none
      real(rkx) , intent(in) :: x
      fsc = (0.23_rkx+x)*4.186e6_rkx
    end function fsc
    pure real(rkx) function fct1(x)
      implicit none
      real(rkx) , intent(in) :: x
      fct1 = wlhf*d_rfour*1.414_rkx/x
    end function fct1

  end subroutine tgrund

end module mod_bats_bndry
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
