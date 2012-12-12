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
!
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_service
  use mod_bats_common
  use mod_bats_leaftemp
  use mod_bats_drag
  use mod_bats_internal
!
  private
!
  public :: bndry
!
  real(rk8) , parameter :: minsigf = 0.001D+00
  real(rk8) , parameter :: lowsice = 1.0D-22
  real(rk8) , parameter :: rainsnowtemp = 2.2D0
  real(rk8) , parameter :: xnu = twopi/secpd
!
  contains
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
!                           tseaice                           root
!                            tgrund                          satur
!                              snow                         lfdrag
!                             water                         condch
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
!
    real(rk8) :: fact , qsatd , rai
    integer(ik4) :: n , i , j
    real(rk8) , parameter :: minwrat = 1.0D-04
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'bndry'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
   
    !=======================================================================
    !   1.   initialize
    !=======================================================================
   
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          htvp(n,j,i) = wlhv
          if ( ( tgrd(n,j,i) < tzero .and. ldmsk1(n,j,i) /= 0 ) .or. &
                sncv(n,j,i) > d_zero ) then
            htvp(n,j,i) = wlhs
          end if
          sdrop(n,j,i) = d_zero
          etrrun(n,j,i) = d_zero
          flnet(n,j,i) = d_zero
          fevpg(n,j,i) = d_zero
          fseng(n,j,i) = d_zero
          vegt(n,j,i) = d_zero
          efpr(n,j,i) = d_zero
          etr(n,j,i) = d_zero
          ! switch between rain and snow /tm is ref. temp set= anemom temp - 2.2
          tm(n,j,i) = sts(n,j,i) - rainsnowtemp
          ! soil moisture ratio (to max) as used in subrouts tgrund,
          ! water, and root (called by lftemp): watu=upper, watr=root,
          ! watt=total
          if ( ldmsk1(n,j,i) /= 0 ) then
            watu(n,j,i) = ssw(n,j,i)/gwmx0(n,j,i)
            watr(n,j,i) = rsw(n,j,i)/gwmx1(n,j,i)
            watt(n,j,i) = tsw(n,j,i)/gwmx2(n,j,i)
            watu(n,j,i) = dmin1(watu(n,j,i),d_one)
            watr(n,j,i) = dmin1(watr(n,j,i),d_one)
            watt(n,j,i) = dmin1(watt(n,j,i),d_one)
            watr(n,j,i) = dmax1(watr(n,j,i),minwrat)
            watu(n,j,i) = dmax1(watu(n,j,i),minwrat)
          end if
        end do
      end do
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
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) <= minsigf .and. ldmsk1(n,j,i) /= 2 ) then
              qsatd = qgrd(n,j,i)*gwet(n,j,i) * &
                      lfta(n,j,i)*(tzero-lftb(n,j,i)) * &
                      (d_one/(tgrd(n,j,i)-lftb(n,j,i)))**2
              rai = cdrx(n,j,i)*vspda(n,j,i)*rhs(n,j,i)
              cgrnds(n,j,i) = rai*cpd
              cgrndl(n,j,i) = rai*qsatd
              cgrnd(n,j,i) = cgrnds(n,j,i) + cgrndl(n,j,i)*htvp(n,j,i)
              ! 3.2  sensible and latent fluxes using soil temperatures
              ! from previous time step
              delq(n,j,i) = (qs(n,j,i)-qgrd(n,j,i))*gwet(n,j,i)
              delt(n,j,i) = sts(n,j,i) - tgrd(n,j,i)
              evpr(n,j,i) = -rai*delq(n,j,i)
              sent(n,j,i) = -cgrnds(n,j,i)*delt(n,j,i)
              ! 3.3  fluxes to subrout tgrund (evap is in kg/m**2/s)
              fseng(n,j,i) = sent(n,j,i)
              fevpg(n,j,i) = evpr(n,j,i)
              ! 3.4  equate canopy to air, for temp, wind over bare grnd;
              ! needed as factors of sigf(=0) in subr water (uaf) and
              ! subr drag (tlef(n,j,i) carried over to next tstep).
              tlef(n,j,i) = sts(n,j,i)
              uaf(n,j,i) = vspda(n,j,i)
            end if
          end if
        end do
      end do
    end do
   
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > minsigf ) then   !  check each point
              !================================================================
              ! 4.   vegetation
              !================================================================
              ! 4.1  add precipitation to leaf water
              ldew(n,j,i) = ldew(n,j,i)+dtbat*sigf(n,j,i)*prcp(n,j,i)
              ldew(n,j,i) = dmax1(ldew(n,j,i),d_zero)
            end if
          end if
        end do
      end do
    end do
!
!   4.2  distribute excess leaf water to soil
    call vcover
    call drip
!
!   4.3  calculate canopy temperature, soil and total fluxes,
!   and leaf water change by evapotranspiration
    call lftemp
!
!   4.4  calculate carbon sources and sinks
!   call co2
!=======================================================================
!   5.   back to any surface but ocean
!=======================================================================
!
!   5.1  over sea ice
    call tseaice
!
!   5.2  over land, calculate soil temp and surface hydrology
    call tgrund
    call snow
    call water
    !
    ! 5.3  over ocean
!   set snow cover to zero in case it was previously sea ice
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) == 0 ) then
            sncv(n,j,i) = d_zero
            snag(n,j,i) = d_zero
          end if
        end do
      end do
    end do

    !=======================================================================
    !   6.   linkage to meteorological model
    !=======================================================================
 
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) == 0 .or. &
               lveg(n,j,i) == 14 .or. lveg(n,j,i) == 15 ) then
            gwet(n,j,i) = d_one
          end if
          ! 6.1  rate of momentum transfer per velocity
          drag(n,j,i) = cdrx(n,j,i)*vspda(n,j,i)*rhs(n,j,i)
          ! 6.3  latent and heat fluxes over ocean, plus a dummy taf
          if ( ldmsk1(n,j,i) == 0 ) then
            tlef(n,j,i) = sts(n,j,i)
            fact = -drag(n,j,i)
            delq(n,j,i) = (qs(n,j,i) - qgrd(n,j,i))*gwet(n,j,i)
            delt(n,j,i) = sts(n,j,i) - tgrd(n,j,i)
            ! evaporation is in kg/m**2/s
            evpr(n,j,i) = fact*delq(n,j,i)
            sent(n,j,i) = fact*cpd*delt(n,j,i)
          end if
          if ( sigf(n,j,i) < minsigf ) taf(n,j,i) = tgrd(n,j,i)
          ! 6.2  parameters for temperature difference at anemometer level
          ! zdelt(i) = zdelt(i)*delt(n,j,i)
          ! 6.4  evaporative flux, accounting for sublimation
          ! evprr(i) = wlhv*(evpr(n,j,i)-fevpg) + htvp(n,j,i)*fevpg
          ! 6.5  nondimensional equivalent bucket capacity for comparisons
          !      with bucket models; usually 1 or less, except where
          !      saturated (then around 2)
          ! rh2ox(i) = (watr-wiltr)/(relfc-wiltr)
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx) 
#endif
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
    integer(ik4) :: n , i , j
#ifdef DEBUG
    character(len=64) :: subroutine_name = 'vcover'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > minsigf ) then
              seasb(n,j,i) = fseas(tgbrd(n,j,i),lveg(n,j,i))
            end if
          end if
        end do
      end do
    end do
   
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > minsigf ) then
              xlai(n,j,i) = xla(lveg(n,j,i))
              xlai(n,j,i) = xlai(n,j,i) + &
                         (xlai0(lveg(n,j,i))-xlai(n,j,i))*(d_one-seasb(n,j,i))
              rlai(n,j,i) = xlai(n,j,i) + sai(lveg(n,j,i))
              xlsai(n,j,i) = xlai(n,j,i) + sai(lveg(n,j,i))
              vegt(n,j,i) = sigf(n,j,i)*xlsai(n,j,i)
            end if
          end if
        end do
      end do
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
    integer(ik4) :: n , i , j
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'drip'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > minsigf ) then
              ! xrun = leaf drip ; sdrop = snow drop off foliage
              etrrun(n,j,i) = d_zero
              xrun(n,j,i) = ldew(n,j,i) - dewmax*vegt(n,j,i)
              sdrop(n,j,i) = d_zero
              ! test on maximum value of dew
              if ( xrun(n,j,i) > d_zero ) then
                etrrun(n,j,i) = xrun(n,j,i)
                ldew(n,j,i) = dewmax*vegt(n,j,i)
              end if
              ! below freezing excess leaf water falls as snow
              if ( (xrun(n,j,i) > d_zero) .and. (tm(n,j,i) < tzero) ) then
                etrrun(n,j,i) = d_zero
                sdrop(n,j,i) = xrun(n,j,i)
              end if
            end if
          end if
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx) 
#endif
  end subroutine drip
! 
!=======================================================================
! TSEAICE
!     Routine provides sensible and latent fluxes
!             and snow melt over ice
!
!     fss = conductive heat flow through ice
!     hrl = latent heat through leads
!     hsl = sensible heat through leads
!     hs  = heat energy balance at surface of ice
!
!=======================================================================
!
  subroutine tseaice
    implicit none
!
    real(rk8) :: bb , fact , fss , hrl , hs , hsl , qgrnd , ratsi ,     &
                rhosw3 , rsd1 , rss , smc4 , smt , tg , tgrnd , wss , wtt
    integer(ik4) :: n , i , j
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'tseaice'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          ! lake model handles this case
          if ( llake .and. iveg1(n,j,i) == 14 ) exit
   
          if ( ldmsk1(n,j,i) == 2 ) then
            ! rhosw = density of snow relative to water
            rhosw3 = rhosw(n,j,i)**3
            ! cice = specific heat of sea-ice per unit volume
            rsd1 = cice*(sfice(n,j,i)*d_r1000)
            if ( sncv(n,j,i) > d_zero ) then
              rss = csnw*(sncv(n,j,i)*d_r1000)
              ratsi = sncv(n,j,i)/(1.4D0*rhosw3*sfice(n,j,i))
              wtt = d_one/(d_one+ratsi)
              wss = (sncv(n,j,i)+2.8D0*rhosw3*sfice(n,j,i)) / &
                    (sncv(n,j,i)+1.4D0*rhosw3*sfice(n,j,i))
              ! include snow heat capacity
              rsd1 = d_half*(wss*rss+wtt*rsd1)
            end if
            tgbrd(n,j,i) = -d_two + tzero
            ! subsurface heat flux through ice
            fss = 7.0D-4*(tgbrd(n,j,i)-tgrd(n,j,i)) * &
                  ch2o*rhosw3/(sncv(n,j,i)+1.4D0*rhosw3*sfice(n,j,i))
            sfice(n,j,i) = sfice(n,j,i) + fss*dtbat/wlhf*1.087D0
            ! set sea ice parameter for melting and return
            if ( sfice(n,j,i) <= d_zero ) then
              sfice(n,j,i) = d_zero
              ldmsk1(n,j,i) = 0
              lveg(n,j,i) = 15
              exit
            end if
            ! assume lead ocean temp is -1.8c
            ! flux of heat and moisture through leads
            ! sat. mixing ratio at t=-1.8c is 3.3e-3
            qice(n,j,i) = 3.3D-3 * stdp/sfcp(n,j,i)
            !
            ! determine effective surface fluxes over ice, allowing for leads;
            ! aarea is set in mod_bats_runparams
            !
            tlef(n,j,i) = sts(n,j,i)
            qgrnd = ((d_one-aarea)*cdr(n,j,i)*qgrd(n,j,i) + &
                    aarea*clead(n,j,i)*qice(n,j,i))/cdrx(n,j,i)
            tgrnd = ((d_one-aarea)*cdr(n,j,i)*tgrd(n,j,i) + &
                    aarea*clead(n,j,i)*(tzero-1.8D0))/cdrx(n,j,i)
            fact = -rhs(n,j,i)*cdrx(n,j,i)*vspda(n,j,i)
            delq(n,j,i) = (qs(n,j,i)-qgrnd)*gwet(n,j,i)
            delt(n,j,i) = sts(n,j,i) - tgrnd
            ! output fluxes, averaged over leads and ice
            evpr(n,j,i) = fact*delq(n,j,i)
            sent(n,j,i) = fact*cpd*delt(n,j,i)
            hrl = rhs(n,j,i)*vspda(n,j,i)*clead(n,j,i) * &
                      (qice(n,j,i)-qs(n,j,i))
            hsl = rhs(n,j,i)*vspda(n,j,i)*clead(n,j,i) * &
                      (tzero-1.8D0-sts(n,j,i))*cpd
            ! get fluxes over ice for sublimation (subrout snow)
            ! and melt (below) calculation
            fseng(n,j,i) = (sent(n,j,i)-aarea*hsl)/(d_one-aarea)
            fevpg(n,j,i) = (evpr(n,j,i)-aarea*hrl)/(d_one-aarea)
            hs = fsw(j,i) - flw(j,i) - fseng(n,j,i) - wlhs*fevpg(n,j,i)
            bb = dtbat*(hs+fss)/rsd1
            ! snow melt
            sm(n,j,i) = d_zero
            if ( tgrd(n,j,i) >= tzero ) sm(n,j,i) = (hs+fss)/wlhf
            if ( sm(n,j,i) <= d_zero ) sm(n,j,i) = d_zero
            smc4 = sm(n,j,i)*dtbat
            if ( sncv(n,j,i) < smc4 ) then
            ! all snow removed, melt ice
              smt = (sncv(n,j,i)/dtbat)
              ! rho(h2o)/rho(ice) = 1.087
              sfice(n,j,i) = sfice(n,j,i) + dtbat*(smt-sm(n,j,i))*1.087D0
              sm(n,j,i) = smt
              tgrd(n,j,i) = tzero
              ! set sea ice parameter for melting and return
              if ( sfice(n,j,i) <= d_zero ) then
                sfice(n,j,i) = d_zero
                ldmsk1(n,j,i) = 0
                lveg(n,j,i) = 15
                exit
              end if
            else
              ! snow or ice with no snow melting
              tg = tgrd(n,j,i) + bb
              if ( tg >= tzero ) tgrd(n,j,i) = tzero
              if ( tg < tzero ) tgrd(n,j,i) = tg
            end if
          end if
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx) 
#endif
  end subroutine tseaice
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
!       rsur = surface runoff
!     trnof(n,j,i) = total runoff mm
!     srnof(n,j,i) = surface runoff mm
!
!     xkmxr and wflux1 determine flow thru upper/root soil interface
!     evmxt, xkmx1, and xkmx2 determine flow thru lower interfaces
!
!     veg type 10 "irrigated crop" is irrigated through reducing
!          the runoff (rsur), i.e., by adding a negative number
!          if the land isn't at least 60% saturated.
!     veg type 13 and 14 are water covered (lake, swamp, rice paddy);
!          negative runoff keeps this land saturated.
!=======================================================================
!
  subroutine water
    implicit none
!
    real(rk8) :: b , bfac , bfac2 , delwat , est0 , evmax , evmxr , &
               evmxt , rap , vakb , wtg2c , xxkb , gwatr , rsubsr , &
               rsubss , xkmx1 , xkmx2 , xkmxr
    integer(ik4) :: n , i , j
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
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) == 1 ) then
            !
            ! 1.1  reduce infiltration for frozen ground
            !
            if ( tgbrd(n,j,i) > tzero ) then
              xkmxr = xkmx(n,j,i)
            else
              xkmxr = d_zero
            end if
            !
            ! 1.11 permafrost or ice sheet
            !
            if ( lveg(n,j,i) == 9 .or. lveg(n,j,i) == 12 ) then
              xkmx1 = d_zero
              xkmx2 = d_zero
            else
              xkmx1 = xkmx(n,j,i)
              xkmx2 = drain
            end if
            !
            ! 1.2  diffusive fluxes
            !
            evmxr = evmx0(n,j,i)*xkmxr/xkmx(n,j,i)
            evmxt = evmx0(n,j,i)*xkmx1/xkmx(n,j,i)
            b = bsw(n,j,i)
            bfac = watr(n,j,i)**(d_three+bfc(n,j,i)) * &
                   watu(n,j,i)**(b-bfc(n,j,i)-d_one)
            bfac2 = watt(n,j,i)**(d_two+bfc(n,j,i)) * &
                    watr(n,j,i)**(b-bfc(n,j,i))
            wfluxc(n,j,i) = evmxr*(depuv(lveg(n,j,i)) / &
                                   deprv(lveg(n,j,i)))**0.4D0*bfac
            wflux1(n,j,i) = wfluxc(n,j,i)*watr(n,j,i)
            wflux2(n,j,i) = evmxt*sqrt(depuv(lveg(n,j,i)) / &
                                       deprv(lveg(n,j,i)))* &
                            bfac2*(watt(n,j,i)-watr(n,j,i))
            !
            ! 1.3  gravitational drainage
            !
            rsubss = xkmxr*watr(n,j,i)**(b+d_half) * &
                                         watu(n,j,i)**(b+2.5D0)
            rsubsr = xkmx1*watt(n,j,i)**(b+d_half) * &
                                         watr(n,j,i)**(b+2.5D0)
            rsubst(n,j,i) = xkmx2*watt(n,j,i)**(d_two*b+d_three)
            !
            ! 1.32 bog or water
            !
            if ( (lveg(n,j,i) >= 13) .and. (lveg(n,j,i) <= 15) ) then
              rsubst(n,j,i) = d_zero
              rsubss = d_zero
              rsubsr = d_zero
            end if
            !
            ! 1.4  fluxes through internal surfaces
            !
            wflux1(n,j,i) = wflux1(n,j,i) - rsubss
            wflux2(n,j,i) = wflux2(n,j,i) - rsubsr
          end if
        end do
      end do
    end do
    !
    ! 1.5  net flux at air interface
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) == 1 ) then
            gwatr = pw(n,j,i) - evapw(n,j,i) + sm(n,j,i) + etrrun(n,j,i)/dtbat
            !
            !=================================================================
            ! 2.   define runoff terms
            !=================================================================
            !
            ! 2.1  surface runoff
            !
            wata(n,j,i) = d_half*(watu(n,j,i)+watr(n,j,i))
            !
            ! 2.11 increase surface runoff over frozen ground
            !
            if ( tgrd(n,j,i) < tzero ) then
              rsur(n,j,i) = dmin1(d_one,wata(n,j,i)**1) * &
                            dmax1(d_zero,gwatr)
            else
              rsur(n,j,i) = dmin1(d_one,wata(n,j,i)**4) * &
                            dmax1(d_zero,gwatr)
            end if
            !
            ! 2.12 irrigate cropland
            !
            if ( lveg(n,j,i) == 10 .and. watr(n,j,i) < relfc(n,j,i) ) then
              rsur(n,j,i) = rsur(n,j,i) + &
                     (rsw(n,j,i)-relfc(n,j,i)*gwmx1(n,j,i))/dtbat
            end if
            ! Imported Lara Kuepper's Irrigated Crop modification from RegCM3
            ! see Kueppers et al. (2008)
            ! if ( lveg(n,j,i) == 10 .and. watr(n,j,i) < relaw(n,j,i) ) then
            !   rsur(n,j,i) = rsur(n,j,i) + &
            !         (rsw(n,j,i)-relaw(n,j,i)*gwmx1(n,j,i))/dtbat
            ! end if
            !
            ! 2.13 saturate swamp or rice paddy
            !
            if ( lveg(n,j,i) == 13 ) then
              ! Graziano. Seems that this runoff is borken at least
              ! at time step zero. Try to mediate using relfc. Is this
              ! correct ? Mhhhh....
              rsur(n,j,i) = rsur(n,j,i) + dmin1(d_zero,(rsw(n,j,i)- &
                                    relfc(n,j,i)*gwmx1(n,j,i))/dtbat)
            end if
            !
            ! 2.2  total runoff
            !
            rnof(n,j,i) = rsur(n,j,i) + rsubst(n,j,i)
            !
            !=================================================================
            !         3.   increment soil moisture
            !=================================================================
            !
            ! 3.1  update top layer with implicit treatment of flux from below
            !
            ssw(n,j,i) = ssw(n,j,i) + dtbat*(gwatr-efpr(n,j,i) * &
                         etr(n,j,i)-rsur(n,j,i)+wflux1(n,j,i))
            ssw(n,j,i) = ssw(n,j,i)/(d_one+wfluxc(n,j,i)*dtbat/gwmx0(n,j,i))
            !
            ! 3.2  update root zone
            !
            rsw(n,j,i) = rsw(n,j,i) + dtbat*(gwatr-etr(n,j,i) - &
                         rsur(n,j,i) + wflux2(n,j,i))
            !
            ! 3.3  update total water
            !
            tsw(n,j,i) = tsw(n,j,i) + &
                           dtbat*(gwatr-etr(n,j,i)-rnof(n,j,i))
          end if
        end do
      end do
    end do
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) == 1 ) then
            !
            !=================================================================
            !  4.   check whether soil water exceeds maximum capacity or
            !       becomes negative (should rarely or never happen)
            !=================================================================
            !
            ! 4.1  surface water assumed to move downward into soil
            !
            if ( ssw(n,j,i) > gwmx0(n,j,i) ) ssw(n,j,i) = gwmx0(n,j,i)
            !
            ! 4.2  excess root layer water assumed to move downward
            !
            if ( rsw(n,j,i) > gwmx1(n,j,i) ) rsw(n,j,i) = gwmx1(n,j,i)
            !
            ! 4.3  excess total water assumed to go to subsurface runoff
            !
            if ( tsw(n,j,i) > gwmx2(n,j,i) ) then
              delwat = tsw(n,j,i) - gwmx2(n,j,i)
              tsw(n,j,i) = gwmx2(n,j,i)
              rsubst(n,j,i) = rsubst(n,j,i) + delwat/dtbat
            end if
            !
            ! 4.4  check for negative water in top layer
            !
            if ( ssw(n,j,i) <= 1.0D-2 ) ssw(n,j,i) = 1.0D-2
            !
            !=================================================================
            !         5.   accumulate leaf interception
            !=================================================================
            !
            ! Graziano : This is NEVER used elsewhere in the code...
            ! ircp = ircp + sigf(n,j,i) * &
            !              (dtbat*prcp(n,j,i)) - (sdrop(n,j,i)+etrrun(n,j,i))
            ! Graziano : This is NEVER used elsewhere in the code...
            !
            !=================================================================
            !         6.   evaluate runoff (incremented in ccm)
            !=================================================================
            !
            ! update total runoff
            !
            rnof(n,j,i)   = rsur(n,j,i) + rsubst(n,j,i)
            ! Graziano: Do not go to mm*day here and then back
            !           to mm in mtrxbats
            trnof(n,j,i)  = rnof(n,j,i)
            srnof(n,j,i) = rsur(n,j,i)
          else                       ! ocean or sea ice
            rnof(n,j,i)   = d_zero
            trnof(n,j,i)  = d_zero
            srnof(n,j,i) = d_zero
          end if
        end do
      end do
    end do
    !
    !=======================================================================
    !   7.   calculate potential evaporation and use mod_to determine
    !        wetness factor, allowing for snow being saturated
    !=======================================================================
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) == 1 ) then
            xxkb = dmin1(rough(lveg(n,j,i)),d_one)
            vakb = (d_one-sigf(n,j,i))*vspda(n,j,i) + sigf(n,j,i) * &
                   (xxkb*uaf(n,j,i)+(d_one-xxkb)*vspda(n,j,i))
            wtg2c = (d_one-sigf(n,j,i))*cdrx(n,j,i)*vakb
            rap = rhs(n,j,i)*(csoilc*uaf(n,j,i)*sigf(n,j,i)*(qgrd(n,j,i) +  &
                  delq(n,j,i)-qs(n,j,i))+wtg2c*(qgrd(n,j,i)-qs(n,j,i)))
            bfac = watr(n,j,i)**(d_three+bfc(n,j,i)) * &
                   watu(n,j,i)**(bsw(n,j,i)-bfc(n,j,i)-d_one)
            est0 = evmx0(n,j,i)*bfac*watu(n,j,i)
            evmax = dmax1(est0,d_zero)
            gwet(n,j,i) = dmin1(d_one,evmax/dmax1(1.0D-14,rap))
            gwet(n,j,i) = scvk(n,j,i) + gwet(n,j,i)*(d_one-scvk(n,j,i))
          end if
        end do
      end do
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
!
    real(rk8) :: age1 , age2 , age3 , arg , arg2 , dela , dela0 , &
                 dels , tage , sold
    integer(ik4) :: n , i , j
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'snow'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
    age3 = 0.3D0
   
    !=======================================================================
    !   1.   partition soil evaporation and precipitation
    !        between water and snow
    !=======================================================================
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            evapw(n,j,i) = fevpg(n,j,i)
            evaps(n,j,i) = scvk(n,j,i)*evapw(n,j,i)
            if ( ldmsk1(n,j,i) == 2 ) then
              evaps(n,j,i) = fevpg(n,j,i)
            end if
            evapw(n,j,i) = (d_one-scvk(n,j,i))*evapw(n,j,i)
            ! tm  is temperature of precipitation
            if ( tm(n,j,i) >= tzero ) then
              pw(n,j,i) = prcp(n,j,i)*(d_one-sigf(n,j,i))
              ps(n,j,i) = d_zero
            else
              ! snowing
              pw(n,j,i) = d_zero
              ps(n,j,i) = prcp(n,j,i)*(d_one-sigf(n,j,i))
            end if
          end if
        end do
      end do
    end do
    !
    !=======================================================================
    !   2.   update snow cover
    !=======================================================================
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            sold = sncv(n,j,i)
            sncv(n,j,i) = sncv(n,j,i) + &
                         dtbat*(ps(n,j,i)-evaps(n,j,i)-sm(n,j,i)) + sdrop(n,j,i)
            if ( sncv(n,j,i) < dlowval ) then
              sncv(n,j,i) = d_zero
              snag(n,j,i) = d_zero
            end if
            ! snow cover except for antarctica
            !==================================================================
            !  3.   increment non-dimensional "age" of snow;
            !       10 mm snow restores surface to that of new snow.
            !==================================================================
            if ( sncv(n,j,i) > d_zero ) then
              arg = 5.0D3*(d_one/tzero-d_one/tgrd(n,j,i))
              age1 = dexp(arg)
              arg2 = dmin1(d_zero,d_10*arg)
              age2 = dexp(arg2)
              tage = age1 + age2 + age3
              dela0 = 1.0D-6*dtbat
              dela = dela0*tage
              dels = d_r10*dmax1(d_zero,sncv(n,j,i)-sold)
              snag(n,j,i) = (snag(n,j,i)+dela)*(d_one-dels)
              if ( snag(n,j,i) < dlowval ) snag(n,j,i) = d_zero
            end if
            ! antarctica
            if ( sncv(n,j,i) > 800.0D0 ) snag(n,j,i) = d_zero
          end if
        end do
      end do
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
    real(rk8) :: bcoefd , bcoefs , c31 , c3t , c41 , c4t , cder , depr , &
             depu , xdt2 , xdtime , dtimea , froze2 , frozen , rscss ,   &
             tbef , tg , tinc , wtas , wtax , wtd , wtds , hs , depann , &
             depdiu , rscsa , rscsd , ska , skd , sks
    real(rk8) :: dtbat2 , rdtbat2 , xlexp , xnua , swtrta , swtrtd
    integer(ik4) :: n , i , j
    real(rk8) , parameter :: xkperi = 1.4D-6
    real(rk8) , parameter :: t3 = 271.0D0 ! permafrost temperature
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'tgrund'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    dtbat2  = dtbat*d_two
    rdtbat2 = d_one/dtbat2

!=======================================================================
!   1.   define thermal conductivity, heat capacity,
!        and other force restore parameters
!=======================================================================

    xnua = xnu/dayspy
    xdtime = dtbat*xnu
    dtimea = dtbat*xnua
    xdt2 = d_half*xdtime
   
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) == 1 ) then
            ! 1.1  frozen ground values using 44 m2/yr for frozen ground
            !      thermal diffusion coefficient, based on the values of
            !      50 and 38 quoted by osterkamp; ice contribution to
            !      specific heat only o.49 that of water
            swtrtd = watu(n,j,i)*porsl(n,j,i)
            if ( tgrd(n,j,i) < tzero ) then
              frozen = 0.85D0*dmin1(d_one,d_rfour*(tzero-tgrd(n,j,i)))
              skd = xkperi
              rscsd = fsc(swtrtd*(d_one-0.51D0*frozen))
            else
              skd = fsk(swtrtd)*texrat(n,j,i)
              rscsd = fsc(swtrtd)
            end if
            swtrta = watr(n,j,i)*porsl(n,j,i)
            if ( tgbrd(n,j,i) < tzero ) then
              froze2 = 0.85D0*dmin1(d_one,d_rfour*(tzero-tgbrd(n,j,i)))
              ska = xkperi
              rscsa = fsc(swtrta*(d_one-0.51D0*froze2))
            else
              ska = fsk(swtrta)*texrat(n,j,i)
              rscsa = fsc(swtrta)
            end if
            ! 1.2  correct for snow cover, if significant
            depdiu = dsqrt(d_two*skd/xnu)
            bcoef(n,j,i) = xdtime*depdiu/(rscsd*skd)
            if ( scrat(n,j,i) > 0.001D0 ) then
              xlexp = -d_two*scrat(n,j,i)/depdiu
              ! Graziano : Limit exponential argument
              if ( xlexp > -25.0D0 ) then
                wtd = dexp(xlexp)
              else
                wtd = d_zero
              end if
              rscss = csnw*rhosw(n,j,i)
              sks = 7.0D-7*cws*rhosw(n,j,i)
              bcoefs = dsqrt(d_two*sks/xnu)/(rscss*sks)
              wtds = (d_one-wtd)*scvk(n,j,i)
              bcoefd = dsqrt(d_two*skd/xnu)/(rscsd*skd)
              bcoef(n,j,i) = xdtime*(wtds*bcoefs+(d_one-wtds)*bcoefd)
              depdiu = wtds*dsqrt(d_two*sks/xnu) + &
                            (d_one-wtds)*depdiu
            end if
            depann = dsqrt(d_two*ska/xnua)
            if ( scrat(n,j,i) > 0.02D0 ) then
              xlexp = -d_two*scrat(n,j,i)/depann
              if ( xlexp > -25.0D0 ) then
                wtax = dexp(xlexp)
              else
                wtax = d_zero
              end if
              wtas = (d_one-wtax)*scvk(n,j,i)
              depann = wtas*dsqrt(d_two*sks/xnua) + &
                            (d_one-wtas)*depann
            end if
            deprat(n,j,i) = depann/depdiu
            !=================================================================
            !         2.   collect force restore terms
            !=================================================================
            cc(n,j,i) = d_one
            fct2(n,j,i) = d_zero
            !
            ! 2.1  add freezing thermal inertia
            if ( (tgrd(n,j,i) < tzero) .and.          &
                 (tgrd(n,j,i) > (tzero-d_four)) .and. &
                 (sfice(n,j,i) < lowsice) ) then
              depu = depuv(lveg(n,j,i))*d_r1000
              cc(n,j,i) = d_one + dmax1(ssw(n,j,i) - &
                        frezu(lveg(n,j,i)),d_zero)*fct1(depu*rscsd)
            end if
            if ( (tgbrd(n,j,i) < tzero) .and.                  &
                 (tgbrd(n,j,i) > (tzero-d_four)) .and.         &
                 (sfice(n,j,i) < lowsice) ) then
              depr = deprv(lveg(n,j,i))*d_r1000
              fct2(n,j,i) = dmax1(rsw(n,j,i)-freza(lveg(n,j,i)),d_zero) * &
                                fct1(depr*rscsa)
            end if
            ! 2.2  large thermal inertial for permanent ice cap
            if ( lveg(n,j,i) == 12 ) fct2(n,j,i) = d_1000*fct2(n,j,i)
            ! 2.3  collect energy flux terms
            rnet(n,j,i) = fsw(j,i) - sigf(n,j,i)*(sabveg(j,i)-flnet(n,j,i)) - &
                     (d_one-sigf(n,j,i))*(flw(j,i)-sigf(n,j,i)*flneto(n,j,i))
            hs = rnet(n,j,i) - fseng(n,j,i) - fevpg(n,j,i)*htvp(n,j,i)
            bb(n,j,i) = bcoef(n,j,i)*hs + xdtime*tgbrd(n,j,i)
            ! 2.4  add in snowmelt (melt enough snow to reach freezing temp)
            sm(n,j,i) = d_zero
            if ( sncv(n,j,i) > d_zero ) then
              cder = bcoef(n,j,i)*cgrnd(n,j,i)
              sm(n,j,i) = (bb(n,j,i) + &
                        (cc(n,j,i)-xdt2+cder)*tgrd(n,j,i) - tzero * &
                        (cc(n,j,i)+xdt2+cder))/(bcoef(n,j,i)*wlhf)
              ! snow melt always between 0 and total snow
              sm(n,j,i) = dmax1(d_zero, &
                                dmin1(sm(n,j,i),sncv(n,j,i)*d_two*rdtbat2))
              bb(n,j,i) = bb(n,j,i) - bcoef(n,j,i)*wlhf*sm(n,j,i)
            end if
          end if
        end do
      end do
    end do
   
    !=======================================================================
    !   3.   update soil temperatures
    !=======================================================================
    !   3.1  update surface soil temperature
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) == 1 ) then
            tbef = tgrd(n,j,i)
            cder = bcoef(n,j,i)*cgrnd(n,j,i)
            tg = (bb(n,j,i)+(cc(n,j,i)-xdt2+cder)*tgrd(n,j,i)) / &
                 (cc(n,j,i)+xdt2+cder)
            tgrd(n,j,i) = tg
            ! 3.2  put brakes on large temperature excursions
            tgrd(n,j,i) = dmin1(tbef+d_10,tgrd(n,j,i))
            tgrd(n,j,i) = dmax1(tbef-d_10,tgrd(n,j,i))
            ! 3.3  correct fluxes to present soil temperature
            tinc = tgrd(n,j,i) - tbef
            sent(n,j,i) = sent(n,j,i) + tinc*cgrnds(n,j,i)
            evpr(n,j,i) = evpr(n,j,i) + tinc*cgrndl(n,j,i)
            fseng(n,j,i) = fseng(n,j,i) + tinc*cgrnds(n,j,i)
            fevpg(n,j,i) = fevpg(n,j,i) + tinc*cgrndl(n,j,i)
            !
            ! 3.5  couple to deep temperature in permafrost
            ! 3.6  update subsoil temperature
            if ( lveg(n,j,i) == 9 .or. lveg(n,j,i) == 12 ) then
              c31 = d_half*dtimea*(d_one+deprat(n,j,i))
              c41 = dtimea*deprat(n,j,i)
              tgbrd(n,j,i) = ((d_one-c31+fct2(n,j,i)) * &
                     tgbrd(n,j,i)+c41*tgrd(n,j,i) + &
                     dtimea*t3)/(d_one+c31+fct2(n,j,i))
            else
              c3t = d_half*dtimea*deprat(n,j,i)
              c4t = dtimea*deprat(n,j,i)
              tgbrd(n,j,i) = ((d_one-c3t+fct2(n,j,i))* &
                    tgbrd(n,j,i)+c4t*tgrd(n,j,i)) / (d_one+c3t+fct2(n,j,i))
            end if
          end if
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx) 
#endif
  contains

    function fsk(x)
      implicit none
      real(rk8) :: fsk
      real(rk8) , intent(in) :: x
      fsk = (2.9D-7*x+4.0D-9)/(((d_one-0.6D0*x)*x+0.09D0)*(0.23D0+x))
    end function fsk
    function fsc(x)
      implicit none
      real(rk8) :: fsc
      real(rk8) , intent(in) :: x
      fsc = (0.23D0+x)*4.186D6
    end function fsc
    function fct1(x)
      implicit none
      real(rk8) :: fct1
      real(rk8) , intent(in) :: x
      fct1 = wlhf*d_rfour*1.414D0/x
    end function fct1
! 
  end subroutine tgrund
!
end module mod_bats_bndry
