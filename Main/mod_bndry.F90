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
 
      module mod_bndry
!
      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_bats
      use mod_leaftemp
      use mod_drag
      use mod_service
!
      private
!
      public :: bndry
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
!
      implicit none
!
      real(8) :: fact , qsatd , rai
      integer :: n , i
      character (len=50) :: subroutine_name='bndry'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
 
!=======================================================================
!l    1.   initialize
!=======================================================================
 
      do i = 2 , iym1
        do n = 1 , nnsg
 
          htvp(n,i) = wlhv
          if ( (tg1d(n,i).lt.tzero .and. ldoc1d(n,i).gt.0.5) .or.       &
             & scv1d(n,i).gt.0. ) htvp(n,i) = wlhs
          sdrop(n,i) = 0.
          etrrun(n,i) = 0.
          flnet(n,i) = 0.
          fevpg(n,i) = 0.
          fseng(n,i) = 0.
          vegt(n,i) = 0.
          efpr(n,i) = 0.
          etr(n,i) = 0.
!         **********            switch between rain and snow /tm is
!         ref. temp set= anemom temp -2.2
          tm(n,i) = ts1d(n,i) - 2.2
 
!*        soil moisture ratio (to max) as used in subrouts tgrund,
!*        water, and root (called by lftemp): watu=upper, watr=root,
 
!         watt=total
          if ( lveg(n,i).ge.1 ) then
            watu(n,i) = ssw1d(n,i)/gwmx0(n,i)
            watr(n,i) = rsw1d(n,i)/gwmx1(n,i)
            watt(n,i) = tsw1d(n,i)/gwmx2(n,i)
            watu(n,i) = dmin1(watu(n,i),1.D0)
            watr(n,i) = dmin1(watr(n,i),1.D0)
            watt(n,i) = dmin1(watt(n,i),1.D0)
            watr(n,i) = dmax1(watr(n,i),1.D-4)
            watu(n,i) = dmax1(watu(n,i),1.D-4)
          end if
 
        end do
      end do
 
!=======================================================================
!l    2.  calculate transfer coefficients at layer 1
!=======================================================================
!     2.1  sigf corrected in drag: remove snow-covered veg
      call drag
 
!     2.2  get saturation vapor pressure of soil surface
      call satur(qg1d,tg1d,p1d)
 
!=======================================================================
!l    3.   bare land
!=======================================================================
!l    3.1  get derivative of fluxes with repect to tg
 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).le.0.001 .and. ldoc1d(n,i).lt.1.5 ) then
              qsatd = qg1d(n,i)*gwet1d(n,i)*lfta(n,i)*(tzero-lftb(n,i)) &
                    & *(1./(tg1d(n,i)-lftb(n,i)))**2
              rai = cdrx(n,i)*vspda(n,i)*rhs1d(n,i)
              cgrnds(n,i) = rai*cpd
              cgrndl(n,i) = rai*qsatd
              cgrnd(n,i) = cgrnds(n,i) + cgrndl(n,i)*htvp(n,i)
 
!l            3.2  sensible and latent fluxes using soil temperatures
!l            from previous time step
              delq1d(n,i) = (qs1d(n,i)-qg1d(n,i))*gwet1d(n,i)
              delt1d(n,i) = ts1d(n,i) - tg1d(n,i)
              evpr1d(n,i) = -rai*delq1d(n,i)
              sent1d(n,i) = -cgrnds(n,i)*delt1d(n,i)
 
!l            3.3  fluxes to subrout tgrund (evap is in kg/m**2/s)
              fseng(n,i) = sent1d(n,i)
              fevpg(n,i) = evpr1d(n,i)
 
!l            3.4  equate canopy to air, for temp, wind over bare grnd;
!l            needed as factors of sigf(=0) in subr water (uaf) and
!l            subr drag (tlef1d(n,i) carried over to next tstep).
              tlef1d(n,i) = ts1d(n,i)
              uaf(n,i) = vspda(n,i)
            end if
          end if
        end do
      end do
 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then   !  check each point
!=======================================================================
!l            4.   vegetation
!=======================================================================
!l            4.1  add precipitation to leaf water
              ldew1d(n,i) = ldew1d(n,i)+dtbat*sigf(n,i)*prcp1d(n,i)
              ldew1d(n,i) = dmax1(ldew1d(n,i),0.D0)
            end if
          end if
        end do
      end do
!
!l    4.2  distribute excess leaf water to soil
      call vcover
      call drip
!
!l    4.3  calculate canopy temperature, soil and total fluxes,
!l    and leaf water change by evapotranspiration
      call lftemp(iemiss)
!
!l    4.4  calculate carbon sources and sinks
!cc   call co2
!=======================================================================
!l    5.   back to any surface but ocean
!=======================================================================
!
!l    5.1  over sea ice
      call tseaice
!
!l    5.2  over land, calculate soil temp and surface hydrology
      call tgrund
      call snow
      call water
!
!l    5.3  over ocean
!l    set snow cover to zero in case it was previously sea ice
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).lt.0.5 ) then
            scv1d(n,i) = 0.
            sag1d(n,i) = 0.
          end if
        end do
      end do
 
!=======================================================================
!l    6.   linkage to meteorological model
!=======================================================================
 
      do i = 2 , iym1
        do n = 1 , nnsg
 
          if ( ldoc1d(n,i).lt.0.5 .or. &
               lveg(n,i).eq.14 .or. lveg(n,i).eq.15 ) gwet1d(n,i) = 1.
 
!l        6.1  rate of momentum transfer per velocity
          drag1d(n,i) = -cdrx(n,i)*vspda(n,i)*rhs1d(n,i)
          drag1d(n,i) = -drag1d(n,i)        ! for coupling with mm4
 
!l        6.3  latent and heat fluxes over ocean, plus a dummy taf
          if ( ldoc1d(n,i).lt.0.5 ) then
            tlef1d(n,i) = ts1d(n,i)
            fact = -rhs1d(n,i)*cdrx(n,i)*vspda(n,i)
            delq1d(n,i) = (qs1d(n,i)-qg1d(n,i))*gwet1d(n,i)
            delt1d(n,i) = ts1d(n,i) - tg1d(n,i)
!           evaporation is in kg/m**2/s
            evpr1d(n,i) = fact*delq1d(n,i)
            sent1d(n,i) = fact*cpd*delt1d(n,i)
          end if
          if ( sigf(n,i).lt.0.001 ) taf1d(n,i) = tg1d(n,i)
 
!l        6.2  parameters for temperature difference at anemometer level
!cc       zdelt(i) = zdelt(i)*delt1d(n,i)
 
!l        6.4  evaporative flux, accounting for sublimation
!cc       evprr(i) = wlhv*(evpr1d(n,i)-fevpg) + htvp(n,i)*fevpg
 
!l        6.5  nondimensional equivalent bucket capacity for comparisons
!l        with bucket models; usually 1 or less, except where
!l        saturated (then around 2)
!cc       rh2ox(i) = (watr-wiltr)/(relfc-wiltr)
 
        end do
      end do
 
      call time_end(subroutine_name,idindx) 
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
!
      integer :: n , i
      character (len=50) :: subroutine_name='vcover'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) seasb(n,i)                        &
               & = dmax1(0.D0,1.-0.0016*dmax1(298.-tgb1d(n,i),0.D0)**2)
          end if
        end do
      end do
 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              xlai(n,i) = xla(lveg(n,i))
              xlai(n,i) = xlai(n,i) + (xlai0(lveg(n,i))-xlai(n,i))      &
                         & *(1.-seasb(n,i))
              rlai(n,i) = xlai(n,i) + sai(lveg(n,i))
              xlsai(n,i) = xlai(n,i) + sai(lveg(n,i))
              vegt(n,i) = sigf(n,i)*xlsai(n,i)
            end if
          end if
        end do
      end do
!
      call time_end(subroutine_name,idindx) 
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
!
      integer :: n , i
      character (len=50) :: subroutine_name='drip'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
!             ***********         xrun = leaf drip ; sdrop = snow drop
!             off foliage
              etrrun(n,i) = 0.
              xrun(n,i) = ldew1d(n,i) - dewmx*vegt(n,i)
              sdrop(n,i) = 0.
 
!             ***********         test on maximum value of dew
              if ( xrun(n,i).gt.0. ) then
                etrrun(n,i) = xrun(n,i)
                ldew1d(n,i) = dewmx*vegt(n,i)
              end if
 
!             ***********         below freezing excess leaf water
!             falls as snow
              if ( (xrun(n,i).gt.0.) .and. (tm(n,i).lt.tzero) ) then
                etrrun(n,i) = 0.
                sdrop(n,i) = xrun(n,i)
              end if
            end if
          end if
        end do
      end do
 
      call time_end(subroutine_name,idindx) 
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
!         sea-ice mask could be reset in here with "imelt" - but not
!                  done at present
!=======================================================================
!
      subroutine tseaice
!
      implicit none
!
      real(8) :: bb , fact , fss , hrl , hs , hsl , qgrnd , ratsi ,     &
               & rhosw3 , rsd1 , rss , smc4 , smt , tg , tgrnd , wss ,  &
               & wtt
      integer :: n , i
      character (len=50) :: subroutine_name='tseaice'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
      do i = 2 , iym1
        do n = 1 , nnsg
 
          ! lake model handles this case
          if ( lakemod.eq.1 .and. oveg(n,i).eq.14 ) exit
 
          if ( ldoc1d(n,i).gt.1.5 ) then
! ******                rhosw = density of snow relative to water
            rhosw3 = rhosw(n,i)**3
            imelt(n,i) = 0
! ******                cice = specific heat of sea-ice per unit volume
            rsd1 = cice*sice1d(n,i)/1000.0
            if ( scv1d(n,i).gt.0.0 ) then
              rss = csnw*scv1d(n,i)/1000.0
              ratsi = scv1d(n,i)/(1.4*rhosw3*sice1d(n,i))
              wtt = 1./(1.+ratsi)
              wss = (scv1d(n,i)+2.8*rhosw3*sice1d(n,i))                 &
                  & /(scv1d(n,i)+1.4*rhosw3*sice1d(n,i))
! ******                include snow heat capacity
              rsd1 = 0.5*(wss*rss+wtt*rsd1)
            end if
            tgb1d(n,i) = -2.0 + tzero
! ******                subsurface heat flux through ice
            fss = 7.E-4*(tgb1d(n,i)-tg1d(n,i))                          &
                & *ch2o*rhosw3/(scv1d(n,i)+1.4*rhosw3*sice1d(n,i))
            sice1d(n,i) = sice1d(n,i) + fss*dtbat/wlhf*1.087
 
! ******                set sea ice parameter for melting and return
            if ( sice1d(n,i).le.0.0 ) then
              imelt(n,i) = 1
              exit
            end if
! ******                assume lead ocean temp is -1.8c
! ******                flux of heat and moisture through leads
! ******                sat. mixing ratio at t=-1.8c is 3.3e-3
            qice(n,i) = 3.3E-3*stdp/p1d(n,i)
!
!  determine effective surface fluxes over ice, allowing for leads;
!  aarea(n,i) is set in subroutine drag.
!
            tlef1d(n,i) = ts1d(n,i)
            qgrnd = ((1.-aarea(n,i))*cdr(n,i)*qg1d(n,i)+aarea(n,i)      &
                  & *clead(n,i)*qice(n,i))/cdrx(n,i)
            tgrnd = ((1.-aarea(n,i))*cdr(n,i)*tg1d(n,i)+aarea(n,i)      &
                  & *clead(n,i)*(tzero-1.8))/cdrx(n,i)
            fact = -rhs1d(n,i)*cdrx(n,i)*vspda(n,i)
            delq1d(n,i) = (qs1d(n,i)-qgrnd)*gwet1d(n,i)
            delt1d(n,i) = ts1d(n,i) - tgrnd
! ******           output fluxes, averaged over leads and ice
            evpr1d(n,i) = fact*delq1d(n,i)
            sent1d(n,i) = fact*cpd*delt1d(n,i)
            hrl = rhs1d(n,i)*vspda(n,i)*clead(n,i)*(qice(n,i)-qs1d(n,i))
            hsl = rhs1d(n,i)*vspda(n,i)*clead(n,i)                      &
                & *(tzero-1.8-ts1d(n,i))*cpd
! ******           get fluxes over ice for sublimation (subrout snow)
! ******              and melt (below) calculation
            fseng(n,i) = (sent1d(n,i)-aarea(n,i)*hsl)/(1.-aarea(n,i))
            fevpg(n,i) = (evpr1d(n,i)-aarea(n,i)*hrl)/(1.-aarea(n,i))
            hs = fsw1d(i) - flw1d(i) - fseng(n,i) - wlhs*fevpg(n,i)
            bb = dtbat*(hs+fss)/rsd1
! ******           snow melt
            sm(n,i) = 0.
            if ( tg1d(n,i).ge.tzero ) sm(n,i) = (hs+fss)/wlhf
            if ( sm(n,i).le.0. ) sm(n,i) = 0.
            smc4 = sm(n,i)*dtbat
            if ( scv1d(n,i).lt.smc4 ) then
! ******                all snow removed, melt ice
              smt = (scv1d(n,i)/dtbat)
! ******                rho(h2o)/rho(ice) = 1.087
              sice1d(n,i) = sice1d(n,i) + dtbat*(smt-sm(n,i))*1.087
              sm(n,i) = smt
              tg1d(n,i) = tzero
! ******                set sea ice parameter for melting and return
              if ( sice1d(n,i).le.0.0 ) then
                imelt(n,i) = 1
                exit
              end if
            else
!  **********             snow or ice with no snow melting
              tg = tg1d(n,i) + bb
              if ( tg.ge.tzero ) tg1d(n,i) = tzero
              if ( tg.lt.tzero ) tg1d(n,i) = tg
            end if
          end if
        end do
      end do
 
      call time_end(subroutine_name,idindx) 
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
!     rno1d(n,i) = total runoff (mm/day)
!     rnos1d(n,i) = surface runoff (mm/day)
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
!
      implicit none
!
      real(8) :: b , bfac , bfac2 , delwat , est0 , evmax , evmxr ,     &
               & evmxt , rap , vakb , wtg2c , xxkb
      real(8) , dimension(nnsg,iym1) :: gwatr , rnof , rsubsr ,         &
           & rsubss , rsubst , rsur , wflux1 , wflux2 , wfluxc , xkmx1 ,&
           & xkmx2 , xkmxr
      integer :: n , i
      character (len=50) :: subroutine_name='water'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
!=======================================================================
!     1.   define soil water fluxes
!=======================================================================
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 .and. ldoc1d(n,i).lt.1.5 ) then
!
!           1.1  reduce infiltration for frozen ground
!
            if ( tgb1d(n,i).gt.tzero ) then
              xkmxr(n,i) = xkmx(n,i)
            else
              xkmxr(n,i) = 0.
            end if
!
!           1.11 permafrost or ice sheet
!
            if ( lveg(n,i).eq.9 .or. lveg(n,i).eq.12 ) then
              xkmx1(n,i) = 0.
              xkmx2(n,i) = 0.
            else
              xkmx1(n,i) = xkmx(n,i)
              xkmx2(n,i) = drain
            end if
!
!           1.2  diffusive fluxes
!
            evmxr = evmx0(n,i)*xkmxr(n,i)/xkmx(n,i)
            evmxt = evmx0(n,i)*xkmx1(n,i)/xkmx(n,i)
            b = bsw(n,i)
            bfac = watr(n,i)**(3.+bfc(n,i))*watu(n,i)**(b-bfc(n,i)-1)
            bfac2 = watt(n,i)**(2.+bfc(n,i))*watr(n,i)**(b-bfc(n,i))
            wfluxc(n,i) = evmxr*(depuv(lveg(n,i))/deprv(lveg(n,i)))     &
                         & **0.4*bfac
            wflux1(n,i) = wfluxc(n,i)*watr(n,i)
            wflux2(n,i) = evmxt*(depuv(lveg(n,i))/deprv(lveg(n,i)))     &
                         & **0.5*bfac2*(watt(n,i)-watr(n,i))
!
!           1.3  gravitational drainage
!
            rsubss(n,i) = xkmxr(n,i)*watr(n,i)**(b+0.5)*watu(n,i)       &
                         & **(b+2.5)
            rsubsr(n,i) = xkmx1(n,i)*watt(n,i)**(b+0.5)*watr(n,i)       &
                         & **(b+2.5)
            rsubst(n,i) = xkmx2(n,i)*watt(n,i)**(2.*b+3.)
!
!           1.32 bog or water
!
            if ( (lveg(n,i).ge.13) .and. (lveg(n,i).le.15) ) then
              rsubst(n,i) = 0.0
              rsubss(n,i) = 0.0
              rsubsr(n,i) = 0.0
            end if
!
!           1.4  fluxes through internal surfaces
!
            wflux1(n,i) = wflux1(n,i) - rsubss(n,i)
            wflux2(n,i) = wflux2(n,i) - rsubsr(n,i)
          end if
        end do
      end do
!
!     1.5  net flux at air interface
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 .and. ldoc1d(n,i).lt.1.5 ) then
            gwatr(n,i) = pw(n,i) - evapw(n,i) + sm(n,i)                 &
                        & + etrrun(n,i)/dtbat
!
!=======================================================================
!           2.   define runoff terms
!=======================================================================
!
!           2.1  surface runoff
!
            wata(n,i) = 0.5*(watu(n,i)+watr(n,i))
!
!           2.11 increase surface runoff over frozen ground
!
            if ( tg1d(n,i).lt.tzero ) then
              rsur(n,i) = dmin1(1.D0,wata(n,i)**1)*                     &
                         & dmax1(0.D0,gwatr(n,i))
            else
              rsur(n,i) = dmin1(1.D0,wata(n,i)**4)                      &
                         & *dmax1(0.D0,gwatr(n,i))
            end if
!
!           2.12 irrigate cropland
!
            if ( lveg(n,i).eq.10 .and. watr(n,i).lt.relfc(n,i) )        &
               & rsur(n,i) = rsur(n,i) + (rsw1d(n,i)-relfc(n,i)*        &
               &             gwmx1(n,i))/dtbat
!
!           2.13 saturate swamp or rice paddy
!
            if ( (lveg(n,i).ge.13) .and. (lveg(n,i).lt.16) )            &
               & rsur(n,i) = rsur(n,i) + dmin1(0.D0,(rsw1d(n,i)-        &
               &             gwmx1(n,i))/dtbat)
!
!           2.2  total runoff
!
            rnof(n,i) = rsur(n,i) + rsubst(n,i)
!
!=======================================================================
!           3.   increment soil moisture
!=======================================================================
!
!           3.1  update top layer with implicit treatment
!           of flux from below
!
            ssw1d(n,i) = ssw1d(n,i) + dtbat*(gwatr(n,i)-efpr(n,i)*      &
                        & etr(n,i)-rsur(n,i)+wflux1(n,i))
            ssw1d(n,i) = ssw1d(n,i)/(1.+wfluxc(n,i)*dtbat/gwmx0(n,i))
!
!           3.2  update root zone
!
            rsw1d(n,i) = rsw1d(n,i) + dtbat*(gwatr(n,i)-etr(n,i)-       &
                        & rsur(n,i)+wflux2(n,i))
!
!           3.3  update total water
!
            tsw1d(n,i) = tsw1d(n,i) + dtbat*(gwatr(n,i)-etr(n,i)-       &
                        & rnof(n,i))
          end if
        end do
      end do
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 .and. ldoc1d(n,i).lt.1.5 ) then
!
!=======================================================================
!           4.   check whether soil water exceeds maximum capacity or
!           becomes negative (should rarely or never happen)
!=======================================================================
!
!           4.1  surface water assumed to move downward into soil
!
            if ( ssw1d(n,i).gt.gwmx0(n,i) ) ssw1d(n,i) = gwmx0(n,i)
!
!           4.2  excess root layer water assumed to move downward
!
            if ( rsw1d(n,i).gt.gwmx1(n,i) ) rsw1d(n,i) = gwmx1(n,i)
!
!           4.3  excess total water assumed to go to subsurface runoff
!
            if ( tsw1d(n,i).gt.gwmx2(n,i) ) then
              delwat = tsw1d(n,i) - gwmx2(n,i)
              tsw1d(n,i) = gwmx2(n,i)
              rsubst(n,i) = rsubst(n,i) + delwat/dtbat
            end if
!
!           4.4  check for negative water in top layer
!
            if ( ssw1d(n,i).le.1.E-2 ) ssw1d(n,i) = 1.E-2
!
!=======================================================================
!           5.   accumulate leaf interception
!=======================================================================
!
            ircp1d(n,i) = ircp1d(n,i) + sigf(n,i)*(dtbat*               &
                         & prcp1d(n,i)) - (sdrop(n,i)+etrrun(n,i))
!
!=======================================================================
!           6.   evaluate runoff (incremented in ccm)
!=======================================================================
!
!*          update total runoff
!
            rnof(n,i) = rsur(n,i) + rsubst(n,i)
            rno1d(n,i) = rnof(n,i)*tau1
            rnos1d(n,i) = rsur(n,i)*tau1
          else                       ! ocean or sea ice
            rnof(n,i) = 0.
            rno1d(n,i) = 0.
            rnos1d(n,i) = 0.
          end if
        end do
      end do
!
!=======================================================================
!     7.   calculate potential evaporation and use mod_to determine
!     wetness factor, allowing for snow being saturated
!=======================================================================
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 .and. ldoc1d(n,i).lt.1.5 ) then
            xxkb = dmin1(rough(lveg(n,i)),1.D0)
            vakb = (1.-sigf(n,i))*vspda(n,i) + sigf(n,i)                &
                 & *(xxkb*uaf(n,i)+(1.-xxkb)*vspda(n,i))
            wtg2c = (1.-sigf(n,i))*cdrx(n,i)*vakb
            rap = rhs1d(n,i)*(csoilc*uaf(n,i)*sigf(n,i)*(qg1d(n,i)+     &
                & delq1d(n,i)-qs1d(n,i))+wtg2c*(qg1d(n,i)-qs1d(n,i)))
            bfac = watr(n,i)**(3.+bfc(n,i))*watu(n,i)                   &
                 & **(bsw(n,i)-bfc(n,i)-1)
            est0 = evmx0(n,i)*bfac*watu(n,i)
            evmax = dmax1(est0,0.D0)
            gwet1d(n,i) = dmin1(1.D0,evmax/dmax1(1.D-14,rap))
            gwet1d(n,i) = scvk(n,i) + gwet1d(n,i)*(1.0-scvk(n,i))
          end if
        end do
      end do
!
      call time_end(subroutine_name,idindx) 
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
!
      implicit none
!
      real(8) :: age1 , age2 , age3 , arg , arg2 , dela , dela0 , dels ,&
               & sge , tage
      integer :: n , i
      real(8) , dimension(nnsg,iym1) :: sold
!
      character (len=50) :: subroutine_name='snow'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
      age3 = 0.3
 
!=======================================================================
!     1.   partition soil evaporation and precipitation
!     between water and snow
!=======================================================================
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
 
            evapw(n,i) = fevpg(n,i)
            evaps(n,i) = scvk(n,i)*evapw(n,i)
            if ( ldoc1d(n,i).gt.1.5 ) evaps(n,i) = fevpg(n,i)
            evapw(n,i) = (1.-scvk(n,i))*evapw(n,i)
!
!           ******                tm  is temperature of precipitation
            if ( tm(n,i).ge.tzero ) then
              pw(n,i) = prcp1d(n,i)*(1.-sigf(n,i))
              ps(n,i) = 0.0
            else
!             ******                snowing
              pw(n,i) = 0.0
              ps(n,i) = prcp1d(n,i)*(1.-sigf(n,i))
            end if
          end if
        end do
      end do
!
!=======================================================================
!     2.   update snow cover
!=======================================================================
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            sold(n,i) = scv1d(n,i)
            scv1d(n,i) = scv1d(n,i) + dtbat                             &
                        & *(ps(n,i)-evaps(n,i)-sm(n,i)) + sdrop(n,i)
            scv1d(n,i) = dmax1(scv1d(n,i),0.D0)
            sag1d(n,i) = dmax1(sag1d(n,i),0.D0)
 
!           ******           snow cover except for antarctica
!=======================================================================
!           3.   increment non-dimensional "age" of snow;
!           10 mm snow restores surface to that of new snow.
!=======================================================================
            if ( scv1d(n,i).gt.0. ) then
              arg = 5.E3*(1./tzero-1./tg1d(n,i))
              age1 = dexp(arg)
              arg2 = dmin1(0.D0,10.*arg)
              age2 = dexp(arg2)
              tage = age1 + age2 + age3
              dela0 = 1.E-6*dtbat
              dela = dela0*tage
              dels = 0.1*dmax1(0.D0,scv1d(n,i)-sold(n,i))
              sge = (sag1d(n,i)+dela)*(1.0-dels)
              sag1d(n,i) = dmax1(0.D0,sge)
            end if
 
!           ******           antarctica
            if ( scv1d(n,i).gt.800. ) sag1d(n,i) = 0.
          end if
        end do
      end do
 
      call time_end(subroutine_name,idindx) 
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
!       tau1 = seconds in a day
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
!
      implicit none
!
      real(8) , dimension(nnsg,iym1) :: bb , bcoef , cc , depann ,      &
           & depdiu , deprat , fct2 , hs , rscsa , rscsd , ska , skd ,  &
           & sks , swtrta , swtrtd
      real(8) :: bcoefd , bcoefs , c31 , c3t , c41 , c4t , cder , depr ,&
             & depu , xdt2 , xdtime , dtimea , froze2 , frozen , rscss ,&
             & t3 , tbef , tg , tinc , wtas , wtax , wtd , wtds ,       &
             & xkperi , xnu , xnua
      real(8) :: dtbat2 , rdtbat2
      integer :: n , i
      character (len=50) :: subroutine_name='tgrund'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
! 
      dtbat2 = dtbat*2.
      rdtbat2 = 1./dtbat2

!=======================================================================
!l    1.   define thermal conductivity, heat capacity,
!l    and other force restore parameters
!=======================================================================
      xnu = 2.*mathpi/tau1
      xnua = xnu/365.
      xdtime = dtbat*xnu
      dtimea = dtbat*xnua
      xdt2 = 0.5*xdtime
      xkperi = 1.4E-6
 
!l    3.4  permafrost temperature
      t3 = 271.
 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 .and. ldoc1d(n,i).lt.1.5 ) then
 
!l          1.1  frozen ground values using 44 m2/yr for frozen ground
!l          thermal diffusion coefficient, based on the values of
!l          50 and 38 quoted by osterkamp; ice contribution to
!l          specific heat only o.49 that of water
 
            swtrtd(n,i) = watu(n,i)*porsl(n,i)
            if ( tg1d(n,i).lt.tzero ) then
              frozen = 0.85*dmin1(1.D0,.25*(tzero-tg1d(n,i)))
              skd(n,i) = xkperi
              rscsd(n,i) = fsc(swtrtd(n,i)*(1.-0.51*frozen))
            else
              skd(n,i) = fsk(swtrtd(n,i))*texrat(n,i)
              rscsd(n,i) = fsc(swtrtd(n,i))
            end if
            swtrta(n,i) = watr(n,i)*porsl(n,i)
            if ( tgb1d(n,i).lt.tzero ) then
              froze2 = 0.85*dmin1(1.D0,.25*(tzero-tgb1d(n,i)))
              ska(n,i) = xkperi
              rscsa(n,i) = fsc(swtrta(n,i)*(1.-0.51*froze2))
            else
              ska(n,i) = fsk(swtrta(n,i))*texrat(n,i)
              rscsa(n,i) = fsc(swtrta(n,i))
            end if
 
!l          1.2  correct for snow cover, if significant
            depdiu(n,i) = dsqrt(2.*skd(n,i)/xnu)
            bcoef(n,i) = xdtime*depdiu(n,i)/(rscsd(n,i)*skd(n,i))
            if ( scv1d(n,i).gt.1.0 ) then
              wtd = dexp(-2.*scrat(n,i)/depdiu(n,i))
              rscss = csnw*rhosw(n,i)
              sks(n,i) = 7.E-7*cws*rhosw(n,i)
              bcoefs = dsqrt(2.*sks(n,i)/xnu)/(rscss*sks(n,i))
              wtds = (1.-wtd)*scvk(n,i)
              bcoefd = dsqrt(2.*skd(n,i)/xnu)/(rscsd(n,i)*skd(n,i))
              bcoef(n,i) = xdtime*(wtds*bcoefs+(1.-wtds)*bcoefd)
              depdiu(n,i) = wtds*dsqrt(2.*sks(n,i)/xnu) + (1.-wtds)     &
                           & *depdiu(n,i)
            end if
            depann(n,i) = dsqrt(2.*ska(n,i)/xnua)
            if ( scv1d(n,i).gt.20. ) then
              wtax = dexp(-2.*scrat(n,i)/depann(n,i))
              wtas = (1.-wtax)*scvk(n,i)
              depann(n,i) = wtas*dsqrt(2.*sks(n,i)/xnua) + (1.-wtas)    &
                           & *depann(n,i)
            end if
            deprat(n,i) = depann(n,i)/depdiu(n,i)
 
!=======================================================================
!l          2.   collect force restore terms
!=======================================================================
            cc(n,i) = 1.0
            fct2(n,i) = 0.
!
!l          2.1  add freezing thermal inertia
            if ( (tg1d(n,i).lt.tzero) .and. (tg1d(n,i).gt.(tzero-4.))   &
               & .and. (sice1d(n,i).le.1.E-22) ) then
              depu = depuv(lveg(n,i))*1.E-3
              cc(n,i) = 1. + dmax1(ssw1d(n,i)-frezu(lveg(n,i)),0.D0)    &
                       & *fct1(depu*rscsd(n,i))
            end if
            if ( (tgb1d(n,i).lt.tzero) .and.                            &
               & (tgb1d(n,i).gt.(tzero-4.)) .and.                       &
               & (sice1d(n,i).le.1.E-22) ) then
              depr = deprv(lveg(n,i))*1.E-3
              fct2(n,i) = dmax1(rsw1d(n,i)-freza(lveg(n,i)),0.D0)       &
                         & *fct1(depr*rscsa(n,i))
            end if
 
!l          2.2  large thermal inertial for permanent ice cap
            if ( lveg(n,i).eq.12 ) fct2(n,i) = 1.E3*fct2(n,i)
 
!l          2.3  collect energy flux terms
            rnet(n,i) = fsw1d(i) - sigf(n,i)*(sabveg(i)-flnet(n,i))     &
                       & - (1.-sigf(n,i))                               &
                       & *(flw1d(i)-sigf(n,i)*flneto(n,i))
            hs(n,i) = rnet(n,i) - fseng(n,i) - fevpg(n,i)*htvp(n,i)
            bb(n,i) = bcoef(n,i)*hs(n,i) + xdtime*tgb1d(n,i)
 
!l          2.4  add in snowmelt (melt enough snow to reach freezing
!           temp)
            sm(n,i) = 0.0
            if ( scv1d(n,i).gt.0.0 ) then
              cder = bcoef(n,i)*cgrnd(n,i)
              sm(n,i) = (bb(n,i)+(cc(n,i)-xdt2+cder)*tg1d(n,i)-tzero    &
                       & *(cc(n,i)+xdt2+cder))/(bcoef(n,i)*wlhf)
!             **********              snow melt always between 0 and
!             total snow
              sm(n,i) = dmax1(0.D0,dmin1(sm(n,i),scv1d(n,i)*2.*         &
                       & rdtbat2))
              bb(n,i) = bb(n,i) - bcoef(n,i)*wlhf*sm(n,i)
            end if
          end if
        end do
      end do
 
!=======================================================================
!l    3.   update soil temperatures
!=======================================================================
!l    3.1  update surface soil temperature
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 .and. ldoc1d(n,i).lt.1.5 ) then
            tbef = tg1d(n,i)
            cder = bcoef(n,i)*cgrnd(n,i)
            tg = (bb(n,i)+(cc(n,i)-xdt2+cder)*tg1d(n,i))/(cc(n,i)+      &
                   & xdt2+cder)
            tg1d(n,i) = tg
 
!l          3.2  put brakes on large temperature excursions
            tg1d(n,i) = dmin1(tbef+10.,tg1d(n,i))
            tg1d(n,i) = dmax1(tbef-10.,tg1d(n,i))
 
!l          3.3  correct fluxes to present soil temperature
            tinc = tg1d(n,i) - tbef
            sent1d(n,i) = sent1d(n,i) + tinc*cgrnds(n,i)
            evpr1d(n,i) = evpr1d(n,i) + tinc*cgrndl(n,i)
            fseng(n,i) = fseng(n,i) + tinc*cgrnds(n,i)
            fevpg(n,i) = fevpg(n,i) + tinc*cgrndl(n,i)
!
!l          3.5  couple to deep temperature in permafrost
!l          3.6  update subsoil temperature
            if ( lveg(n,i).eq.9 .or. lveg(n,i).eq.12 ) then
              c31 = 0.5*dtimea*(1.+deprat(n,i))
              c41 = dtimea*deprat(n,i)
              tgb1d(n,i) = ((1.-c31+fct2(n,i))*tgb1d(n,i)+c41*tg1d(n,i)+&
                    & dtimea*t3)/(1.+c31+fct2(n,i))
            else
              c3t = 0.5*dtimea*deprat(n,i)
              c4t = dtimea*deprat(n,i)
              tgb1d(n,i) = ((1.-c3t+fct2(n,i))*tgb1d(n,i)+c4t*tg1d(n,i))&
                     & /(1.+c3t+fct2(n,i))
            end if
          end if
        end do
      end do
      call time_end(subroutine_name,idindx) 

      contains

      function fsk(x)
        implicit none
        real(8) :: fsk
        real(8) , intent(in) :: x
        fsk = (2.9D-7*x+4.D-9)/(((1.0D0-0.6D0*x)*x+0.09D0)*(0.23D0+x))
      end function fsk
      function fsc(x)
        implicit none
        real(8) :: fsc
        real(8) , intent(in) :: x
        fsc = (0.23D0+x)*4.186D6
      end function fsc
      function fct1(x)
        implicit none
        real(8) :: fct1
        real(8) , intent(in) :: x
        fct1 = wlhf*0.25D0*1.414D0/x
      end function fct1
! 
      end subroutine tgrund
!
      end module mod_bndry
