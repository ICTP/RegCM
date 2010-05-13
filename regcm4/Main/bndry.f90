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
 
      subroutine bndry(iemiss)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     this is the main routine when interfacing with a
!     meteorological model
!
!
!                f l o w   d i a g r a m   f o r   b n d r y
!
!                bndry ===>    drag ===> dragdn ===> depth
!                             satur
!                            vcover
!                              drip
!                            lftemp ======================> stomat
!                               co2 ===> carbon             frawat
!                            tseice                           root
!                            tgrund                          satur
!                              snow                         lfdrag
!                             water                         condch
!                                                           condcq
!                                                            deriv
!
!
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
      use mod_dynparam
      use mod_param1 , only : dtbat
      use mod_bats , only : htvp , sdrop , etrrun , flnet , fevpg ,     &
                  & fseng , vegt , efpr , etr , ts1d , ssw1d , tm ,     &
                  & watu , watr , watt , gwmx0 , gwmx1 , gwmx2 ,        &
                  & cgrnds , cgrndl , cgrnd , tg1d , delq1d , delt1d ,  &
                  & evpr1d , sent1d , tlef1d , uaf , vspda , ldew1d ,   &
                  & scv1d , sag1d , gwet1d , drag1d , rhs1d , qs1d ,    &
                  & taf1d , p1d , ldoc1d , lveg , rsw1d , tsw1d , qg1d ,&
                  & sigf , cdrx , prcp1d , vcover , snow , water
      use mod_constants , only : tzero , wlhv , wlhs , cpd
      use mod_ictp01
      implicit none
!
! Dummy arguments
!
      integer :: iemiss
!
! Local variables
!
      real(8) :: fact , qsatd , rai
      integer :: n , i
 
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
              qsatd = qg1d(n,i)*gwet1d(n,i)*a(n,i)*(tzero-b(n,i))       &
                    & *(1./(tg1d(n,i)-b(n,i)))**2
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
      call tseice
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
 
          if ( ldoc1d(n,i).lt.0.5 .or. lveg(n,i).eq.14 .or. lveg(n,i)   &
             & .eq.15 ) gwet1d(n,i) = 1.
 
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
 
      end subroutine bndry
