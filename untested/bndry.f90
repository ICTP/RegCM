!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
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
      use regcm_param
      use bats
      use ictp01
      implicit none
!
! Dummy arguments
!
      integer :: iemiss
!
! Local variables
!
      real(8) :: fact , qsatd , rai
      integer :: n , np
 
!=======================================================================
!l    1.   initialize
!=======================================================================
 
      do np = np1 , npts
        do n = 1 , nnsg
 
          htvp(n,np) = c(125)
          if ( (tg1d(n,np).lt.c(67) .and. ldoc1d(n,np).gt.0.5) .or.     &
             & scv1d(n,np).gt.0. ) htvp(n,np) = c(126)
          sdrop(n,np) = 0.
          etrrun(n,np) = 0.
          flnet(n,np) = 0.
          fevpg(n,np) = 0.
          fseng(n,np) = 0.
          vegt(n,np) = 0.
          efpr(n,np) = 0.
          etr(n,np) = 0.
!         **********            switch between rain and snow /tm is
!         ref. temp set= anemom temp -2.2
          tm(n,np) = ts1d(n,np) - 2.2
 
!*        soil moisture ratio (to max) as used in subrouts tgrund,
!*        water, and root (called by lftemp): watu=upper, watr=root,
 
!         watt=total
          if ( lveg(n,np).ge.1 ) then
            watu(n,np) = ssw1d(n,np)/gwmx0(n,np)
            watr(n,np) = rsw1d(n,np)/gwmx1(n,np)
            watt(n,np) = tsw1d(n,np)/gwmx2(n,np)
            watu(n,np) = dmin1(watu(n,np),1.D0)
            watr(n,np) = dmin1(watr(n,np),1.D0)
            watt(n,np) = dmin1(watt(n,np),1.D0)
            watr(n,np) = dmax1(watr(n,np),1.D-4)
            watu(n,np) = dmax1(watu(n,np),1.D-4)
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
 
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).le.0.001 .and. ldoc1d(n,np).lt.1.5 ) then
              qsatd = qg1d(n,np)*gwet1d(n,np)*a(n,np)*(c(67)-b(n,np))   &
                    & *(1./(tg1d(n,np)-b(n,np)))**2
              rai = cdrx(n,np)*vspda(n,np)*rhs1d(n,np)
              cgrnds(n,np) = rai*c(58)
              cgrndl(n,np) = rai*qsatd
              cgrnd(n,np) = cgrnds(n,np) + cgrndl(n,np)*htvp(n,np)
 
!l            3.2  sensible and latent fluxes using soil temperatures
!l            from previous time step
              delq1d(n,np) = (qs1d(n,np)-qg1d(n,np))*gwet1d(n,np)
              delt1d(n,np) = ts1d(n,np) - tg1d(n,np)
              evpr1d(n,np) = -rai*delq1d(n,np)
              sent1d(n,np) = -cgrnds(n,np)*delt1d(n,np)
 
!l            3.3  fluxes to subrout tgrund (evap is in kg/m**2/s)
              fseng(n,np) = sent1d(n,np)
              fevpg(n,np) = evpr1d(n,np)
 
!l            3.4  equate canopy to air, for temp, wind over bare grnd;
!l            needed as factors of sigf(=0) in subr water (uaf) and
!l            subr drag (tlef1d(n,np) carried over to next tstep).
              tlef1d(n,np) = ts1d(n,np)
              uaf(n,np) = vspda(n,np)
            end if
          end if
        end do
      end do
 
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) then   !  check each point
!=======================================================================
!l            4.   vegetation
!=======================================================================
!l            4.1  add precipitation to leaf water
              ldew1d(n,np) = ldew1d(n,np) + c(4)*sigf(n,np)*prcp1d(n,np)
              ldew1d(n,np) = dmax1(ldew1d(n,np),0.D0)
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
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).lt.0.5 ) then
            scv1d(n,np) = 0.
            sag1d(n,np) = 0.
          end if
        end do
      end do
 
!=======================================================================
!l    6.   linkage to meteorological model
!=======================================================================
 
      do np = np1 , npts
        do n = 1 , nnsg
 
          if ( ldoc1d(n,np).lt.0.5 .or. lveg(n,np).eq.14 .or. lveg(n,np)&
             & .eq.15 ) gwet1d(n,np) = 1.
 
!l        6.1  rate of momentum transfer per velocity
          drag1d(n,np) = -cdrx(n,np)*vspda(n,np)*rhs1d(n,np)
          drag1d(n,np) = -drag1d(n,np)        ! for coupling with mm4
 
!l        6.3  latent and heat fluxes over ocean, plus a dummy taf
          if ( ldoc1d(n,np).lt.0.5 ) then
            tlef1d(n,np) = ts1d(n,np)
            fact = -rhs1d(n,np)*cdrx(n,np)*vspda(n,np)
            delq1d(n,np) = (qs1d(n,np)-qg1d(n,np))*gwet1d(n,np)
            delt1d(n,np) = ts1d(n,np) - tg1d(n,np)
!           evaporation is in kg/m**2/s
            evpr1d(n,np) = fact*delq1d(n,np)
            sent1d(n,np) = fact*c(58)*delt1d(n,np)
          end if
          if ( sigf(n,np).lt.0.001 ) taf1d(n,np) = tg1d(n,np)
 
!l        6.2  parameters for temperature difference at anemometer level
!cc       zdelt(i) = zdelt(i)*delt1d(n,np)
 
!l        6.4  evaporative flux, accounting for sublimation
!cc       evprr(i) = c(125)*(evpr1d(n,np)-fevpg) + htvp(n,np)*fevpg
 
!l        6.5  nondimensional equivalent bucket capacity for comparisons
!l        with bucket models; usually 1 or less, except where
!l        saturated (then around 2)
!cc       rh2ox(i) = (watr-wiltr)/(relfc-wiltr)
 
        end do
      end do
 
      end subroutine bndry
