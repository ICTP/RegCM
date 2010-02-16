C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine bndry(iemiss)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     this is the main routine when interfacing with a
c     meteorological model
c
c
c                f l o w   d i a g r a m   f o r   b n d r y
c
c                bndry ===>    drag ===> dragdn ===> depth
c                             satur
c                            vcover
c                              drip
c                            lftemp ======================> stomat
c                               co2 ===> carbon             frawat
c                            tseice                           root
c                            tgrund                          satur
c                              snow                         lfdrag
c                             water                         condch
c                                                           condcq
c                                                            deriv
c
c
c  **  type1  = crop
c  **  type2  = short grass
c  **  type3  = evergreen needle leaf tree
c  **  type4  = deciduous needle leaf tree
c  **  type5  = deciduous broadleaf tree
c  **  type6  = evergreen brodaleaf tree
c  **  type7  = tall grass
c  **  type8  = desert
c  **  type9  = tundra
c  **  type10 = irrig crop
c  **  type11 = semi-desert
c  **  type12 = ice
c  **  type13 = bog or marsh
c  **  type14 = inland water
c  **  type15 = sea
c  **  type16 = evgr shrub
c  **  type17 = decid shrub
c  **  type18 = mixed tree
c
c  ** note: water and soil parameters are in mm
c
      implicit none
#     include "parameter.inc"

#     include "include/bats.cb"
      integer iemiss

      real(kind=8)  a,b
      common/ab/ a(NNSG,nbmax), b(NNSG,nbmax)

      integer n,np
      real(kind=8)  qsatd,rai,fact

c=======================================================================
cl            1.   initialize
c=======================================================================

      do 90 np=np1,npts
      do n=1,NNSG

       htvp(n,np) = c(125)
       if ((tg1d(n,np) .lt. c(67) .and. ldoc1d(n,np) .gt. 0.5)
     &    .or. scv1d(n,np) .gt. 0.) htvp(n,np) = c(126)
       sdrop(n,np) = 0.
       etrrun(n,np) = 0.
       flnet(n,np) = 0.
       fevpg(n,np) = 0.
       fseng(n,np) = 0.
       vegt(n,np) = 0.
       efpr(n,np) = 0.
       etr(n,np) = 0.
c **********            switch between rain and snow /tm is ref. temp
c                         set= anemom temp -2.2
       tm(n,np) = ts1d(n,np) - 2.2

c* soil moisture ratio (to max) as used in subrouts tgrund, water,
c*   and root (called by lftemp): watu=upper, watr=root, watt=total

       if (lveg(n,np).ge.1) then
        watu(n,np) = ssw1d(n,np)/gwmx0(n,np)
        watr(n,np) = rsw1d(n,np)/gwmx1(n,np)
        watt(n,np) = tsw1d(n,np)/gwmx2(n,np)
        watu(n,np) = dmin1(watu(n,np),1.d0)
        watr(n,np) = dmin1(watr(n,np),1.d0)
        watt(n,np) = dmin1(watt(n,np),1.d0)
        watr(n,np) = dmax1(watr(n,np),1.d-4)
        watu(n,np) = dmax1(watu(n,np),1.d-4)
       end if

      end do
90    continue

c=======================================================================
cl            2.  calculate transfer coefficients at layer 1
c=======================================================================
c             2.1  sigf corrected in drag: remove snow-covered veg
      call drag

c             2.2  get saturation vapor pressure of soil surface
      call satur( qg1d, tg1d, p1d )

c=======================================================================
cl       3.   bare land
c=======================================================================
cl       3.1  get derivative of fluxes with repect to tg

      do 91 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).le.0.001 .and. ldoc1d(n,np).lt.1.5) then
            qsatd = qg1d(n,np)*gwet1d(n,np)*a(n,np)*(c(67)-b(n,np))
     a            * (1./(tg1d(n,np)-b(n,np)))**2
            rai = cdrx(n,np) * vspda(n,np) * rhs1d(n,np)
            cgrnds(n,np) = rai*c(58)
            cgrndl(n,np) = rai*qsatd
            cgrnd(n,np)  = cgrnds(n,np) + cgrndl(n,np)*htvp(n,np)

cl       3.2  sensible and latent fluxes using soil temperatures
cl            from previous time step
            delq1d(n,np) = (qs1d(n,np)-qg1d(n,np))*gwet1d(n,np)
            delt1d(n,np) = ts1d(n,np)-tg1d(n,np)
            evpr1d(n,np) = -rai*delq1d(n,np)
            sent1d(n,np) = -cgrnds(n,np)*delt1d(n,np)

cl       3.3  fluxes to subrout tgrund (evap is in kg/m**2/s)
            fseng(n,np) = sent1d(n,np)
            fevpg(n,np) = evpr1d(n,np)

cl       3.4  equate canopy to air, for temp, wind over bare grnd;
cl                 needed as factors of sigf(=0) in subr water (uaf) and
cl                 subr drag (tlef1d(n,np) carried over to next tstep).
            tlef1d(n,np) = ts1d(n,np)
            uaf(n,np) = vspda(n,np)
        end if
       end if
      end do
91    continue

      do 92 np = np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if(sigf(n,np).gt.0.001) then          !  check each point
c=======================================================================
cl       4.   vegetation
c=======================================================================
cl       4.1  add precipitation to leaf water
            ldew1d(n,np) = ldew1d(n,np)+c(4)*sigf(n,np)*prcp1d(n,np)
            ldew1d(n,np) = dmax1( ldew1d(n,np), 0.d0 )
        end if
       end if
      end do
92    continue
c
cl       4.2  distribute excess leaf water to soil
      call vcover
      call drip
c
cl       4.3  calculate canopy temperature, soil and total fluxes,
cl            and leaf water change by evapotranspiration
      call lftemp(iemiss)
c
cl       4.4  calculate carbon sources and sinks
ccc   call co2
c=======================================================================
cl       5.   back to any surface but ocean
c=======================================================================
c
cl       5.1  over sea ice
      call tseice
c
cl       5.2  over land, calculate soil temp and surface hydrology
      call tgrund
      call snow
      call water
c
cl       5.3  over ocean
cl            set snow cover to zero in case it was previously sea ice
      do 93 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).lt.0.5) then
         scv1d(n,np) = 0.
         sag1d(n,np) = 0.
       end if
      end do
93    continue

c=======================================================================
cl       6.   linkage to meteorological model
c=======================================================================

      do 94 np=np1,npts
      do n=1,NNSG

       if (ldoc1d(n,np).lt.0.5.or.lveg(n,np).eq.14.or.lveg(n,np).eq.15)
     a         gwet1d(n,np) = 1.

cl       6.1  rate of momentum transfer per velocity
       drag1d(n,np) = -cdrx(n,np)*vspda(n,np)*rhs1d(n,np)
       drag1d(n,np) = -drag1d(n,np)           ! for coupling with mm4

cl       6.3  latent and heat fluxes over ocean, plus a dummy taf
       if (ldoc1d(n,np).lt.0.5) then
         tlef1d(n,np) = ts1d(n,np)
         fact = -rhs1d(n,np)*cdrx(n,np)*vspda(n,np)
         delq1d(n,np) = (qs1d(n,np)-qg1d(n,np))*gwet1d(n,np)
         delt1d(n,np) = ts1d(n,np)-tg1d(n,np)
c               evaporation is in kg/m**2/s
         evpr1d(n,np) = fact*delq1d(n,np)
         sent1d(n,np) = fact*c(58)*delt1d(n,np)
       end if
       if (sigf(n,np).lt.0.001) taf1d(n,np) = tg1d(n,np)

cl       6.2  parameters for temperature difference at anemometer level
ccc    zdelt(i) = zdelt(i)*delt1d(n,np)

cl       6.4  evaporative flux, accounting for sublimation
ccc    evprr(i) = c(125)*(evpr1d(n,np)-fevpg) + htvp(n,np)*fevpg

cl       6.5  nondimensional equivalent bucket capacity for comparisons
cl            with bucket models; usually 1 or less, except where
cl            saturated (then around 2)
ccc    rh2ox(i) = (watr-wiltr)/(relfc-wiltr)

      end do
94    continue

      return
      end subroutine bndry
