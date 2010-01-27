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

      subroutine vecbats(jslc)
c=======================================================================
cl  based on: bats version 1e          copyright 18 august 1989
c=======================================================================
c      written by robert e. dickinson and patrick j. kennedy with
c contributions from klaus bluemel, filippo giorgi, and ann henderson-
c sellers.  version 1b (1 dec 1988) greatly modifies versions 1 and 1a,
c for example:  general cleanup and stabilization of iterations;
c three-layer soil moisture; improvement and bug fixes by bluemel
c including an improved canopy model; improvement in soil temperature
c force-restore method for overlying snow; stomatal resistance
c dependence on solar radiation made dependent on lai; inclusion of a
c non-zero displacement level; elimination of soil water movement
c in frozen ground; fixing of the carbon routine;
c and derivation of leaf iteration scheme.
c
c      version 1c (1 march 1989) modifies version 1b by:
c adding a dependence of stomatal resistance on vapor pressure deficit,
c allowing for solar zenith dependence of incident light in
c calculating dependence of stomatal resistance on light;
c further streamlining the iterative calculation of leaf temperature;
c and fixing a bug in the soil moisture calculation of bats1b.
c
c      version 1d (1 june 1989) fixes two bugs introduced into 1c,
c in subroutines lftemp and tgrund, and includes minor changes
c suggested by k. bluemel.  it redefines output fluxes of sensible and
c latent heat to agree exactly with those used for soil energy balance
c to ensure conservation of energy.  some other cleanup is done also.
c
c      verson 1e (18 august 1989) makes changes consistent with 1d to
c insure water conservation.  latent heat of sublimation now applies
c to all soil terms where snow-covered or below freezing.  two minor
c bugs from misplaced lines of code are fixed.  leaf temperature
c calculation reformulated to use average longwave rather than over
c bare soil as done previously, to better match ccm1 tie-in.
c
c      vector bats (1989) was developed from bats1e by gary bates.
c it is designed such that one call to bats perfoms calculations at
c a specified arbitrary number of gridpoints instead of just one point
c as in the standard version of bats1e.  this allows a more efficient
c vectorization of the code.  vector bats was coupled to the mm4
c in februrary 1992 by gary bates.
c
c
c                  main drive program
c
c the core model begins with subroutine bndry, which is driven by
c this routine.
c
c *********************************************************************
c
c
c                      * ts (lowest model layer "midpt" at 75m)
c
c
c                         i-----------i
c                         i         ======= taf (temp air in leaves)
c                         i     tlef  i
c                         i-----------i
c                              i  i
c          ta (anemom)         i  i
c          is actually a ss    i  i
c             i   i temp (4')  i  i
c  -----------i---i------------i--i------------------------------------
c                                                  tg
c -------------------------------------------------------------------
c                                                    tgb
c ********************************************************************
c  **  type1-crop
c  **  type2-short grass                  ****************************
c  **  type3-evergreen needle leaf tree   * soils parameters are a fn of
c  **  type4-deciduous ""  """ "" "" "    * soil colour -- from 1(light)
c  **  type5- "" """"  broadleaf tree     *                  to 8(dark)
c  **  type6- evergreen   """  """"       * soil texture -- from 1(sand)
c  **  type7-tall grass                   *     thru 6(loam) to 12(clay)
c  **  type8-desert                       *
c  **  type9-tundra                       * soil drainage - 7(free),
c  **  type10-irrig crop                  *      8(poor), 9(impeded)
c  **  type11-semi-desert                 *
c  **  type12-ice                         *  (drainage to be used yet)
c  **  type13-bog or marsh                *
c  **  type14-(inland water)              *
c  **  type15-(sea)                       *****************************
c  **  type16-evgr shrub
c  **  type17-decid shrub
c  **  type18-mixed tree
c
c
c
c  flow from this driver and order of subroutines is:
c
c   bconst
c   initb
c  bmarch   ==> soilbc
c    outb       albedo
c               zenith
c                bndry   ==>   drag  ==> dragdn  ==> depth
c               bufout        satur
c                            vcover
c                              drip
c                            lftemp   =====================> stomat
c                              drip                          frawat
c                               co2 ===> carbon             frawat
c                            tseice                           root
c                            tgrund                          satur
c                              snow                         lfdrag
c                             water                         condch
c                                                           condcq
c                                                            deriv
c
c=======================================================================
c  si version  - water fluxes are generally calculated in kg/m**2/s.
c  1000 kg/m**2/s = 1 m/s  - converted to energy units for display:
c                            1 kg/m**2/s = 2.5 x 10.e6  watts/m**2.
c  note also  1 kg/m**2/s = 1 mm/m**2/s so fluxes can be thought of
c                            in mm/s.
c=======================================================================

      implicit none
      integer jslc
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#     include "include/bats.cb"
#     include "include/date.cb"
#     include "include/main.cb"
#     include "include/bats.cb2"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif

      real(kind=8)  rs
      real(kind=8)  ra
      common /ictp01/ rs(nbmax), ra(nbmax)
c---------------------------------------------------------------------
      real(kind=8)  rategb,fc1,rm,rp

c ** note: constants related to the time step are defined here instead
c          of at initialization to allow for the possibility
c          of variable time step
#ifdef MPP1
      if (jslc .le. jbegin) then
#else
      if (jslc .le. 2) then
#endif
         c(4) = dtbat
         c(6) = dtbat * 2.
         c(7) = 1./c(6)
         rategb=0.2/tau1
         fc1=c(4)*rategb*0.5
         rm=1.-fc1
         rp=1.+fc1
      end if

      call interf(1,jslc)

c ****** calculate surface fluxes and hydrology budgets
      call soilbc
      call bndry(iemiss)

c ****************************************************************

c ******    call hostetler lake model hourly at lake points
      if(lakemod.eq.1)then
        if( (jyear.eq.jyear0.and.ktau.eq.0) .or.
     &       mod(ktau+1,klake).eq.0 ) then
          call lakedrv(jslc)
        end if
      end if

c ****************************************************************

      return
      end subroutine vecbats
