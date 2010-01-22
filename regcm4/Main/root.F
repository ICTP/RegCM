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

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine root
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     this subroutine provides root function in terms of maximum
c     transpiration rate plants can sustain depending on soil moisture.
c
c     trsmx0 is a prescribed constant (kg/m**2/s).
c     trsmx is the maximum transpiration rate,
c        including a low temperature correction (=seasb)
c        and a correction for fractional vegetation (=sigf).
c
c     rotf is ratio of moisture extracton from top to total when
c           fully saturated
c     rootf is ratio of roots in upper soil layer
c                    to roots in root soil layer
c     bsw is the b param in clapp and hornberger
c
c     "wlt  " are ratios factors controlling the saturation
c                 cf wilting (see ewing paper)
c     wlttb (total) & wltub (upper) become 1 at the wilting point
c     (eqn 14 in ewing paper) n.b. etrc=etrmx in ewing paper
c
c     etrc= max poss transpiration given the soil moisture distributions
c     efpr = the relative contribution of upper soil layer to
c     evapotranspiration - need soil moist. budget (subrout water)
c
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/bats.cb"

      real(kind=8)  etrc
      Common /abd/ etrc(NNSG,nbmax)
c
      real(kind=8)  rootf(20)
      data rootf/.30,.80,.67,.67,.50,.80,.80,.90,.90,.30,.80,9*.50/
      integer n,np
      real(kind=8)  trsmx,rotf,bneg,wmli,wlttb,wltub
c
      do 90 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
C        trsmx = trsmx0*sigf(n,np)*seasb(n,np)
         trsmx = trsmx0*sigf(n,np)
         rotf = rootf(lveg(n,np))
         bneg = -bsw(n,np)
         wmli = 1. / (wiltr(n,np)**bneg-1.)
         wlttb = (watr(n,np)**bneg-1.)*wmli
         wltub = (watu(n,np)**bneg-1.)*wmli
         wlttb = dmin1(wlttb,1.d0)
         wltub = dmin1(wltub,1.d0)
         etrc(n,np) = trsmx * (1.-(1.-rotf)*wlttb-rotf*wltub)
         efpr(n,np) = trsmx * rotf * (1.-wltub)
         if(etrc(n,np).lt.1.e-12) then
          etrc(n,np) = 1.e-12
          efpr(n,np) = 1.0
         else
          efpr(n,np) = efpr(n,np) / etrc(n,np)
         endif
        endif
       endif
      enddo
90    continue
c
      return
      end
