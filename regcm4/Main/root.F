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

      subroutine root

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C     this subroutine provides root function in terms of maximum
C     transpiration rate plants can sustain depending on soil moisture.
C
C     trsmx0 is a prescribed constant (kg/m**2/s).
C     trsmx is the maximum transpiration rate,
C        including a low temperature correction (=seasb)
C        and a correction for fractional vegetation (=sigf).
C
C     rotf is ratio of moisture extracton from top to total when
C           fully saturated
C     rootf is ratio of roots in upper soil layer
C                    to roots in root soil layer
C     bsw is the b param in clapp and hornberger
C
C     "wlt  " are ratios factors controlling the saturation
C                 cf wilting (see ewing paper)
C     wlttb (total) & wltub (upper) become 1 at the wilting point
C     (eqn 14 in ewing paper) n.b. etrc=etrmx in ewing paper
C
C     etrc= max poss transpiration given the soil moisture distributions
C     efpr = the relative contribution of upper soil layer to
C     evapotranspiration - need soil moist. budget (subrout water)
C
      implicit none
#     include "parameter.inc"
#     include "include/bats.cb"

      real(kind=8)  etrc
      Common /abd/ etrc(NNSG,nbmax)
C
      real(kind=8)  rootf(20)
      data rootf/.30,.80,.67,.67,.50,.80,.80,.90,.90,.30,.80,9*.50/

      real(kind=8)  trsmx,rotf,bneg,wmli,wlttb,wltub
      integer n, np
C
      do np = np1, npts
        do n = 1, NNSG
          if (ldoc1d(n,np).gt.0.5) then
            if (sigf(n,np).gt.0.001) then
C             trsmx = trsmx0*sigf(n,np)*seasb(n,np)
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
              if (etrc(n,np).lt.1.e-12) then
                etrc(n,np) = 1.e-12
                efpr(n,np) = 1.0
              else
                efpr(n,np) = efpr(n,np) / etrc(n,np)
              end if
            end if
          end if
        end do
      end do
C
      return
      end subroutine root
