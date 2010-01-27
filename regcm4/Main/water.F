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
C
      subroutine water
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C              update soil moisture and runoff
C
C     new algorithms for three soil layers (dickinson & kennedy 8-88)
C     calculate fluxes through air, surface layer, and root layer faces
C
C          b = b of clapp and hornberger
C       est0 = soil water flux, out top
C      gwatr = net input of water to the soil surface
C       ircp = leaf interception
C     wflux1 = soil water flux, 10 cm
C     wflux2 = soil water flux, 1 m
C     rsubss = soil water flux by grav. drainage thru 10 cm interface
C     rsubsr = soil water flux by grav. drainage thru 1 m interface
C     rsubst = soil water flux by grav. drainage thru 10 m interface
C       rsur = surface runoff
C     rno1d(n,np) = total runoff (mm/day)
C     rnos1d(n,np) = surface runoff (mm/day)
C
C     xkmxr and wflux1 determine flow thru upper/root soil interface
C     evmxt, xkmx1, and xkmx2 determine flow thru lower interfaces
C
C     veg type 10 "irrigated crop" is irrigated through reducing
C          the runoff (rsur), i.e., by adding a negative number
C          if the land isn't at least 60% saturated.
C     veg type 13 and 14 are water covered (lake, swamp, rice paddy);
C          negative runoff keeps this land saturated.
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      implicit none
C
#     include "parameter.inc"
#     include "include/bats.cb"
C
      real(kind=8) rsur(nnsg,nbmax)
      real(kind=8) rnof(nnsg,nbmax)
      real(kind=8) gwatr(nnsg,nbmax)
      real(kind=8) wfluxc(nnsg,nbmax)
      real(kind=8) wflux1(nnsg,nbmax)
      real(kind=8) wflux2(nnsg,nbmax)
      real(kind=8) rsubsr(nnsg,nbmax)
      real(kind=8) rsubss(nnsg,nbmax)
      real(kind=8) rsubst(nnsg,nbmax)
      real(kind=8) xkmxr(nnsg,nbmax)
      real(kind=8) xkmx1(nnsg,nbmax)
      real(kind=8) xkmx2(nnsg,nbmax)
C
      integer n , np
      real(kind=8) evmxr , evmxt , b , bfac , bfac2 , delwat , xxkb , 
     &             vakb , wtg2c
      real(kind=8) rap , est0 , evmax
C
C***********************************************************************
C

C=======================================================================
C     1.   define soil water fluxes
C=======================================================================
C
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np)>0.5 .and. ldoc1d(n,np)<1.5 ) then
C
C           1.1  reduce infiltration for frozen ground
C
            if ( tgb1d(n,np)>c(67) ) then
              xkmxr(n,np) = xkmx(n,np)
            else
              xkmxr(n,np) = 0.
            end if
C
C           1.11 permafrost or ice sheet
C
            if ( lveg(n,np)==9 .or. lveg(n,np)==12 ) then
              xkmx1(n,np) = 0.
              xkmx2(n,np) = 0.
            else
              xkmx1(n,np) = xkmx(n,np)
              xkmx2(n,np) = drain
            end if
C
C           1.2  diffusive fluxes
C
            evmxr = evmx0(n,np)*xkmxr(n,np)/xkmx(n,np)
            evmxt = evmx0(n,np)*xkmx1(n,np)/xkmx(n,np)
            b = bsw(n,np)
            bfac = watr(n,np)**(3.+bfc(n,np))*watu(n,np)
     &             **(b-bfc(n,np)-1)
            bfac2 = watt(n,np)**(2.+bfc(n,np))*watr(n,np)**(b-bfc(n,np))
            wfluxc(n,np) = evmxr*(depuv(lveg(n,np))/deprv(lveg(n,np)))
     &                     **0.4*bfac
            wflux1(n,np) = wfluxc(n,np)*watr(n,np)
            wflux2(n,np) = evmxt*(depuv(lveg(n,np))/deprv(lveg(n,np)))
     &                     **0.5*bfac2*(watt(n,np)-watr(n,np))
C
C           1.3  gravitational drainage
C
            rsubss(n,np) = xkmxr(n,np)*watr(n,np)**(b+0.5)*watu(n,np)
     &                     **(b+2.5)
            rsubsr(n,np) = xkmx1(n,np)*watt(n,np)**(b+0.5)*watr(n,np)
     &                     **(b+2.5)
            rsubst(n,np) = xkmx2(n,np)*watt(n,np)**(2.*b+3.)
C
C           1.32 bog or water
C
            if ( (lveg(n,np)>=13) .and. (lveg(n,np)<=15) ) then
              rsubst(n,np) = 0.0
              rsubss(n,np) = 0.0
              rsubsr(n,np) = 0.0
            end if
C
C           1.4  fluxes through internal surfaces
C
            wflux1(n,np) = wflux1(n,np) - rsubss(n,np)
            wflux2(n,np) = wflux2(n,np) - rsubsr(n,np)
          end if
        end do
      end do
C
C     1.5  net flux at air interface
C
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np)>0.5 .and. ldoc1d(n,np)<1.5 ) then
            gwatr(n,np) = pw(n,np) - evapw(n,np) + sm(n,np)
     &                    + etrrun(n,np)/c(4)
C
C=======================================================================
C           2.   define runoff terms
C=======================================================================
C
C           2.1  surface runoff
C
            wata(n,np) = 0.5*(watu(n,np)+watr(n,np))
C
C           2.11 increase surface runoff over frozen ground
C
            if ( tg1d(n,np)<c(67) ) then
              rsur(n,np) = dmin1(1.D0,wata(n,np)**1)
     &                     *dmax1(0.D0,gwatr(n,np))
            else
              rsur(n,np) = dmin1(1.D0,wata(n,np)**4)
     &                     *dmax1(0.D0,gwatr(n,np))
            end if
C
C           2.12 irrigate cropland
C
            if ( lveg(n,np)==10 .and. watr(n,np)<relfc(n,np) )
     &           rsur(n,np) = rsur(n,np)
     &                        + (rsw1d(n,np)-relfc(n,np)*gwmx1(n,np))
     &                        /c(4)
C
C           2.13 saturate swamp or rice paddy
C
            if ( (lveg(n,np)>=13) .and. (lveg(n,np)<16) ) rsur(n,np)
     &           = rsur(n,np)
     &           + dmin1(0.D0,(rsw1d(n,np)-gwmx1(n,np))/c(4))
C
C           2.2  total runoff
C
            rnof(n,np) = rsur(n,np) + rsubst(n,np)
C
C=======================================================================
C           3.   increment soil moisture
C=======================================================================
C
C           3.1  update top layer with implicit treatment
C           of flux from below
C
            ssw1d(n,np) = ssw1d(n,np) + c(4)
     &                    *(gwatr(n,np)-efpr(n,np)*etr(n,np)-rsur(n,np)
     &                    +wflux1(n,np))
            ssw1d(n,np) = ssw1d(n,np)/(1.+wfluxc(n,np)*c(4)/gwmx0(n,np))
C
C           3.2  update root zone
C
            rsw1d(n,np) = rsw1d(n,np) + c(4)
     &                    *(gwatr(n,np)-etr(n,np)-rsur(n,np)
     &                    +wflux2(n,np))
C
C           3.3  update total water
C
            tsw1d(n,np) = tsw1d(n,np) + c(4)
     &                    *(gwatr(n,np)-etr(n,np)-rnof(n,np))
          end if
        end do
      end do
C
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np)>0.5 .and. ldoc1d(n,np)<1.5 ) then
C
C=======================================================================
C           4.   check whether soil water exceeds maximum capacity or
C           becomes negative (should rarely or never happen)
C=======================================================================
C
C           4.1  surface water assumed to move downward into soil
C
            if ( ssw1d(n,np)>gwmx0(n,np) ) ssw1d(n,np) = gwmx0(n,np)
C
C           4.2  excess root layer water assumed to move downward
C
            if ( rsw1d(n,np)>gwmx1(n,np) ) rsw1d(n,np) = gwmx1(n,np)
C
C           4.3  excess total water assumed to go to subsurface runoff
C
            if ( tsw1d(n,np)>gwmx2(n,np) ) then
              delwat = tsw1d(n,np) - gwmx2(n,np)
              tsw1d(n,np) = gwmx2(n,np)
              rsubst(n,np) = rsubst(n,np) + delwat/c(4)
            end if
C
C           4.4  check for negative water in top layer
C
            if ( ssw1d(n,np)<=1.E-2 ) ssw1d(n,np) = 1.E-2
C
C=======================================================================
C           5.   accumulate leaf interception
C=======================================================================
C
            ircp1d(n,np) = ircp1d(n,np) + sigf(n,np)*(c(4)*prcp1d(n,np))
     &                     - (sdrop(n,np)+etrrun(n,np))
C
C=======================================================================
C           6.   evaluate runoff (incremented in ccm)
C=======================================================================
C
C*          update total runoff
C
            rnof(n,np) = rsur(n,np) + rsubst(n,np)
            rno1d(n,np) = rnof(n,np)*tau1
            rnos1d(n,np) = rsur(n,np)*tau1
          else                       ! ocean or sea ice
            rnof(n,np) = 0.
            rno1d(n,np) = 0.
            rnos1d(n,np) = 0.
          end if
        end do
      end do
C
C=======================================================================
C     7.   calculate potential evaporation and use to determine
C     wetness factor, allowing for snow being saturated
C=======================================================================
C
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np)>0.5 .and. ldoc1d(n,np)<1.5 ) then
            xxkb = dmin1(rough(lveg(n,np)),1.D0)
            vakb = (1.-sigf(n,np))*vspda(n,np) + sigf(n,np)
     &             *(xxkb*uaf(n,np)+(1.-xxkb)*vspda(n,np))
            wtg2c = (1.-sigf(n,np))*cdrx(n,np)*vakb
            rap = rhs1d(n,np)
     &            *(csoilc*uaf(n,np)*sigf(n,np)*(qg1d(n,np)+delq1d(n,np)
     &            -qs1d(n,np))+wtg2c*(qg1d(n,np)-qs1d(n,np)))
            bfac = watr(n,np)**(3.+bfc(n,np))*watu(n,np)
     &             **(bsw(n,np)-bfc(n,np)-1)
            est0 = evmx0(n,np)*bfac*watu(n,np)
            evmax = dmax1(est0,0.D0)
            gwet1d(n,np) = dmin1(1.D0,evmax/dmax1(1.D-14,rap))
            gwet1d(n,np) = scvk(n,np) + gwet1d(n,np)*(1.0-scvk(n,np))
          end if
        end do
      end do
C
      end subroutine water
