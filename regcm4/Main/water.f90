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
!
      subroutine water
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!              update soil moisture and runoff
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
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      use mod_regcm_param
      use mod_param1 , only : dtbat
      use mod_bats
      use mod_constants , only : drain , tau1 , csoilc , tzero
      implicit none
!
! Local variables
!
      real(8) :: b , bfac , bfac2 , delwat , est0 , evmax , evmxr ,     &
               & evmxt , rap , vakb , wtg2c , xxkb
      real(8) , dimension(nnsg,iym1) :: gwatr , rnof , rsubsr ,         &
           & rsubss , rsubst , rsur , wflux1 , wflux2 , wfluxc , xkmx1 ,&
           & xkmx2 , xkmxr
      integer :: n , i
!
!***********************************************************************
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
      end subroutine water
