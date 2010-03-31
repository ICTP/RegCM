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
 
      subroutine root
 
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     this subroutine provides root function in terms of maximum
!     transpiration rate plants can sustain depending on soil moisture.
!
!     trsmx0 is a prescribed constant (kg/m**2/s).
!     trsmx is the maximum transpiration rate,
!        including a low temperature correction (=seasb)
!        and a correction for fractional vegetation (=sigf).
!
!     rotf is ratio of moisture extracton from top to total when
!           fully saturated
!     rootf is ratio of roots in upper soil layer
!                    to roots in root soil layer
!     bsw is the b param in clapp and hornberger
!
!     "wlt  " are ratios factors controlling the saturation
!                 cf wilting (see ewing paper)
!     wlttb (total) & wltub (upper) become 1 at the wilting point
!     (eqn 14 in ewing paper) n.b. etrc=etrmx in ewing paper
!
!     etrc= max poss transpiration given the soil moisture distributions
!     efpr = the relative contribution of upper soil layer to
!     evapotranspiration - need soil moist. budget (subrout water)
!
      use mod_regcm_param
      use mod_bats
      use mod_ictp01
      use mod_constants , only : trsmx0
      implicit none
!
! Local variables
!
      real(8) :: bneg , rotf , trsmx , wlttb , wltub , wmli
      integer :: n , i
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
!             trsmx = trsmx0*sigf(n,i)*seasb(n,i)
              trsmx = trsmx0*sigf(n,i)
              rotf = rootf(lveg(n,i))
              bneg = -bsw(n,i)
              wmli = 1./(wiltr(n,i)**bneg-1.)
              wlttb = (watr(n,i)**bneg-1.)*wmli
              wltub = (watu(n,i)**bneg-1.)*wmli
              wlttb = dmin1(wlttb,1.D0)
              wltub = dmin1(wltub,1.D0)
              etrc(n,i) = trsmx*(1.-(1.-rotf)*wlttb-rotf*wltub)
              efpr(n,i) = trsmx*rotf*(1.-wltub)
              if ( etrc(n,i).lt.1.E-12 ) then
                etrc(n,i) = 1.E-12
                efpr(n,i) = 1.0
              else
                efpr(n,i) = efpr(n,i)/etrc(n,i)
              end if
            end if
          end if
        end do
      end do
!
      end subroutine root
