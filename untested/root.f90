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
      integer :: n , np
      real(8) , dimension(20) :: rootf
!
      data rootf/.30 , .80 , .67 , .67 , .50 , .80 , .80 , .90 , .90 ,  &
         & .30 , .80 , 9*.50/
!
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) then
!             trsmx = trsmx0*sigf(n,np)*seasb(n,np)
              trsmx = trsmx0*sigf(n,np)
              rotf = rootf(lveg(n,np))
              bneg = -bsw(n,np)
              wmli = 1./(wiltr(n,np)**bneg-1.)
              wlttb = (watr(n,np)**bneg-1.)*wmli
              wltub = (watu(n,np)**bneg-1.)*wmli
              wlttb = dmin1(wlttb,1.D0)
              wltub = dmin1(wltub,1.D0)
              etrc(n,np) = trsmx*(1.-(1.-rotf)*wlttb-rotf*wltub)
              efpr(n,np) = trsmx*rotf*(1.-wltub)
              if ( etrc(n,np).lt.1.E-12 ) then
                etrc(n,np) = 1.E-12
                efpr(n,np) = 1.0
              else
                efpr(n,np) = efpr(n,np)/etrc(n,np)
              end if
            end if
          end if
        end do
      end do
!
      end subroutine root
