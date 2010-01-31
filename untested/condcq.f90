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
 
      subroutine condcq

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     dimensional and non-dimensional latent heat conductances
!               for canopy and soil flux calculations
!
!     latent fluxes differ from sensible due to stomatal resistance
!
      use regcm_param
      use bats
      implicit none
!
! Local variables
!
      integer :: n , np
!
!     symbols used for weights are:   wt : weight
!     a : air
!     l : leaf
!     i : inverse
!     s : sum
!     h : sensible heat
!     q : water vapor
!     0 : normalized (sums to one)
!     g : ground
!
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) then
              rgr(n,np) = gwet1d(n,np)
              wtlq(n,np) = wtlh(n,np)*rpp(n,np)
              wtgq(n,np) = wtg(n,np)*rgr(n,np)
              wtsqi(n,np) = 1./(wta(n,np)+wtlq(n,np)+wtgq(n,np))
              wtgq0(n,np) = wtgq(n,np)*wtsqi(n,np)
              wtlq0(n,np) = wtlq(n,np)*wtsqi(n,np)
              wtglq(n,np) = wtgq0(n,np) + wtlq0(n,np)
              wtaq0(n,np) = 1. - wtglq(n,np)
              wtgaq(n,np) = wtaq0(n,np) + wtgq0(n,np)
            end if
          end if
        end do
      end do
 
      end subroutine condcq
