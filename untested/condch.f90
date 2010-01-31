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
 
      subroutine condch

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     dimensional and non-dimensional sensible heat conductances
!               for canopy and soil flux calculations
!
      use regcm_param
      use bats
      implicit none
!
! Local variables
!
      integer :: n , np
!
!     csoilc = constant drag coefficient for soil under canopy;
!     prescribed in subroutine bconst.
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
              uaf(n,np) = vspda(n,np)*dsqrt(cdr(n,np))
              cf(n,np) = 0.01*sqrtdi(lveg(n,np))/dsqrt(uaf(n,np))
              wta(n,np) = sigf(n,np)*cdr(n,np)*vspda(n,np)
              wtlh(n,np) = cf(n,np)*uaf(n,np)*vegt(n,np)
              wtg(n,np) = csoilc*uaf(n,np)*sigf(n,np)
              wtshi(n,np) = 1./(wta(n,np)+wtlh(n,np)+wtg(n,np))
              wtl0(n,np) = wtlh(n,np)*wtshi(n,np)
              wtg0(n,np) = wtg(n,np)*wtshi(n,np)
              wtgl(n,np) = wtl0(n,np) + wtg0(n,np)
              wta0(n,np) = 1. - wtgl(n,np)
              wtga(n,np) = wta0(n,np) + wtg0(n,np)
            end if
          end if
        end do
      end do
!
      end subroutine condch
