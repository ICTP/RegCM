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
      use mod_regcm_param
      use mod_bats , only : vegt , vspda , cdr , lveg , sigf , uaf ,    &
                    & cf , wta , wtlh , wtg , wtshi , wtl0 , wtg0 ,     &
                    & wtgl , wta0 , wtga , ldoc1d , sqrtdi
      use mod_constants , only : csoilc
      implicit none
!
! Local variables
!
      integer :: n , i
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
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              uaf(n,i) = vspda(n,i)*dsqrt(cdr(n,i))
              cf(n,i) = 0.01*sqrtdi(lveg(n,i))/dsqrt(uaf(n,i))
              wta(n,i) = sigf(n,i)*cdr(n,i)*vspda(n,i)
              wtlh(n,i) = cf(n,i)*uaf(n,i)*vegt(n,i)
              wtg(n,i) = csoilc*uaf(n,i)*sigf(n,i)
              wtshi(n,i) = 1./(wta(n,i)+wtlh(n,i)+wtg(n,i))
              wtl0(n,i) = wtlh(n,i)*wtshi(n,i)
              wtg0(n,i) = wtg(n,i)*wtshi(n,i)
              wtgl(n,i) = wtl0(n,i) + wtg0(n,i)
              wta0(n,i) = 1. - wtgl(n,i)
              wtga(n,i) = wta0(n,i) + wtg0(n,i)
            end if
          end if
        end do
      end do
!
      end subroutine condch
