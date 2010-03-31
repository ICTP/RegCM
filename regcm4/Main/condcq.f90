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
      use mod_regcm_param
      use mod_bats , only : wtlq , gwet1d , wtlq , wtlh , rpp , sigf ,  &
                   & ldoc1d , wtgaq , wtaq0 , wtglq , wtlq0 , wtgq0 ,   &
                   & wtsqi , wtgq , rgr , wtg , wta
      implicit none
!
! Local variables
!
      integer :: n , i
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
              rgr(n,i) = gwet1d(n,i)
              wtlq(n,i) = wtlh(n,i)*rpp(n,i)
              wtgq(n,i) = wtg(n,i)*rgr(n,i)
              wtsqi(n,i) = 1./(wta(n,i)+wtlq(n,i)+wtgq(n,i))
              wtgq0(n,i) = wtgq(n,i)*wtsqi(n,i)
              wtlq0(n,i) = wtlq(n,i)*wtsqi(n,i)
              wtglq(n,i) = wtgq0(n,i) + wtlq0(n,i)
              wtaq0(n,i) = 1. - wtglq(n,i)
              wtgaq(n,i) = wtaq0(n,i) + wtgq0(n,i)
            end if
          end if
        end do
      end do
 
      end subroutine condcq
