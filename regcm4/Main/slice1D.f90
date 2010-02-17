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
 
      subroutine slice1D(j)
 
      use mod_regcm_param
      use mod_main
      use mod_pbldim
      use mod_bats
      use mod_constants , only : rgti
      implicit none
!
! Dummy arguments
!
      integer :: j
      intent (in) j
!
! Local variables
!
      real(8) :: amxtem , sfac
      integer :: n , i
! 
!     ********* For albedov
!     **********            This is taken from subroutine interf so that
!     **********            radiation can be called in tend (not
!     vecbats).
      do i = 2 , ixm1
        do n = 1 , nnsg
          ldoc1d(n,i) = ocld2d(n,i,j)
          sice1d(n,i) = sice2d(n,i,j)
          tgb1d(n,i) = tgb2d(n,i,j)
          ssw1d(n,i) = ssw2d(n,i,j)
          lveg(n,i) = nint(veg2d1(n,i,j))
          amxtem = dmax1(298.-tgb1d(n,i),0.D0)
          sfac = 1. - dmax1(0.D0,1.-0.0016*amxtem**2)
          if ( lveg(n,i).eq.0 ) then
            veg1d(n,i) = 0.
          else
            veg1d(n,i) = vegc(lveg(n,i)) - seasf(lveg(n,i))*sfac
          end if
          ts1d(n,i) = thx3d(i,kx,j)-6.5E-3*rgti*(ht1(n,i,j)-ht(i,j))
          scv1d(n,i) = scv2d(n,i,j)
          sag1d(n,i) = sag2d(n,i,j)
        end do
      end do
 
      end subroutine slice1D
