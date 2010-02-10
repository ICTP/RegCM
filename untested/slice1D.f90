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
 
      use regcm_param
      use param3
      use main
      use pbldim
      use mod_bats
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
      integer :: n , np
! 
!     ********* For albedov
!     **********            This is taken from subroutine interf so that
!     **********            radiation can be called in tend (not
!     vecbats).
      npts = nbmax
      do np = np1 , npts
        do n = 1 , nnsg
          ldoc1d(n,np) = ocld2d(n,np,j)
          sice1d(n,np) = sice2d(n,np,j)
          tgb1d(n,np) = tgb2d(n,np,j)
          ssw1d(n,np) = ssw2d(n,np,j)
          lveg(n,np) = nint(veg2d1(n,np,j))
          amxtem = dmax1(298.-tgb1d(n,np),0.D0)
          sfac = 1. - dmax1(0.D0,1.-0.0016*amxtem**2)
          if ( lveg(n,np).eq.0 ) then
            veg1d(n,np) = 0.
          else
            veg1d(n,np) = vegc(lveg(n,np)) - seasf(lveg(n,np))*sfac
          end if
          ts1d(n,np) = thx3d(np,kx,j) - 6.5E-3/g*(ht1(n,np,j)-ht(np,j))
          scv1d(n,np) = scv2d(n,np,j)
          sag1d(n,np) = sag2d(n,np,j)
        end do
      end do
 
      end subroutine slice1D
