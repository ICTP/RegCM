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

      subroutine maskme(landmask,vals,vmisdat,nlon,nlat,nlev,ntim)
      implicit none
!
! Dummy arguments
!
      integer :: nlat , nlev , nlon , ntim
      real(4) :: vmisdat
      real(4) , dimension(nlon,nlat) :: landmask
      real(4) , dimension(nlon,nlat,nlev,ntim) :: vals
      intent (in) landmask , nlat , nlev , nlon , ntim , vmisdat
      intent (inout) vals
!
! Local variables
!
      integer :: i , j , k , l
!
!     ** Variables Passed in
!     ** Local variables
 
      do l = 1 , ntim
        do k = 1 , nlev
          do j = 1 , nlat
            do i = 1 , nlon
              if ( landmask(i,j)>0.5 ) then
                vals(i,j,k,l) = vals(i,j,k,l)
              else
                vals(i,j,k,l) = vmisdat
              end if
            end do
          end do
        end do
      end do
 
      end subroutine maskme
