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

      subroutine top2btm(x,nlon1,nlat1,nlev1)
      implicit none
!
! Dummy arguments
!
      integer :: nlat1 , nlev1 , nlon1
      real , dimension(nlon1,nlat1,nlev1) :: x
      intent (in) nlat1 , nlev1 , nlon1
      intent (inout) x
!
! Local variables
!
      integer :: i , j , k , kr
      real , dimension(nlev1) :: work
!
      do i = 1 , nlon1
        do j = 1 , nlat1
          do k = 1 , nlev1
            work(k) = x(i,j,k)
          end do
          do k = 1 , nlev1
            kr = nlev1 - k + 1
            x(i,j,k) = work(kr)
          end do
        end do
      end do
      end subroutine top2btm
