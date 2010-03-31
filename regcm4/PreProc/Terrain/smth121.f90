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

      subroutine smth121(htgrid,iy,jx,hscr1)
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx
      real(4) , dimension(iy,jx) :: hscr1 , htgrid
      intent (in) iy , jx
      intent (inout) hscr1 , htgrid
!
! Local variables
!
      integer :: i , j
!
!     PURPOSE :  PERFORMS THE 1-2-1 SMOOTHING TO REMOVE PRIMARILY THE
!     2DX WAVES FROM THE FIELDS htgrid
!
      do j = 1 , jx
        do i = 1 , iy
          hscr1(i,j) = htgrid(i,j)
        end do
      end do
      do i = 1 , iy
        do j = 2 , jx - 1
          if ( (htgrid(i,j)<=-.1) .or. (htgrid(i,j)>0.) ) hscr1(i,j)    &
             & = .25*(2.*htgrid(i,j)+htgrid(i,j+1)+htgrid(i,j-1))
        end do
      end do
      do j = 1 , jx
        do i = 2 , iy - 1
          if ( (hscr1(i,j)<=-.1) .or. (hscr1(i,j)>0.) ) htgrid(i,j)     &
             & = .25*(2.*hscr1(i,j)+hscr1(i+1,j)+hscr1(i-1,j))
        end do
      end do
      end subroutine smth121
