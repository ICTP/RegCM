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

      subroutine smthtr(slab1,is1,is2)
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: nocean = 20000
!
! Dummy arguments
!
      integer :: is1 , is2
      real(4) , dimension(is1,is2) :: slab1
      intent (inout) slab1
!
! Local variables
!
      integer :: i , iflg , j , k , n , n1 , npass
      integer , dimension(nocean) :: ii , jj
      character(5) :: point
!
!     smooth terrain arrays
!
      n = 1
      do i = 1 , is1
        do j = 1 , is2
          if ( slab1(i,j)<0.0 ) then
            ii(n) = i
            jj(n) = j
            slab1(i,j) = 0.0
            n = n + 1
          end if
        end do
      end do
      n1 = n - 1
      print 99001 , n1
      if ( n>nocean ) print 99002
      point = 'cross'
      npass = 10
      iflg = 0          ! 0 = smoothing only at boundary
      call smther(slab1,is1,is2,npass,point,iflg)
!     npass = 2000
!     iflg = 0          ! 1 = extensive smoothing over south boundary
!     call smther( slab1, is1, is2, npass, point, iflg )
      do i = 1 , is1
        do j = 1 , is2
          slab1(i,j) = slab1(i,j)
          if ( slab1(i,j)<0.0 ) slab1(i,j) = 0.0
        end do
      end do
      do k = 1 , n - 1
        i = ii(k)
        j = jj(k)
        slab1(i,j) = -0.001
      end do
99001 format (5x,'  there are a total of ',i5,' points of ocean')
99002 format ('0',2x,'dimension exceeded in subr smthtr')
      end subroutine smthtr
