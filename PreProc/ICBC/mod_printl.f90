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

      module mod_printl

      contains

      subroutine printl(a,n1,n2)
      implicit none
!
      integer :: n1 , n2
      real(4) , dimension(n1,n2) :: a
      intent (in) a , n1 , n2
!
      integer :: inc1 , inc2 , k , l
!
      inc2 = -1
      inc1 = 1
 
      do l = n2 , 1 , inc2
        print 99001 , l , (a(k,l),k=1,17,inc1)
      end do
99001 format (1x,i3,2x,17F7.2)
 
      end subroutine printl

      end module mod_printl
