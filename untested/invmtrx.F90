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
 
      subroutine invmtrx(a,na,v,nv,n,d,ip,ier,work)
      use mod_message
      implicit none
      integer na , nv , n , ier , info , ip(n)
      real(kind=8) a(n,n) , v(n,n) , work(n) , d(2)
      integer i , j , job
!
!     08/23/91 version 1.0
!     12/10/92 updated to correct bugs
!     note: different from cray routine invmtx
!     uses subroutines sgefa/sgedi from library linpack
!     see dick valent, (SCD, consulting) if problems
!
      if ( n.ne.na .or. n.ne.nv ) call fatal(__FILE__,__LINE__,         &
          &'valent invmtx: equate n, na, nv')
!
      do j = 1 , n
        do i = 1 , n
          v(i,j) = a(i,j)
        end do
      end do
      call sgefa(v,n,n,ip,info)
      if ( info.ne.0 ) then
        write (aline,*) 'sgefa info = ' , info
        call say
        call fatal(__FILE__,__LINE__,'sgefa error')
      end if
      job = 11
      call sgedi(v,n,n,ip,d,work,job)
      ier = info
      end subroutine invmtrx
