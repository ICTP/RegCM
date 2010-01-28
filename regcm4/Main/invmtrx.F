C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine invmtrx(a, na, v, nv, n, d, ip, ier, work)
      implicit none
#     include "include/message.cb"
      integer na, nv, n, ier, info, ip(n)
      real(kind=8)  a(n,n), v(n,n), work(n), d(2)
      integer i,j,job
c
c 08/23/91 version 1.0
c 12/10/92 updated to correct bugs
c note: different from cray routine invmtx
c uses subroutines sgefa/sgedi from library linpack
c see dick valent, (SCD, consulting) if problems
c
      if (n .ne. na .or. n .ne. nv) then
        call fatal(__FILE__,__LINE__, 'valent invmtx: equate n, na, nv')
      end if
c
      do 10 j=1,n
      do 10 i=1,n
        v(i,j) = a(i,j)
   10 continue
      call sgefa(v,n,n,ip,info)
      if (info .ne. 0) then
        write (aline, *) 'sgefa info = ', info
        call say
        call fatal(__FILE__,__LINE__, 'sgefa error')
      end if
      job = 11
      call sgedi(v,n,n,ip,d,work,job)
      ier = info
      return
      end subroutine invmtrx
