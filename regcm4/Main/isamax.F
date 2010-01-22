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

      integer function isamax(n,sx,incx)
      implicit none
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      integer incx,n
      real(kind=8)  sx(1+(n-1)*incx),smax
      integer i,ix
c
      isamax = 0
      if( n .lt. 1 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = dabs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = dabs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = dabs(sx(1))
      do 30 i = 2,n
         if(dabs(sx(i)).le.smax) go to 30
         isamax = i
         smax = dabs(sx(i))
   30 continue
      return
      end
