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

      subroutine vprntm(a, n1, n2, nam)
      implicit none
      integer n1,n2
      character nam*8
      real(kind=8)  a(n1,n2)
      integer k,l
c
      print 10, nam
      do 1 k=1,n1
    1 print 11,k,(a(k,l),l=1,n2)
   10 format(1h1,a8,/)
   11 format(1x,i3,5x,1p,11g11.3,1x,/,9x,1p,11g11.3)
      return
      end
