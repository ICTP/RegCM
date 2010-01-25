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

      subroutine vnorml(z,s,nk,nk1)
      implicit none
      integer nk,nk1
      real(kind=8)  z(nk,nk),s(nk1)
      integer k,l,kmax
      real(kind=8)  zmax,v,a
c
c  this routine normalizes the columns of z such that the component
c  with the largest absolute value is positive, and the sum of the
c  mass-weighted squares equals one.
c
      do 3 l=1,nk
      zmax=-1.
      v=0.
c
      do 1 k=1,nk
      a=dabs(z(k,l))
      if (a.gt.zmax) then
         zmax=a
         kmax=k
      end if
      v=(s(k+1)-s(k))*a*a + v
    1 continue
c
      a=(z(kmax,l)/zmax)/dsqrt(v)
      do 2 k=1,nk
    2 z(k,l)=a*z(k,l)
    3 continue
c
      return
      end
