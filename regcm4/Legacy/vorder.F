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

      subroutine vorder(z,hbar,wz,wh,nk)
      implicit none
      integer nk
      real(kind=8)  z(nk,nk),hbar(nk),wz(nk,nk),wh(nk,2)
      integer k,l,kmax
      real(kind=8)  hmax
c
c  this routine orders the components of hbar so they are largest to
c  smallest valued.  the columns of z are reorded so that they
c  correspond to the same (but reordered) components of hbar.
c
      do 1 k=1,nk
      wh(k,1)=hbar(k)
      wh(k,2)=0.
      do 1 l=1,nk
    1 wz(k,l)=z(k,l)
c
      do 4 l=1,nk
      hmax=-1.d100
      do 2 k=1,nk
      if ((wh(k,2).eq.0.).and.(wh(k,1).gt.hmax)) then
         hmax=wh(k,1)
         kmax=k
      end if
    2 continue
c
      hbar(l)=hmax
      wh(kmax,2)=1.
      do 3 k=1,nk
    3 z(k,l)=wz(k,kmax)
    4 continue
c
      return
      end subroutine vorder
