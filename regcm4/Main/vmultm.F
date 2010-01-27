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

      subroutine vmultm(a1,a2,a3,nk)
      implicit none
      integer nk
      real(kind=8)  a1(nk,nk),a2(nk,nk),a3(nk,nk)
      integer k,l,mm
      real(kind=8)  c
c
      do 10 l=1,nk
      do 10 k=1,nk
      a1(k,l)=0.
      do 10 mm=1,nk
   10 a1(k,l)=a1(k,l)+a2(k,mm)*a3(mm,l)
      return
c
      entry vaddms(a1,a2,a3,nk)
      do 20 l=1,nk
      do 20 k=1,nk
   20 a1(k,l)=a2(k,l)+a3(k,l)
      return
c
      entry vsubtm(a1,a2,a3,nk)
      do 30 l=1,nk
      do 30 k=1,nk
   30 a1(k,l)=a2(k,l)-a3(k,l)
      return
c
      entry vmultc(a1,a2,nk,c)
      do 40 l=1,nk
      do 40 k=1,nk
   40 a1(k,l)=a2(k,l)*c
      return
      end subroutine vmultm
