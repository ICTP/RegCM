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

        subroutine minimi(array,ix,kx,ks,kend,kt,istart,iend)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer ix,kx,kend,istart,iend
      real(kind=8)  array(ix,kx)
      integer kt(ix),ks(ix)
      integer i,k
      real(kind=8)  x
c
      do 200 i=istart,iend
       kt(i)=ks(i)
       x=array(i,ks(i))
       do 100 k=ks(i)+1,kend
        if(array(i,k).lt.x)then
         x=array(i,k)
         kt(i)=k
        end if
100    continue
200   continue
c
      return
      end subroutine minimi
