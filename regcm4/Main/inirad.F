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

      subroutine inirad

      implicit none
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/rad.cb"

      integer i,j,k
c
c  compute ozone mixing ratio distribution
c
      if (jyear.eq.jyear0.and.ktau.eq.0) then
       do 10 k=1,kx
#ifdef MPP1
       do 10 j=1,jendl
#else
       do 10 j=1,jx-1
#endif
       do 10 i=1,ix-1
        heatrt(i,k,j)=0.
        o3prof(i,k,j)=0.
10     continue
#ifdef MPP1
       do 11 j=1,jendl
#else
       do 11 j=1,jx-1
#endif
       do 11 i=1,ix-1
        o3prof(i,kxp1,j)=0.
11     continue
       call o3data
#ifdef MPP1
       if (myid .eq. 0) then
#endif
       print *,'ozone profiles'
       do 34 k=1,kx+1
        write(6,1212)o3prof(3,k,3)
34     continue
1212   format (1x,7e12.4)
#ifdef MPP1
        end if
#endif
      end if

      return
      end subroutine inirad
