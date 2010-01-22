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

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine satur( qsat, t, p )
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     ****  calculates saturation vapor pressure (eg)
c           and qsat = saturated specific humidity (dimensionless)
c
c           uses tetens formula (1930) (ref. riegel,1974,jam,p606
c                                                 equation 1)
c
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/bats.cb"

      real(kind=8)  a,b
      common/ab/ a(NNSG,nbmax), b(NNSG,nbmax)
      real(kind=8)  qsat(NNSG,nbmax), t(NNSG,nbmax), p(NNSG,nbmax)
      integer n,np

      do 90 np=np1,npts
      do n=1,NNSG
       if(t(n,np)-c(67) .le. 0) then
         goto 1
       else
         goto 2
       end if
1      a(n,np) = c(70)
       b(n,np) = c(71)
       go to 3
2      a(n,np) = c(72)
       b(n,np) = c(73)
3      continue
       eg(n,np)   = c(74) * dexp( a(n,np)*(t(n,np)-c(67)) 
     &                                   /(t(n,np)-b(n,np)) )
       qsat(n,np) = c(75)*eg(n,np) / (p(n,np)-c(76)*eg(n,np))
      enddo
90    continue

      return
      end
