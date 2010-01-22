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

      subroutine zenitm(coszrs,ivmx,jslc)
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c   this subroutine calculates the cosine of the solar zenith angle
c   for all longitude points of the mm42d domain. it needs as inputs
c   the longitude and latitude of the points, the initial date of the
c   simulation and the gmt. all these quantities are specified
c   in the initialization procedure of mm42d
c
      implicit none
      integer ivmx,jslc
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/param1.cb"
#     include "include/param3.cb"
#     include "include/main.cb"
#     include "include/date.cb"
c
      real(kind=8)  coszrs(ix)
      integer ill
      real(kind=8)  xt24,tlocap,omega,xxlat
c
      xt24 = dmod(lhour*60.+xtime,1440.d0)
      do 10 ill = 1,ivmx
       tlocap = xt24/60. + xlong(ill,jslc)/15.
       tlocap = dmod(tlocap+24., 24.d0)
       omega  = 15.*(tlocap-12.)*degrad
       xxlat  = xlat(ill,jslc)*degrad
c      coszrs = cosine of solar zenith angle
       coszrs(ill)=dsin(declin)*dsin(xxlat) +
     1            dcos(declin)*dcos(xxlat)*dcos(omega)
       coszrs(ill) = dmax1(0.d0,coszrs(ill))
10    continue

      return
      end
