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
C
      subroutine zenitm(coszrs,ivmx,jslc)
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C   this subroutine calculates the cosine of the solar zenith angle
C   for all longitude points of the mm42d domain. it needs as inputs
C   the longitude and latitude of the points, the initial date of the
C   simulation and the gmt. all these quantities are specified
C   in the initialization procedure of RegCM
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      implicit none
C
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param3.cb"
#     include "include/main.cb"
#     include "include/date.cb"
C
      integer ivmx , jslc
      real(kind=8) coszrs(ix)
C
      integer ill
      real(kind=8) xt24 , tlocap , omega , xxlat
C
C***********************************************************************
C
      xt24 = dmod(lhour*60.+xtime,1440.D0)
      do ill = 1 , ivmx
        tlocap = xt24/60. + xlong(ill,jslc)/15.
        tlocap = dmod(tlocap+24.,24.D0)
        omega = 15.*(tlocap-12.)*degrad
        xxlat = xlat(ill,jslc)*degrad
C       coszrs = cosine of solar zenith angle
        coszrs(ill) = dsin(declin)*dsin(xxlat) + dcos(declin)
     &                *dcos(xxlat)*dcos(omega)
        coszrs(ill) = dmax1(0.D0,coszrs(ill))
      end do
C
      end subroutine zenitm
