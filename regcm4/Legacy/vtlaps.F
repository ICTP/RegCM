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

      subroutine vtlaps(t,sigma,r,pt,pd,nk)
      implicit none
      integer nk
      real(kind=8)  r,pt,pd
      real(kind=8)  rate,g,tstrat,zstrat,t0,p0
      parameter (rate=.0065, g=9.8, tstrat=218.15, zstrat=10769
     a,          t0=288.15, p0=101.325)
      real(kind=8)  t(nk),sigma(nk)
      real(kind=8)  fac,p,z
      integer k
c
c  this routine computes the temperature corresponding to a u. s.
c  standard atmosphere (see text by hess). units of p are cb.
c
      fac=r*rate/g
      do 1 k=1,nk
      p=sigma(k)*pd+pt
      t(k)=t0*((p/p0)**fac)
      z=(t0-t(k))/rate
      if (z.gt.zstrat) t(k)=tstrat
    1 continue
c
      return
      end subroutine vtlaps
