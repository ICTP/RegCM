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


c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine deriv
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     derivatives of energy fluxes with respect to leaf temperature for
c     newton-raphson calculation of leaf temperature.
c     input: rs,ra,cdrd,rppq,efe.    output: qsatld,dcd.
c
c     approximate by derivatives of cdr and ef.  many weaker
c     dependences on leaf temperature are omitted, as convergence
c     rate is not affected.
c
      implicit none
#     include "parameter.inc"
#     include "include/bats.cb"

      real(kind=8)  a,b
      common/ab/ a(NNSG,nbmax), b(NNSG,nbmax)
      real(kind=8)  cdrd
      Common/aba/ cdrd(NNSG,nbmax)
      real(kind=8)  rppq, efe, qsatld, dcd
      Common/abc/ rppq(NNSG,nbmax), efe(NNSG,nbmax)
     &         ,qsatld(NNSG,nbmax), dcd(NNSG,nbmax)
      real(kind=8)  rs
      real(kind=8)  ra
      common /ictp01/ rs(NNSG,nbmax), ra(NNSG,nbmax)
      integer n,np
      real(kind=8)  dne,xkb,hfl

      do 90 np =1, npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
         dne  = 1. / (tlef1d(n,np)-b(n,np))
         qsatld(n,np) = qsatl(n,np)*a(n,np)*(c(67)-b(n,np))*dne**2
         xkb = cdrd(n,np) / cdr(n,np)
         hfl = df(n,np) * ( wtga(n,np)*tlef1d(n,np)
     a               - wtg0(n,np)*tg1d(n,np) - wta0(n,np)*ts1d(n,np) )
         dcd(n,np) = cn1(n,np)*rppq(n,np)*wtgaq(n,np)*qsatld(n,np)
     a             + (1.-wtgaq(n,np))*efe(n,np)*xkb
     b             + (1.- wtga(n,np))*hfl    *xkb
         dcd(n,np) = dmax1(dcd(n,np),0.d0)
         dcd(n,np) = dmin1(dcd(n,np),500.d0)
        end if
       end if
      enddo
90    continue

      return
      end
