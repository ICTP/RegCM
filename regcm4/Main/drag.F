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
      subroutine drag
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     *** determines surface transfer coeffs. at anemom. level from
c     *** lowest model level based on monin-obukov theory using
c     *** deardorff parameterization in terms of bulk richardson no.

c     ****  a.  calculates neutral drag coefficient (cdrn) as a fn of
c     ****             underlying surface

c     ****  b.  modifies cdrn as fn of bulk rich. no. of surface layer

      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/bats.cb"

      real(kind=8)  rib(NNSG,nbmax),ribl(NNSG,nbmax),ribn(NNSG,nbmax)
      real(kind=8)  cdrmin(NNSG,nbmax)
      integer n,np
      real(kind=8)  zatild,dthdz,u1,u2

ccc   nlat  = 40       ! ccm specific code
ccc   nlat2 = nlat/2

c=======================================================================
cl       1.   get neutral drag coefficient
c=======================================================================
      call dragdn

      do 90 np=np1,npts
      do n=1,NNSG
c=======================================================================
cl       2.   compute stability as bulk rich. no. = rin/rid =
cl                    ri(numerator)/ri(denominator)
c=======================================================================
       if (lveg(n,np) .ne. 0) then
         zatild=(z1(n,np)-displa(lveg(n,np)))*sigf(n,np)
     &                              +z1(n,np)*(1.-sigf(n,np))
       else
         zatild = z1(n,np)
       endif
       ribn(n,np)   = zatild * c(54)
     &           * (ts1d(n,np)-sigf(n,np)*taf1d(n,np)-(1.-sigf(n,np))
     &           *  tg1d(n,np)) /  ts1d(n,np)
c=======================================================================
cl       2.1  compute the bulk richardson number;
cl            first get avg winds to use for ri number by summing the
cl            squares of horiz., vertical, and convective velocities
c=======================================================================
       if(ribn(n,np).le.0.) then
         dthdz=(1.-sigf(n,np))*tg1d(n,np)+sigf(n,np)*taf1d(n,np)
     &        -ts1d(n,np)
         u1 = c(90) + 2.*dsqrt(dthdz)
         ribd(n,np) = us1d(np)**2 + vs1d(np)**2 + u1**2
       else
         u2 = c(90)
         ribd(n,np) = us1d(np)**2 + vs1d(np)**2 + u2**2
       endif
       vspda(n,np) = dsqrt(ribd(n,np))
       if (vspda(n,np).lt.1.) then
         vspda(n,np)=1.
         ribd(n,np)=1.
       endif
       rib(n,np) = ribn(n,np) / ribd(n,np)
c=======================================================================
cl       3.   obtain drag coefficient as product of neutral value
cl                     and stability correction
c=======================================================================
c     ****   -0.4 < rib < 0.2   (deardorff, jgr, 1968, 2549-2557)
       if (rib(n,np) .lt. 0.) then
         cdr(n,np) = cdrn(n,np)*(1.0+24.5*dsqrt(-cdrn(n,np)*rib(n,np)))
       else
         cdr(n,np) = cdrn(n,np)/(1.0+11.5*rib(n,np))
       endif
cl       3.1  apply lower limit to drag coefficient value
       cdrmin(n,np) = dmax1(0.25*cdrn(n,np),6.d-4)
       if(cdr(n,np).lt.cdrmin(n,np)) cdr(n,np) = cdrmin(n,np)
       cdrx(n,np) = cdr(n,np)

      enddo
90    continue

c=======================================================================
cl       4.   obtain drag coefficient over sea ice as weighted average
cl                          over ice and leads
cl            warning! the lat test below (4.1-4.3) is model dependent!
c=======================================================================

cl       4.1  test if northern or southern hemisphere
      do 91 np=np1,npts
      do n=1,NNSG
       if(ldoc1d(n,np).gt.1.5) then         ! check each point
        aarea(n,np) = 0.02                      !  4.3  arctic
ccc     if(lat(np).eq.    1) aarea(np) = 0.005   ! ccm specific code
ccc     if(lat(np).eq.    2) aarea(np) = 0.01
ccc     if(lat(np).ge.nlat2) aarea(np) = 0.04     !  4.2  antarctic
       endif
      enddo
91    continue

cl       4.4  neutral cd over lead water
      do 92 np=np1,npts
      do n=1,NNSG
       if(ldoc1d(n,np).gt.1.5) then         !  check each point
        cdrn(n,np) = (vonkar/dlog(z1(n,np)/zoce))**2

cl       4.5  drag coefficient over leads
        ribl(n,np)    = (1.-271.5/ts1d(n,np))*z1(n,np)*c(54)/ribd(n,np)
        if(ribl(n,np).ge.0) then
         clead(n,np)=cdrn(n,np)/(1.+11.5*ribl(n,np))
        else
         clead(n,np)=cdrn(n,np)*(1.+24.5*dsqrt(-cdrn(n,np)*ribl(n,np)))
        endif

cl       4.6  calculate weighted avg of ice and lead drag coefficients
        cdrx(n,np) = (1.-aarea(n,np))*cdr(n,np)+aarea(n,np)*clead(n,np)
       endif
      enddo
92    continue

      return
      end
