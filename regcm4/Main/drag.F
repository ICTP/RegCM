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

      subroutine drag

C     *** determines surface transfer coeffs. at anemom. level from
C     *** lowest model level based on monin-obukov theory using
C     *** deardorff parameterization in terms of bulk richardson no.

C     ****  a.  calculates neutral drag coefficient (cdrn) as a fn of
C     ****             underlying surface

C     ****  b.  modifies cdrn as fn of bulk rich. no. of surface layer

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
      real(kind=8)  zatild,dthdz,u1,u2

      integer n,np

Ccc   nlat  = 40       ! ccm specific code
Ccc   nlat2 = nlat/2

C=======================================================================
C        1.   get neutral drag coefficient
C=======================================================================
      call dragdn

      do np = np1, npts
        do n = 1, NNSG
C=======================================================================
C        2.   compute stability as bulk rich. no. = rin/rid =
C                     ri(numerator)/ri(denominator)
C=======================================================================
         if (lveg(n,np) .ne. 0) then
           zatild=(z1(n,np)-displa(lveg(n,np)))*sigf(n,np)
     &                     +z1(n,np)*(1.-sigf(n,np))
         else
           zatild = z1(n,np)
         end if
         ribn(n,np)   = zatild * c(54)
     &           * (ts1d(n,np)-sigf(n,np)*taf1d(n,np)-(1.-sigf(n,np))
     &           *  tg1d(n,np)) /  ts1d(n,np)
C=======================================================================
C        2.1  compute the bulk richardson number;
C             first get avg winds to use for ri number by summing the
C             squares of horiz., vertical, and convective velocities
C=======================================================================
         if (ribn(n,np).le.0.) then
           dthdz=(1.-sigf(n,np))*tg1d(n,np)+sigf(n,np)*taf1d(n,np)
     &        -ts1d(n,np)
           u1 = c(90) + 2.*dsqrt(dthdz)
           ribd(n,np) = us1d(np)**2 + vs1d(np)**2 + u1**2
         else
           u2 = c(90)
           ribd(n,np) = us1d(np)**2 + vs1d(np)**2 + u2**2
         end if
         vspda(n,np) = dsqrt(ribd(n,np))
         if (vspda(n,np).lt.1.) then
           vspda(n,np)=1.
           ribd(n,np)=1.
         end if
         rib(n,np) = ribn(n,np) / ribd(n,np)
C=======================================================================
C        3.   obtain drag coefficient as product of neutral value
C                      and stability correction
C=======================================================================
C     ****   -0.4 < rib < 0.2   (deardorff, jgr, 1968, 2549-2557)
         if (rib(n,np) .lt. 0.) then
           cdr(n,np) = cdrn(n,np)*(1.0+24.5*
     $                 dsqrt(-cdrn(n,np)*rib(n,np)))
         else
           cdr(n,np) = cdrn(n,np)/(1.0+11.5*rib(n,np))
         end if
C        3.1  apply lower limit to drag coefficient value
          cdrmin(n,np) = dmax1(0.25*cdrn(n,np),6.d-4)
          if (cdr(n,np).lt.cdrmin(n,np)) cdr(n,np) = cdrmin(n,np)
          cdrx(n,np) = cdr(n,np)

        end do
      end do

C=======================================================================
C        4.   obtain drag coefficient over sea ice as weighted average
C                           over ice and leads
C             warning! the lat test below (4.1-4.3) is model dependent!
C=======================================================================

C        4.1  test if northern or southern hemisphere
      do np = np1, npts
        do n = 1, NNSG
          if (ldoc1d(n,np).gt.1.5) then             ! check each point
            aarea(n,np) = 0.02                      !  4.3  arctic
Ccc         if(lat(np).eq.    1) aarea(np) = 0.005  ! ccm specific code
Ccc         if(lat(np).eq.    2) aarea(np) = 0.01
Ccc         if(lat(np).ge.nlat2) aarea(np) = 0.04   !  4.2  antarctic
          end if
        end do
      end do

C        4.4  neutral cd over lead water
      do np = np1, npts
        do n = 1, NNSG
          if (ldoc1d(n,np).gt.1.5) then         !  check each point
            cdrn(n,np) = (vonkar/dlog(z1(n,np)/zoce))**2

C        4.5  drag coefficient over leads
            ribl(n,np) = (1.-271.5/ts1d(n,np))*z1(n,np)*c(54)/ribd(n,np)
            if (ribl(n,np).ge.0) then
              clead(n,np)=cdrn(n,np)/(1.+11.5*ribl(n,np))
            else
              clead(n,np) = cdrn(n,np)*(1.+24.5*
     $                    dsqrt(-cdrn(n,np)*ribl(n,np)))
            end if

C        4.6  calculate weighted avg of ice and lead drag coefficients
            cdrx(n,np) = (1.-aarea(n,np))*cdr(n,np)+
     $                   aarea(n,np)*clead(n,np)
          end if
        end do
      end do

      return
      end subroutine drag
