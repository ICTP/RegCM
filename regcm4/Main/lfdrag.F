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

      subroutine lfdrag
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     recalculate stability dependent drag coefficient for vegetation,
c     given the neutral drag coefficient.
c
      implicit none
#     include "parameter.inc"
#     include "include/bats.cb"

      real(kind=8)  a,b
      common/ab/ a(nbmax), b(nbmax)
      real(kind=8)  cdrd
      Common/aba/ cdrd(NNSG,nbmax)
      real(kind=8)  dlstaf(NNSG,nbmax), rib1(NNSG,nbmax)
      real(kind=8)  cdrmin(NNSG,nbmax), rib(NNSG,nbmax)
      integer n,np
      real(kind=8)  tkb,dthdz,u1,u2,zatild,sqrtf,ribi

      do 90 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
         tkb    = wta0(n,np)*ts1d(n,np) + wtl0(n,np)*tlef1d(n,np)
     a                             + wtg0(n,np)*tg1d(n,np)
         dlstaf(n,np)=ts1d(n,np)-sigf(n,np)*tkb
     &             -(1.-sigf(n,np))*tg1d(n,np)
         if(dlstaf(n,np).le.0) then
          dthdz = (1.-sigf(n,np))*tg1d(n,np)+sigf(n,np)*tkb-ts1d(n,np)
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
         zatild=(z1(n,np)-displa(lveg(n,np)))*sigf(n,np)
     &                              +z1(n,np)*(1.-sigf(n,np))
         rib1(n,np) = c(54)*zatild / (ribd(n,np)*ts1d(n,np))
         rib(n,np)  = rib1(n,np)*dlstaf(n,np)
         if (rib(n,np) .lt. 0.) then
          cdr(n,np) = cdrn(n,np)*(1.+24.5*dsqrt(-cdrn(n,np)*rib(n,np)))
          sqrtf = dmin1(dsqrt(-cdrn(n,np)/rib(n,np)),11.5d0/12.25d0)
          cdrd(n,np) = cdrn(n,np)*12.25*wtl0(n,np)*rib1(n,np)
     &                           *sigf(n,np)*sqrtf
         else
          ribi = 1./(1.+11.5*rib(n,np))
          cdr(n,np) = cdrn(n,np)*ribi
          cdrd(n,np)= cdr(n,np)*ribi*11.5*rib1(n,np)*wtl0(n,np)
     &                         *sigf(n,np)
          cdrmin(n,np) = dmax1(0.25*cdrn(n,np),6.d-4)
         end if
         if( (rib(n,np).ge.0.) ) then
           if( (cdr(n,np).lt.cdrmin(n,np)) ) then
             cdr(n,np)  = cdrmin(n,np)
             cdrd(n,np) = 0.
           end if
         end if
        end if
       end if
      end do
90    continue

      return
      end subroutine lfdrag
